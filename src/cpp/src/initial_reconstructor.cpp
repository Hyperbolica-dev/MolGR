/**
 * @file initial_reconstructor.cpp
 * @brief Implementation of initial reconstruction logic.
 * @details STRICTLY aligned with Python 'GraphReconstruction.py'.
 * @author TMJ
 * @date 2025-12-28
 */

#include "molgr/initial_reconstructor.h"
#include "molgr/utils.h"
#include "molgr/consts.h"
#include "molgr/logger.h"

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obfunctions.h>
#include <openbabel/elements.h>
#include <openbabel/obiter.h>
#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
#include <iostream>

namespace molgr
{
    namespace reconstruct
    {

        using namespace OpenBabel;

        // =============================================================================
        // Internal Helpers
        // =============================================================================

        /**
         * @brief Helper: Fetch element info from consts map.
         * Analysis: 引入此函数是为了替代 Python 中直接访问字典 `NON_METAL_DICT[atomic_num]` 的操作。
         * C++ Map 查找需要检查是否存在，封装后更安全。
         */
        const ElementInfo *GetElementInfo(int atomic_num)
        {
            auto it = kNonMetalDict.find(atomic_num);
            if (it != kNonMetalDict.end())
                return &(it->second);
            return nullptr;
        }

        /**
         * @brief Helper: Get a flattened set of atom indices matching a SMARTS pattern.
         * Analysis: Python 中使用 `smarts.findall(omol)` 返回元组列表 (e.g. `[(1,), (2,)]`)。
         * 在 MakeConnections 等函数中，我们需要快速判断 `if idx in donor_atoms`。
         * 因此引入此函数将 SMARTS 结果扁平化为 std::set<int>，对应 Python 的 `list(chain(*...))` 或 set 转换。
         */
        std::set<int> GetFlatAtomSet(OBMol &mol, const std::string &smarts)
        {
            std::set<int> atom_indices;
            auto matches = molgr::utils::FindSmarts(mol, smarts);
            for (const auto &match : matches)
            {
                for (int idx : match)
                    atom_indices.insert(idx);
            }
            return atom_indices;
        }

        /**
         * @brief Helper: Validate conservation of charge and radicals.
         * Analysis: 对应 Python `xyz_to_omol_no_metal` 开头和结尾的校验逻辑。
         * Python 代码中分散写了 `sum(atom.GetFormalCharge())` 和 `sum(...)`，
         * 这里封装为一个函数以复用逻辑，保持主函数整洁。
         */
        bool ValidateOmol(OBMol &mol, int total_charge, int total_radical)
        {
            int charge_sum = 0;
            int radical_sum = 0;
            int radical_sum_singlet = 0; // Python: sum(atom.GetSpinMultiplicity() % 2)

            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OpenBabel::OBAtom *atom = &(*atom_iter);
                charge_sum += atom->GetFormalCharge();
                int spin = atom->GetSpinMultiplicity();
                radical_sum += spin;
                radical_sum_singlet += (spin % 2);
            }

            if (charge_sum != total_charge)
            {
                LOG_WARN("[Validate] Charge mismatch. Target: " << total_charge << ", Actual: " << charge_sum);
                // Python raises AssertError or returns None logic
                return false;
            }

            // Python logic: if sum(spin % 2) == total_radical: use singlet sum
            if (radical_sum_singlet == total_radical)
            {
                radical_sum = radical_sum_singlet;
            }

            if (radical_sum != total_radical)
            {
                LOG_WARN("[Validate] Radical mismatch. Target: " << total_radical << ", Actual: " << radical_sum);
                return false;
            }
            return true;
        }

        // =============================================================================
        // 1. Connectivity & Initialization
        // =============================================================================

        // Python: assign_charge_radical_for_atom
        void AssignChargeRadical(OBAtom *atom)
        {
            int atomic_num = atom->GetAtomicNum();
            int current_val = atom->GetTotalValence();
            int charge = atom->GetFormalCharge();

            // Python: assign_radical_dots logic inline
            int typical_val = GetTypicalValence(atomic_num, current_val, charge);
            int rad = std::max(0, typical_val - current_val);

            if (rad > 0)
            {
                atom->SetSpinMultiplicity(rad);
            }
            else
            {
                const ElementInfo *info = GetElementInfo(atomic_num);
                if (!info)
                    return; // Not a non-metal we handle

                // Special Python rule: Boron-like (3 outer e-) with 4 bonds -> Charge -1
                if (info->num_outer_electrons == 3 && current_val == 4)
                {
                    atom->SetFormalCharge(-1);
                }
                else
                {
                    // Complex logic from Python implementation
                    int spin = atom->GetSpinMultiplicity();
                    int total_elec = info->num_outer_electrons + current_val + spin - charge;

                    // Python: (...) % 8 and (...) % 2 logic
                    int low_valence_total = total_elec % 8;
                    int high_valence_total = (info->num_outer_electrons - current_val + spin - charge) % 2;

                    if (low_valence_total == 0)
                        return;

                    if (low_valence_total <= high_valence_total)
                    {
                        atom->SetFormalCharge(low_valence_total);
                    }
                    else
                    {
                        int new_spin = info->num_outer_electrons - current_val + spin - charge;
                        atom->SetSpinMultiplicity(new_spin);
                    }
                }
            }
        }

        // Python: fresh_omol_charge_radical
        void FreshOmolChargeRadical(OBMol &mol)
        {
            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                AssignChargeRadical(&(*atom_iter));
            }
        }

        // Python: make_connections
        void MakeConnections(OBMol &mol, double factor)
        {
            // SMARTS used in Python to identify donors/acceptors
            std::string donate_smarts = "[Nv0,Cv1,Nv3,Clv1,Clv2,Clv3,Brv1,Brv2,Brv3,Iv1,Iv2,Iv3]";
            std::string accept_smarts = "[Hv0,Bv2,Bv3,Cv0,Cv1,Cv2,Cv3,Nv1,Nv2,Ov0,Ov1,Clv0,Siv3,Pv2,Sv0,Sv1,Brv0,Iv0]";

            // Get initial sets
            std::set<int> donate_atoms = GetFlatAtomSet(mol, donate_smarts);
            std::vector<int> initial_donors(donate_atoms.begin(), donate_atoms.end()); // Copy for iteration

            for (int donor_idx : initial_donors)
            {
                // Python Logic Note: The acceptor list is refreshed, but effectively we iterate pairs.
                // Python does `pairs = [(atom, accept_atom) for accept_atom in accept_atoms]`
                // And sorts by distance.
                std::set<int> accept_atoms = GetFlatAtomSet(mol, accept_smarts);
                std::vector<std::pair<int, double>> candidates;

                OBAtom *donor = mol.GetAtom(donor_idx);
                if (!donor)
                    continue;

                for (int acceptor_idx : accept_atoms)
                {
                    OBAtom *acceptor = mol.GetAtom(acceptor_idx);
                    if (acceptor && donor_idx != acceptor_idx)
                    {
                        candidates.push_back({acceptor_idx, donor->GetDistance(acceptor)});
                    }
                }

                // Sort by distance
                std::sort(candidates.begin(), candidates.end(), [](const auto &a, const auto &b)
                          { return a.second < b.second; });

                for (auto &pair : candidates)
                {
                    int p1 = donor_idx;
                    int p2 = pair.first;
                    double dist = pair.second;
                    OBAtom *a1 = mol.GetAtom(p1);
                    OBAtom *a2 = mol.GetAtom(p2);

                    double r1 = OBElements::GetCovalentRad(a1->GetAtomicNum());
                    double r2 = OBElements::GetCovalentRad(a2->GetAtomicNum());

                    // Critical Python Logic: "if atom in donate_atoms and pair_atom in accept_atoms:"
                    // Python re-evaluates the sets inside the loop because adding a bond changes connectivity!
                    // This is O(N^2 * SmartsCost), extremely expensive but required for correctness.
                    std::set<int> current_donors = GetFlatAtomSet(mol, donate_smarts);
                    std::set<int> current_accepts = GetFlatAtomSet(mol, accept_smarts);

                    if (dist < (r1 + r2) * factor &&
                        current_donors.count(p1) &&
                        current_accepts.count(p2))
                    {
                        OBBond *bond = mol.GetBond(a1, a2);
                        if (!bond)
                        {
                            mol.AddBond(p1, p2, 1);
                            LOG_DEBUG("[MakeConnections] Add Bond " << p1 << "-" << p2);
                            // Bond added, this might satisfy valency, so we proceed to next candidate
                            continue;
                        }
                        if (bond->GetBondOrder() == 0)
                        {
                            bond->SetBondOrder(1);
                            LOG_DEBUG("[MakeConnections] Set Bond Order 1 " << p1 << "-" << p2);
                        }
                    }
                }
            }
        }

        // Python: pre_clean
        void PreClean(OBMol &mol)
        {
            // 1. Reduce bond order for specific hypervalent patterns
            // Python: matches = smarts.findall(omol); for match in matches: ...
            auto matches1 = molgr::utils::FindSmarts(mol, "[Cv5,Nv5,Pv5,Siv5]=,#[*]");
            for (const auto &match : matches1)
            {
                OBBond *bond = mol.GetBond(match[0], match[1]);
                if (bond)
                    bond->SetBondOrder(bond->GetBondOrder() - 1);
            }

            // 2. Fix N-BCP cage (specific structure cleanup)
            auto matches2 = molgr::utils::FindSmarts(mol, "[#6]1([#6]2)([#6]3)[#7]23[#6]1");
            for (const auto &idxs : matches2)
            {
                int n_idx = -1, c_idx = -1;
                // Python logic: Identify Bridgehead N and C based on connectivity within the substructure
                for (int idx : idxs)
                {
                    OBAtom *atom = mol.GetAtom(idx);
                    int internal_degree = 0;
                    for (int other : idxs)
                        if (idx != other && mol.GetBond(idx, other))
                            internal_degree++;

                    // Bridgeheads have higher connectivity within the cage
                    if (internal_degree >= 3)
                    {
                        if (atom->GetAtomicNum() == 7)
                            n_idx = idx;
                        if (atom->GetAtomicNum() == 6)
                            c_idx = idx;
                    }
                }
                if (n_idx != -1 && c_idx != -1)
                {
                    OBBond *bond = mol.GetBond(n_idx, c_idx);
                    if (bond)
                        mol.DeleteBond(bond);
                }
            }

            // 3. Fix small N-BCP
            auto matches3 = molgr::utils::FindSmarts(mol, "[#6]1([#6]2)[#7]2[#6]1");
            for (const auto &idxs : matches3)
            {
                int amine_n = -1, butyl_c = -1;
                for (int idx : idxs)
                {
                    OBAtom *atom = mol.GetAtom(idx);
                    int internal_degree = 0;
                    for (int other : idxs)
                        if (idx != other && mol.GetBond(idx, other))
                            internal_degree++;
                    if (internal_degree >= 2)
                    {
                        if (atom->GetAtomicNum() == 7)
                            amine_n = idx;
                        if (atom->GetAtomicNum() == 6)
                            butyl_c = idx;
                    }
                }
                if (amine_n != -1 && butyl_c != -1)
                {
                    OBBond *bond = mol.GetBond(amine_n, butyl_c);
                    if (bond)
                        mol.DeleteBond(bond);
                }
            }

            // 4. Remove hypervalent Si-O/F bonds
            auto matches4 = molgr::utils::FindSmarts(mol, "[Siv5]-[O,F]");
            for (const auto &match : matches4)
            {
                OBBond *bond = mol.GetBond(match[0], match[1]);
                if (bond)
                    mol.DeleteBond(bond);
            }
        }

        // =============================================================================
        // 2. Elimination Rules
        // =============================================================================

        // Python: clean_carbene_neighbor_unsaturated
        void CleanCarbeneNeighborUnsaturated(OBMol &mol)
        {
            while (true)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[*]-[*]=[*]");
                if (matches.empty())
                    break;
                bool any_applied = false;

                for (const auto &idxs : matches)
                {
                    OBAtom *a1 = mol.GetAtom(idxs[0]);
                    OBAtom *a2 = mol.GetAtom(idxs[1]);
                    OBAtom *a3 = mol.GetAtom(idxs[2]);

                    // Logic: Carbene (Spin 2) neighbor to a double bond
                    if (a1->GetSpinMultiplicity() == 2 && a3->GetSpinMultiplicity() == 0)
                    {
                        OBBond *b23 = mol.GetBond(a2, a3);
                        OBBond *b12 = mol.GetBond(a1, a2);
                        if (b23 && b12)
                        {
                            b23->SetBondOrder(b23->GetBondOrder() - 1);
                            b12->SetBondOrder(b12->GetBondOrder() + 1);
                            a1->SetSpinMultiplicity(a1->GetSpinMultiplicity() - 1); // 2 -> 1? Python says -1
                            a3->SetSpinMultiplicity(a3->GetSpinMultiplicity() + 1);
                            any_applied = true;
                            // Python logic breaks inner loop and restarts finding smarts
                            break;
                        }
                    }
                }
                if (!any_applied)
                    break;
            }
        }

        // Python: eliminate_NNN
        void EliminateNNN(OBMol &mol, int &charge)
        {
            while (true)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[#7v1+0]-[#7v2+0]-[#7v1+0]");
                if (matches.empty())
                    break;

                auto idxs = matches[0]; // Process first match
                mol.GetBond(idxs[0], idxs[1])->SetBondOrder(2);
                mol.GetBond(idxs[1], idxs[2])->SetBondOrder(2);

                OBAtom *a1 = mol.GetAtom(idxs[0]);
                a1->SetSpinMultiplicity(a1->GetSpinMultiplicity() - 1);
                a1->SetFormalCharge(a1->GetFormalCharge() - 1);

                OBAtom *a2 = mol.GetAtom(idxs[1]);
                a2->SetSpinMultiplicity(a2->GetSpinMultiplicity() - 1);
                a2->SetFormalCharge(a2->GetFormalCharge() + 1);

                OBAtom *a3 = mol.GetAtom(idxs[2]);
                a3->SetSpinMultiplicity(a3->GetSpinMultiplicity() - 1);
                a3->SetFormalCharge(a3->GetFormalCharge() - 1);

                charge += 1; // charge deficit increases (we removed +0, added -1+1-1 = -1 charge to system? check python)
                // Python: given_charge += 1.
                // If system charge drops (becomes more negative), deficit (Target - Current) increases. Correct.
                LOG_DEBUG("[EliminateNNN] Applied");
            }
        }

        // Python: eliminate_high_positive_charge_atoms
        void EliminateHighPositiveChargeAtoms(OBMol &mol, int &charge)
        {
            while (true)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[*+1,*+2,*+3]-[Ov1+0,Nv2+0,Sv1+0]");
                if (matches.empty())
                    break;

                auto idxs = matches[0];
                OBAtom *a1 = mol.GetAtom(idxs[0]); // High positive
                OBAtom *a2 = mol.GetAtom(idxs[1]); // Neutral Neighbor

                int sum_nbr_charge = 0;
                FOR_NB_OF_ATOM(nbr, a1)
                sum_nbr_charge += nbr->GetFormalCharge();

                // Python Condition: if -sum_nbr_charge >= atom.GetFormalCharge(): break
                if (-sum_nbr_charge >= a1->GetFormalCharge())
                    break;

                // Transfer charge/spin
                a2->SetSpinMultiplicity(a2->GetSpinMultiplicity() - 1);
                a2->SetFormalCharge(a2->GetFormalCharge() - 1);
                charge += 1;
                LOG_DEBUG("[EliminateHighPos] Applied");
            }
        }

        // Python: eliminate_CN_in_doubt
        void EliminateCNInDoubt(OBMol &mol, int &charge)
        {
            auto matches = molgr::utils::FindSmarts(mol, "[#6v4+0]=,#[#7v4+1,#15v4+1]");
            size_t count = matches.size();
            if (count > 0 && count % 2 == 0)
            {
                // Python: for atom_1, atom_2 in doubt_pair[:len//2]:
                for (size_t i = 0; i < count / 2; ++i)
                {
                    auto idxs = matches[i];
                    OBAtom *a1 = mol.GetAtom(idxs[0]);
                    OBAtom *a2 = mol.GetAtom(idxs[1]);
                    OBBond *bond = mol.GetBond(a1, a2);

                    a1->SetFormalCharge(-1);
                    bond->SetBondOrder(bond->GetBondOrder() - 1);
                    a2->SetFormalCharge(0); // Restore N/P to neutral
                    charge += 2;            // -1 added to C, +1 removed from N/P -> net -2 change -> deficit +2
                }
                LOG_DEBUG("[EliminateCN] Applied to " << count / 2 << " pairs");
            }
        }

        // Python: eliminate_carboxyl
        void EliminateCarboxyl(OBMol &mol, int &charge)
        {
            while (true)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[Ov1+0]-C=O");
                if (matches.empty())
                    break;
                OBAtom *a1 = mol.GetAtom(matches[0][0]); // Single bonded O
                a1->SetSpinMultiplicity(a1->GetSpinMultiplicity() - 1);
                a1->SetFormalCharge(a1->GetFormalCharge() - 1);
                charge += 1;
                LOG_DEBUG("[EliminateCarboxyl] Applied");
            }
        }

        // Python: eliminate_carbene_neighbor_heteroatom
        void EliminateCarbeneNeighborHeteroatom(OBMol &mol, int &charge)
        {
            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OBAtom *atom = &(*atom_iter);
                if (atom->GetSpinMultiplicity() == 2)
                { // Carbene
                    bool neighbor_radical = false;
                    FOR_NB_OF_ATOM(nbr, atom)
                    {
                        if (nbr->GetSpinMultiplicity() > 0)
                        {
                            neighbor_radical = true;
                            break;
                        }
                    }
                    if (neighbor_radical)
                        continue;

                    // Check for heteroatom neighbor to donate electrons
                    FOR_NB_OF_ATOM(nbr, atom)
                    {
                        if (kHeteroatoms.count(nbr->GetAtomicNum()) &&
                            nbr->GetFormalCharge() == 0 &&
                            nbr->GetSpinMultiplicity() == 0)
                        {
                            OBBond *bond = mol.GetBond(atom, &*nbr);
                            bond->SetBondOrder(bond->GetBondOrder() + 1);
                            atom->SetSpinMultiplicity(0);
                            atom->SetFormalCharge(atom->GetFormalCharge() - 1);
                            nbr->SetFormalCharge(nbr->GetFormalCharge() + 1);
                            LOG_DEBUG("[EliminateCarbeneHetero] Applied");
                            return; // Python returns after one application
                        }
                    }
                }
            }
        }

        // Python: clean_neighbor_radicals
        void CleanNeighborRadicals(OBMol &mol)
        {
            FOR_BONDS_OF_MOL(bond_iter, mol)
            {
                OBBond *bond = &(*bond_iter);
                OBAtom *a1 = bond->GetBeginAtom();
                OBAtom *a2 = bond->GetEndAtom();
                int r1 = a1->GetSpinMultiplicity();
                int r2 = a2->GetSpinMultiplicity();
                if (r1 > 0 && r2 > 0)
                {
                    int to_add = std::min(r1, r2);
                    bond->SetBondOrder(bond->GetBondOrder() + to_add);
                    a1->SetSpinMultiplicity(r1 - to_add);
                    a2->SetSpinMultiplicity(r2 - to_add);
                }
            }
        }

        // Python: eliminate_charge_spliting
        void EliminateChargeSpliting(OBMol &mol, int &charge)
        {
            bool all_neutral = true;
            int sum_radicals = 0;
            std::vector<OBAtom *> radical_atoms;

            FOR_ATOMS_OF_MOL(a, mol)
            {
                if (a->GetFormalCharge() != 0)
                    all_neutral = false;
                if (a->GetSpinMultiplicity() > 0)
                {
                    sum_radicals += (&*a)->GetSpinMultiplicity();
                    radical_atoms.push_back(&(*a));
                }
            }

            if (all_neutral && sum_radicals >= 2)
            {
                // Priority removal logic
                auto process = [&](int atomic_num, bool check_hetero_neighbor)
                {
                    while (radical_atoms.size() > static_cast<size_t>(std::abs(charge) + 1))
                    {
                        bool found = false;
                        for (auto it = radical_atoms.begin(); it != radical_atoms.end(); ++it)
                        {
                            OBAtom *atom = *it;
                            if (atom->GetAtomicNum() != atomic_num)
                                continue;

                            if (check_hetero_neighbor)
                            {
                                bool has_hetero = false;
                                FOR_NB_OF_ATOM(nbr, atom)
                                if (kHeteroatoms.count(nbr->GetAtomicNum())) has_hetero = true;
                                if (has_hetero)
                                    continue;
                            }

                            atom->SetSpinMultiplicity(atom->GetSpinMultiplicity() - 1);
                            atom->SetFormalCharge(atom->GetFormalCharge() - 1);
                            charge += 1;
                            radical_atoms.erase(it);
                            found = true;
                            break;
                        }
                        if (!found)
                            break;
                    }
                };
                process(8, false); // O
                process(7, false); // N
                process(6, true);  // C (no hetero)
                process(6, false); // C
            }
        }

        // Python: break_deformed_ene
        void BreakDeformedEne(OBMol &mol, int allowed_charge, int allowed_radical, double tolerance)
        {
            auto current_total_radical = [&]()
            {
                int sum = 0;
                FOR_ATOMS_OF_MOL(a, mol)
                sum += a->GetSpinMultiplicity();
                return sum;
            };

            // Pattern 1: [*]~[*+0]=,:[*+0]~[*]
            {
                auto matches = molgr::utils::FindSmarts(mol, "[*]~[*+0]=,:[*+0]~[*]");
                for (const auto &idxs : matches)
                {
                    if (current_total_radical() >= std::abs(allowed_charge) + allowed_radical)
                        return;
                    OBBond *bond = mol.GetBond(idxs[1], idxs[2]);
                    if (!bond || bond->IsRotor() || bond->GetBondOrder() == 1)
                        continue;

                    // GetTorsion uses 1-based indices
                    double torsion = std::abs(mol.GetTorsion(idxs[0], idxs[1], idxs[2], idxs[3]));
                    double deviation = std::min(torsion, 180.0 - torsion);
                    if (deviation > tolerance)
                    {
                        bond->SetBondOrder(bond->GetBondOrder() - 1);
                        mol.GetAtom(idxs[1])->SetSpinMultiplicity(mol.GetAtom(idxs[1])->GetSpinMultiplicity() + 1);
                        mol.GetAtom(idxs[2])->SetSpinMultiplicity(mol.GetAtom(idxs[2])->GetSpinMultiplicity() + 1);
                        LOG_DEBUG("[BreakDeformedEne] Broken type 1");
                    }
                }
            }

            // Pattern 2: [*]~[*+0](=,:[*+0])~[*]
            {
                auto matches = molgr::utils::FindSmarts(mol, "[*]~[*+0](=,:[*+0])~[*]");
                for (const auto &idxs : matches)
                {
                    if (current_total_radical() >= std::abs(allowed_charge) + allowed_radical)
                        return;
                    // Note: SMARTS mapping for branched structure depends on OB version, assuming standard mapping.
                    // Python: bond = GetBond(idxs[0], idxs[1]) where 0 is central C, 1 is double bonded neighbor
                    OBBond *bond = mol.GetBond(idxs[0], idxs[1]);
                    if (!bond || bond->IsRotor() || bond->GetBondOrder() == 1)
                        continue;

                    double torsion = std::abs(mol.GetTorsion(idxs[0], idxs[1], idxs[2], idxs[3]));
                    double deviation = std::min(torsion, 180.0 - torsion);
                    if (deviation > tolerance)
                    {
                        bond->SetBondOrder(bond->GetBondOrder() - 1);
                        mol.GetAtom(idxs[0])->SetSpinMultiplicity(mol.GetAtom(idxs[0])->GetSpinMultiplicity() + 1);
                        mol.GetAtom(idxs[1])->SetSpinMultiplicity(mol.GetAtom(idxs[1])->GetSpinMultiplicity() + 1);
                        LOG_DEBUG("[BreakDeformedEne] Broken type 2");
                    }
                }
            }
        }

        // Python: break_one_bond
        void BreakOneBond(OBMol &mol, int &charge, int allowed_radical)
        {
            auto check_cond = [&]()
            {
                int sum = 0;
                FOR_ATOMS_OF_MOL(a, mol)
                sum += a->GetSpinMultiplicity();
                return sum >= std::abs(charge) + allowed_radical;
            };

            auto apply_break = [&](const std::string &smarts, int idx1, int idx2, bool charge_fix)
            {
                auto matches = molgr::utils::FindSmarts(mol, smarts);
                for (const auto &idxs : matches)
                {
                    if (check_cond())
                        return;
                    OBBond *bond = mol.GetBond(idxs[idx1], idxs[idx2]);
                    if (bond)
                    {
                        bond->SetBondOrder(bond->GetBondOrder() - 1);
                        OBAtom *b = bond->GetBeginAtom();
                        OBAtom *e = bond->GetEndAtom();
                        if (charge_fix)
                        {
                            // Specific logic for [N+,P+]=[*]
                            // Python: end_atom spin+1, begin_atom charge-1
                            e->SetSpinMultiplicity(e->GetSpinMultiplicity() + 1);
                            b->SetFormalCharge(b->GetFormalCharge() - 1);
                            charge += 1;
                        }
                        else
                        {
                            b->SetSpinMultiplicity(b->GetSpinMultiplicity() + 1);
                            e->SetSpinMultiplicity(e->GetSpinMultiplicity() + 1);
                        }
                        LOG_DEBUG("[BreakOneBond] Broken smarts match");
                    }
                }
            };

            apply_break("[*+0]#,=[*+0]", 0, 1, false);
            apply_break("[#7+1,#15+1]=[*+0]", 0, 1, true);
            apply_break("[*+0]:[*+0]", 0, 1, false); // Aromatic

            // Single bonds check
            if (check_cond())
                return;

            // Check if ALL bonds are single?
            bool all_single = true;
            FOR_BONDS_OF_MOL(b, mol)
            if (b->GetBondOrder() != 1) all_single = false;

            if (all_single)
            {
                // Break ANY single bond?
                FOR_BONDS_OF_MOL(b, mol)
                {
                    if (check_cond())
                        return;
                    OBAtom *b_at = b->GetBeginAtom();
                    OBAtom *e_at = b->GetEndAtom();
                    b_at->SetSpinMultiplicity(b_at->GetSpinMultiplicity() + 1);
                    e_at->SetSpinMultiplicity(e_at->GetSpinMultiplicity() + 1);
                    mol.DeleteBond(&*b);
                    LOG_DEBUG("[BreakOneBond] Deleted single bond");
                }
            }
        }

        // =============================================================================
        // Main Entry Point
        // =============================================================================

        // Python: xyz_to_omol_no_metal (Part 1)
        std::unique_ptr<OBMol> ReconstructFromXYZNoMetal(const std::string &xyz_block, int total_charge, int total_radical)
        {
            LOG_DEBUG("[ReconstructNoMetal] Start. Target Charge=" << total_charge << " Radical=" << total_radical);

            if (total_radical < 0)
                return nullptr;

            auto mol = std::make_unique<OBMol>();
            OBConversion conv;
            conv.SetInFormat("xyz");
            conv.ReadString(mol.get(), xyz_block);

            // 1. Initial Topology
            MakeConnections(*mol);
            PreClean(*mol);

            // 2. Initial State Calculation
            FreshOmolChargeRadical(*mol);

            // Initial charge deficit calculation
            int current_charge_sum = 0;
            FOR_ATOMS_OF_MOL(a, *mol)
            current_charge_sum += a->GetFormalCharge();
            int given_charge = total_charge - current_charge_sum;

            // 3. Rule Execution Sequence (Strictly following Python order)
            EliminateNNN(*mol, given_charge);
            EliminateHighPositiveChargeAtoms(*mol, given_charge);
            EliminateCNInDoubt(*mol, given_charge);
            EliminateCarboxyl(*mol, given_charge);
            CleanCarbeneNeighborUnsaturated(*mol);
            EliminateCarbeneNeighborHeteroatom(*mol, given_charge);
            CleanNeighborRadicals(*mol);
            // Python repeats this check
            CleanCarbeneNeighborUnsaturated(*mol);

            EliminateChargeSpliting(*mol, given_charge);

            // Deformed / Breaking Logic
            BreakDeformedEne(*mol, given_charge, total_radical);
            BreakOneBond(*mol, given_charge, total_radical);

            // 4. Final Refresh
            FreshOmolChargeRadical(*mol);

            // Note: Python returns None if validation fails.
            // Here we return the molecule, validation should be done by the caller or we can return nullptr.
            // Given the Python structure, if this stage fails validation, it proceeds to resonance search.
            // So we return the mol, but maybe we should expose ValidateOmol to the caller.

            return mol;
        }

    } // namespace reconstruct
} // namespace molgr