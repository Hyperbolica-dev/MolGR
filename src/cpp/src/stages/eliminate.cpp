#include "molgr/stages/eliminate.h"

#include "molgr/stages/internal_helpers.h"

#include "molgr/utils/consts.h"
#include "molgr/utils/logger.h"
#include "molgr/utils/smarts.h"
#include "molgr/utils/utils.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include "molgr/compat/openbabel_iter.h"

#include <algorithm>
#include <set>
#include <tuple>
#include <vector>

namespace molgr
{
    namespace reconstruct
    {
        using namespace OpenBabel;

        struct ChargeAssignmentAction
        {
            int atom_idx;
            int formal_charge;
            int spin_consumed;
            int charge_delta;
            std::vector<int> score_key;
        };

        struct NegativeChargeAssignmentPattern
        {
            molgr::smarts::PatternId pattern_id;
            int tier;
            int target_idx;
            bool requires_negative_deficit;
        };

        const std::vector<NegativeChargeAssignmentPattern> &NegativeChargeAssignmentPatterns()
        {
            static const std::vector<NegativeChargeAssignmentPattern> patterns = {
                {molgr::smarts::PatternId::ELIM_NEGATIVE_F, 10, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_O, 20, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_O_1, 21, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_CL, 30, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_N, 40, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_N_1, 41, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_N_2, 42, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_BR, 50, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_I, 60, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_S, 70, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_S_1, 71, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_SE, 80, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_SE_1, 81, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_P, 90, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_P_1, 91, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_P_2, 92, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_B, 95, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_B_1, 96, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_B_2, 97, 0, false},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_C_V3, 100, 0, true},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_C_LOW, 110, 0, true},
                {molgr::smarts::PatternId::ELIM_NEGATIVE_H, 120, 0, true},
            };
            return patterns;
        }

        int NegativeChargeAssignmentAmount(OBAtom *atom, int charge)
        {
            if (atom == nullptr)
            {
                return 0;
            }
            return std::min(atom->GetSpinMultiplicity(), std::max(1, std::abs(charge)));
        }

        bool ApplyChargeAssignmentAction(OBMol &mol, const ChargeAssignmentAction &action)
        {
            OBAtom *atom = mol.GetAtom(action.atom_idx);
            if (atom == nullptr || action.spin_consumed <= 0)
            {
                return false;
            }
            if (atom->GetSpinMultiplicity() < action.spin_consumed)
            {
                return false;
            }
            atom->SetSpinMultiplicity(atom->GetSpinMultiplicity() - action.spin_consumed);
            atom->SetFormalCharge(action.formal_charge);
            return true;
        }

        void AppendPositiveChargeAssignmentAction(
            std::vector<ChargeAssignmentAction> &actions,
            std::set<std::tuple<int, int, int>> &seen,
            OBAtom *atom,
            int charge,
            int tier,
            int match_order,
            int amount)
        {
            if (atom == nullptr || amount <= 0)
            {
                return;
            }
            const int atom_idx = atom->GetIdx();
            const auto seen_key = std::make_tuple(atom_idx, tier, amount);
            if (seen.count(seen_key) != 0)
            {
                return;
            }
            seen.insert(seen_key);

            const int charge_after = charge - amount;
            const int atomic_num = atom->GetAtomicNum();
            actions.push_back(ChargeAssignmentAction{
                atom_idx,
                amount,
                amount,
                -amount,
                {
                    tier,
                    std::abs(charge_after),
                    std::max(charge_after, 0),
                    atomic_num,
                    atom_idx,
                    match_order,
                }});
        }

        void AppendNegativeChargeAssignmentAction(
            std::vector<ChargeAssignmentAction> &actions,
            std::set<std::tuple<int, int, int>> &seen,
            OBAtom *atom,
            int charge,
            int tier,
            int match_order,
            int amount)
        {
            if (atom == nullptr || amount <= 0)
            {
                return;
            }
            const int atom_idx = atom->GetIdx();
            const auto seen_key = std::make_tuple(atom_idx, tier, amount);
            if (seen.count(seen_key) != 0)
            {
                return;
            }
            seen.insert(seen_key);

            const int charge_after = charge + amount;
            const int atomic_num = atom->GetAtomicNum();
            actions.push_back(ChargeAssignmentAction{
                atom_idx,
                -amount,
                amount,
                amount,
                {
                    tier,
                    std::abs(charge_after),
                    std::max(charge_after, 0),
                    atomic_num,
                    atom_idx,
                    match_order,
                }});
        }

        std::vector<ChargeAssignmentAction> PositiveChargeAssignmentActions(
            OBMol &mol,
            int charge)
        {
            std::vector<ChargeAssignmentAction> actions;
            if (charge <= 0)
            {
                return actions;
            }

            std::set<std::tuple<int, int, int>> seen;
            auto n_matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::ELIM_POSITIVE_N);
            for (std::size_t match_order = 0; match_order < n_matches.size(); ++match_order)
            {
                const auto &idxs = n_matches[match_order];
                if (idxs.size() < 2)
                {
                    continue;
                }
                OBAtom *atom = mol.GetAtom(idxs[1]);
                if (atom != nullptr && atom->GetFormalCharge() == 0 && atom->GetSpinMultiplicity() >= 1)
                {
                    AppendPositiveChargeAssignmentAction(
                        actions,
                        seen,
                        atom,
                        charge,
                        0,
                        static_cast<int>(match_order),
                        1);
                }
            }

            auto c_h_matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::ELIM_POSITIVE_C_H);
            for (std::size_t match_order = 0; match_order < c_h_matches.size(); ++match_order)
            {
                const auto &idxs = c_h_matches[match_order];
                if (idxs.empty())
                {
                    continue;
                }
                OBAtom *atom = mol.GetAtom(idxs[0]);
                if (atom != nullptr && atom->GetFormalCharge() == 0 && atom->GetSpinMultiplicity() >= 1)
                {
                    AppendPositiveChargeAssignmentAction(
                        actions,
                        seen,
                        atom,
                        charge,
                        10,
                        static_cast<int>(match_order),
                        1);
                }
            }

            int match_order = 0;
            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OBAtom *atom = &(*atom_iter);
                if (atom->GetFormalCharge() == 0 && atom->GetSpinMultiplicity() >= 1)
                {
                    AppendPositiveChargeAssignmentAction(
                        actions,
                        seen,
                        atom,
                        charge,
                        100,
                        match_order,
                        std::min(atom->GetSpinMultiplicity(), charge));
                }
                ++match_order;
            }

            std::sort(
                actions.begin(),
                actions.end(),
                [](const auto &left, const auto &right)
                {
                    return left.score_key < right.score_key;
                });
            return actions;
        }

        std::vector<ChargeAssignmentAction> NegativeChargeAssignmentActions(
            OBMol &mol,
            int charge)
        {
            std::vector<ChargeAssignmentAction> actions;
            if (charge > 0)
            {
                return actions;
            }

            std::set<std::tuple<int, int, int>> seen;

            for (const auto &pattern : NegativeChargeAssignmentPatterns())
            {
                if (pattern.requires_negative_deficit && charge >= 0)
                {
                    continue;
                }
                auto matches = molgr::smarts::FindAll(mol, pattern.pattern_id);
                for (std::size_t match_order = 0; match_order < matches.size(); ++match_order)
                {
                    const auto &idxs = matches[match_order];
                    if (idxs.size() <= static_cast<std::size_t>(pattern.target_idx))
                    {
                        continue;
                    }
                    OBAtom *atom = mol.GetAtom(idxs[pattern.target_idx]);
                    if (atom != nullptr && atom->GetFormalCharge() == 0 && atom->GetSpinMultiplicity() >= 1)
                    {
                        AppendNegativeChargeAssignmentAction(
                            actions,
                            seen,
                            atom,
                            charge,
                            pattern.tier,
                            static_cast<int>(match_order),
                            NegativeChargeAssignmentAmount(atom, charge));
                    }
                }
            }

            if (charge < 0)
            {
                int match_order = 0;
                FOR_ATOMS_OF_MOL(atom_iter, mol)
                {
                    OBAtom *atom = &(*atom_iter);
                    if (atom->GetFormalCharge() == 0 && atom->GetSpinMultiplicity() >= 1)
                    {
                        AppendNegativeChargeAssignmentAction(
                            actions,
                            seen,
                            atom,
                            charge,
                            1000,
                            match_order,
                            NegativeChargeAssignmentAmount(atom, charge));
                    }
                    ++match_order;
                }
            }

            std::sort(
                actions.begin(),
                actions.end(),
                [](const auto &left, const auto &right)
                {
                    return left.score_key < right.score_key;
                });
            return actions;
        }

        bool EliminateNNN(OBMol &mol, int &charge, bool positive)
        {
            bool hit = false;
            if (!positive)
            {
                while (true)
                {
                    auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::ELIM_NNN_NEGATIVE);
                    if (matches.empty())
                    {
                        break;
                    }

                    const auto idxs = matches[0];
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    OBAtom *a1 = mol.GetAtom(idxs[0]);
                    OBAtom *a2 = mol.GetAtom(idxs[1]);
                    OBAtom *a3 = mol.GetAtom(idxs[2]);
                    if (!bond1 || !bond2 || !a1 || !a2 || !a3)
                    {
                        break;
                    }

                    bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                    bond2->SetBondOrder(bond2->GetBondOrder() + 1);

                    a1->SetSpinMultiplicity(a1->GetSpinMultiplicity() - 2);
                    a1->SetFormalCharge(a1->GetFormalCharge() - 1);

                    a2->SetSpinMultiplicity(a2->GetSpinMultiplicity() - 1);
                    a2->SetFormalCharge(a2->GetFormalCharge() + 1);

                    a3->SetSpinMultiplicity(a3->GetSpinMultiplicity() - 2);
                    a3->SetFormalCharge(a3->GetFormalCharge() - 1);

                    charge += 1;
                    hit = true;
                    LOG_DEBUG("[EliminateNNN] Applied negative branch");
                }
                return hit;
            }

            while (true)
            {
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::ELIM_NNN_POSITIVE);
                if (matches.empty())
                {
                    break;
                }

                const auto idxs = matches[0];
                OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                OBAtom *a1 = mol.GetAtom(idxs[0]);
                OBAtom *a2 = mol.GetAtom(idxs[1]);
                if (!bond1 || !a1 || !a2)
                {
                    break;
                }

                bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                a1->SetFormalCharge(a1->GetFormalCharge() + 1);
                a2->SetSpinMultiplicity(a2->GetSpinMultiplicity() - 1);
                charge -= 1;
                hit = true;
                LOG_DEBUG("[EliminateNNN] Applied positive branch");
            }
            return hit;
        }

        bool EliminateHighPositiveChargeAtoms(OBMol &mol, int &charge)
        {
            bool hit = false;
            auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::ELIM_HIGH_POSITIVE);
            while (!matches.empty())
            {
                auto idxs = matches.front();
                matches.erase(matches.begin());
                OBAtom *a1 = mol.GetAtom(idxs[0]);
                OBAtom *a2 = mol.GetAtom(idxs[1]);
                if (!a1 || !a2)
                {
                    continue;
                }

                int sum_nbr_charge = 0;
                FOR_NB_OF_ATOM(nbr, a1)
                sum_nbr_charge += nbr->GetFormalCharge();

                if (-sum_nbr_charge > a1->GetFormalCharge() || a2->GetSpinMultiplicity() != 1)
                {
                    continue;
                }

                a2->SetSpinMultiplicity(a2->GetSpinMultiplicity() - 1);
                a2->SetFormalCharge(a2->GetFormalCharge() - 1);
                charge += 1;
                hit = true;
                LOG_DEBUG("[EliminateHighPos] Applied");
            }
            return hit;
        }

        bool EliminateCNInDoubt(OBMol &mol, int &charge)
        {
            bool hit = false;
            auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::ELIM_CN_IN_DOUBT);
            size_t count = matches.size();
            std::vector<int> atom_indices;
            atom_indices.reserve(count * 2);
            for (const auto &idxs : matches)
            {
                atom_indices.push_back(idxs[0]);
                atom_indices.push_back(idxs[1]);
            }
            std::sort(atom_indices.begin(), atom_indices.end());
            if (std::unique(atom_indices.begin(), atom_indices.end()) != atom_indices.end())
            {
                return false;
            }
            if (count > 0 && count % 2 == 0)
            {
                for (size_t i = 0; i < count / 2; ++i)
                {
                    auto idxs = matches[i];
                    OBAtom *a1 = mol.GetAtom(idxs[0]);
                    OBAtom *a2 = mol.GetAtom(idxs[1]);
                    OBBond *bond = mol.GetBond(a1, a2);

                    a1->SetFormalCharge(-1);
                    bond->SetBondOrder(bond->GetBondOrder() - 1);
                    a2->SetFormalCharge(0);
                    charge += 2;
                    hit = true;
                }
                LOG_DEBUG("[EliminateCN] Applied to " << count / 2 << " pairs");
            }
            return hit;
        }

        bool EliminateCarboxyl(OBMol &mol, int &charge)
        {
            bool hit = false;
            while (true)
            {
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::ELIM_CARBOXYL);
                if (matches.empty())
                    break;
                OBAtom *a1 = mol.GetAtom(matches[0][0]);
                a1->SetSpinMultiplicity(a1->GetSpinMultiplicity() - 1);
                a1->SetFormalCharge(a1->GetFormalCharge() - 1);
                charge += 1;
                hit = true;
                LOG_DEBUG("[EliminateCarboxyl] Applied");
            }
            return hit;
        }

        bool EliminateCarbeneNeighborHeteroatom(OBMol &mol, int &charge)
        {
            const auto possible_carbene_atom = [](const OBAtom &atom) -> bool
            {
                const int atomic_num = atom.GetAtomicNum();
                if (kHeteroatoms.count(atomic_num) == 0 && atomic_num != 6)
                {
                    return false;
                }
                if (atom.GetSpinMultiplicity() != 2)
                {
                    return false;
                }
                return atomic_num != 8 && atomic_num != 9 && atomic_num != 16 &&
                       atomic_num != 17 && atomic_num != 35 && atomic_num != 53;
            };

            bool hit = false;
            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OBAtom *atom = &(*atom_iter);
                if (possible_carbene_atom(*atom))
                {
                    bool neighbor_radical = false;
                    FOR_NB_OF_ATOM(nbr, atom)
                    {
                        if (nbr->GetSpinMultiplicity() == 1)
                        {
                            neighbor_radical = true;
                            break;
                        }
                    }
                    if (neighbor_radical)
                    {
                        continue;
                    }

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
                            hit = true;
                            LOG_DEBUG("[EliminateCarbeneHetero] Applied");
                            break;
                        }
                    }
                }
            }
            return hit;
        }

        bool EliminateChargeSpliting(OBMol &mol, int &charge)
        {
            bool hit = false;
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
                    if (a->GetSpinMultiplicity() == 1)
                    {
                        radical_atoms.push_back(&(*a));
                    }
                }
            }

            if (all_neutral && sum_radicals >= 2)
            {
                int total_radicals = 0;
                for (OBAtom *atom : radical_atoms)
                {
                    total_radicals += atom->GetSpinMultiplicity();
                }
                auto process = [&](int atomic_num, bool check_hetero_neighbor)
                {
                    while (total_radicals > std::abs(charge))
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

                            const int charge_delta = atom->GetSpinMultiplicity();
                            atom->SetFormalCharge(atom->GetFormalCharge() - charge_delta);
                            atom->SetSpinMultiplicity(atom->GetSpinMultiplicity() - charge_delta);
                            charge += charge_delta;
                            total_radicals -= charge_delta;
                            radical_atoms.erase(it);
                            hit = true;
                            found = true;
                            break;
                        }
                        if (!found)
                            break;
                    }
                };
                for (const int atomic_num : {8, 9, 17, 35, 53})
                {
                    process(atomic_num, false);
                }
                process(16, false);
                process(7, false);
                process(6, true);
                process(6, false);
            }
            return hit;
        }

        bool Eliminate13Dipole(OBMol &mol, int &charge)
        {
            bool hit = false;
            auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::ELIM_1_3_DIPOLE);
            while (!matches.empty())
            {
                auto idxs = matches.front();
                matches.erase(matches.begin());

                OBAtom *atom2 = mol.GetAtom(idxs[1]);
                OBAtom *atom3 = mol.GetAtom(idxs[2]);
                OBBond *bond23 = mol.GetBond(idxs[1], idxs[2]);
                if (!atom2 || !atom3 || !bond23)
                {
                    continue;
                }

                const ElementInfo *info = GetElementInfo(atom2->GetAtomicNum());
                if (info != nullptr &&
                    atom3->GetSpinMultiplicity() > 0 &&
                    (info->num_outer_electrons + atom2->GetTotalValence()) == 8)
                {
                    atom2->SetFormalCharge(atom2->GetFormalCharge() + 1);
                    bond23->SetBondOrder(static_cast<int>(bond23->GetBondOrder() + 1));
                    atom3->SetSpinMultiplicity(atom3->GetSpinMultiplicity() - 1);
                    charge -= 1;
                    hit = true;
                }
            }
            return hit;
        }

        bool EliminateCPLikeRadicalAnion(OBMol &mol, int &charge)
        {
            bool hit = false;
            while (charge < 0)
            {
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::ELIM_NEGATIVE_CP);
                if (matches.empty())
                {
                    break;
                }

                auto idxs = matches.front();
                if (idxs.size() < 5)
                {
                    break;
                }
                OBAtom *atom = mol.GetAtom(idxs[4]);
                if (atom == nullptr || atom->GetFormalCharge() != 0)
                {
                    break;
                }

                const int to_add = atom->GetSpinMultiplicity();
                if (to_add <= 0)
                {
                    break;
                }
                atom->SetSpinMultiplicity(atom->GetSpinMultiplicity() - to_add);
                atom->SetFormalCharge(-to_add);
                charge += to_add;
                hit = true;
            }
            return hit;
        }

        bool EliminatePositiveCharges(OBMol &mol, int &charge)
        {
            bool hit = false;
            while (charge > 0)
            {
                const auto actions = PositiveChargeAssignmentActions(mol, charge);
                if (actions.empty())
                {
                    break;
                }
                const auto &action = actions.front();
                if (!ApplyChargeAssignmentAction(mol, action))
                {
                    break;
                }
                charge += action.charge_delta;
                hit = true;
            }
            return hit;
        }

        bool EliminateNegativeCharges(OBMol &mol, int &charge)
        {
            bool hit = false;
            while (charge <= 0)
            {
                const auto actions = NegativeChargeAssignmentActions(mol, charge);
                if (actions.empty())
                {
                    break;
                }
                const auto &action = actions.front();
                if (!ApplyChargeAssignmentAction(mol, action))
                {
                    break;
                }
                charge += action.charge_delta;
                hit = true;
            }
            return hit;
        }
    }
}
