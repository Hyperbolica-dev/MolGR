#include "molgr/stages/preprocess.h"

#include "molgr/stages/internal_helpers.h"

#include "molgr/utils/logger.h"
#include "molgr/utils/electrons.h"
#include "molgr/utils/smarts.h"
#include "molgr/utils/utils.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/elements.h>
#include "molgr/compat/openbabel_iter.h"

#include <algorithm>

namespace molgr
{
    namespace reconstruct
    {
        using namespace OpenBabel;

        std::vector<int> GetFlatAtomList(OBMol &mol, molgr::smarts::PatternId pattern_id)
        {
            std::vector<int> atom_indices;
            auto matches = molgr::smarts::FindAll(mol, pattern_id);
            for (const auto &match : matches)
            {
                for (int idx : match)
                    atom_indices.push_back(idx);
            }
            return atom_indices;
        }

        // Electron bookkeeping: fixed atoms contribute only real unpaired
        // electrons to the radical budget. Resolve deferred centers deterministically
        // as triplet (2,0) or singlet (0,1); active lone pairs never count as spin.
        // Reject an impossible budget before mutating any center.
        bool ValidateOmol(OBMol &mol, int total_charge, int total_radical, bool emit_warnings)
        {
            int charge_sum = 0;
            int radical_sum = 0;
            std::vector<OpenBabel::OBAtom *> unresolved_atoms;

            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OpenBabel::OBAtom *atom = &(*atom_iter);
                if (molgr::utils::HasUnresolvedTwoElectronCenter(*atom))
                {
                    unresolved_atoms.push_back(atom);
                }
                charge_sum += atom->GetFormalCharge();
                if (!molgr::utils::HasUnresolvedTwoElectronCenter(*atom))
                {
                    radical_sum += molgr::utils::GetUnpairedElectronCount(*atom);
                }
            }

            // Resolve the specific open-shell anion form ``C-`` + one radical
            // only when both global budgets require it. Generic charge stages
            // intentionally never partially consume an unresolved marker.
            if (total_charge - charge_sum == -1 &&
                !unresolved_atoms.empty() &&
                total_radical - radical_sum == 1)
            {
                auto *atom = unresolved_atoms.front();
                molgr::utils::SetUnresolvedTwoElectronCenter(*atom, false);
                molgr::utils::SetUnpairedElectronCount(*atom, 1);
                molgr::utils::SetLonePairCount(*atom, 0);
                atom->SetFormalCharge(-1);
                --charge_sum;
                ++radical_sum;
                unresolved_atoms.erase(unresolved_atoms.begin());
            }

            if (charge_sum != total_charge)
            {
                if (emit_warnings)
                {
                    LOG_WARN("[Validate] Charge mismatch. Target: " << total_charge << ", Actual: " << charge_sum);
                }
                return false;
            }

            const int required_unpaired_electrons = total_radical - radical_sum;
            if (required_unpaired_electrons < 0 ||
                required_unpaired_electrons % 2 != 0 ||
                required_unpaired_electrons >
                    2 * static_cast<int>(unresolved_atoms.size()))
            {
                if (emit_warnings)
                {
                    LOG_WARN("[Validate] Radical budget cannot resolve "
                             << unresolved_atoms.size()
                             << " two-electron centers. Target: " << total_radical
                             << ", Fixed: " << radical_sum);
                }
                return false;
            }

            const int triplet_center_count = required_unpaired_electrons / 2;
            for (std::size_t index = 0; index < unresolved_atoms.size(); ++index)
            {
                auto &atom = *unresolved_atoms[index];
                molgr::utils::SetUnresolvedTwoElectronCenter(atom, false);
                if (static_cast<int>(index) < triplet_center_count)
                {
                    molgr::utils::SetUnpairedElectronCount(atom, 2);
                    molgr::utils::SetLonePairCount(atom, 0);
                }
                else
                {
                    molgr::utils::SetUnpairedElectronCount(atom, 0);
                    molgr::utils::SetLonePairCount(atom, 1);
                }
            }
            return true;
        }

        bool MakeConnections(OBMol &mol, double extra_tolerance_angstrom)
        {
            bool hit = false;
            while (true)
            {
                std::vector<int> donate_atoms = GetFlatAtomList(mol, molgr::smarts::PatternId::PREPROCESS_DONATE);
                std::vector<int> accept_atoms = GetFlatAtomList(mol, molgr::smarts::PatternId::PREPROCESS_ACCEPT);
                if (donate_atoms.empty() || accept_atoms.empty())
                {
                    break;
                }

                bool changed = false;
                for (int donor_idx : donate_atoms)
                {
                    OBAtom *donor = mol.GetAtom(donor_idx);
                    if (!donor)
                    {
                        continue;
                    }

                    std::vector<std::pair<int, int>> pairs;
                    for (int acceptor_idx : accept_atoms)
                    {
                        OBAtom *acceptor = mol.GetAtom(acceptor_idx);
                        if (acceptor && acceptor_idx != donor_idx)
                        {
                            pairs.push_back({donor_idx, acceptor_idx});
                        }
                    }

                    std::sort(pairs.begin(), pairs.end(), [&mol](const auto &a, const auto &b)
                              {
                        OBAtom *a1 = mol.GetAtom(a.first);
                        OBAtom *a2 = mol.GetAtom(a.second);
                        OBAtom *b1 = mol.GetAtom(b.first);
                        OBAtom *b2 = mol.GetAtom(b.second);
                        return a1->GetDistance(a2) < b1->GetDistance(b2); });

                    for (const auto &pair : pairs)
                    {
                        int p1 = pair.first;
                        int p2 = pair.second;
                        OBAtom *a1 = mol.GetAtom(p1);
                        OBAtom *a2 = mol.GetAtom(p2);
                        if (!a1 || !a2)
                        {
                            continue;
                        }
                        double dist = a1->GetDistance(a2);

                        double r1 = OBElements::GetCovalentRad(a1->GetAtomicNum());
                        double r2 = OBElements::GetCovalentRad(a2->GetAtomicNum());

                        if (dist >= (r1 + r2) + extra_tolerance_angstrom)
                        {
                            continue;
                        }

                        OBBond *bond = mol.GetBond(a1, a2);
                        if (!bond)
                        {
                            mol.AddBond(p1, p2, 1);
                            hit = true;
                            changed = true;
                            LOG_DEBUG("[MakeConnections] Add Bond " << p1 << "-" << p2);
                            break;
                        }
                        if (bond->GetBondOrder() == 0)
                        {
                            bond->SetBondOrder(1);
                            hit = true;
                            changed = true;
                            LOG_DEBUG("[MakeConnections] Set Bond Order 1 " << p1 << "-" << p2);
                            break;
                        }
                    }
                    if (changed)
                    {
                        break;
                    }
                }
                if (!changed)
                {
                    break;
                }
            }
            return hit;
        }

        bool PreClean(OBMol &mol)
        {
            bool hit = false;
            while (true)
            {
                auto matches1 = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::PRE_CLEAN_HYPERVALENT);
                if (matches1.empty())
                {
                    break;
                }
                const auto &match = matches1.front();
                OBBond *bond = mol.GetBond(match[0], match[1]);
                if (bond)
                {
                    bond->SetBondOrder(bond->GetBondOrder() - 1);
                    hit = true;
                }
            }

            while (true)
            {
                auto matches2 = molgr::smarts::FindAll(
                    mol,
                    molgr::smarts::PatternId::PRE_CLEAN_HYPER_PI_BOND);
                if (matches2.empty())
                {
                    break;
                }
                const auto &match = matches2.front();
                OBBond *bond = mol.GetBond(match[0], match[1]);
                if (bond)
                {
                    bond->SetBondOrder(1);
                    hit = true;
                }
            }

            while (true)
            {
                auto matches3 = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::PRE_CLEAN_BCP_RING_5);
                if (matches3.empty())
                {
                    break;
                }
                const auto &idxs = matches3.front();
                int n_idx = -1, c_idx = -1;
                for (int idx : idxs)
                {
                    OBAtom *atom = mol.GetAtom(idx);
                    int internal_degree = 0;
                    for (int other : idxs)
                        if (idx != other && mol.GetBond(idx, other))
                            internal_degree++;

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
                    {
                        mol.DeleteBond(bond);
                        hit = true;
                    }
                }
            }

            while (true)
            {
                auto matches4 = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::PRE_CLEAN_BCP_RING_4);
                if (matches4.empty())
                {
                    break;
                }
                const auto &idxs = matches4.front();
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
                    {
                        mol.DeleteBond(bond);
                        hit = true;
                    }
                }
            }

            while (true)
            {
                auto matches5 = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::PRE_CLEAN_SI_O_F);
                if (matches5.empty())
                {
                    break;
                }
                const auto &match = matches5.front();
                OBBond *bond = mol.GetBond(match[0], match[1]);
                if (bond)
                {
                    mol.DeleteBond(bond);
                    hit = true;
                }
            }
            return hit;
        }
    }
}
