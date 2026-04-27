#include "molgr/stages/eliminate.h"

#include "molgr/stages/internal_helpers.h"

#include "molgr/utils/consts.h"
#include "molgr/utils/logger.h"
#include "molgr/utils/smarts.h"
#include "molgr/utils/utils.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>

#include <algorithm>
#include <vector>

namespace molgr
{
    namespace reconstruct
    {
        using namespace OpenBabel;

        bool EliminateNNN(OBMol &mol, int &charge, bool positive)
        {
            bool hit = false;
            if (!positive)
            {
                while (true)
                {
                    auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::ELIM_NNN_NEGATIVE);
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
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::ELIM_NNN_POSITIVE);
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
            while (true)
            {
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::ELIM_HIGH_POSITIVE);
                if (matches.empty())
                    break;

                auto idxs = matches[0];
                OBAtom *a1 = mol.GetAtom(idxs[0]);
                OBAtom *a2 = mol.GetAtom(idxs[1]);

                int sum_nbr_charge = 0;
                FOR_NB_OF_ATOM(nbr, a1)
                sum_nbr_charge += nbr->GetFormalCharge();

                if (-sum_nbr_charge >= a1->GetFormalCharge())
                    break;

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
            auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::ELIM_CN_IN_DOUBT);
            size_t count = matches.size();
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
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::ELIM_CARBOXYL);
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
            bool hit = false;
            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OBAtom *atom = &(*atom_iter);
                if (atom->GetSpinMultiplicity() == 2)
                {
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
                    {
                        return hit;
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
                    radical_atoms.push_back(&(*a));
                }
            }

            if (all_neutral && sum_radicals >= 2)
            {
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
                process(7, false);
                process(6, true);
                process(6, false);
            }
            return hit;
        }

        bool Eliminate13Dipole(OBMol &mol, int &charge)
        {
            bool hit = false;
            auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::ELIM_1_3_DIPOLE);
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

        bool EliminatePositiveCharges(OBMol &mol, int &charge)
        {
            bool hit = false;
            while (charge > 0)
            {
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::ELIM_POSITIVE_N);
                if (matches.empty())
                {
                    break;
                }

                auto idxs = matches.front();
                OBAtom *atom = mol.GetAtom(idxs[1]);
                if (!atom)
                {
                    break;
                }

                atom->SetSpinMultiplicity(atom->GetSpinMultiplicity() - 1);
                atom->SetFormalCharge(1);
                charge -= 1;
                hit = true;
            }

            while (charge > 0)
            {
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::ELIM_POSITIVE_C_H);
                if (matches.empty())
                {
                    break;
                }

                auto idxs = matches.front();
                OBAtom *atom = mol.GetAtom(idxs[0]);
                if (!atom)
                {
                    break;
                }

                atom->SetSpinMultiplicity(atom->GetSpinMultiplicity() - 1);
                atom->SetFormalCharge(1);
                charge -= 1;
                hit = true;
            }

            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OBAtom *atom = &(*atom_iter);
                if (charge <= 0)
                {
                    break;
                }
                if (atom->GetSpinMultiplicity() >= 1 && atom->GetFormalCharge() == 0)
                {
                    int to_add = std::min(atom->GetSpinMultiplicity(), charge);
                    atom->SetSpinMultiplicity(atom->GetSpinMultiplicity() - to_add);
                    atom->SetFormalCharge(to_add);
                    charge -= to_add;
                    hit = true;
                }
            }
            return hit;
        }

        bool EliminateNegativeCharges(OBMol &mol, int &charge)
        {
            bool hit = false;
            const std::vector<int> heteroatom_priority = {9, 8, 17, 7, 35, 53, 16, 34, 15};

            std::vector<std::pair<OBAtom *, size_t>> possible_heteroatoms;
            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OBAtom *atom = &(*atom_iter);
                if (atom->GetFormalCharge() != 0 || atom->GetSpinMultiplicity() < 1)
                {
                    continue;
                }

                auto pos = std::find(heteroatom_priority.begin(), heteroatom_priority.end(), atom->GetAtomicNum());
                if (pos != heteroatom_priority.end())
                {
                    possible_heteroatoms.push_back({atom, static_cast<size_t>(pos - heteroatom_priority.begin())});
                }
            }

            std::sort(possible_heteroatoms.begin(), possible_heteroatoms.end(),
                      [](const auto &a, const auto &b)
                      { return a.second < b.second; });

            for (const auto &entry : possible_heteroatoms)
            {
                if (charge >= 0)
                {
                    break;
                }
                OBAtom *atom = entry.first;
                int to_add = std::min(atom->GetSpinMultiplicity(), std::abs(charge));
                atom->SetSpinMultiplicity(atom->GetSpinMultiplicity() - to_add);
                atom->SetFormalCharge(-to_add);
                charge += to_add;
                if (to_add > 0)
                {
                    hit = true;
                }
            }

            while (charge < 0)
            {
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::ELIM_NEGATIVE_C_V3);
                if (matches.empty())
                {
                    break;
                }

                auto idxs = matches.front();
                OBAtom *atom = mol.GetAtom(idxs[0]);
                if (!atom)
                {
                    break;
                }

                int to_add = std::min(atom->GetSpinMultiplicity(), std::abs(charge));
                atom->SetSpinMultiplicity(atom->GetSpinMultiplicity() - to_add);
                atom->SetFormalCharge(-to_add);
                charge += to_add;
                if (to_add > 0)
                {
                    hit = true;
                }
            }

            while (charge < 0)
            {
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::ELIM_NEGATIVE_H);
                if (matches.empty())
                {
                    break;
                }

                auto idxs = matches.front();
                OBAtom *atom = mol.GetAtom(idxs[0]);
                if (!atom)
                {
                    break;
                }

                int to_add = std::min(atom->GetSpinMultiplicity(), std::abs(charge));
                atom->SetSpinMultiplicity(atom->GetSpinMultiplicity() - to_add);
                atom->SetFormalCharge(-to_add);
                charge += to_add;
                if (to_add > 0)
                {
                    hit = true;
                }
            }

            while (charge < 0)
            {
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::ELIM_NEGATIVE_C_LOW);
                if (matches.empty())
                {
                    break;
                }

                auto idxs = matches.front();
                OBAtom *atom = mol.GetAtom(idxs[0]);
                if (!atom)
                {
                    break;
                }

                int to_add = std::min(atom->GetSpinMultiplicity(), std::abs(charge));
                atom->SetSpinMultiplicity(atom->GetSpinMultiplicity() - to_add);
                atom->SetFormalCharge(-to_add);
                charge += to_add;
                if (to_add > 0)
                {
                    hit = true;
                }
            }

            while (charge < 0)
            {
                bool updated = false;
                FOR_ATOMS_OF_MOL(atom_iter, mol)
                {
                    OBAtom *atom = &(*atom_iter);
                    if (atom->GetSpinMultiplicity() < 1 || atom->GetFormalCharge() != 0)
                    {
                        continue;
                    }

                    const int to_add = std::min(atom->GetSpinMultiplicity(), std::abs(charge));
                    atom->SetSpinMultiplicity(atom->GetSpinMultiplicity() - to_add);
                    atom->SetFormalCharge(-to_add);
                    charge += to_add;
                    if (to_add > 0)
                    {
                        hit = true;
                        updated = true;
                    }
                    if (charge >= 0)
                    {
                        break;
                    }
                }
                if (!updated)
                {
                    break;
                }
            }
            return hit;
        }
    }
}
