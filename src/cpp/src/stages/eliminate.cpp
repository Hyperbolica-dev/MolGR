#include "molgr/stages/eliminate.h"

#include "molgr/stages/internal_helpers.h"

#include "molgr/utils/consts.h"
#include "molgr/utils/logger.h"
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

                    if (a1->GetSpinMultiplicity() == 2 && a3->GetSpinMultiplicity() == 0)
                    {
                        OBBond *b23 = mol.GetBond(a2, a3);
                        OBBond *b12 = mol.GetBond(a1, a2);
                        if (b23 && b12)
                        {
                            b23->SetBondOrder(b23->GetBondOrder() - 1);
                            b12->SetBondOrder(b12->GetBondOrder() + 1);
                            a1->SetSpinMultiplicity(a1->GetSpinMultiplicity() - 1);
                            a3->SetSpinMultiplicity(a3->GetSpinMultiplicity() + 1);
                            any_applied = true;
                            break;
                        }
                    }
                }
                if (!any_applied)
                    break;
            }
        }

        void EliminateNNN(OBMol &mol, int &charge)
        {
            while (true)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[#7v1+0]-[#7v2+0]-[#7v1+0]");
                if (matches.empty())
                    break;

                auto idxs = matches[0];
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

                charge += 1;
                LOG_DEBUG("[EliminateNNN] Applied");
            }
        }

        void EliminateHighPositiveChargeAtoms(OBMol &mol, int &charge)
        {
            while (true)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[*+1,*+2,*+3]-[Ov1+0,Nv2+0,Sv1+0]");
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
                LOG_DEBUG("[EliminateHighPos] Applied");
            }
        }

        void EliminateCNInDoubt(OBMol &mol, int &charge)
        {
            auto matches = molgr::utils::FindSmarts(mol, "[#6v4+0]=,#[#7v4+1,#15v4+1]");
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
                }
                LOG_DEBUG("[EliminateCN] Applied to " << count / 2 << " pairs");
            }
        }

        void EliminateCarboxyl(OBMol &mol, int &charge)
        {
            while (true)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[Ov1+0]-C=O");
                if (matches.empty())
                    break;
                OBAtom *a1 = mol.GetAtom(matches[0][0]);
                a1->SetSpinMultiplicity(a1->GetSpinMultiplicity() - 1);
                a1->SetFormalCharge(a1->GetFormalCharge() - 1);
                charge += 1;
                LOG_DEBUG("[EliminateCarboxyl] Applied");
            }
        }

        void EliminateCarbeneNeighborHeteroatom(OBMol &mol, int &charge)
        {
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
                        continue;

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
                            return;
                        }
                    }
                }
            }
        }

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
                process(8, false);
                process(7, false);
                process(6, true);
                process(6, false);
            }
        }

        void Eliminate13Dipole(OBMol &mol, int &charge)
        {
            auto matches = molgr::utils::FindSmarts(mol, "[*-1]-,=[N+0,O+0]-,=[*]");
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
                }
            }
        }

        void EliminatePositiveCharges(OBMol &mol, int &charge)
        {
            while (charge > 0)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[Nv3+0]=[Nv2+0]");
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
            }

            while (charge > 0)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[#6v3+0,#6v2+0,#1v0+0]");
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
                }
            }
        }

        void EliminateNegativeCharges(OBMol &mol, int &charge)
        {
            const std::vector<int> heteroatom_priority = {9, 8, 17, 7, 35, 54, 16, 34, 15};

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
            }

            while (charge < 0)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[#6v3+0]");
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
            }

            while (charge < 0)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[#1v0+0]");
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
            }

            while (charge < 0)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[#6v2+0,#6v1+0,#6v0+0]");
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
            }
        }
    }
}
