#include "molgr/stages/clean.h"

#include "molgr/stages/fresh.h"
#include "molgr/stages/internal_helpers.h"

#include "molgr/utils/utils.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>

#include <algorithm>

namespace molgr
{
    namespace reconstruct
    {
        using namespace OpenBabel;

        namespace
        {
            int TotalValenceForRoomCheck(const OBAtom *atom)
            {
                if (atom == nullptr)
                {
                    return 0;
                }
                return atom->GetExplicitValence() + atom->GetImplicitHCount();
            }

            void CleanResonances0(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[*-]-[*]=[*]~[*+]");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom4 = mol.GetAtom(idxs[3]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    OBBond *bond3 = mol.GetBond(idxs[2], idxs[3]);
                    if (!atom1 || !atom4 || !bond1 || !bond2 || !bond3)
                    {
                        continue;
                    }

                    const ElementInfo *info1 = GetElementInfo(atom1->GetAtomicNum());
                    const ElementInfo *info4 = GetElementInfo(atom4->GetAtomicNum());
                    if (info1 != nullptr && info4 != nullptr &&
                        info1->default_valence > atom1->GetTotalValence() &&
                        info4->default_valence > atom4->GetTotalValence())
                    {
                        bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                        bond2->SetBondOrder(bond2->GetBondOrder() - 1);
                        bond3->SetBondOrder(bond3->GetBondOrder() + 1);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() + 1);
                        atom4->SetFormalCharge(atom4->GetFormalCharge() - 1);
                    }
                }
            }

            void CleanResonances1(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[*-]=[*+]=[*+0]");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom3 = mol.GetAtom(idxs[2]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    if (!atom1 || !atom3 || !bond1 || !bond2)
                    {
                        continue;
                    }

                    bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                    bond2->SetBondOrder(bond2->GetBondOrder() - 1);
                    atom1->SetFormalCharge(atom1->GetFormalCharge() + 1);
                    atom3->SetFormalCharge(atom3->GetFormalCharge() - 1);
                }
            }

            void CleanResonances2(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[#8]=[#6](-[!-])-[*]=[*]-[#7-,#6-]");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom6 = mol.GetAtom(idxs[5]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[3]);
                    OBBond *bond3 = mol.GetBond(idxs[3], idxs[4]);
                    OBBond *bond4 = mol.GetBond(idxs[4], idxs[5]);
                    if (!atom1 || !atom6 || !bond1 || !bond2 || !bond3 || !bond4)
                    {
                        continue;
                    }

                    if (bond4->GetBondOrder() == 1 &&
                        bond3->GetBondOrder() == 2 &&
                        bond2->GetBondOrder() == 1 &&
                        bond1->GetBondOrder() == 2 &&
                        atom1->GetFormalCharge() == 0 &&
                        atom6->GetFormalCharge() == -1)
                    {
                        bond4->SetBondOrder(bond4->GetBondOrder() + 1);
                        bond3->SetBondOrder(bond3->GetBondOrder() - 1);
                        bond2->SetBondOrder(bond2->GetBondOrder() + 1);
                        bond1->SetBondOrder(bond1->GetBondOrder() - 1);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() - 1);
                        atom6->SetFormalCharge(atom6->GetFormalCharge() + 1);
                    }
                }
            }

            void CleanResonances3(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[#7v2+]=[*]-[*]=[*]-[#8-]");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom5 = mol.GetAtom(idxs[4]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    OBBond *bond3 = mol.GetBond(idxs[2], idxs[3]);
                    OBBond *bond4 = mol.GetBond(idxs[3], idxs[4]);
                    if (!atom1 || !atom5 || !bond1 || !bond2 || !bond3 || !bond4)
                    {
                        continue;
                    }

                    bond1->SetBondOrder(bond1->GetBondOrder() - 1);
                    bond2->SetBondOrder(bond2->GetBondOrder() + 1);
                    bond3->SetBondOrder(bond3->GetBondOrder() - 1);
                    bond4->SetBondOrder(bond4->GetBondOrder() + 1);
                    atom1->SetFormalCharge(atom1->GetFormalCharge() - 1);
                    atom5->SetFormalCharge(atom5->GetFormalCharge() + 1);
                    FreshOmolChargeRadical(mol);
                }
            }

            void CleanResonances4(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[#7+,#8+]=[*]-[#6-,#7-,#8-]");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom3 = mol.GetAtom(idxs[2]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    if (!atom1 || !atom3 || !bond1 || !bond2)
                    {
                        continue;
                    }

                    bond2->SetBondOrder(bond2->GetBondOrder() + 1);
                    bond1->SetBondOrder(bond1->GetBondOrder() - 1);
                    atom1->SetFormalCharge(atom1->GetFormalCharge() - 1);
                    atom3->SetFormalCharge(atom3->GetFormalCharge() + 1);
                }
            }

            void CleanResonances5(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[#7+0,#8+0,#16+0]=[*+0]-[#6-,#7-]");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom3 = mol.GetAtom(idxs[2]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    if (!atom1 || !atom3 || !bond1 || !bond2)
                    {
                        continue;
                    }

                    if (bond1->GetBondOrder() == 2 &&
                        bond2->GetBondOrder() == 1 &&
                        atom3->GetFormalCharge() == -1 &&
                        atom1->GetFormalCharge() == 0)
                    {
                        bond2->SetBondOrder(bond2->GetBondOrder() + 1);
                        bond1->SetBondOrder(bond1->GetBondOrder() - 1);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() - 1);
                        atom3->SetFormalCharge(atom3->GetFormalCharge() + 1);
                    }
                }
            }

            void CleanResonances6(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[#6]=[#6]=[#6-,#7-]");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom3 = mol.GetAtom(idxs[2]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    if (!atom1 || !atom3 || !bond1 || !bond2)
                    {
                        continue;
                    }

                    bond2->SetBondOrder(bond2->GetBondOrder() + 1);
                    bond1->SetBondOrder(bond1->GetBondOrder() - 1);
                    atom1->SetFormalCharge(atom1->GetFormalCharge() - 1);
                    atom3->SetFormalCharge(atom3->GetFormalCharge() + 1);
                }
            }

            void CleanResonances7(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[*-]1-[*](=[*])-[*]=[*]-[*]=[*]1");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom3 = mol.GetAtom(idxs[2]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    if (!atom1 || !atom3 || !bond1 || !bond2)
                    {
                        continue;
                    }

                    bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                    bond2->SetBondOrder(bond2->GetBondOrder() - 1);
                    atom1->SetFormalCharge(atom1->GetFormalCharge() + 1);
                    atom3->SetFormalCharge(atom3->GetFormalCharge() - 1);
                }
            }

            void CleanResonances8(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[*-]1-[*]=[*]-[*](=[*])-[*]=[*]1");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom5 = mol.GetAtom(idxs[4]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    OBBond *bond3 = mol.GetBond(idxs[2], idxs[3]);
                    OBBond *bond4 = mol.GetBond(idxs[3], idxs[4]);
                    if (!atom1 || !atom5 || !bond1 || !bond2 || !bond3 || !bond4)
                    {
                        continue;
                    }

                    bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                    bond2->SetBondOrder(bond2->GetBondOrder() - 1);
                    bond3->SetBondOrder(bond3->GetBondOrder() + 1);
                    bond4->SetBondOrder(bond4->GetBondOrder() - 1);
                    atom1->SetFormalCharge(atom1->GetFormalCharge() + 1);
                    atom5->SetFormalCharge(atom5->GetFormalCharge() - 1);
                }
            }

            void CleanResonances9(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[*+,*+2,*+3]-,=[*-,*-2,*-3]");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom2 = mol.GetAtom(idxs[1]);
                    OBBond *bond = mol.GetBond(idxs[0], idxs[1]);
                    if (!atom1 || !atom2 || !bond)
                    {
                        continue;
                    }

                    const ElementInfo *info1 = GetElementInfo(atom1->GetAtomicNum());
                    const ElementInfo *info2 = GetElementInfo(atom2->GetAtomicNum());
                    if (info1 == nullptr || info2 == nullptr)
                    {
                        continue;
                    }

                    int room1 = info1->default_valence - atom1->GetTotalValence();
                    int room2 = info2->default_valence - atom2->GetTotalValence();
                    if (room1 >= 1 && room2 >= 1)
                    {
                        int bond_to_add = std::min(room1, room2);
                        bond->SetBondOrder(bond->GetBondOrder() + bond_to_add);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() - bond_to_add);
                        atom2->SetFormalCharge(atom2->GetFormalCharge() + bond_to_add);
                    }
                }
            }

            void CleanResonances10(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[*]-[*]=,#[*]-[*]");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom4 = mol.GetAtom(idxs[3]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    OBBond *bond3 = mol.GetBond(idxs[2], idxs[3]);
                    if (!atom1 || !atom4 || !bond1 || !bond2 || !bond3)
                    {
                        continue;
                    }

                    if (atom1->GetSpinMultiplicity() == 1 && atom4->GetSpinMultiplicity() == 1)
                    {
                        bond2->SetBondOrder(bond2->GetBondOrder() - 1);
                        bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                        bond3->SetBondOrder(bond3->GetBondOrder() + 1);
                        atom1->SetSpinMultiplicity(atom1->GetSpinMultiplicity() - 1);
                        atom4->SetSpinMultiplicity(atom4->GetSpinMultiplicity() - 1);
                    }
                }
            }

            void CleanResonances11(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[#7v3+0,#8v2+0,#16v2+0]-,=,:[*+1]");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom2 = mol.GetAtom(idxs[1]);
                    OBBond *bond = mol.GetBond(idxs[0], idxs[1]);
                    if (!atom1 || !atom2 || !bond)
                    {
                        continue;
                    }

                    const ElementInfo *info2 = GetElementInfo(atom2->GetAtomicNum());
                    const int atom2_room = info2 == nullptr
                                               ? -1
                                               : info2->default_valence - TotalValenceForRoomCheck(atom2);
                    if (info2 != nullptr &&
                        atom2_room >= 1 &&
                        (bond->GetBondOrder() == 1 || bond->GetBondOrder() == 2))
                    {
                        bond->SetBondOrder(bond->GetBondOrder() + 1);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() + 1);
                        atom2->SetFormalCharge(atom2->GetFormalCharge() - 1);
                    }
                }
            }

            void CleanResonances12(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[#7v3+0,#8v2+0,#16v2+0]-,:[*]=,:[*]-,:[*+1]");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom4 = mol.GetAtom(idxs[3]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    OBBond *bond3 = mol.GetBond(idxs[2], idxs[3]);
                    if (!atom1 || !atom4 || !bond1 || !bond2 || !bond3)
                    {
                        continue;
                    }

                    const ElementInfo *info4 = GetElementInfo(atom4->GetAtomicNum());
                    const int atom4_room = info4 == nullptr
                                               ? -1
                                               : info4->default_valence - TotalValenceForRoomCheck(atom4);
                    if (info4 != nullptr &&
                        atom4_room >= 1 &&
                        bond1->GetBondOrder() == 1 &&
                        bond2->GetBondOrder() == 2 &&
                        bond3->GetBondOrder() == 1)
                    {
                        bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                        bond2->SetBondOrder(bond2->GetBondOrder() - 1);
                        bond3->SetBondOrder(bond3->GetBondOrder() + 1);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() + 1);
                        atom4->SetFormalCharge(atom4->GetFormalCharge() - 1);
                    }
                }
            }

            void CleanResonances13(OBMol &mol)
            {
                auto matches = molgr::utils::FindSmarts(mol, "[*-]:[*]=[#7+0,#8+0]");
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom3 = mol.GetAtom(idxs[2]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    if (!atom1 || !atom3 || !bond1 || !bond2)
                    {
                        continue;
                    }

                    const ElementInfo *info1 = GetElementInfo(atom1->GetAtomicNum());
                    if (info1 != nullptr &&
                        info1->default_valence - atom1->GetTotalValence() >= 1 &&
                        bond1->GetBondOrder() == 1)
                    {
                        bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                        bond2->SetBondOrder(bond2->GetBondOrder() - 1);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() + 1);
                        atom3->SetFormalCharge(atom3->GetFormalCharge() - 1);
                    }
                }
            }
        }

        void CleanResonances(OBMol &mol)
        {
            CleanResonances0(mol);
            CleanResonances1(mol);
            CleanResonances2(mol);
            CleanResonances3(mol);
            CleanResonances4(mol);
            CleanResonances9(mol);
            CleanResonances5(mol);
            CleanResonances6(mol);
            CleanResonances7(mol);
            CleanResonances8(mol);
            CleanResonances9(mol);
            CleanResonances10(mol);
            CleanResonances11(mol);
            CleanResonances12(mol);
            CleanResonances13(mol);
        }
    }
}
