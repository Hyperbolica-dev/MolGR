#include "molgr/stages/clean.h"

#include "molgr/stages/fresh.h"
#include "molgr/stages/internal_helpers.h"

#include "molgr/utils/smarts.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/kekulize.h>
#include <openbabel/obiter.h>

#include <openbabel/obconversion.h>
#include <algorithm>

namespace molgr
{
    namespace reconstruct
    {
        using namespace OpenBabel;

        bool CleanCarbeneNeighborUnsaturated(OBMol &mol)
        {
            bool hit = false;
            while (true)
            {
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_CARBENE_NEIGHBOR_UNSAT);
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
                            hit = true;
                            any_applied = true;
                            break;
                        }
                    }
                }
                if (!any_applied)
                    break;
            }
            return hit;
        }

        bool CleanNeighborRadicals(OBMol &mol)
        {
            bool hit = false;
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
                    hit = true;
                }
            }
            return hit;
        }

        namespace
        {
            bool CleanResonances0(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_0);
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
                        hit = true;
                    }
                }
                return hit;
            }

            bool CleanResonances1(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_1);
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
                    hit = true;
                }
                return hit;
            }

            bool CleanResonances2(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_2);
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
                        hit = true;
                    }
                }
                return hit;
            }

            bool CleanResonances3(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_3);
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
                    hit = true;
                    hit = FreshOmolChargeRadical(mol) || hit;
                }
                return hit;
            }

            bool CleanResonances4(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_4);
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
                    hit = true;
                }
                return hit;
            }

            bool CleanResonances5(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_5);
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
                        hit = true;
                    }
                }
                return hit;
            }

            bool CleanResonances6(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_6);
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
                    hit = true;
                }
                return hit;
            }

            bool CleanResonances7(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_7);
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
                    if (bond1->GetBondOrder() == 1 &&
                        bond2->GetBondOrder() == 2 &&
                        atom1->GetFormalCharge() == -1 &&
                        atom3->GetFormalCharge() == 0)
                    {
                        bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                        bond2->SetBondOrder(bond2->GetBondOrder() - 1);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() + 1);
                        atom3->SetFormalCharge(atom3->GetFormalCharge() - 1);
                        hit = true;
                    }
                }
                return hit;
            }

            bool CleanResonances8(OBMol &mol)
            {
                bool hit = false;
                mol.SetAromaticPerceived(false);
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_8);
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
                    if (bond1->GetBondOrder() == 1 &&
                        bond2->GetBondOrder() == 2 &&
                        bond3->GetBondOrder() == 1 &&
                        bond4->GetBondOrder() == 2 &&
                        atom1->GetFormalCharge() == -1 &&
                        atom5->GetFormalCharge() == 0)
                    {
                        bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                        bond2->SetBondOrder(bond2->GetBondOrder() - 1);
                        bond3->SetBondOrder(bond3->GetBondOrder() + 1);
                        bond4->SetBondOrder(bond4->GetBondOrder() - 1);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() + 1);
                        atom5->SetFormalCharge(atom5->GetFormalCharge() - 1);
                        hit = true;
                    }
                }
                return hit;
            }

            bool CleanResonances9(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_9);
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
                        hit = true;
                    }
                }
                return hit;
            }

            bool CleanResonances10(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_10);
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
                        hit = true;
                    }
                }
                return hit;
            }

            bool CleanResonances11(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_11);
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
                                               : info2->default_valence - atom2->GetTotalValence();
                    if (info2 != nullptr &&
                        atom2_room >= 1 &&
                        (bond->GetBondOrder() == 1 || bond->GetBondOrder() == 2))
                    {
                        bond->SetBondOrder(bond->GetBondOrder() + 1);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() + 1);
                        atom2->SetFormalCharge(atom2->GetFormalCharge() - 1);
                        hit = true;
                    }
                }
                return hit;
            }

            bool CleanResonances12(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_12);
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
                                               : info4->default_valence - atom4->GetTotalValence();
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
                        hit = true;
                    }
                }
                return hit;
            }

            bool CleanResonances13(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::Match(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_13);
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
                        hit = true;
                    }
                }
                return hit;
            }
        }

        bool CleanResonances(OBMol &mol)
        {
            bool hit = false;
            hit = CleanResonances11(mol) || hit;
            hit = CleanResonances0(mol) || hit;
            hit = CleanResonances1(mol) || hit;
            hit = CleanResonances2(mol) || hit;
            hit = CleanResonances3(mol) || hit;
            hit = CleanResonances4(mol) || hit;
            hit = CleanResonances9(mol) || hit;
            hit = CleanResonances5(mol) || hit;
            hit = CleanResonances6(mol) || hit;
            hit = CleanResonances7(mol) || hit;
            hit = CleanResonances8(mol) || hit;
            hit = CleanResonances9(mol) || hit;
            hit = CleanResonances10(mol) || hit;
            hit = CleanResonances12(mol) || hit;
            hit = CleanResonances13(mol) || hit;
            return hit;
        }
    }
}
