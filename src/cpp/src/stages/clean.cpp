#include "molgr/stages/clean.h"

#include "molgr/stages/fresh.h"
#include "molgr/stages/internal_helpers.h"

#include "molgr/utils/electrons.h"
#include "molgr/utils/smarts.h"
#include "molgr/vendor/openbabel_threading.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/kekulize.h>
#include "molgr/compat/openbabel_iter.h"

#include <algorithm>
#include <cstdlib>

namespace molgr
{
    namespace reconstruct
    {
        using namespace OpenBabel;

        // Electron bookkeeping: move an unresolved A center through A-B=C,
        // consume its marker, and produce one real unpaired electron at A and C.
        // Active lone pairs at C are preserved.
        bool CleanCarbeneNeighborUnsaturated(OBMol &mol)
        {
            bool hit = false;
            while (true)
            {
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_CARBENE_NEIGHBOR_UNSAT);
                if (matches.empty())
                    break;
                bool any_applied = false;

                for (const auto &idxs : matches)
                {
                    OBAtom *a1 = mol.GetAtom(idxs[0]);
                    OBAtom *a2 = mol.GetAtom(idxs[1]);
                    OBAtom *a3 = mol.GetAtom(idxs[2]);

                    OBBond *b23 = mol.GetBond(a2, a3);
                    OBBond *b12 = mol.GetBond(a1, a2);
                    if (molgr::utils::HasUnresolvedTwoElectronCenter(*a1) &&
                        molgr::utils::GetUnpairedElectronCount(*a3) == 0 &&
                        !molgr::utils::HasUnresolvedTwoElectronCenter(*a3) &&
                        b23 && b12 && b12->GetBondOrder() == 1 &&
                        b23->GetBondOrder() == 2)
                    {
                        b23->SetBondOrder(b23->GetBondOrder() - 1);
                        b12->SetBondOrder(b12->GetBondOrder() + 1);
                        molgr::utils::SetUnresolvedTwoElectronCenter(*a1, false);
                        molgr::utils::SetLonePairCount(*a1, 0);
                        molgr::utils::SetUnpairedElectronCount(*a1, 1);
                        molgr::utils::SetUnpairedElectronCount(*a3, molgr::utils::GetUnpairedElectronCount(*a3) + 1);
                        hit = true;
                        any_applied = true;
                        break;
                    }
                }
                if (!any_applied)
                    break;
            }
            return hit;
        }

        // Convert two excess terminal radicals into one neutral 1,3-dipole.
        // The target radical count and pending absolute charge reserve are global
        // constraints; each accepted fragment consumes one radical at each end.
        bool CleanPossible13Dipole(
            OBMol &mol,
            int given_charge,
            int total_radical_electrons)
        {
            const auto available_unpaired_electrons = [&]()
            {
                int count = 0;
                FOR_ATOMS_OF_MOL(atom_iter, mol)
                {
                    count += molgr::utils::GetUnpairedElectronCount(*atom_iter);
                }
                return count - total_radical_electrons - std::abs(given_charge);
            };

            bool hit = false;
            const auto matches = molgr::smarts::FindAll(
                mol,
                molgr::smarts::PatternId::CLEAN_POSSIBLE_1_3_DIPOLE);
            for (const auto &idxs : matches)
            {
                if (available_unpaired_electrons() < 2)
                {
                    break;
                }
                if (idxs.size() != 3)
                {
                    continue;
                }

                OBAtom *atom1 = mol.GetAtom(idxs[0]);
                OBAtom *atom2 = mol.GetAtom(idxs[1]);
                OBAtom *atom3 = mol.GetAtom(idxs[2]);
                OBBond *bond12 = mol.GetBond(idxs[0], idxs[1]);
                OBBond *bond23 = mol.GetBond(idxs[1], idxs[2]);
                if (!atom1 || !atom2 || !atom3 || !bond12 || !bond23)
                {
                    continue;
                }
                if (atom1->GetFormalCharge() != 0 ||
                    atom2->GetFormalCharge() != 0 ||
                    atom3->GetFormalCharge() != 0 ||
                    molgr::utils::GetUnpairedElectronCount(*atom1) < 1 ||
                    molgr::utils::GetUnpairedElectronCount(*atom3) < 1 ||
                    molgr::utils::HasUnresolvedTwoElectronCenter(*atom1) ||
                    molgr::utils::HasUnresolvedTwoElectronCenter(*atom3))
                {
                    continue;
                }

                const int degree1 = static_cast<int>(atom1->GetExplicitDegree());
                const int degree3 = static_cast<int>(atom3->GetExplicitDegree());
                if (degree1 <= 0 || degree3 <= 0)
                {
                    continue;
                }
                const int valence1 = static_cast<int>(atom1->GetExplicitValence());
                const int valence3 = static_cast<int>(atom3->GetExplicitValence());
                const int left_average_key = valence1 * degree3;
                const int right_average_key = valence3 * degree1;
                const bool add_left_bond =
                    left_average_key < right_average_key ||
                    (left_average_key == right_average_key && atom1->GetIdx() < atom3->GetIdx());
                OBBond *bond_to_increase = add_left_bond ? bond12 : bond23;
                OBAtom *negative_atom = add_left_bond ? atom3 : atom1;
                if (bond_to_increase->GetBondOrder() >= 3)
                {
                    continue;
                }

                bond_to_increase->SetBondOrder(bond_to_increase->GetBondOrder() + 1);
                atom2->SetFormalCharge(1);
                negative_atom->SetFormalCharge(-1);
                molgr::utils::SetUnpairedElectronCount(
                    *atom1,
                    molgr::utils::GetUnpairedElectronCount(*atom1) - 1);
                molgr::utils::SetUnpairedElectronCount(
                    *atom3,
                    molgr::utils::GetUnpairedElectronCount(*atom3) - 1);
                hit = true;
            }
            return hit;
        }

        // Electron bookkeeping: resolve adjacent radical-compatible states.
        // Each endpoint contributes either its real unpaired electrons or the
        // two electrons represented by an unresolved center. A matched pair
        // consumes one electron at each endpoint and raises the bond order by
        // one. The operation is limited by the excess-electron budget, so an
        // unresolved/unresolved pair may stop at a double bond and leave real
        // unpaired electrons to be classified later.
        bool CleanNeighborRadicals(
            OBMol &mol,
            int given_charge,
            int total_radical_electrons)
        {
            const auto available_electron_budget = [&]()
            {
                int count = 0;
                FOR_ATOMS_OF_MOL(atom_iter, mol)
                {
                    count += molgr::utils::GetUnpairedElectronCount(*atom_iter);
                    if (molgr::utils::HasUnresolvedTwoElectronCenter(*atom_iter))
                    {
                        count += 2;
                    }
                }
                return count - total_radical_electrons - std::abs(given_charge);
            };

            bool hit = false;
            FOR_BONDS_OF_MOL(bond_iter, mol)
            {
                OBBond *bond = &(*bond_iter);
                OBAtom *a1 = bond->GetBeginAtom();
                OBAtom *a2 = bond->GetEndAtom();
                int r1 = molgr::utils::GetUnpairedElectronCount(*a1);
                int r2 = molgr::utils::GetUnpairedElectronCount(*a2);
                const bool u1 = molgr::utils::HasUnresolvedTwoElectronCenter(*a1);
                const bool u2 = molgr::utils::HasUnresolvedTwoElectronCenter(*a2);
                const int capacity1 = u1 ? 2 : r1;
                const int capacity2 = u2 ? 2 : r2;
                if (capacity1 <= 0 || capacity2 <= 0)
                {
                    continue;
                }

                const int available = available_electron_budget();
                if (available < 2)
                {
                    continue;
                }
                const int to_add = std::min({capacity1, capacity2, available / 2});
                if (to_add <= 0)
                {
                    continue;
                }

                bond->SetBondOrder(bond->GetBondOrder() + to_add);
                const auto consume_endpoint = [&](OBAtom *atom, bool unresolved, int capacity)
                {
                    const int remaining = capacity - to_add;
                    if (unresolved)
                    {
                        molgr::utils::SetUnresolvedTwoElectronCenter(*atom, false);
                        molgr::utils::SetLonePairCount(*atom, 0);
                        molgr::utils::SetUnpairedElectronCount(*atom, remaining);
                        // Do not reclassify a partial unresolved center: its
                        // post-bond deficit would immediately recreate the
                        // unresolved marker.
                        if (remaining == 0)
                        {
                            AssignChargeRadicalForAtom(*atom);
                        }
                    }
                    else
                    {
                        molgr::utils::SetUnpairedElectronCount(*atom, remaining);
                        AssignChargeRadicalForAtom(*atom);
                    }
                };
                consume_endpoint(a1, u1, capacity1);
                consume_endpoint(a2, u2, capacity2);
                hit = true;
            }
            return hit;
        }

        // Electron bookkeeping: resolve separated terminal radical-compatible
        // states across A-B=C-D by shifting the middle pi bond outward.
        bool Clean14Radicals(
            OBMol &mol,
            int given_charge,
            int total_radical_electrons)
        {
            const auto available_unpaired_electrons = [&]()
            {
                int count = 0;
                FOR_ATOMS_OF_MOL(atom_iter, mol)
                {
                    count += molgr::utils::GetUnpairedElectronCount(*atom_iter);
                }
                return count - total_radical_electrons - std::abs(given_charge);
            };

            bool hit = false;
            const auto matches = molgr::smarts::FindAll(
                mol,
                molgr::smarts::PatternId::CLEAN_1_4_RADICALS);
            for (const auto &idxs : matches)
            {
                if (idxs.size() != 4)
                {
                    continue;
                }
                OBAtom *a1 = mol.GetAtom(idxs[0]);
                OBAtom *a2 = mol.GetAtom(idxs[1]);
                OBAtom *a3 = mol.GetAtom(idxs[2]);
                OBAtom *a4 = mol.GetAtom(idxs[3]);
                OBBond *b1 = mol.GetBond(a1, a2);
                OBBond *b2 = mol.GetBond(a2, a3);
                OBBond *b3 = mol.GetBond(a3, a4);
                if (!a1 || !a2 || !a3 || !a4 || !b1 || !b2 || !b3 ||
                    b1->GetBondOrder() != 1 || b2->GetBondOrder() != 2 ||
                    b3->GetBondOrder() != 1)
                {
                    continue;
                }

                const int r1 = molgr::utils::GetUnpairedElectronCount(*a1);
                const int r4 = molgr::utils::GetUnpairedElectronCount(*a4);
                const bool u1 = molgr::utils::HasUnresolvedTwoElectronCenter(*a1);
                const bool u4 = molgr::utils::HasUnresolvedTwoElectronCenter(*a4);
                if ((r1 <= 0 && !u1) || (r4 <= 0 && !u4))
                {
                    continue;
                }
                if (static_cast<int>(r1 > 0) + static_cast<int>(r4 > 0) == 2 &&
                    available_unpaired_electrons() < 2)
                {
                    continue;
                }

                b1->SetBondOrder(2);
                b2->SetBondOrder(1);
                b3->SetBondOrder(2);
                if (r1 > 0)
                {
                    molgr::utils::SetUnpairedElectronCount(*a1, r1 - 1);
                }
                if (r4 > 0)
                {
                    molgr::utils::SetUnpairedElectronCount(*a4, r4 - 1);
                }
                if (u1)
                {
                    molgr::utils::SetUnresolvedTwoElectronCenter(*a1, false);
                    molgr::utils::SetLonePairCount(*a1, 0);
                }
                if (u4)
                {
                    molgr::utils::SetUnresolvedTwoElectronCenter(*a4, false);
                    molgr::utils::SetLonePairCount(*a4, 0);
                }
                AssignChargeRadicalForAtom(*a1);
                AssignChargeRadicalForAtom(*a4);
                hit = true;
            }
            return hit;
        }

        // Electron bookkeeping: resolve terminal radical-compatible states
        // across A-B=C-D=E-F by shifting both middle pi bonds outward.
        bool Clean16Radicals(
            OBMol &mol,
            int given_charge,
            int total_radical_electrons)
        {
            const auto available_unpaired_electrons = [&]()
            {
                int count = 0;
                FOR_ATOMS_OF_MOL(atom_iter, mol)
                {
                    count += molgr::utils::GetUnpairedElectronCount(*atom_iter);
                }
                return count - total_radical_electrons - std::abs(given_charge);
            };

            bool hit = false;
            const auto matches = molgr::smarts::FindAll(
                mol,
                molgr::smarts::PatternId::CLEAN_1_6_RADICALS);
            for (const auto &idxs : matches)
            {
                if (idxs.size() != 6)
                {
                    continue;
                }
                OBAtom *a1 = mol.GetAtom(idxs[0]);
                OBAtom *a2 = mol.GetAtom(idxs[1]);
                OBAtom *a3 = mol.GetAtom(idxs[2]);
                OBAtom *a4 = mol.GetAtom(idxs[3]);
                OBAtom *a5 = mol.GetAtom(idxs[4]);
                OBAtom *a6 = mol.GetAtom(idxs[5]);
                OBBond *b1 = mol.GetBond(a1, a2);
                OBBond *b2 = mol.GetBond(a2, a3);
                OBBond *b3 = mol.GetBond(a3, a4);
                OBBond *b4 = mol.GetBond(a4, a5);
                OBBond *b5 = mol.GetBond(a5, a6);
                if (!a1 || !a2 || !a3 || !a4 || !a5 || !a6 ||
                    !b1 || !b2 || !b3 || !b4 || !b5 ||
                    b1->GetBondOrder() != 1 || b2->GetBondOrder() != 2 ||
                    b3->GetBondOrder() != 1 || b4->GetBondOrder() != 2 ||
                    b5->GetBondOrder() != 1)
                {
                    continue;
                }

                const int r1 = molgr::utils::GetUnpairedElectronCount(*a1);
                const int r6 = molgr::utils::GetUnpairedElectronCount(*a6);
                const bool u1 = molgr::utils::HasUnresolvedTwoElectronCenter(*a1);
                const bool u6 = molgr::utils::HasUnresolvedTwoElectronCenter(*a6);
                if ((r1 <= 0 && !u1) || (r6 <= 0 && !u6))
                {
                    continue;
                }
                if (static_cast<int>(r1 > 0) + static_cast<int>(r6 > 0) == 2 &&
                    available_unpaired_electrons() < 2)
                {
                    continue;
                }

                b1->SetBondOrder(2);
                b2->SetBondOrder(1);
                b3->SetBondOrder(2);
                b4->SetBondOrder(1);
                b5->SetBondOrder(2);
                if (r1 > 0)
                {
                    molgr::utils::SetUnpairedElectronCount(*a1, r1 - 1);
                }
                if (r6 > 0)
                {
                    molgr::utils::SetUnpairedElectronCount(*a6, r6 - 1);
                }
                if (u1)
                {
                    molgr::utils::SetUnresolvedTwoElectronCenter(*a1, false);
                    molgr::utils::SetLonePairCount(*a1, 0);
                }
                if (u6)
                {
                    molgr::utils::SetUnresolvedTwoElectronCenter(*a6, false);
                    molgr::utils::SetLonePairCount(*a6, 0);
                }
                AssignChargeRadicalForAtom(*a1);
                AssignChargeRadicalForAtom(*a6);
                hit = true;
            }
            return hit;
        }

        namespace
        {
            bool CleanResonances0(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_0);
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
                    if (atom1->GetFormalCharge() == -1 &&
                        atom4->GetFormalCharge() == 1 &&
                        bond1->GetBondOrder() == 1 &&
                        bond2->GetBondOrder() == 2 &&
                        info1 != nullptr && info4 != nullptr &&
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
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_1);
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom2 = mol.GetAtom(idxs[1]);
                    OBAtom *atom3 = mol.GetAtom(idxs[2]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    if (!atom1 || !atom2 || !atom3 || !bond1 || !bond2)
                    {
                        continue;
                    }

                    if (atom1->GetFormalCharge() == -1 &&
                        atom2->GetFormalCharge() == 1 &&
                        atom3->GetFormalCharge() == 0 &&
                        bond1->GetBondOrder() == 2 &&
                        bond2->GetBondOrder() == 2)
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

            bool CleanResonances2(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_2);
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

                    const ElementInfo *info6 = GetElementInfo(atom6->GetAtomicNum());
                    const int atom6_room = info6 == nullptr
                                               ? -1
                                               : info6->default_valence -
                                                     (atom6->GetTotalValence() + 1);

                    if (bond4->GetBondOrder() == 1 &&
                        bond3->GetBondOrder() == 2 &&
                        bond2->GetBondOrder() == 1 &&
                        bond1->GetBondOrder() == 2 &&
                        atom1->GetFormalCharge() == 0 &&
                        atom6->GetFormalCharge() == -1 &&
                        atom6_room >= 0)
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

            // Electron bookkeeping: rule 3 changes total valence only at its two
            // charged endpoints; internal atoms exchange bond-order units with net
            // zero change. Refresh endpoints only and preserve unrelated labels.
            bool CleanResonances3(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_3);
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

                    if (atom1->GetFormalCharge() == 1 &&
                        atom5->GetFormalCharge() == -1 &&
                        bond1->GetBondOrder() == 2 &&
                        bond2->GetBondOrder() == 1 &&
                        bond3->GetBondOrder() == 2 &&
                        bond4->GetBondOrder() == 1)
                    {
                        bond1->SetBondOrder(bond1->GetBondOrder() - 1);
                        bond2->SetBondOrder(bond2->GetBondOrder() + 1);
                        bond3->SetBondOrder(bond3->GetBondOrder() - 1);
                        bond4->SetBondOrder(bond4->GetBondOrder() + 1);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() - 1);
                        atom5->SetFormalCharge(atom5->GetFormalCharge() + 1);
                        AssignChargeRadicalForAtom(*atom1);
                        AssignChargeRadicalForAtom(*atom5);
                        hit = true;
                    }
                }
                return hit;
            }

            bool CleanResonances4(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_4);
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

                    if (atom1->GetFormalCharge() == 1 &&
                        atom3->GetFormalCharge() == -1 &&
                        bond1->GetBondOrder() == 2 &&
                        bond2->GetBondOrder() == 1)
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

            bool CleanResonances5(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_5);
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
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_6);
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

                    if (atom1->GetFormalCharge() == 0 &&
                        atom3->GetFormalCharge() == -1 &&
                        bond1->GetBondOrder() == 2 &&
                        bond2->GetBondOrder() == 2)
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

            bool CleanResonances7(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_7);
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
                    const int atom1_room = info1 == nullptr
                                               ? -1
                                               : info1->default_valence -
                                                     (atom1->GetTotalValence() + 1);
                    if (bond1->GetBondOrder() == 1 &&
                        bond2->GetBondOrder() == 2 &&
                        atom1->GetFormalCharge() == -1 &&
                        atom3->GetFormalCharge() == 0 &&
                        atom1_room >= 0)
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

            // Electron bookkeeping: shift charge/pi bonds and locally rebuild both
            // endpoints so newly exposed unpaired/lone-pair/unresolved states are
            // not hidden by stale fields.
            bool CleanResonances8(OBMol &mol)
            {
                bool hit = false;
                molgr::vendor::openbabel_threading::SetAromaticPerceived(mol, false);
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_8);
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
                    if (atom1->GetFormalCharge() == -1 &&
                        atom5->GetFormalCharge() == 0 &&
                        bond1->GetBondOrder() == 1 &&
                        bond2->GetBondOrder() == 2 &&
                        bond3->GetBondOrder() == 1 &&
                        bond4->GetBondOrder() == 2)
                    {
                        bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                        bond2->SetBondOrder(bond2->GetBondOrder() - 1);
                        bond3->SetBondOrder(bond3->GetBondOrder() + 1);
                        bond4->SetBondOrder(bond4->GetBondOrder() - 1);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() + 1);
                        atom5->SetFormalCharge(atom5->GetFormalCharge() - 1);
                        AssignChargeRadicalForAtom(*atom1);
                        AssignChargeRadicalForAtom(*atom5);
                        hit = true;
                    }
                }
                return hit;
            }

            bool CleanResonances9(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_9);
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
                    if (atom1->GetFormalCharge() > 0 &&
                        atom2->GetFormalCharge() < 0 &&
                        (bond->GetBondOrder() == 1 || bond->GetBondOrder() == 2) &&
                        room1 >= 1 && room2 >= 1)
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

            // Electron bookkeeping: require and consume one real unpaired electron
            // at each terminal to create two new bond-order units.
            bool CleanResonances10(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_10);
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

                    if (molgr::utils::GetUnpairedElectronCount(*atom1) == 1 &&
                        molgr::utils::GetUnpairedElectronCount(*atom4) == 1 &&
                        bond1->GetBondOrder() == 1 &&
                        (bond2->GetBondOrder() == 2 || bond2->GetBondOrder() == 3) &&
                        bond3->GetBondOrder() == 1)
                    {
                        bond2->SetBondOrder(bond2->GetBondOrder() - 1);
                        bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                        bond3->SetBondOrder(bond3->GetBondOrder() + 1);
                        molgr::utils::SetUnpairedElectronCount(*atom1, molgr::utils::GetUnpairedElectronCount(*atom1) - 1);
                        molgr::utils::SetUnpairedElectronCount(*atom4, molgr::utils::GetUnpairedElectronCount(*atom4) - 1);
                        hit = true;
                    }
                }
                return hit;
            }

            bool CleanResonances11(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_11);
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
                        atom1->GetFormalCharge() == 0 &&
                        atom2->GetFormalCharge() == 1 &&
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
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_12);
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
                        atom1->GetFormalCharge() == 0 &&
                        atom4->GetFormalCharge() == 1 &&
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
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_13);
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
                        atom1->GetFormalCharge() == -1 &&
                        atom3->GetFormalCharge() == 0 &&
                        info1->default_valence - atom1->GetTotalValence() >= 1 &&
                        bond1->GetBondOrder() == 1 &&
                        bond2->GetBondOrder() == 2)
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

            // Electron bookkeeping: lower the charged triple bond, neutralize both
            // endpoints, and rebuild their explicit electron classifications.
            bool CleanResonances14Impl(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_14);
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

                    if (atom1->GetFormalCharge() == -1 &&
                        atom2->GetFormalCharge() == 1 &&
                        bond->GetBondOrder() == 3)
                    {
                        bond->SetBondOrder(bond->GetBondOrder() - 1);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() + 1);
                        atom2->SetFormalCharge(atom2->GetFormalCharge() - 1);
                        AssignChargeRadicalForAtom(*atom1);
                        AssignChargeRadicalForAtom(*atom2);
                        hit = true;
                    }
                }
                return hit;
            }

            // Electron bookkeeping: neutralize the conjugated ion pair, then
            // rebuild both endpoints to expose any newly created explicit or
            // unresolved electron state.
            bool CleanResonances16Impl(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_16);
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

                    if (atom1->GetFormalCharge() == -1 &&
                        molgr::utils::GetUnpairedElectronCount(*atom1) == 0 &&
                        molgr::utils::GetLonePairCount(*atom1) == 0 &&
                        atom5->GetFormalCharge() == 1 &&
                        molgr::utils::GetUnpairedElectronCount(*atom5) == 0 &&
                        molgr::utils::GetLonePairCount(*atom5) == 0 &&
                        bond1->GetBondOrder() == 1 &&
                        bond2->GetBondOrder() == 2 &&
                        bond3->GetBondOrder() == 1 &&
                        bond4->GetBondOrder() == 2)
                    {
                        bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                        bond2->SetBondOrder(bond2->GetBondOrder() - 1);
                        bond3->SetBondOrder(bond3->GetBondOrder() + 1);
                        bond4->SetBondOrder(bond4->GetBondOrder() - 1);
                        atom1->SetFormalCharge(atom1->GetFormalCharge() + 1);
                        atom5->SetFormalCharge(atom5->GetFormalCharge() - 1);
                        AssignChargeRadicalForAtom(*atom1);
                        AssignChargeRadicalForAtom(*atom5);
                        hit = true;
                    }
                }
                return hit;
            }

            // Shift A(-)-B=C=D to A=B-C(-)=D only when B, C, and D are
            // ring atoms, replacing a ring allene with conjugated pi bonds.
            bool CleanResonances17Impl(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_17);
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());

                    OBAtom *atom1 = mol.GetAtom(idxs[0]);
                    OBAtom *atom3 = mol.GetAtom(idxs[2]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    OBBond *bond3 = mol.GetBond(idxs[2], idxs[3]);
                    if (!atom1 || !atom3 || !bond1 || !bond2 || !bond3)
                    {
                        continue;
                    }
                    const ElementInfo *info1 = GetElementInfo(atom1->GetAtomicNum());
                    const int atom1_room = info1 == nullptr
                                               ? -1
                                               : info1->default_valence -
                                                     (atom1->GetTotalValence() + 1);
                    if (atom1->GetFormalCharge() == -1 &&
                        atom3->GetFormalCharge() == 0 &&
                        bond1->GetBondOrder() == 1 &&
                        bond2->GetBondOrder() == 2 &&
                        bond3->GetBondOrder() == 2 &&
                        atom1_room >= 0)
                    {
                        bond1->SetBondOrder(2);
                        bond2->SetBondOrder(1);
                        atom1->SetFormalCharge(0);
                        atom3->SetFormalCharge(-1);
                        hit = true;
                    }
                }
                return hit;
            }

            // Electron bookkeeping: resolve an unresolved terminal diazene as
            // the charge-separated azide resonance form.
            bool CleanResonances18Impl(OBMol &mol)
            {
                bool hit = false;
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::CLEAN_RESONANCE_18);
                while (!matches.empty())
                {
                    auto idxs = matches.front();
                    matches.erase(matches.begin());
                    if (idxs.size() != 4)
                    {
                        continue;
                    }

                    OBAtom *atom2 = mol.GetAtom(idxs[1]);
                    OBAtom *atom3 = mol.GetAtom(idxs[2]);
                    OBAtom *atom4 = mol.GetAtom(idxs[3]);
                    OBBond *bond1 = mol.GetBond(idxs[0], idxs[1]);
                    OBBond *bond2 = mol.GetBond(idxs[1], idxs[2]);
                    OBBond *bond3 = mol.GetBond(idxs[2], idxs[3]);
                    if (!atom2 || !atom3 || !atom4 || !bond1 || !bond2 || !bond3)
                    {
                        continue;
                    }
                    if (atom2->GetFormalCharge() != 0 ||
                        atom3->GetFormalCharge() != 0 ||
                        atom4->GetFormalCharge() != 0 ||
                        molgr::utils::GetUnpairedElectronCount(*atom2) != 0 ||
                        molgr::utils::GetUnpairedElectronCount(*atom3) != 0 ||
                        molgr::utils::GetUnpairedElectronCount(*atom4) != 0 ||
                        molgr::utils::HasUnresolvedTwoElectronCenter(*atom2) ||
                        molgr::utils::HasUnresolvedTwoElectronCenter(*atom3) ||
                        !molgr::utils::HasUnresolvedTwoElectronCenter(*atom4) ||
                        bond1->GetBondOrder() != 1 ||
                        bond2->GetBondOrder() != 2 ||
                        bond3->GetBondOrder() != 1)
                    {
                        continue;
                    }

                    bond3->SetBondOrder(bond3->GetBondOrder() + 1);
                    atom3->SetFormalCharge(1);
                    atom4->SetFormalCharge(-1);
                    molgr::utils::SetUnresolvedTwoElectronCenter(*atom4, false);
                    molgr::utils::SetLonePairCount(*atom4, 0);
                    molgr::utils::SetUnpairedElectronCount(*atom4, 0);
                    hit = true;
                }
                return hit;
            }
        }

        bool CleanResonances14(OBMol &mol)
        {
            return CleanResonances14Impl(mol);
        }

        bool CleanResonances16(OBMol &mol)
        {
            return CleanResonances16Impl(mol);
        }

        bool CleanResonances17(OBMol &mol)
        {
            return CleanResonances17Impl(mol);
        }

        bool CleanResonances18(OBMol &mol)
        {
            return CleanResonances18Impl(mol);
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
            hit = CleanResonances14(mol) || hit;
            hit = CleanResonances16(mol) || hit;
            hit = CleanResonances17(mol) || hit;
            hit = CleanResonances18(mol) || hit;
            return hit;
        }
    }
}
