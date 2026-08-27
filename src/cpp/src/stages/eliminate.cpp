#include "molgr/stages/eliminate.h"

#include "molgr/stages/internal_helpers.h"

#include "molgr/utils/consts.h"
#include "molgr/utils/electrons.h"
#include "molgr/utils/logger.h"
#include "molgr/utils/smarts.h"
#include "molgr/utils/utils.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/elements.h>
#include "molgr/compat/openbabel_iter.h"

#include <algorithm>
#include <map>
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
            bool consume_unresolved_center;
            int charge_delta;
            std::vector<int> score_key;
            int bond_idx = -1;
            int charge_atom_idx = -1;
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

        // Electron bookkeeping: choose how many real unpaired electrons to bundle
        // into one multivalent anion action; never count lone-pair/marker state.
        int NegativeChargeAssignmentAmount(OBAtom *atom, int charge)
        {
            if (atom == nullptr)
            {
                return 0;
            }
            return std::min(molgr::utils::GetUnpairedElectronCount(*atom), std::max(1, std::abs(charge)));
        }

        // Electron bookkeeping: radical actions consume real unpaired electrons.
        // An unresolved action atomically consumes one pure deferred 2e marker as
        // +/-2 without first creating a diradical. Active lone pairs are excluded.
        bool ApplyChargeAssignmentAction(OBMol &mol, const ChargeAssignmentAction &action)
        {
            OBAtom *atom = mol.GetAtom(action.atom_idx);
            if (atom == nullptr)
            {
                return false;
            }
            if (action.bond_idx >= 0 && action.charge_atom_idx >= 0)
            {
                OBAtom *charge_atom = mol.GetAtom(action.charge_atom_idx);
                OBBond *bond = mol.GetBond(atom, charge_atom);
                if (charge_atom == nullptr || bond == nullptr ||
                    molgr::utils::HasUnresolvedTwoElectronCenter(*atom) ||
                    molgr::utils::GetUnpairedElectronCount(*atom) < action.spin_consumed)
                {
                    return false;
                }
                molgr::utils::SetUnpairedElectronCount(
                    *atom,
                    molgr::utils::GetUnpairedElectronCount(*atom) - action.spin_consumed);
                bond->SetBondOrder(bond->GetBondOrder() + 1);
                charge_atom->SetFormalCharge(charge_atom->GetFormalCharge() + 1);
                return true;
            }
            if (atom->GetFormalCharge() != 0)
            {
                return false;
            }
            if (action.consume_unresolved_center)
            {
                if (action.spin_consumed != 0 ||
                    std::abs(action.formal_charge) != 2 ||
                    !molgr::utils::HasUnresolvedTwoElectronCenter(*atom) ||
                    molgr::utils::GetUnpairedElectronCount(*atom) != 0 ||
                    molgr::utils::GetLonePairCount(*atom) != 0)
                {
                    return false;
                }
                molgr::utils::SetUnresolvedTwoElectronCenter(*atom, false);
            }
            else
            {
                if (action.spin_consumed <= 0 ||
                    molgr::utils::HasUnresolvedTwoElectronCenter(*atom) ||
                    molgr::utils::GetUnpairedElectronCount(*atom) < action.spin_consumed)
                {
                    return false;
                }
                molgr::utils::SetUnpairedElectronCount(
                    *atom,
                    molgr::utils::GetUnpairedElectronCount(*atom) - action.spin_consumed);
            }
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
            if (atom == nullptr || amount <= 0 ||
                molgr::utils::HasUnresolvedTwoElectronCenter(*atom))
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
                false,
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
            if (atom == nullptr || amount <= 0 ||
                molgr::utils::HasUnresolvedTwoElectronCenter(*atom))
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
                false,
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

        // Consume one pure unresolved center directly as a +/-2 formal-charge
        // action. A deficit smaller than two is rejected by the caller.
        void AppendUnresolvedChargeAssignmentAction(
            std::vector<ChargeAssignmentAction> &actions,
            OBAtom *atom,
            int charge,
            int formal_charge,
            int tier,
            int match_order)
        {
            if (atom == nullptr ||
                std::abs(formal_charge) != 2 ||
                atom->GetFormalCharge() != 0 ||
                !molgr::utils::HasUnresolvedTwoElectronCenter(*atom) ||
                molgr::utils::GetUnpairedElectronCount(*atom) != 0 ||
                molgr::utils::GetLonePairCount(*atom) != 0)
            {
                return;
            }
            const int charge_after = charge - formal_charge;
            const int atom_idx = static_cast<int>(atom->GetIdx());
            const int atomic_num = static_cast<int>(atom->GetAtomicNum());
            actions.push_back(ChargeAssignmentAction{
                atom_idx,
                formal_charge,
                0,
                true,
                -formal_charge,
                {
                    tier,
                    std::abs(charge_after),
                    std::max(charge_after, 0),
                    atomic_num,
                    atom_idx,
                    match_order,
                }});
        }

        // Electron bookkeeping: rank cation actions from real radicals or, with a
        // remaining deficit of at least two, one pure unresolved 2e center.
        std::vector<ChargeAssignmentAction> PositiveChargeAssignmentActions(
            OBMol &mol,
            int charge,
            int total_radical_electrons = 0)
        {
            std::vector<ChargeAssignmentAction> actions;
            std::set<std::tuple<int, int, int>> seen;
            int real_radicals = 0;
            int unresolved_electrons = 0;
            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OBAtom *atom = &(*atom_iter);
                real_radicals += molgr::utils::GetUnpairedElectronCount(*atom);
                if (molgr::utils::HasUnresolvedTwoElectronCenter(*atom))
                {
                    unresolved_electrons += 1;
                }
            }
            const bool tier5_allowed =
                std::abs(charge - 1) <= real_radicals + unresolved_electrons * 2 - 1;
            auto n_matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::ELIM_POSITIVE_N);
            for (std::size_t match_order = 0; match_order < n_matches.size(); ++match_order)
            {
                const auto &idxs = n_matches[match_order];
                if (idxs.size() < 2)
                {
                    continue;
                }
                OBAtom *atom = mol.GetAtom(idxs[1]);
                if (charge > 0 && atom != nullptr && atom->GetFormalCharge() == 0 && molgr::utils::GetUnpairedElectronCount(*atom) >= 1)
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
                if (charge > 0 && atom != nullptr && atom->GetFormalCharge() == 0 && molgr::utils::GetUnpairedElectronCount(*atom) >= 1)
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

            auto dipole_matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::ELIM_POSITIVE_DIPOLE);
            for (std::size_t match_order = 0; match_order < dipole_matches.size(); ++match_order)
            {
                const auto &idxs = dipole_matches[match_order];
                if (idxs.size() < 2)
                {
                    continue;
                }
                OBAtom *radical_atom = mol.GetAtom(idxs[0]);
                OBAtom *charge_atom = mol.GetAtom(idxs[1]);
                OBBond *bond = mol.GetBond(radical_atom, charge_atom);
                if (tier5_allowed && radical_atom != nullptr && charge_atom != nullptr && bond != nullptr &&
                    radical_atom->GetFormalCharge() == 0 &&
                    molgr::utils::GetUnpairedElectronCount(*radical_atom) == 1 &&
                    !molgr::utils::HasUnresolvedTwoElectronCenter(*radical_atom))
                {
                    const int charge_after = charge - 1;
                    actions.push_back(ChargeAssignmentAction{
                        static_cast<int>(radical_atom->GetIdx()),
                        0,
                        1,
                        false,
                        -1,
                        {5, std::abs(charge_after), std::max(charge_after, 0), static_cast<int>(radical_atom->GetIdx()), static_cast<int>(match_order)},
                        static_cast<int>(bond->GetIdx()),
                        static_cast<int>(charge_atom->GetIdx()),
                    });
                }
            }

            int match_order = 0;
            FOR_ATOMS_OF_MOL(atom_iter, mol)
            {
                OBAtom *atom = &(*atom_iter);
                if (charge > 0 && atom->GetFormalCharge() == 0 && molgr::utils::GetUnpairedElectronCount(*atom) >= 1)
                {
                    AppendPositiveChargeAssignmentAction(
                        actions,
                        seen,
                        atom,
                        charge,
                        100,
                        match_order,
                        std::min(molgr::utils::GetUnpairedElectronCount(*atom), charge));
                }
                if (charge >= 2)
                {
                    AppendUnresolvedChargeAssignmentAction(
                        actions, atom, charge, 2, 100, match_order);
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

        // Electron bookkeeping: rank anion actions from real radicals or, with a
        // remaining deficit of at least two, one pure unresolved 2e center.
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
                    if (atom != nullptr && atom->GetFormalCharge() == 0 && molgr::utils::GetUnpairedElectronCount(*atom) >= 1)
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
                    if (atom->GetFormalCharge() == 0 && molgr::utils::GetUnpairedElectronCount(*atom) >= 1)
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
                    if (charge <= -2)
                    {
                        AppendUnresolvedChargeAssignmentAction(
                            actions, atom, charge, -2, 1000, match_order);
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

        // Electron bookkeeping: the negative N-N-N closure consumes one classified
        // two-electron state at each terminal N and one real unpaired electron at
        // the middle N. The positive closure consumes one middle-N monoradical.
        // All electronic preconditions are checked before bonds or charges change.
        bool EliminateNNN(OBMol &mol, int &charge, bool positive)
        {
            // Accept exactly one pure unresolved, triplet, or active-singlet
            // two-electron state. Mixed marker/occupancy states are invalid here.
            const auto has_consumable_two_electron_center = [](const OBAtom &atom)
            {
                const bool unresolved =
                    molgr::utils::HasUnresolvedTwoElectronCenter(atom);
                const int unpaired = molgr::utils::GetUnpairedElectronCount(atom);
                const int lone_pairs = molgr::utils::GetLonePairCount(atom);
                return (unresolved && unpaired == 0 && lone_pairs == 0) ||
                       (!unresolved && unpaired == 2 && lone_pairs == 0) ||
                       (!unresolved && unpaired == 0 && lone_pairs >= 1);
            };
            // Consume the prevalidated state from its own field; do not infer a
            // lone pair from unpaired-electron parity or vice versa.
            const auto consume_two_electron_center = [](OBAtom &atom)
            {
                if (molgr::utils::HasUnresolvedTwoElectronCenter(atom))
                {
                    molgr::utils::SetUnresolvedTwoElectronCenter(atom, false);
                }
                else if (molgr::utils::GetUnpairedElectronCount(atom) == 2)
                {
                    molgr::utils::SetUnpairedElectronCount(atom, 0);
                }
                else if (molgr::utils::GetLonePairCount(atom) >= 1)
                {
                    molgr::utils::SetLonePairCount(
                        atom,
                        molgr::utils::GetLonePairCount(atom) - 1);
                }
            };
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
                    if (!has_consumable_two_electron_center(*a1) ||
                        molgr::utils::GetUnpairedElectronCount(*a2) != 1 ||
                        !has_consumable_two_electron_center(*a3))
                    {
                        break;
                    }

                    bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                    bond2->SetBondOrder(bond2->GetBondOrder() + 1);

                    consume_two_electron_center(*a1);
                    a1->SetFormalCharge(a1->GetFormalCharge() - 1);

                    molgr::utils::SetUnpairedElectronCount(*a2, molgr::utils::GetUnpairedElectronCount(*a2) - 1);
                    a2->SetFormalCharge(a2->GetFormalCharge() + 1);

                    consume_two_electron_center(*a3);
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
                    if (molgr::utils::GetUnpairedElectronCount(*a2) != 1)
                    {
                        break;
                    }

                bond1->SetBondOrder(bond1->GetBondOrder() + 1);
                a1->SetFormalCharge(a1->GetFormalCharge() + 1);
                molgr::utils::SetUnpairedElectronCount(*a2, molgr::utils::GetUnpairedElectronCount(*a2) - 1);
                charge -= 1;
                hit = true;
                LOG_DEBUG("[EliminateNNN] Applied positive branch");
            }
            return hit;
        }

        // Neutralize each positive center from eligible neighboring
        // monoradicals in descending electronegativity order. Stop once the
        // center charge is balanced by neighboring negative formal charges.
        bool EliminateHighPositiveChargeAtoms(OBMol &mol, int &charge)
        {
            bool hit = false;
            const auto matches = molgr::smarts::FindAll(
                mol,
                molgr::smarts::PatternId::ELIM_HIGH_POSITIVE);
            std::map<int, std::set<int>> neighbors_by_center;
            for (const auto &idxs : matches)
            {
                if (idxs.size() == 2)
                {
                    neighbors_by_center[idxs[0]].insert(idxs[1]);
                }
            }

            for (const auto &[center_idx, matched_neighbor_indices] : neighbors_by_center)
            {
                OBAtom *center = mol.GetAtom(center_idx);
                if (center == nullptr)
                {
                    continue;
                }
                int remaining_positive_charge = center->GetFormalCharge();
                FOR_NB_OF_ATOM(nbr, center)
                {
                    remaining_positive_charge += std::min(0, nbr->GetFormalCharge());
                }
                if (remaining_positive_charge <= 0)
                {
                    continue;
                }

                std::vector<int> neighbor_indices(
                    matched_neighbor_indices.begin(),
                    matched_neighbor_indices.end());
                std::sort(
                    neighbor_indices.begin(),
                    neighbor_indices.end(),
                    [&mol](int left_idx, int right_idx)
                    {
                        const OBAtom *left = mol.GetAtom(left_idx);
                        const OBAtom *right = mol.GetAtom(right_idx);
                        const double left_electronegativity =
                            OBElements::GetElectroNeg(left->GetAtomicNum());
                        const double right_electronegativity =
                            OBElements::GetElectroNeg(right->GetAtomicNum());
                        if (left_electronegativity != right_electronegativity)
                        {
                            return left_electronegativity > right_electronegativity;
                        }
                        return left_idx < right_idx;
                    });

                for (int neighbor_idx : neighbor_indices)
                {
                    if (remaining_positive_charge <= 0)
                    {
                        break;
                    }
                    OBAtom *neighbor = mol.GetAtom(neighbor_idx);
                    if (neighbor == nullptr)
                    {
                        continue;
                    }

                    bool neighbor_has_pending_electrons = false;
                    FOR_NB_OF_ATOM(adjacent, neighbor)
                    {
                        if (adjacent->GetIdx() == center_idx)
                        {
                            continue;
                        }
                        if (molgr::utils::GetUnpairedElectronCount(*adjacent) > 0 ||
                            molgr::utils::HasUnresolvedTwoElectronCenter(*adjacent))
                        {
                            neighbor_has_pending_electrons = true;
                            break;
                        }
                    }

                    if (neighbor->GetFormalCharge() != 0 ||
                        molgr::utils::GetUnpairedElectronCount(*neighbor) != 1 ||
                        molgr::utils::HasUnresolvedTwoElectronCenter(*neighbor) ||
                        neighbor_has_pending_electrons)
                    {
                        continue;
                    }

                    molgr::utils::SetUnpairedElectronCount(*neighbor, 0);
                    neighbor->SetFormalCharge(-1);
                    ++charge;
                    --remaining_positive_charge;
                    hit = true;
                    LOG_DEBUG("[EliminateHighPos] Applied");
                }
            }
            return hit;
        }

        // Electron bookkeeping: localize exactly one oxygen monoradical as O-.
        // A missing or multi-electron radical label rejects the transformation.
        bool EliminateCarboxyl(OBMol &mol, int &charge)
        {
            bool hit = false;
            while (true)
            {
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::ELIM_CARBOXYL);
                if (matches.empty())
                    break;
                OBAtom *a1 = mol.GetAtom(matches[0][0]);
                if (a1 == nullptr ||
                    molgr::utils::GetUnpairedElectronCount(*a1) != 1)
                {
                    break;
                }
                molgr::utils::SetUnpairedElectronCount(*a1, 0);
                a1->SetFormalCharge(a1->GetFormalCharge() - 1);
                charge += 1;
                hit = true;
                LOG_DEBUG("[EliminateCarboxyl] Applied");
            }
            return hit;
        }

        // Electron bookkeeping: close an unresolved two-electron center with a
        // closed-shell heteroatom donor. The center marker becomes center-/donor+;
        // consume one explicitly tracked active donor lone pair when present.
        // Real unpaired electrons cannot act as the donor in this rule.
        bool EliminateCarbeneNeighborHeteroatom(OBMol &mol, int &charge)
        {
            const auto possible_carbene_atom = [](const OBAtom &atom) -> bool
            {
                return molgr::utils::HasUnresolvedTwoElectronCenter(atom);
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
                        if (molgr::utils::GetUnpairedElectronCount(*nbr) == 1)
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
                            molgr::utils::GetUnpairedElectronCount(*nbr) == 0 &&
                            !molgr::utils::HasUnresolvedTwoElectronCenter(*nbr))
                        {
                            OBBond *bond = mol.GetBond(atom, &*nbr);
                            bond->SetBondOrder(bond->GetBondOrder() + 1);
                            molgr::utils::SetLonePairCount(*atom, 0);
                            molgr::utils::SetUnpairedElectronCount(*atom, 0);
                            molgr::utils::SetUnresolvedTwoElectronCenter(*atom, false);
                            atom->SetFormalCharge(atom->GetFormalCharge() - 1);
                            if (molgr::utils::GetLonePairCount(*nbr) > 0)
                            {
                                molgr::utils::SetLonePairCount(
                                    *nbr,
                                    molgr::utils::GetLonePairCount(*nbr) - 1);
                            }
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

        // Electron bookkeeping: consume exactly one terminal unpaired electron
        // to raise the adjacent bond order and shift one formal-charge unit.
        bool Eliminate13DipolePostive(OBMol &mol, int &charge)
        {
            bool hit = false;
            auto matches = molgr::smarts::FindAll(
                mol,
                molgr::smarts::PatternId::ELIM_1_3_DIPOLE_POSTIVE);
            while (charge >= 0 && !matches.empty())
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
                    molgr::utils::GetUnpairedElectronCount(*atom3) == 1 &&
                    (info->num_outer_electrons + atom2->GetTotalValence()) == 8)
                {
                    atom2->SetFormalCharge(atom2->GetFormalCharge() + 1);
                    bond23->SetBondOrder(static_cast<int>(bond23->GetBondOrder() + 1));
                    molgr::utils::SetUnpairedElectronCount(*atom3, molgr::utils::GetUnpairedElectronCount(*atom3) - 1);
                    charge -= 1;
                    hit = true;
                }
            }
            return hit;
        }

        // Electron bookkeeping: consume a CP-like monoradical only when it is
        // not reserved for the target radical count or pending charge budget.
        bool EliminatePossibleCPLikeRadicalAnion(
            OBMol &mol,
            int &charge,
            int total_radical_electrons)
        {
            const auto available_unpaired_electrons = [&]()
            {
                int count = 0;
                FOR_ATOMS_OF_MOL(atom_iter, mol)
                {
                    count += molgr::utils::GetUnpairedElectronCount(*atom_iter);
                }
                return count - total_radical_electrons - std::abs(charge);
            };

            bool hit = false;
            while (true)
            {
                auto matches = molgr::smarts::FindAll(mol, molgr::smarts::PatternId::ELIM_NEGATIVE_CP);
                if (matches.empty() || available_unpaired_electrons() < 1)
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

                if (molgr::utils::GetUnpairedElectronCount(*atom) != 1)
                {
                    break;
                }
                molgr::utils::SetUnpairedElectronCount(*atom, 0);
                atom->SetFormalCharge(-1);
                charge += 1;
                hit = true;
            }
            return hit;
        }

        // Electron bookkeeping: encode real radical actions or a pure unresolved
        // 2e center as positive charge. The latter requires at least two remaining
        // charge units, becomes +2, and clears its marker atomically.
        bool EliminatePositiveChargesWithTarget(OBMol &mol, int &charge, int total_radical_electrons)
        {
            bool hit = false;
            while (true)
            {
                const auto actions = PositiveChargeAssignmentActions(mol, charge, total_radical_electrons);
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

        bool EliminatePositiveCharges(OBMol &mol, int &charge)
        {
            return EliminatePositiveChargesWithTarget(mol, charge, 0);
        }

        // Electron bookkeeping: encode real radical actions or a pure unresolved
        // 2e center as negative charge. The latter requires at least two remaining
        // charge units and becomes -2; zero-budget separation remains radical-only.
        bool EliminateNegativeCharges(OBMol &mol, int &charge, int total_radical_electrons)
        {
            bool hit = false;
            while (charge <= 0)
            {
                if (total_radical_electrons >= 0 && charge == 0)
                {
                    int real_unpaired = 0;
                    FOR_ATOMS_OF_MOL(atom_iter, mol)
                    {
                        real_unpaired += molgr::utils::GetUnpairedElectronCount(*atom_iter);
                    }
                    if (real_unpaired <= total_radical_electrons)
                    {
                        break;
                    }
                }

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
