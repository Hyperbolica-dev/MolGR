#include "molgr/pipeline/resonance.h"

#include "molgr/state.h"
#include "molgr/utils/consts.h"
#include "molgr/utils/scoring.h"
#include "molgr/utils/smarts.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

#include <algorithm>
#include <deque>
#include <memory>
#include <queue>
#include <set>
#include <sstream>

namespace
{
    constexpr double kDirectGainConjugationScoreWeight = 4.0;
    constexpr double kDirectGainDeviationScoreWeight = 10.0;
    constexpr double kDirectGainBranchBoundStepSlack = 10.0;

    using BondOrderOverrides = std::map<std::pair<int, int>, int>;

    struct IndexedResonanceTraversalMove
    {
        std::tuple<int, int, int> idxs;
        molgr::reconstruct::ResonanceStateKey next_state_key;
    };

    struct QueueEntry
    {
        int discrepancy = 0;
        int depth = 0;
        int order = 0;
        std::shared_ptr<OpenBabel::OBMol> omol;
        molgr::reconstruct::ResonanceStateKey state_key;
    };

    struct QueueEntryGreater
    {
        bool operator()(const QueueEntry &lhs, const QueueEntry &rhs) const
        {
            return std::tie(lhs.discrepancy, lhs.depth, lhs.order) >
                   std::tie(rhs.discrepancy, rhs.depth, rhs.order);
        }
    };

    std::pair<int, int> BondPair(int atom_idx_a, int atom_idx_b)
    {
        if (atom_idx_a <= atom_idx_b)
        {
            return {atom_idx_a, atom_idx_b};
        }
        return {atom_idx_b, atom_idx_a};
    }

    int GetBondOrderWithOverrides(
        const OpenBabel::OBBond &bond,
        const BondOrderOverrides &bond_order_overrides)
    {
        const int begin_idx = bond.GetBeginAtom()->GetIdx();
        const int end_idx = bond.GetEndAtom()->GetIdx();
        const auto it = bond_order_overrides.find(BondPair(begin_idx, end_idx));
        if (it != bond_order_overrides.end())
        {
            return it->second;
        }
        return bond.GetBondOrder();
    }

    int ConjugatedBondKind(
        const OpenBabel::OBBond &bond,
        const BondOrderOverrides &bond_order_overrides)
    {
        if (bond.IsAromatic())
        {
            return 2;
        }
        const int bond_order = GetBondOrderWithOverrides(bond, bond_order_overrides);
        if (bond_order <= 0)
        {
            return 0;
        }
        if (bond_order == 1)
        {
            return 1;
        }
        return 2;
    }

    double CalculateRadicalPenalty(const OpenBabel::OBAtom &atom)
    {
        const int radical_num = atom.GetSpinMultiplicity();
        if (radical_num == 0)
        {
            return 0.0;
        }
        if (molgr::kHeteroatoms.find(atom.GetAtomicNum()) != molgr::kHeteroatoms.end())
        {
            return radical_num * 10.0;
        }
        return (3.0 - static_cast<double>(atom.GetHvyDegree())) * 1.5;
    }

    int EstimateRadicalConjugationSize(
        const OpenBabel::OBMol &mol,
        int radical_atom_idx,
        const BondOrderOverrides &bond_order_overrides = {})
    {
        std::set<std::pair<int, int>> visited_states{{radical_atom_idx, 0}};
        std::set<int> visited_atoms{radical_atom_idx};
        std::vector<std::pair<int, int>> frontier{{radical_atom_idx, 0}};

        while (!frontier.empty())
        {
            const auto [atom_idx, previous_bond_kind] = frontier.back();
            frontier.pop_back();
            OpenBabel::OBAtom *atom = const_cast<OpenBabel::OBMol &>(mol).GetAtom(atom_idx);
            if (atom == nullptr)
            {
                continue;
            }

            FOR_BONDS_OF_ATOM(bond_iter, atom)
            {
                const int bond_kind = ConjugatedBondKind(*bond_iter, bond_order_overrides);
                if (bond_kind == 0)
                {
                    continue;
                }
                if (previous_bond_kind != 0 && bond_kind == previous_bond_kind)
                {
                    continue;
                }

                const int begin_idx = bond_iter->GetBeginAtom()->GetIdx();
                const int end_idx = bond_iter->GetEndAtom()->GetIdx();
                const int neighbor_idx = begin_idx == atom_idx ? end_idx : begin_idx;
                const std::pair<int, int> state{neighbor_idx, bond_kind};
                if (visited_states.find(state) != visited_states.end())
                {
                    continue;
                }
                visited_states.insert(state);
                visited_atoms.insert(neighbor_idx);
                frontier.push_back(state);
            }
        }

        return static_cast<int>(visited_atoms.size());
    }

    double CalculateAllDoubleBondCarbonBonus(
        const OpenBabel::OBMol &mol,
        int atom_idx,
        const BondOrderOverrides &bond_order_overrides = {})
    {
        OpenBabel::OBAtom *atom = const_cast<OpenBabel::OBMol &>(mol).GetAtom(atom_idx);
        if (atom == nullptr || atom->GetAtomicNum() != 6)
        {
            return 0.0;
        }

        bool saw_bond = false;
        FOR_BONDS_OF_ATOM(bond_iter, atom)
        {
            saw_bond = true;
            if (GetBondOrderWithOverrides(*bond_iter, bond_order_overrides) != 2)
            {
                return 0.0;
            }
        }
        return saw_bond ? 5.0 : 0.0;
    }

    molgr::reconstruct::DirectGainMetrics ComputeDirectGainResonanceMetrics(
        const OpenBabel::OBMol &mol,
        const std::tuple<int, int, int> &move_path)
    {
        const int old_radical_idx = std::get<0>(move_path);
        const int center_idx = std::get<1>(move_path);
        const int new_radical_idx = std::get<2>(move_path);

        OpenBabel::OBBond *bond_old_center =
            const_cast<OpenBabel::OBMol &>(mol).GetBond(old_radical_idx, center_idx);
        OpenBabel::OBBond *bond_center_new =
            const_cast<OpenBabel::OBMol &>(mol).GetBond(center_idx, new_radical_idx);
        OpenBabel::OBAtom *old_atom = const_cast<OpenBabel::OBMol &>(mol).GetAtom(old_radical_idx);
        OpenBabel::OBAtom *new_atom = const_cast<OpenBabel::OBMol &>(mol).GetAtom(new_radical_idx);
        if (bond_old_center == nullptr || bond_center_new == nullptr || old_atom == nullptr || new_atom == nullptr)
        {
            return {0.0, 0.0, 0.0, 0.0};
        }

        const BondOrderOverrides bond_order_overrides = {
            {BondPair(old_radical_idx, center_idx), bond_old_center->GetBondOrder() + 1},
            {BondPair(center_idx, new_radical_idx), bond_center_new->GetBondOrder() - 1},
        };

        const double conjugation_gain =
            static_cast<double>(EstimateRadicalConjugationSize(mol, new_radical_idx, bond_order_overrides) -
                                EstimateRadicalConjugationSize(mol, old_radical_idx));
        const double deviation_gain =
            molgr::scoring::GetDeviationScore(mol, old_atom) -
            molgr::scoring::GetDeviationScore(mol, new_atom);
        const double radical_penalty_gain =
            CalculateRadicalPenalty(*old_atom) - CalculateRadicalPenalty(*new_atom);
        const double double_bond_bonus_gain =
            CalculateAllDoubleBondCarbonBonus(mol, new_radical_idx, bond_order_overrides) -
            CalculateAllDoubleBondCarbonBonus(mol, old_radical_idx);
        return {
            conjugation_gain,
            deviation_gain,
            radical_penalty_gain,
            double_bond_bonus_gain,
        };
    }

    bool HasPositiveDirectGain(const molgr::reconstruct::DirectGainMetrics &metrics)
    {
        return std::any_of(metrics.begin(), metrics.end(), [](double value) { return value > 0.0; });
    }

    bool DirectGainMoveLess(
        const std::pair<IndexedResonanceTraversalMove, molgr::reconstruct::DirectGainMetrics> &lhs,
        const std::pair<IndexedResonanceTraversalMove, molgr::reconstruct::DirectGainMetrics> &rhs)
    {
        const auto lhs_key = std::make_tuple(
            HasPositiveDirectGain(lhs.second) ? 0 : 1,
            -lhs.second[0],
            -lhs.second[1],
            -lhs.second[2],
            -lhs.second[3],
            lhs.first.idxs);
        const auto rhs_key = std::make_tuple(
            HasPositiveDirectGain(rhs.second) ? 0 : 1,
            -rhs.second[0],
            -rhs.second[1],
            -rhs.second[2],
            -rhs.second[3],
            rhs.first.idxs);
        return lhs_key < rhs_key;
    }

    std::vector<std::pair<IndexedResonanceTraversalMove, molgr::reconstruct::DirectGainMetrics>>
    OrderDirectGainMoves(
        const OpenBabel::OBMol &mol,
        const std::vector<IndexedResonanceTraversalMove> &moves,
        bool positive_only)
    {
        std::vector<std::pair<IndexedResonanceTraversalMove, molgr::reconstruct::DirectGainMetrics>>
            ordered_moves;
        for (const auto &move : moves)
        {
            const auto metrics = ComputeDirectGainResonanceMetrics(mol, move.idxs);
            if (positive_only && !HasPositiveDirectGain(metrics))
            {
                continue;
            }
            ordered_moves.emplace_back(move, metrics);
        }
        std::sort(ordered_moves.begin(), ordered_moves.end(), DirectGainMoveLess);
        return ordered_moves;
    }

    molgr::reconstruct::DirectGainMetrics AddDirectGainMetrics(
        const molgr::reconstruct::DirectGainMetrics &left,
        const molgr::reconstruct::DirectGainMetrics &right)
    {
        molgr::reconstruct::DirectGainMetrics result{};
        for (std::size_t idx = 0; idx < result.size(); ++idx)
        {
            result[idx] = left[idx] + right[idx];
        }
        return result;
    }

    molgr::reconstruct::DirectGainMetrics ComponentwiseMaxDirectGainMetrics(
        const molgr::reconstruct::DirectGainMetrics &left,
        const molgr::reconstruct::DirectGainMetrics &right)
    {
        molgr::reconstruct::DirectGainMetrics result{};
        for (std::size_t idx = 0; idx < result.size(); ++idx)
        {
            result[idx] = std::max(left[idx], right[idx]);
        }
        return result;
    }

    std::pair<molgr::reconstruct::ResonanceStateKey, molgr::reconstruct::ResonanceBondIndexMap>
    BuildResonanceSearchContext(const OpenBabel::OBMol &mol)
    {
        const auto state_key = molgr::reconstruct::BuildResonanceStateKey(mol);
        molgr::reconstruct::ResonanceBondIndexMap bond_index_map;
        for (std::size_t idx = 0; idx < state_key.bond_keys.size(); ++idx)
        {
            const auto &bond_key = state_key.bond_keys[idx];
            bond_index_map[{
                std::get<0>(bond_key),
                std::get<1>(bond_key),
            }] = idx;
        }
        return {state_key, bond_index_map};
    }

    molgr::reconstruct::ResonanceBondIndexMap BuildBondIndexMapFromStateKey(
        const molgr::reconstruct::ResonanceStateKey &state_key)
    {
        molgr::reconstruct::ResonanceBondIndexMap bond_index_map;
        for (std::size_t idx = 0; idx < state_key.bond_keys.size(); ++idx)
        {
            const auto &bond_key = state_key.bond_keys[idx];
            bond_index_map[{
                std::get<0>(bond_key),
                std::get<1>(bond_key),
            }] = idx;
        }
        return bond_index_map;
    }

    molgr::reconstruct::ResonanceStateKey IncrementResonanceStateKey(
        const molgr::reconstruct::ResonanceStateKey &state_key,
        const molgr::reconstruct::ResonanceBondIndexMap &bond_index_map,
        const std::tuple<int, int, int> &idxs)
    {
        auto next_state_key = state_key;

        const std::size_t atom1_idx = static_cast<std::size_t>(std::get<0>(idxs) - 1);
        const std::size_t atom3_idx = static_cast<std::size_t>(std::get<2>(idxs) - 1);

        auto atom1 = next_state_key.atom_keys[atom1_idx];
        auto atom3 = next_state_key.atom_keys[atom3_idx];
        std::get<2>(atom1) -= 1;
        std::get<2>(atom3) += 1;
        next_state_key.atom_keys[atom1_idx] = atom1;
        next_state_key.atom_keys[atom3_idx] = atom3;

        const auto bond1_it =
            bond_index_map.find(BondPair(std::get<0>(idxs), std::get<1>(idxs)));
        const auto bond2_it =
            bond_index_map.find(BondPair(std::get<1>(idxs), std::get<2>(idxs)));
        if (bond1_it == bond_index_map.end() || bond2_it == bond_index_map.end())
        {
            return next_state_key;
        }

        auto bond1 = next_state_key.bond_keys[bond1_it->second];
        auto bond2 = next_state_key.bond_keys[bond2_it->second];
        std::get<2>(bond1) += 1;
        std::get<2>(bond2) -= 1;
        next_state_key.bond_keys[bond1_it->second] = bond1;
        next_state_key.bond_keys[bond2_it->second] = bond2;
        return next_state_key;
    }

    std::vector<IndexedResonanceTraversalMove> EnumerateOneStepResonanceMoves(
        const OpenBabel::OBMol &mol,
        const molgr::reconstruct::ResonanceStateKey &state_key,
        const molgr::reconstruct::ResonanceBondIndexMap &bond_index_map)
    {
        OpenBabel::OBMol &query_mol = const_cast<OpenBabel::OBMol &>(mol);
        const auto matches =
            molgr::smarts::Match(query_mol, molgr::smarts::PatternId::RESONANCE_ONE_STEP);
        std::vector<IndexedResonanceTraversalMove> result;
        result.reserve(matches.size());

        for (const auto &idxs_vec : matches)
        {
            if (idxs_vec.size() < 3)
            {
                continue;
            }

            const std::tuple<int, int, int> idxs{idxs_vec[0], idxs_vec[1], idxs_vec[2]};
            OpenBabel::OBAtom *atom1 = query_mol.GetAtom(std::get<0>(idxs));
            OpenBabel::OBAtom *atom3 = query_mol.GetAtom(std::get<2>(idxs));
            OpenBabel::OBBond *bond1 = query_mol.GetBond(std::get<0>(idxs), std::get<1>(idxs));
            OpenBabel::OBBond *bond2 = query_mol.GetBond(std::get<1>(idxs), std::get<2>(idxs));
            if (atom1 == nullptr || atom3 == nullptr || bond1 == nullptr || bond2 == nullptr)
            {
                continue;
            }

            if (atom1->GetSpinMultiplicity() == 1 &&
                atom3->GetSpinMultiplicity() == 0 &&
                bond1->GetBondOrder() <= 2 &&
                bond2->GetBondOrder() >= 2)
            {
                result.push_back(
                    IndexedResonanceTraversalMove{
                        idxs,
                        IncrementResonanceStateKey(state_key, bond_index_map, idxs),
                    });
            }
        }

        return result;
    }

    OpenBabel::OBMol MaterializeOneStepResonance(
        const OpenBabel::OBMol &mol,
        const std::tuple<int, int, int> &idxs)
    {
        OpenBabel::OBMol new_mol(mol);
        OpenBabel::OBAtom *atom1 = new_mol.GetAtom(std::get<0>(idxs));
        OpenBabel::OBAtom *atom3 = new_mol.GetAtom(std::get<2>(idxs));
        OpenBabel::OBBond *bond1 = new_mol.GetBond(std::get<0>(idxs), std::get<1>(idxs));
        OpenBabel::OBBond *bond2 = new_mol.GetBond(std::get<1>(idxs), std::get<2>(idxs));
        if (atom1 == nullptr || atom3 == nullptr || bond1 == nullptr || bond2 == nullptr)
        {
            return new_mol;
        }

        bond1->SetBondOrder(bond1->GetBondOrder() + 1);
        bond2->SetBondOrder(bond2->GetBondOrder() - 1);
        atom1->SetSpinMultiplicity(atom1->GetSpinMultiplicity() - 1);
        atom3->SetSpinMultiplicity(atom3->GetSpinMultiplicity() + 1);
        return new_mol;
    }

    std::shared_ptr<OpenBabel::OBMol> MaterializeOneStepResonancePtr(
        const OpenBabel::OBMol &mol,
        const std::tuple<int, int, int> &idxs)
    {
        auto new_mol = std::make_shared<OpenBabel::OBMol>(mol);
        OpenBabel::OBAtom *atom1 = new_mol->GetAtom(std::get<0>(idxs));
        OpenBabel::OBAtom *atom3 = new_mol->GetAtom(std::get<2>(idxs));
        OpenBabel::OBBond *bond1 = new_mol->GetBond(std::get<0>(idxs), std::get<1>(idxs));
        OpenBabel::OBBond *bond2 = new_mol->GetBond(std::get<1>(idxs), std::get<2>(idxs));
        if (atom1 == nullptr || atom3 == nullptr || bond1 == nullptr || bond2 == nullptr)
        {
            return new_mol;
        }

        bond1->SetBondOrder(bond1->GetBondOrder() + 1);
        bond2->SetBondOrder(bond2->GetBondOrder() - 1);
        atom1->SetSpinMultiplicity(atom1->GetSpinMultiplicity() - 1);
        atom3->SetSpinMultiplicity(atom3->GetSpinMultiplicity() + 1);
        return new_mol;
    }

    std::vector<IndexedResonanceTraversalMove> SelectLimitedDiscrepancyMoves(
        const OpenBabel::OBMol &mol,
        const std::vector<IndexedResonanceTraversalMove> &moves,
        const molgr::reconstruct::LimitedDiscrepancyTraversalConfig &config)
    {
        const auto positive_moves = OrderDirectGainMoves(mol, moves, true);
        if (!positive_moves.empty())
        {
            std::vector<IndexedResonanceTraversalMove> selected_moves;
            selected_moves.reserve(positive_moves.size());
            for (const auto &entry : positive_moves)
            {
                selected_moves.push_back(entry.first);
            }
            return selected_moves;
        }

        if (config.fallback_to_full_frontier)
        {
            const auto ordered_moves = OrderDirectGainMoves(mol, moves, false);
            std::vector<IndexedResonanceTraversalMove> selected_moves;
            selected_moves.reserve(ordered_moves.size());
            for (const auto &entry : ordered_moves)
            {
                selected_moves.push_back(entry.first);
            }
            return selected_moves;
        }

        return {};
    }

    molgr::reconstruct::DirectGainMetrics ComputeNStepDirectGainUpperBound(
        const OpenBabel::OBMol &mol,
        const molgr::reconstruct::ResonanceStateKey &state_key,
        int remaining_steps,
        molgr::reconstruct::DirectGainBoundCache &cache)
    {
        if (remaining_steps <= 0)
        {
            return {0.0, 0.0, 0.0, 0.0};
        }

        const molgr::reconstruct::DirectGainBoundCacheKey cache_key{state_key, remaining_steps};
        const auto cache_it = cache.find(cache_key);
        if (cache_it != cache.end())
        {
            return cache_it->second;
        }

        molgr::reconstruct::DirectGainMetrics upper_bound{0.0, 0.0, 0.0, 0.0};
        const auto bond_index_map = BuildBondIndexMapFromStateKey(state_key);
        const auto moves = EnumerateOneStepResonanceMoves(mol, state_key, bond_index_map);
        for (const auto &move : moves)
        {
            const auto direct_metrics = ComputeDirectGainResonanceMetrics(mol, move.idxs);
            if (!HasPositiveDirectGain(direct_metrics))
            {
                continue;
            }

            auto total_metrics = direct_metrics;
            if (remaining_steps > 1)
            {
                auto next_omol = MaterializeOneStepResonance(mol, move.idxs);
                const auto future_upper_bound = ComputeNStepDirectGainUpperBound(
                    next_omol,
                    move.next_state_key,
                    remaining_steps - 1,
                    cache);
                total_metrics = AddDirectGainMetrics(direct_metrics, future_upper_bound);
            }

            upper_bound = ComponentwiseMaxDirectGainMetrics(upper_bound, total_metrics);
        }

        cache[cache_key] = upper_bound;
        return upper_bound;
    }

    double EstimateDirectGainScoreImprovementUpperBound(
        const molgr::reconstruct::DirectGainMetrics &metrics,
        int remaining_steps)
    {
        return std::max(
            0.0,
            metrics[0] * kDirectGainConjugationScoreWeight +
                metrics[1] * kDirectGainDeviationScoreWeight +
                metrics[2] +
                metrics[3] +
                static_cast<double>(std::max(remaining_steps, 0)) * kDirectGainBranchBoundStepSlack);
    }
}

namespace molgr
{
    namespace reconstruct
    {
        bool ResonanceStateKey::operator==(const ResonanceStateKey &other) const
        {
            return atom_keys == other.atom_keys && bond_keys == other.bond_keys;
        }

        bool ResonanceStateKey::operator<(const ResonanceStateKey &other) const
        {
            return std::tie(atom_keys, bond_keys) < std::tie(other.atom_keys, other.bond_keys);
        }

        std::string SmilesFirstToken(const OpenBabel::OBMol &mol)
        {
            thread_local OpenBabel::OBConversion conv;
            thread_local bool initialized = false;
            if (!initialized)
            {
                conv.SetOutFormat("smi");
                initialized = true;
            }
            OpenBabel::OBMol temp(mol);
            const std::string smi = conv.WriteString(&temp, true);
            std::istringstream iss(smi);
            std::string token;
            iss >> token;
            return token;
        }

        ResonanceStateKey BuildResonanceStateKey(const OpenBabel::OBMol &mol)
        {
            ResonanceStateKey state_key;
            OpenBabel::OBMol &mutable_mol = const_cast<OpenBabel::OBMol &>(mol);

            state_key.atom_keys.reserve(static_cast<std::size_t>(mutable_mol.NumAtoms()));
            FOR_ATOMS_OF_MOL(atom_iter, mutable_mol)
            {
                const OpenBabel::OBAtom &atom = *atom_iter;
                state_key.atom_keys.emplace_back(
                    atom.GetAtomicNum(),
                    atom.GetFormalCharge(),
                    atom.GetSpinMultiplicity(),
                    atom.IsAromatic());
            }

            state_key.bond_keys.reserve(static_cast<std::size_t>(mutable_mol.NumBonds()));
            FOR_BONDS_OF_MOL(bond_iter, mutable_mol)
            {
                int begin_idx = bond_iter->GetBeginAtom()->GetIdx();
                int end_idx = bond_iter->GetEndAtom()->GetIdx();
                if (begin_idx > end_idx)
                {
                    std::swap(begin_idx, end_idx);
                }
                state_key.bond_keys.emplace_back(
                    begin_idx,
                    end_idx,
                    bond_iter->GetBondOrder(),
                    bond_iter->IsAromatic());
            }
            return state_key;
        }

        ProcessedResonanceKey BuildProcessedResonanceKey(const OpenBabel::OBMol &mol)
        {
            return BuildResonanceStateKey(mol);
        }

        std::vector<OpenBabel::OBMol> GetOneStepResonance(const OpenBabel::OBMol &mol)
        {
            const auto [state_key, bond_index_map] = BuildResonanceSearchContext(mol);
            const auto moves = EnumerateOneStepResonanceMoves(mol, state_key, bond_index_map);
            std::vector<OpenBabel::OBMol> result;
            result.reserve(moves.size());
            for (const auto &move : moves)
            {
                result.push_back(MaterializeOneStepResonance(mol, move.idxs));
            }
            return result;
        }

        void WalkRadicalResonances(
            const OpenBabel::OBMol &mol,
            int max_depth,
            ResonanceVisitCallback visit)
        {
            if (!visit)
            {
                visit = [](const ResonanceSearchNode &) { return true; };
            }

            const auto [root_key, bond_index_map] = BuildResonanceSearchContext(mol);
            std::set<ResonanceStateKey> seen{root_key};
            std::deque<std::tuple<std::shared_ptr<OpenBabel::OBMol>, ResonanceStateKey, int>> frontier{
                {std::make_shared<OpenBabel::OBMol>(mol), root_key, 0}};

            while (!frontier.empty())
            {
                auto [current, current_key, depth] = std::move(frontier.front());
                frontier.pop_front();

                const bool should_expand = visit(ResonanceSearchNode{
                    *current,
                    current_key,
                    depth,
                });
                if (depth >= max_depth || !should_expand)
                {
                    continue;
                }

                const auto moves = EnumerateOneStepResonanceMoves(*current, current_key, bond_index_map);
                for (const auto &move : moves)
                {
                    if (seen.find(move.next_state_key) != seen.end())
                    {
                        continue;
                    }
                    seen.insert(move.next_state_key);
                    frontier.emplace_back(
                        MaterializeOneStepResonancePtr(*current, move.idxs),
                        move.next_state_key,
                        depth + 1);
                }
            }
        }

        void WalkRadicalResonancesLimitedDiscrepancy(
            const OpenBabel::OBMol &mol,
            int max_depth,
            ResonanceVisitCallback visit,
            const LimitedDiscrepancyTraversalConfig &config)
        {
            if (!visit)
            {
                visit = [](const ResonanceSearchNode &) { return true; };
            }

            const auto [root_key, bond_index_map] = BuildResonanceSearchContext(mol);
            std::map<ResonanceStateKey, int> best_discrepancy_by_state{{root_key, 0}};
            std::set<ResonanceStateKey> emitted_states;
            std::priority_queue<QueueEntry, std::vector<QueueEntry>, QueueEntryGreater> frontier;
            frontier.push(QueueEntry{0, 0, 0, std::make_shared<OpenBabel::OBMol>(mol), root_key});

            int push_order = 0;
            while (!frontier.empty())
            {
                QueueEntry current_entry = std::move(frontier.top());
                frontier.pop();

                const auto best_it = best_discrepancy_by_state.find(current_entry.state_key);
                if (best_it == best_discrepancy_by_state.end() ||
                    current_entry.discrepancy != best_it->second)
                {
                    continue;
                }

                bool should_expand = true;
                if (emitted_states.find(current_entry.state_key) == emitted_states.end())
                {
                    emitted_states.insert(current_entry.state_key);
                    should_expand = visit(ResonanceSearchNode{
                        *current_entry.omol,
                        current_entry.state_key,
                        current_entry.depth,
                    });
                }

                if (current_entry.depth >= max_depth || !should_expand)
                {
                    continue;
                }

                const auto moves =
                    EnumerateOneStepResonanceMoves(*current_entry.omol, current_entry.state_key, bond_index_map);
                const auto selected_moves =
                    SelectLimitedDiscrepancyMoves(*current_entry.omol, moves, config);
                for (std::size_t move_rank = 0; move_rank < selected_moves.size(); ++move_rank)
                {
                    const int next_discrepancy =
                        current_entry.discrepancy + static_cast<int>(move_rank);
                    if (next_discrepancy > config.max_discrepancy)
                    {
                        break;
                    }

                    const auto &move = selected_moves[move_rank];
                    const auto best_known_it = best_discrepancy_by_state.find(move.next_state_key);
                    if (best_known_it != best_discrepancy_by_state.end() &&
                        best_known_it->second <= next_discrepancy)
                    {
                        continue;
                    }

                    best_discrepancy_by_state[move.next_state_key] = next_discrepancy;
                    frontier.push(QueueEntry{
                        next_discrepancy,
                        current_entry.depth + 1,
                        ++push_order,
                        MaterializeOneStepResonancePtr(*current_entry.omol, move.idxs),
                        move.next_state_key,
                    });
                }
            }
        }

        std::vector<OpenBabel::OBMol> GetRadicalResonances(const OpenBabel::OBMol &mol, int max_depth)
        {
            std::vector<OpenBabel::OBMol> resonances;
            WalkRadicalResonances(
                mol,
                max_depth,
                [&resonances](const ResonanceSearchNode &node)
                {
                    resonances.push_back(node.omol);
                    return true;
                });
            return resonances;
        }

        std::vector<OpenBabel::OBMol> GetRadicalResonancesLimitedDiscrepancy(
            const OpenBabel::OBMol &mol,
            int max_depth,
            const LimitedDiscrepancyTraversalConfig &config)
        {
            std::vector<OpenBabel::OBMol> resonances;
            WalkRadicalResonancesLimitedDiscrepancy(
                mol,
                max_depth,
                [&resonances](const ResonanceSearchNode &node)
                {
                    resonances.push_back(node.omol);
                    return true;
                },
                config);
            return resonances;
        }

        double EstimateRemainingResonanceScoreImprovementUpperBound(
            const OpenBabel::OBMol &mol,
            const ResonanceStateKey &state_key,
            int remaining_steps,
            DirectGainBoundCache *cache)
        {
            if (remaining_steps <= 0)
            {
                return 0.0;
            }

            DirectGainBoundCache local_cache;
            DirectGainBoundCache &cache_ref = cache == nullptr ? local_cache : *cache;
            const auto optimistic_metrics =
                ComputeNStepDirectGainUpperBound(mol, state_key, remaining_steps, cache_ref);
            return EstimateDirectGainScoreImprovementUpperBound(
                optimistic_metrics,
                remaining_steps);
        }

        std::tuple<OpenBabel::OBMol, int, bool> ProcessResonanceDetailed(
            const OpenBabel::OBMol &mol,
            int charge)
        {
            auto machine = molgr::state::OmolStateMachine(
                std::make_shared<OpenBabel::OBMol>(mol),
                charge);
            bool hit = false;
            hit = machine.RunOmolChargeStage(std::nullopt, Eliminate13Dipole) || hit;
            hit = machine.RunOmolChargeStage(std::nullopt, EliminatePositiveCharges) || hit;
            hit = machine.RunOmolChargeStage(std::nullopt, EliminateNegativeCharges) || hit;
            hit = machine.RunOmolStage(std::nullopt, CleanNeighborRadicals) || hit;
            hit = machine.RunOmolStage(std::nullopt, CleanResonances) || hit;
            return std::make_tuple(OpenBabel::OBMol(*machine.omol), machine.given_charge, hit);
        }

        std::pair<OpenBabel::OBMol, int> ProcessResonance(const OpenBabel::OBMol &mol, int charge)
        {
            auto processed = ProcessResonanceDetailed(mol, charge);
            return {std::get<0>(processed), std::get<1>(processed)};
        }
    }
}
