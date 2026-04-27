#include "molgr/utils/resonance.h"

#include "molgr/stages/clean.h"
#include "molgr/stages/eliminate.h"
#include "molgr/state.h"
#include "molgr/utils/consts.h"
#include "molgr/utils/force_field.h"
#include "molgr/utils/lru_cache.h"
#include "molgr/utils/scoring.h"
#include "molgr/utils/smarts.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>

#include <algorithm>
#include <cstddef>
#include <exception>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

namespace
{
    constexpr std::size_t kDefaultResonanceMoveScoreCacheMaxSize = 4096;
    constexpr double kDirectGainConjugationScoreWeight = 4.0;
    constexpr double kDirectGainDeviationScoreWeight = 10.0;
    constexpr double kDirectGainBranchBoundStepSlack = 10.0;

    using BondOrderOverrides = std::map<std::pair<int, int>, int>;

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

    molgr::resonance::DirectGainMetrics ComputeDirectGainResonanceMetrics(
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
        if (bond_old_center == nullptr || bond_center_new == nullptr ||
            old_atom == nullptr || new_atom == nullptr)
        {
            return {0.0, 0.0, 0.0, 0.0};
        }

        const BondOrderOverrides bond_order_overrides = {
            {BondPair(old_radical_idx, center_idx), bond_old_center->GetBondOrder() + 1},
            {BondPair(center_idx, new_radical_idx), bond_center_new->GetBondOrder() - 1},
        };

        const double conjugation_gain =
            static_cast<double>(
                EstimateRadicalConjugationSize(mol, new_radical_idx, bond_order_overrides) -
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

    bool HasPositiveDirectGain(const molgr::resonance::DirectGainMetrics &metrics)
    {
        return std::any_of(metrics.begin(), metrics.end(), [](double value) { return value > 0.0; });
    }

    molgr::resonance::DirectGainMetrics AddDirectGainMetrics(
        const molgr::resonance::DirectGainMetrics &left,
        const molgr::resonance::DirectGainMetrics &right)
    {
        molgr::resonance::DirectGainMetrics result{};
        for (std::size_t idx = 0; idx < result.size(); ++idx)
        {
            result[idx] = left[idx] + right[idx];
        }
        return result;
    }

    molgr::resonance::DirectGainMetrics ComponentwiseMaxDirectGainMetrics(
        const molgr::resonance::DirectGainMetrics &left,
        const molgr::resonance::DirectGainMetrics &right)
    {
        molgr::resonance::DirectGainMetrics result{};
        for (std::size_t idx = 0; idx < result.size(); ++idx)
        {
            result[idx] = std::max(left[idx], right[idx]);
        }
        return result;
    }

    molgr::resonance::DirectGainMetrics ComputeNStepDirectGainUpperBound(
        const OpenBabel::OBMol &mol,
        const molgr::resonance::ResonanceStateKey &state_key,
        int remaining_steps,
        molgr::resonance::DirectGainBoundCache &cache)
    {
        if (remaining_steps <= 0)
        {
            return {0.0, 0.0, 0.0, 0.0};
        }

        const molgr::resonance::DirectGainBoundCacheKey cache_key{state_key, remaining_steps};
        const auto cache_it = cache.find(cache_key);
        if (cache_it != cache.end())
        {
            return cache_it->second;
        }

        molgr::resonance::DirectGainMetrics upper_bound{0.0, 0.0, 0.0, 0.0};
        const auto bond_index_map = molgr::resonance::BuildBondIndexMapFromStateKey(state_key);
        const auto moves =
            molgr::resonance::EnumerateOneStepResonanceMoves(mol, state_key, bond_index_map);
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
                auto next_omol = molgr::resonance::MaterializeOneStepResonance(mol, move.idxs);
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
        const molgr::resonance::DirectGainMetrics &metrics,
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

    double DirectGainMoveScore(const molgr::resonance::DirectGainMetrics &metrics)
    {
        return -(
            metrics[0] * kDirectGainConjugationScoreWeight +
            metrics[1] * kDirectGainDeviationScoreWeight +
            metrics[2] +
            metrics[3]);
    }

    void AppendStringVector(std::string &out, const std::vector<std::string> &values)
    {
        out.push_back('[');
        for (std::size_t idx = 0; idx < values.size(); ++idx)
        {
            if (idx > 0)
            {
                out.push_back(',');
            }
            out += values[idx];
        }
        out.push_back(']');
    }

    std::string BuildSelectionForceFieldCacheConfigKey(
        const molgr::config::MolGRConfig &config)
    {
        std::string key;
        key.reserve(192);
        key += "selection=";
        key += config.force_field.selection_force_field;
        key += ";auto_metal_free=";
        AppendStringVector(key, config.force_field.auto_force_fields_metal_free);
        key += ";auto_with_metals=";
        AppendStringVector(key, config.force_field.auto_force_fields_with_metals);
        return key;
    }

    std::string BuildResonanceMoveScoreCacheKey(
        const OpenBabel::OBMol &mol,
        const std::tuple<int, int, int> &idxs,
        const molgr::config::MolGRConfig &config)
    {
        std::string key = molgr::scoring::BuildScoreKey(mol);
        key += "|move=";
        key += std::to_string(std::get<0>(idxs));
        key.push_back(',');
        key += std::to_string(std::get<1>(idxs));
        key.push_back(',');
        key += std::to_string(std::get<2>(idxs));
        key.push_back('|');
        key += BuildSelectionForceFieldCacheConfigKey(config);
        return key;
    }

    std::string SerializeResonanceStateKey(
        const molgr::resonance::ResonanceStateKey &state_key)
    {
        std::string key;
        key.reserve(
            32 +
            state_key.atom_keys.size() * 16 +
            state_key.bond_keys.size() * 18);

        key += "A";
        key += std::to_string(state_key.atom_keys.size());
        key.push_back(':');
        for (const auto &atom_key : state_key.atom_keys)
        {
            key += std::to_string(std::get<0>(atom_key));
            key.push_back(',');
            key += std::to_string(std::get<1>(atom_key));
            key.push_back(',');
            key += std::to_string(std::get<2>(atom_key));
            key.push_back(',');
            key.push_back(std::get<3>(atom_key) ? '1' : '0');
            key.push_back(';');
        }

        key += "|B";
        key += std::to_string(state_key.bond_keys.size());
        key.push_back(':');
        for (const auto &bond_key : state_key.bond_keys)
        {
            key += std::to_string(std::get<0>(bond_key));
            key.push_back(',');
            key += std::to_string(std::get<1>(bond_key));
            key.push_back(',');
            key += std::to_string(std::get<2>(bond_key));
            key.push_back(',');
            key.push_back(std::get<3>(bond_key) ? '1' : '0');
            key.push_back(';');
        }
        return key;
    }

    molgr::utils::StringLruCache<double> &ResonanceMoveScoreCache()
    {
        static molgr::utils::StringLruCache<double> cache(
            kDefaultResonanceMoveScoreCacheMaxSize);
        return cache;
    }
}

namespace molgr
{
    namespace resonance
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
            thread_local OpenBabel::OBConversion *conv = nullptr;
            if (conv == nullptr)
            {
                conv = new OpenBabel::OBConversion();
                conv->SetOutFormat("smi");
            }
            OpenBabel::OBMol temp(mol);
            const std::string smi = conv->WriteString(&temp, true);
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
            return SerializeResonanceStateKey(BuildResonanceStateKey(mol));
        }

        std::pair<ResonanceStateKey, ResonanceBondIndexMap> BuildResonanceSearchContext(
            const OpenBabel::OBMol &mol)
        {
            const auto state_key = BuildResonanceStateKey(mol);
            return {state_key, BuildBondIndexMapFromStateKey(state_key)};
        }

        ResonanceBondIndexMap BuildBondIndexMapFromStateKey(
            const ResonanceStateKey &state_key)
        {
            ResonanceBondIndexMap bond_index_map;
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

        ResonanceStateKey IncrementResonanceStateKey(
            const ResonanceStateKey &state_key,
            const ResonanceBondIndexMap &bond_index_map,
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
            const ResonanceStateKey &state_key,
            const ResonanceBondIndexMap &bond_index_map)
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

        double ScoreOneStepResonanceWithForceField(
            const OpenBabel::OBMol &mol,
            const std::tuple<int, int, int> &idxs,
            const molgr::config::MolGRConfig &config)
        {
            const std::string cache_key = BuildResonanceMoveScoreCacheKey(mol, idxs, config);
            double cached_score = 0.0;
            if (ResonanceMoveScoreCache().Get(cache_key, cached_score))
            {
                return cached_score;
            }

            double score = std::numeric_limits<double>::infinity();
            const OpenBabel::OBMol moved_mol = MaterializeOneStepResonance(mol, idxs);
            try
            {
                score = molgr::scoring::SelectionForceFieldEnergy(moved_mol, config);
            }
            catch (const std::exception &)
            {
                score = std::numeric_limits<double>::infinity();
            }
            ResonanceMoveScoreCache().Put(cache_key, score);
            return score;
        }

        std::tuple<std::size_t, std::size_t, std::size_t> ResonanceMoveScoreCacheInfo()
        {
            return ResonanceMoveScoreCache().Info();
        }

        void ResonanceMoveScoreCacheClear()
        {
            ResonanceMoveScoreCache().Clear();
        }

        std::vector<std::pair<IndexedResonanceTraversalMove, double>> OrderForceFieldMoves(
            const OpenBabel::OBMol &mol,
            const std::vector<IndexedResonanceTraversalMove> &moves,
            const molgr::config::MolGRConfig &config)
        {
            std::vector<std::pair<IndexedResonanceTraversalMove, double>> ordered_moves;
            ordered_moves.reserve(moves.size());
            for (const auto &move : moves)
            {
                ordered_moves.emplace_back(
                    move,
                    ScoreOneStepResonanceWithForceField(mol, move.idxs, config));
            }
            std::sort(
                ordered_moves.begin(),
                ordered_moves.end(),
                [](const auto &lhs, const auto &rhs)
                {
                    return std::tie(lhs.second, lhs.first.idxs) <
                           std::tie(rhs.second, rhs.first.idxs);
                });
            return ordered_moves;
        }

        std::vector<std::pair<IndexedResonanceTraversalMove, double>> OrderDirectGainMoves(
            const OpenBabel::OBMol &mol,
            const std::vector<IndexedResonanceTraversalMove> &moves)
        {
            std::vector<std::pair<IndexedResonanceTraversalMove, double>> ordered_moves;
            ordered_moves.reserve(moves.size());
            for (const auto &move : moves)
            {
                ordered_moves.emplace_back(
                    move,
                    DirectGainMoveScore(ComputeDirectGainResonanceMetrics(mol, move.idxs)));
            }
            std::sort(
                ordered_moves.begin(),
                ordered_moves.end(),
                [](const auto &lhs, const auto &rhs)
                {
                    return std::tie(lhs.second, lhs.first.idxs) <
                           std::tie(rhs.second, rhs.first.idxs);
                });
            return ordered_moves;
        }

        std::vector<IndexedResonanceTraversalMove> SelectLimitedDiscrepancyMoves(
            const OpenBabel::OBMol &mol,
            const std::vector<IndexedResonanceTraversalMove> &moves,
            const LimitedDiscrepancyTraversalConfig &traversal_config,
            const molgr::config::MolGRConfig &config)
        {
            (void)traversal_config;
            if (config.resonance.traversal_score == "input_order")
            {
                return moves;
            }

            std::vector<std::pair<IndexedResonanceTraversalMove, double>> ordered_moves;
            if (config.resonance.traversal_score == "direct_gain")
            {
                ordered_moves = OrderDirectGainMoves(mol, moves);
            }
            else if (config.resonance.traversal_score == "force_field")
            {
                ordered_moves = OrderForceFieldMoves(mol, moves, config);
            }
            else
            {
                throw std::invalid_argument(
                    "Unsupported resonance traversal_score. Expected 'force_field', "
                    "'direct_gain', or 'input_order'.");
            }

            std::vector<IndexedResonanceTraversalMove> selected_moves;
            selected_moves.reserve(ordered_moves.size());
            for (const auto &entry : ordered_moves)
            {
                selected_moves.push_back(entry.first);
            }
            return selected_moves;
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
            hit = machine.RunOmolChargeStage(std::nullopt, molgr::reconstruct::Eliminate13Dipole) || hit;
            hit = machine.RunOmolChargeStage(std::nullopt, molgr::reconstruct::EliminatePositiveCharges) || hit;
            hit = machine.RunOmolChargeStage(std::nullopt, molgr::reconstruct::EliminateNegativeCharges) || hit;
            hit = machine.RunOmolStage(std::nullopt, molgr::reconstruct::CleanNeighborRadicals) || hit;
            hit = machine.RunOmolStage(std::nullopt, molgr::reconstruct::CleanResonances) || hit;
            return std::make_tuple(OpenBabel::OBMol(*machine.omol), machine.given_charge, hit);
        }

        std::pair<OpenBabel::OBMol, int> ProcessResonance(const OpenBabel::OBMol &mol, int charge)
        {
            auto processed = ProcessResonanceDetailed(mol, charge);
            return {std::get<0>(processed), std::get<1>(processed)};
        }
    }
}
