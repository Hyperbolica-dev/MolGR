/**
 * @file metal_handler.cpp
 * @brief Implementation of metal handling logic.
 * @author TMJ
 * @date 2025-12-28
 */

#include "molgr/pipeline/reconstruct_with_metals.h"
#include "molgr/pipeline/reconstruct_without_metals.h"
#include "molgr/state.h"
#include "molgr/utils/conversions.h"
#include "molgr/utils/consts.h"
#include "molgr/utils/logger.h"
#include "molgr/utils/scoring.h"
#include "molgr/utils/utils.h"

#include <openbabel/obconversion.h>
#include <openbabel/elements.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <future>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <sstream>
#include <thread>
#include <tuple>
#include <vector>

namespace molgr
{
    namespace metal
    {
        using namespace OpenBabel;

        OpenBabel::OBConversion &ThreadLocalXyzOutConversion()
        {
            thread_local OpenBabel::OBConversion conv;
            thread_local bool initialized = false;
            if (!initialized)
            {
                conv.SetOutFormat("xyz");
                initialized = true;
            }
            return conv;
        }

        OpenBabel::OBConversion &ThreadLocalXyzInConversion()
        {
            thread_local OpenBabel::OBConversion conv;
            thread_local bool initialized = false;
            if (!initialized)
            {
                conv.SetInFormat("xyz");
                initialized = true;
            }
            return conv;
        }

        void ReinsertMetalStates(OBMol &mol, const std::vector<MetalAtomPosition> &metals)
        {
            mol.BeginModify();

            int num_organic = mol.NumAtoms();
            int num_metals = static_cast<int>(metals.size());
            int total_atoms = num_organic + num_metals;

            for (const auto &m : metals)
            {
                OBAtom *atom = mol.NewAtom();
                atom->SetAtomicNum(m.element_idx);
                atom->SetFormalCharge(m.valence);
                atom->SetSpinMultiplicity(m.radical_num);
                atom->SetVector(m.position_x, m.position_y, m.position_z);
            }

            std::vector<int> new_order(total_atoms, 0);
            bool error_flag = false;

            for (int i = 0; i < num_metals; ++i)
            {
                int current_idx = num_organic + 1 + i;
                int target_slot = metals[i].idx - 1;

                if (target_slot >= 0 && target_slot < total_atoms)
                {
                    if (new_order[target_slot] != 0)
                    {
                        LOG_ERROR("[ReinsertMetalStates] Index collision at slot " << target_slot);
                        error_flag = true;
                    }
                    new_order[target_slot] = current_idx;
                }
                else
                {
                    LOG_ERROR("[ReinsertMetalStates] Original index out of bounds: " << metals[i].idx);
                    error_flag = true;
                }
            }

            int current_organic_idx = 1;
            for (int i = 0; i < total_atoms; ++i)
            {
                if (new_order[i] == 0)
                {
                    if (current_organic_idx <= num_organic)
                    {
                        new_order[i] = current_organic_idx;
                        current_organic_idx++;
                    }
                    else
                    {
                        LOG_ERROR("[ReinsertMetalStates] Not enough organic atoms to fill slots.");
                        error_flag = true;
                    }
                }
            }

            for (int idx : new_order)
            {
                if (idx == 0)
                {
                    LOG_ERROR("[ReinsertMetalStates] Invalid 0 index in renumber map. Aborting renumber.");
                    error_flag = true;
                    break;
                }
            }

            if (!error_flag)
            {
                mol.RenumberAtoms(new_order);
            }

            mol.EndModify();
        }

    } // namespace metal
} // namespace molgr

namespace molgr
{
    namespace pipeline
    {
        namespace reconstruct_with_metals
        {
            namespace
            {
                constexpr int kDefaultMaxMixedValenceSpread = 3;
                constexpr int kDefaultMaxAssignmentsPerTarget = 64;
                constexpr std::size_t kCandidateScoreParallelThreshold = 32;

                using ValenceBoundsEntry = std::tuple<std::string, int, int>;
                using ValenceBoundsKey = std::vector<ValenceBoundsEntry>;

                struct PartialMetalAssignment
                {
                    std::vector<molgr::metal::MetalAtomPosition> metal_states;
                    int total_metal_charge = 0;
                    int total_metal_radicals = 0;
                    double metal_assignment_rank = 0.0;
                    ValenceBoundsKey valence_bounds;
                    int order = 0;
                };

                struct FrontierKey
                {
                    int charge = 0;
                    int radicals = 0;
                    ValenceBoundsKey valence_bounds;

                    bool operator<(const FrontierKey &other) const
                    {
                        return std::tie(charge, radicals, valence_bounds) <
                               std::tie(other.charge, other.radicals, other.valence_bounds);
                    }
                };

                using PartialAssignmentFrontier = std::map<FrontierKey, std::vector<PartialMetalAssignment>>;
                using ChargeGroupedAssignments =
                    std::map<int, std::vector<std::pair<ValenceBoundsKey, std::vector<PartialMetalAssignment>>>>;
                using RadicalBucketIndex = std::map<int, ChargeGroupedAssignments>;
                using TargetBucket = std::pair<int, int>;
                using TargetCandidateBuckets =
                    std::map<TargetBucket, std::vector<molgr::state::MetalCandidateState>>;

                struct TargetBucketTask
                {
                    TargetBucket target;
                    std::vector<molgr::state::MetalCandidateState> candidates;
                };

                struct PreparedTargetBucket
                {
                    TargetBucket target;
                    std::shared_ptr<molgr::state::ReconstructionState> no_metal_state;
                };

                struct CandidateScoreJob
                {
                    std::size_t bucket_index = 0;
                    std::size_t candidate_index = 0;
                };

                std::size_t HardwareParallelism()
                {
                    const unsigned int concurrency = std::thread::hardware_concurrency();
                    return std::max<std::size_t>(
                        1,
                        concurrency == 0 ? 1 : static_cast<std::size_t>(concurrency));
                }

                template <typename Func>
                void ParallelForIndices(std::size_t count, std::size_t worker_count, Func &&func)
                {
                    if (count == 0)
                    {
                        return;
                    }

                    worker_count = std::max<std::size_t>(1, std::min(worker_count, count));
                    if (worker_count == 1)
                    {
                        for (std::size_t idx = 0; idx < count; ++idx)
                        {
                            func(idx);
                        }
                        return;
                    }

                    std::atomic<std::size_t> next_index{0};
                    std::exception_ptr first_exception;
                    std::mutex exception_mutex;

                    const auto worker = [&]()
                    {
                        try
                        {
                            while (true)
                            {
                                const std::size_t idx = next_index.fetch_add(1);
                                if (idx >= count)
                                {
                                    break;
                                }
                                func(idx);
                            }
                        }
                        catch (...)
                        {
                            std::lock_guard<std::mutex> lock(exception_mutex);
                            if (!first_exception)
                            {
                                first_exception = std::current_exception();
                                next_index.store(count);
                            }
                        }
                    };

                    std::vector<std::thread> workers;
                    workers.reserve(worker_count - 1);
                    for (std::size_t worker_idx = 1; worker_idx < worker_count; ++worker_idx)
                    {
                        workers.emplace_back(worker);
                    }
                    worker();
                    for (auto &thread : workers)
                    {
                        thread.join();
                    }
                    if (first_exception)
                    {
                        std::rethrow_exception(first_exception);
                    }
                }

                int CandidateCombinationIndex(const molgr::state::MetalCandidateState &candidate)
                {
                    const auto metadata_it = candidate.metadata.find("combination_index");
                    if (metadata_it == candidate.metadata.end())
                    {
                        return std::numeric_limits<int>::max();
                    }
                    if (const auto *value = std::get_if<int>(&metadata_it->second))
                    {
                        return *value;
                    }
                    return std::numeric_limits<int>::max();
                }

                bool CandidateScoreLess(
                    const molgr::state::MetalCandidateState &lhs,
                    const molgr::state::MetalCandidateState &rhs)
                {
                    const double lhs_score = lhs.score.has_value() ? *lhs.score : lhs.CombinedScore();
                    const double rhs_score = rhs.score.has_value() ? *rhs.score : rhs.CombinedScore();
                    return std::make_tuple(lhs_score, CandidateCombinationIndex(lhs)) <
                           std::make_tuple(rhs_score, CandidateCombinationIndex(rhs));
                }

                molgr::state::MetalPreparationState PrepareMetalState(
                    const std::string &xyz_block,
                    int total_charge,
                    int total_radical_electrons)
                {
                    OpenBabel::OBMol mol;
                    OpenBabel::OBConversion &conv = molgr::metal::ThreadLocalXyzInConversion();
                    if (!conv.ReadString(&mol, xyz_block))
                    {
                        return {};
                    }

                    std::vector<OpenBabel::OBAtom *> removable_metal_atoms;
                    std::vector<std::vector<molgr::metal::MetalAtomPosition>> available_states;
                    FOR_ATOMS_OF_MOL(atom_iter, mol)
                    {
                        OpenBabel::OBAtom *atom = &(*atom_iter);
                        if (!atom->IsMetal())
                        {
                            continue;
                        }
                        removable_metal_atoms.push_back(atom);
                        available_states.push_back(build_metal_states(*atom));
                    }

                    for (OpenBabel::OBAtom *atom : removable_metal_atoms)
                    {
                        mol.DeleteAtom(atom);
                    }

                    molgr::state::MetalPreparationState state;
                    state.no_metal_xyz_block = molgr::metal::ThreadLocalXyzOutConversion().WriteString(&mol);
                    state.available_valence_radical_states = std::move(available_states);
                    state.total_charge = total_charge;
                    state.total_radical_electrons = total_radical_electrons;
                    state.phase_history = {
                        "read_xyz",
                        "build_metal_state_options",
                        "remove_metal_atoms",
                        "serialize_no_metal_xyz",
                    };
                    state.metadata["metal_atom_count"] = static_cast<int>(removable_metal_atoms.size());
                    return state;
                }

                double MetalStateAssignmentPenalty(const molgr::metal::MetalAtomPosition &metal_state)
                {
                    double penalty = 0.0;
                    if (metal_state.valence <= 0)
                    {
                        penalty += 10.0 * std::max(std::abs(metal_state.valence), 1);
                    }

                    const auto prior_it = kMetalValencePrior.find(metal_state.symbol);
                    const auto minor_it = kMetalValenceMinor.find(metal_state.symbol);
                    const auto contains = [](const std::vector<int> &values, int target)
                    {
                        return std::find(values.begin(), values.end(), target) != values.end();
                    };

                    const bool in_prior =
                        prior_it != kMetalValencePrior.end() && contains(prior_it->second, metal_state.valence);
                    const bool in_minor =
                        minor_it != kMetalValenceMinor.end() && contains(minor_it->second, metal_state.valence);

                    if (!in_prior)
                    {
                        penalty += in_minor ? 10.0 : 20.0;
                    }
                    return penalty;
                }

                std::optional<ValenceBoundsKey> UpdateValenceBounds(
                    const ValenceBoundsKey &bounds,
                    const molgr::metal::MetalAtomPosition &metal_state,
                    const std::optional<int> &max_mixed_valence_spread)
                {
                    if (!max_mixed_valence_spread.has_value() || *max_mixed_valence_spread < 0)
                    {
                        return ValenceBoundsKey{};
                    }

                    ValenceBoundsKey updated_bounds = bounds;
                    for (auto &entry : updated_bounds)
                    {
                        if (std::get<0>(entry) != metal_state.symbol)
                        {
                            continue;
                        }

                        const int next_lower = std::min(std::get<1>(entry), metal_state.valence);
                        const int next_upper = std::max(std::get<2>(entry), metal_state.valence);
                        if (next_upper - next_lower > *max_mixed_valence_spread)
                        {
                            return std::nullopt;
                        }
                        entry = std::make_tuple(metal_state.symbol, next_lower, next_upper);
                        return updated_bounds;
                    }

                    updated_bounds.emplace_back(
                        metal_state.symbol,
                        metal_state.valence,
                        metal_state.valence);
                    std::sort(updated_bounds.begin(), updated_bounds.end());
                    return updated_bounds;
                }

                std::vector<PartialMetalAssignment> TrimPartialAssignments(
                    std::vector<PartialMetalAssignment> entries,
                    int max_assignments_per_target)
                {
                    const int limit = std::max(1, max_assignments_per_target);
                    std::sort(
                        entries.begin(),
                        entries.end(),
                        [](const PartialMetalAssignment &lhs, const PartialMetalAssignment &rhs)
                        {
                            return std::tie(lhs.metal_assignment_rank, lhs.order) <
                                   std::tie(rhs.metal_assignment_rank, rhs.order);
                        });
                    if (static_cast<int>(entries.size()) > limit)
                    {
                        entries.resize(limit);
                    }
                    return entries;
                }

                PartialAssignmentFrontier EnumeratePartialAssignmentFrontier(
                    const std::vector<std::vector<molgr::metal::MetalAtomPosition>> &available_valence_radical_states,
                    const std::optional<int> &max_mixed_valence_spread,
                    const std::optional<int> &max_total_metal_radicals,
                    int max_assignments_per_state)
                {
                    PartialAssignmentFrontier partial_assignments;
                    partial_assignments[FrontierKey{}] = {
                        PartialMetalAssignment{{}, 0, 0, 0.0, {}, 0}};
                    int next_order = 1;

                    for (const auto &metal_state_options : available_valence_radical_states)
                    {
                        PartialAssignmentFrontier next_partial_assignments;
                        for (const auto &frontier_entry : partial_assignments)
                        {
                            for (const auto &entry : frontier_entry.second)
                            {
                                for (const auto &metal_state : metal_state_options)
                                {
                                    const int next_total_metal_radicals =
                                        entry.total_metal_radicals + metal_state.radical_num;
                                    if (max_total_metal_radicals.has_value() &&
                                        next_total_metal_radicals > *max_total_metal_radicals)
                                    {
                                        continue;
                                    }

                                    const auto next_valence_bounds = UpdateValenceBounds(
                                        entry.valence_bounds,
                                        metal_state,
                                        max_mixed_valence_spread);
                                    if (!next_valence_bounds.has_value())
                                    {
                                        continue;
                                    }

                                    PartialMetalAssignment next_entry;
                                    next_entry.metal_states = entry.metal_states;
                                    next_entry.metal_states.push_back(metal_state);
                                    next_entry.total_metal_charge = entry.total_metal_charge + metal_state.valence;
                                    next_entry.total_metal_radicals = next_total_metal_radicals;
                                    next_entry.metal_assignment_rank =
                                        entry.metal_assignment_rank + MetalStateAssignmentPenalty(metal_state);
                                    next_entry.valence_bounds = *next_valence_bounds;
                                    next_entry.order = next_order++;

                                    FrontierKey next_key{
                                        next_entry.total_metal_charge,
                                        next_entry.total_metal_radicals,
                                        next_entry.valence_bounds};
                                    next_partial_assignments[next_key].push_back(std::move(next_entry));
                                }
                            }
                        }

                        PartialAssignmentFrontier trimmed_frontier;
                        for (auto &entry : next_partial_assignments)
                        {
                            trimmed_frontier.emplace(
                                entry.first,
                                TrimPartialAssignments(
                                    std::move(entry.second),
                                    max_assignments_per_state));
                        }
                        partial_assignments = std::move(trimmed_frontier);
                        if (partial_assignments.empty())
                        {
                            break;
                        }
                    }

                    return partial_assignments;
                }

                std::optional<ValenceBoundsKey> MergeValenceBounds(
                    const ValenceBoundsKey &left_bounds,
                    const ValenceBoundsKey &right_bounds,
                    const std::optional<int> &max_mixed_valence_spread)
                {
                    if (!max_mixed_valence_spread.has_value() || *max_mixed_valence_spread < 0)
                    {
                        return ValenceBoundsKey{};
                    }

                    std::map<std::string, std::pair<int, int>> merged_bounds;
                    for (const auto &entry : left_bounds)
                    {
                        merged_bounds[std::get<0>(entry)] = {std::get<1>(entry), std::get<2>(entry)};
                    }

                    for (const auto &entry : right_bounds)
                    {
                        const std::string &symbol = std::get<0>(entry);
                        const int lower = std::get<1>(entry);
                        const int upper = std::get<2>(entry);
                        const auto it = merged_bounds.find(symbol);
                        if (it == merged_bounds.end())
                        {
                            merged_bounds.emplace(symbol, std::make_pair(lower, upper));
                            continue;
                        }

                        const int next_lower = std::min(it->second.first, lower);
                        const int next_upper = std::max(it->second.second, upper);
                        if (next_upper - next_lower > *max_mixed_valence_spread)
                        {
                            return std::nullopt;
                        }
                        it->second = {next_lower, next_upper};
                    }

                    ValenceBoundsKey result;
                    result.reserve(merged_bounds.size());
                    for (const auto &entry : merged_bounds)
                    {
                        result.emplace_back(entry.first, entry.second.first, entry.second.second);
                    }
                    return result;
                }

                RadicalBucketIndex BucketPartialAssignmentsByChargeRadicals(
                    const PartialAssignmentFrontier &frontier)
                {
                    RadicalBucketIndex bucket_index;
                    for (const auto &frontier_entry : frontier)
                    {
                        bucket_index[frontier_entry.first.radicals][frontier_entry.first.charge].push_back(
                            std::make_pair(frontier_entry.first.valence_bounds, frontier_entry.second));
                    }
                    return bucket_index;
                }

                std::map<TargetBucket, std::vector<PartialMetalAssignment>> CombinePartialAssignmentFrontiers(
                    const PartialAssignmentFrontier &left_frontier,
                    const PartialAssignmentFrontier &right_frontier,
                    int total_charge,
                    int total_radical_electrons,
                    const std::optional<int> &max_mixed_valence_spread,
                    const std::optional<int> &max_total_metal_radicals,
                    int max_assignments_per_target)
                {
                    std::map<TargetBucket, std::vector<PartialMetalAssignment>> grouped_entries;
                    int next_order = 0;
                    const int trim_trigger = std::max(1, max_assignments_per_target) * 4;
                    const RadicalBucketIndex left_bucket_index =
                        BucketPartialAssignmentsByChargeRadicals(left_frontier);
                    const RadicalBucketIndex right_bucket_index =
                        BucketPartialAssignmentsByChargeRadicals(right_frontier);

                    int max_combined_metal_radicals = total_radical_electrons;
                    if (max_total_metal_radicals.has_value())
                    {
                        max_combined_metal_radicals =
                            std::min(max_combined_metal_radicals, *max_total_metal_radicals);
                    }

                    for (const auto &left_radical_bucket : left_bucket_index)
                    {
                        const int left_radicals = left_radical_bucket.first;
                        const int max_right_radicals = max_combined_metal_radicals - left_radicals;
                        if (max_right_radicals < 0)
                        {
                            continue;
                        }

                        for (auto right_it = right_bucket_index.begin();
                             right_it != right_bucket_index.end() && right_it->first <= max_right_radicals;
                             ++right_it)
                        {
                            const int right_radicals = right_it->first;
                            const int total_metal_radicals = left_radicals + right_radicals;
                            const int target_radicals = total_radical_electrons - total_metal_radicals;

                            for (const auto &left_charge_group : left_radical_bucket.second)
                            {
                                const int left_charge = left_charge_group.first;
                                for (const auto &right_charge_group : right_it->second)
                                {
                                    const int right_charge = right_charge_group.first;
                                    const TargetBucket target{
                                        total_charge - (left_charge + right_charge),
                                        target_radicals};
                                    auto &bucket = grouped_entries[target];
                                    const std::size_t initial_bucket_size = bucket.size();

                                    for (const auto &left_bounds_group : left_charge_group.second)
                                    {
                                        for (const auto &right_bounds_group : right_charge_group.second)
                                        {
                                            const auto merged_bounds = MergeValenceBounds(
                                                left_bounds_group.first,
                                                right_bounds_group.first,
                                                max_mixed_valence_spread);
                                            if (!merged_bounds.has_value())
                                            {
                                                continue;
                                            }

                                            for (const auto &left_entry : left_bounds_group.second)
                                            {
                                                for (const auto &right_entry : right_bounds_group.second)
                                                {
                                                    PartialMetalAssignment combined_entry;
                                                    combined_entry.metal_states = left_entry.metal_states;
                                                    combined_entry.metal_states.insert(
                                                        combined_entry.metal_states.end(),
                                                        right_entry.metal_states.begin(),
                                                        right_entry.metal_states.end());
                                                    combined_entry.total_metal_charge =
                                                        left_entry.total_metal_charge +
                                                        right_entry.total_metal_charge;
                                                    combined_entry.total_metal_radicals =
                                                        left_entry.total_metal_radicals +
                                                        right_entry.total_metal_radicals;
                                                    combined_entry.metal_assignment_rank =
                                                        left_entry.metal_assignment_rank +
                                                        right_entry.metal_assignment_rank;
                                                    combined_entry.valence_bounds = *merged_bounds;
                                                    combined_entry.order = next_order++;
                                                    bucket.push_back(std::move(combined_entry));
                                                }
                                            }
                                        }
                                    }

                                    if (bucket.size() == initial_bucket_size)
                                    {
                                        if (initial_bucket_size == 0)
                                        {
                                            grouped_entries.erase(target);
                                        }
                                        continue;
                                    }

                                    if (static_cast<int>(bucket.size()) > trim_trigger)
                                    {
                                        bucket = TrimPartialAssignments(
                                            std::move(bucket),
                                            max_assignments_per_target);
                                    }
                                }
                            }
                        }
                    }

                    for (auto &entry : grouped_entries)
                    {
                        entry.second = TrimPartialAssignments(
                            std::move(entry.second),
                            max_assignments_per_target);
                    }
                    return grouped_entries;
                }

                TargetCandidateBuckets GroupCandidatesByTargetDp(
                    const std::vector<std::string> &base_phase_history,
                    const std::vector<std::vector<molgr::metal::MetalAtomPosition>> &available_valence_radical_states,
                    int total_charge,
                    int total_radical_electrons)
                {
                    const std::size_t split_index = available_valence_radical_states.size() / 2;
                    auto left_options = std::vector<std::vector<molgr::metal::MetalAtomPosition>>(
                        available_valence_radical_states.begin(),
                        available_valence_radical_states.begin() + static_cast<std::ptrdiff_t>(split_index));
                    auto right_options = std::vector<std::vector<molgr::metal::MetalAtomPosition>>(
                        available_valence_radical_states.begin() + static_cast<std::ptrdiff_t>(split_index),
                        available_valence_radical_states.end());

                    PartialAssignmentFrontier left_frontier;
                    PartialAssignmentFrontier right_frontier;
                    const bool parallelize_frontiers =
                        HardwareParallelism() > 1 &&
                        !left_options.empty() &&
                        !right_options.empty();
                    if (parallelize_frontiers)
                    {
                        auto left_future = std::async(
                            std::launch::async,
                            [left_options = std::move(left_options), total_radical_electrons]()
                            {
                                return EnumeratePartialAssignmentFrontier(
                                    left_options,
                                    kDefaultMaxMixedValenceSpread,
                                    total_radical_electrons,
                                    kDefaultMaxAssignmentsPerTarget);
                            });
                        right_frontier = EnumeratePartialAssignmentFrontier(
                            right_options,
                            kDefaultMaxMixedValenceSpread,
                            total_radical_electrons,
                            kDefaultMaxAssignmentsPerTarget);
                        left_frontier = left_future.get();
                    }
                    else
                    {
                        left_frontier = EnumeratePartialAssignmentFrontier(
                            left_options,
                            kDefaultMaxMixedValenceSpread,
                            total_radical_electrons,
                            kDefaultMaxAssignmentsPerTarget);
                        right_frontier = EnumeratePartialAssignmentFrontier(
                            right_options,
                            kDefaultMaxMixedValenceSpread,
                            total_radical_electrons,
                            kDefaultMaxAssignmentsPerTarget);
                    }
                    if (left_frontier.empty() || right_frontier.empty())
                    {
                        return {};
                    }

                    const auto grouped_entries = CombinePartialAssignmentFrontiers(
                        left_frontier,
                        right_frontier,
                        total_charge,
                        total_radical_electrons,
                        kDefaultMaxMixedValenceSpread,
                        total_radical_electrons,
                        kDefaultMaxAssignmentsPerTarget);

                    TargetCandidateBuckets grouped_candidates;
                    int combination_index = 0;
                    for (const auto &target_entry : grouped_entries)
                    {
                        auto &bucket = grouped_candidates[target_entry.first];
                        for (const auto &entry : target_entry.second)
                        {
                            molgr::state::MetalCandidateStateMachine machine(
                                entry.metal_states,
                                target_entry.first.first,
                                target_entry.first.second,
                                base_phase_history,
                                {{"combination_index", combination_index}});
                            machine.Annotate("enumerate_metal_combination");
                            machine.Annotate("reconstruct_no_metal_candidate");
                            machine.metadata["metal_assignment_rank"] = entry.metal_assignment_rank;
                            bucket.push_back(machine.Freeze());
                            ++combination_index;
                        }
                    }
                    return grouped_candidates;
                }

                std::vector<TargetBucketTask> BuildTargetBucketTasks(TargetCandidateBuckets &&grouped_candidates)
                {
                    std::vector<TargetBucketTask> tasks;
                    tasks.reserve(grouped_candidates.size());
                    for (auto &entry : grouped_candidates)
                    {
                        tasks.push_back(TargetBucketTask{
                            entry.first,
                            std::move(entry.second),
                        });
                    }
                    return tasks;
                }

                molgr::state::MetalCandidateState ScoreCandidateWithNoMetalState(
                    const molgr::state::MetalCandidateState &candidate,
                    const std::shared_ptr<molgr::state::ReconstructionState> &no_metal_state)
                {
                    auto machine = molgr::state::MetalCandidateStateMachine::FromCandidateState(candidate);
                    machine.SetNoMetalState("reconstruct_no_metal", no_metal_state);
                    machine.Annotate("score_candidate");
                    auto scored_candidate = machine.Freeze();
                    const double score = scored_candidate.CombinedScore();
                    scored_candidate.score = score;
                    scored_candidate.metadata["score"] = score;
                    return scored_candidate;
                }
            }

            std::set<int> get_possible_metal_radicals(const std::string &metal_symbol, int valence)
            {
                return molgr::GetPossibleMetalRadicals(metal_symbol, valence);
            }

            std::vector<molgr::metal::MetalAtomPosition> build_metal_states(const OpenBabel::OBAtom &obatom)
            {
                const int atomic_num = static_cast<int>(obatom.GetAtomicNum());
                const std::string symbol = OpenBabel::OBElements::GetSymbol(atomic_num);

                const auto default_state = [&]()
                {
                    return molgr::metal::MetalAtomPosition{
                        static_cast<int>(obatom.GetIdx()),
                        symbol,
                        atomic_num,
                        0,
                        0,
                        obatom.GetX(),
                        obatom.GetY(),
                        obatom.GetZ()};
                };

                std::vector<int> valences;
                std::set<int> seen_valences;
                const auto add_valences = [&](const std::vector<int> &source)
                {
                    for (const int valence : source)
                    {
                        if (seen_valences.insert(valence).second)
                        {
                            valences.push_back(valence);
                        }
                    }
                };

                if (kMetalValencePrior.count(symbol))
                {
                    add_valences(kMetalValencePrior.at(symbol));
                }
                if (kMetalValenceMinor.count(symbol))
                {
                    add_valences(kMetalValenceMinor.at(symbol));
                }
                if (valences.empty())
                {
                    valences.push_back(0);
                }

                if (!kMetalFDSP.count(symbol))
                {
                    return {default_state()};
                }

                std::vector<molgr::metal::MetalAtomPosition> states;
                for (const int valence : valences)
                {
                    const auto radicals = get_possible_metal_radicals(symbol, valence);
                    for (const int radical_num : radicals)
                    {
                        states.push_back(
                            molgr::metal::MetalAtomPosition{
                                static_cast<int>(obatom.GetIdx()),
                                symbol,
                                atomic_num,
                                valence,
                                radical_num,
                                obatom.GetX(),
                                obatom.GetY(),
                                obatom.GetZ()});
                    }
                }

                if (states.empty())
                {
                    return {default_state()};
                }
                return states;
            }

            void combine_metal_with_omol(
                OpenBabel::OBMol &mol,
                const std::vector<molgr::metal::MetalAtomPosition> &metals)
            {
                molgr::metal::ReinsertMetalStates(mol, metals);
            }

            std::unique_ptr<molgr::utils::MoleculeData> Xyz2OmolMolData(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons)
            {
                molgr::pipeline::perf::RunTimingScope timing_scope;
                auto &timing_reducer = timing_scope.Reducer();
                if (total_radical_electrons < 0)
                {
                    return nullptr;
                }

                const auto metal_enum_started = std::chrono::steady_clock::now();
                const auto base_state = PrepareMetalState(
                    xyz_block,
                    total_charge,
                    total_radical_electrons);
                if (base_state.phase_history.empty())
                {
                    return nullptr;
                }

                auto grouped_candidates = GroupCandidatesByTargetDp(
                    base_state.phase_history,
                    base_state.available_valence_radical_states,
                    total_charge,
                    total_radical_electrons);
                const auto metal_enum_now = std::chrono::steady_clock::now();
                const double metal_enum_ms =
                    std::chrono::duration<double, std::milli>(metal_enum_now - metal_enum_started).count();
                timing_reducer.AddMetalEnumerationCombinationMs(metal_enum_ms);

                auto target_bucket_tasks = BuildTargetBucketTasks(std::move(grouped_candidates));
                if (target_bucket_tasks.empty())
                {
                    return nullptr;
                }

                std::vector<std::optional<PreparedTargetBucket>> prepared_buckets(target_bucket_tasks.size());
                const std::size_t bucket_parallelism =
                    target_bucket_tasks.size() > 1
                        ? std::min(HardwareParallelism(), target_bucket_tasks.size())
                        : 1;
                ParallelForIndices(
                    target_bucket_tasks.size(),
                    bucket_parallelism,
                    [&](std::size_t bucket_index)
                    {
                        const auto &task = target_bucket_tasks[bucket_index];
                        auto no_metal_state =
                            molgr::pipeline::reconstruct_without_metals::XyzToOmolNoMetalState(
                                base_state.no_metal_xyz_block,
                                task.target.first,
                                task.target.second,
                                &timing_reducer);
                        if (!no_metal_state.has_value())
                        {
                            return;
                        }

                        prepared_buckets[bucket_index] = PreparedTargetBucket{
                            task.target,
                            std::make_shared<molgr::state::ReconstructionState>(std::move(*no_metal_state)),
                        };
                    });

                std::vector<CandidateScoreJob> score_jobs;
                std::vector<std::vector<molgr::state::MetalCandidateState>> scored_buckets(
                    target_bucket_tasks.size());
                for (std::size_t bucket_index = 0; bucket_index < target_bucket_tasks.size(); ++bucket_index)
                {
                    if (!prepared_buckets[bucket_index].has_value())
                    {
                        continue;
                    }
                    auto &scored_bucket = scored_buckets[bucket_index];
                    scored_bucket.resize(target_bucket_tasks[bucket_index].candidates.size());
                    for (std::size_t candidate_index = 0;
                         candidate_index < target_bucket_tasks[bucket_index].candidates.size();
                         ++candidate_index)
                    {
                        score_jobs.push_back(CandidateScoreJob{
                            bucket_index,
                            candidate_index,
                        });
                    }
                }

                if (score_jobs.empty())
                {
                    return nullptr;
                }

                const std::size_t score_parallelism =
                    score_jobs.size() >= kCandidateScoreParallelThreshold
                        ? std::min(HardwareParallelism(), score_jobs.size())
                        : 1;
                ParallelForIndices(
                    score_jobs.size(),
                    score_parallelism,
                    [&](std::size_t job_index)
                    {
                        const auto &job = score_jobs[job_index];
                        const auto &prepared_bucket = *prepared_buckets[job.bucket_index];
                        scored_buckets[job.bucket_index][job.candidate_index] =
                            ScoreCandidateWithNoMetalState(
                                target_bucket_tasks[job.bucket_index].candidates[job.candidate_index],
                                prepared_bucket.no_metal_state);
                    });

                std::vector<molgr::state::MetalCandidateState> possible_candidates;
                possible_candidates.reserve(score_jobs.size());
                for (std::size_t bucket_index = 0; bucket_index < scored_buckets.size(); ++bucket_index)
                {
                    auto &scored_bucket = scored_buckets[bucket_index];
                    possible_candidates.insert(
                        possible_candidates.end(),
                        std::make_move_iterator(scored_bucket.begin()),
                        std::make_move_iterator(scored_bucket.end()));
                }

                auto best_it = std::min_element(
                    possible_candidates.begin(),
                    possible_candidates.end(),
                    CandidateScoreLess);

                if (best_it == possible_candidates.end())
                {
                    return nullptr;
                }

                auto best_candidate = *best_it;
                if (!best_candidate.combined_omol)
                {
                    best_candidate.MaterializeCombinedOmol(
                        [](const OpenBabel::OBMol &no_metal_omol,
                           const std::vector<molgr::metal::MetalAtomPosition> &metal_states)
                        {
                            auto combined_omol = std::make_shared<OpenBabel::OBMol>(no_metal_omol);
                            combine_metal_with_omol(*combined_omol, metal_states);
                            return combined_omol;
                        });
                    auto winner_machine =
                        molgr::state::MetalCandidateStateMachine::FromCandidateState(best_candidate);
                    winner_machine.Annotate("combine_metal_with_omol");
                    best_candidate = winner_machine.Freeze();
                }

                auto winner_machine =
                    molgr::state::MetalCandidateStateMachine::FromCandidateState(best_candidate);
                winner_machine.Annotate("select_best_candidate");
                best_candidate = winner_machine.Freeze();

                if (!best_candidate.combined_omol)
                {
                    return nullptr;
                }
                return std::make_unique<molgr::utils::MoleculeData>(
                    molgr::utils::MoleculeDataFromOBMol(*best_candidate.combined_omol));
            }
        }
    }
}
