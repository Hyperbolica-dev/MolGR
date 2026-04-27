#pragma once

#include "molgr/context.h"
#include "molgr/state.h"
#include "molgr/types.h"

#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace molgr
{
    namespace metal
    {
        namespace search
        {
            using ValenceBoundsEntry = std::tuple<std::string, int, int>;
            using ValenceBoundsKey = std::vector<ValenceBoundsEntry>;
            using MetalStateChoice = std::vector<molgr::MetalAtomPosition>;
            using MetalStateChoiceGroup = std::vector<MetalStateChoice>;
            using MetalStateSearchLayer = std::vector<MetalStateChoiceGroup>;

            struct PartialMetalAssignment
            {
                std::vector<molgr::MetalAtomPosition> metal_states;
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
            using TargetBucket = molgr::context::TargetBucket;
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

            double MetalStateAssignmentPenalty(const molgr::metal::MetalAtomPosition &metal_state);

            std::optional<ValenceBoundsKey> UpdateValenceBounds(
                const ValenceBoundsKey &bounds,
                const molgr::metal::MetalAtomPosition &metal_state,
                const std::optional<int> &max_mixed_valence_spread);

            std::vector<PartialMetalAssignment> TrimPartialAssignments(
                std::vector<PartialMetalAssignment> entries,
                int max_assignments_per_target);

            std::vector<MetalStateChoiceGroup> BuildMetalStateSearchGroups(
                const std::vector<std::vector<molgr::metal::MetalAtomPosition>> &available_valence_radical_states,
                const molgr::config::MolGRConfig &config);

            std::vector<MetalStateSearchLayer> BuildLayeredMetalStateSearchGroups(
                const std::vector<MetalStateChoiceGroup> &available_state_search_groups,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config);

            PartialAssignmentFrontier EnumeratePartialAssignmentFrontier(
                const std::vector<std::vector<molgr::metal::MetalAtomPosition>> &available_valence_radical_states,
                const std::optional<int> &max_mixed_valence_spread,
                const std::optional<int> &max_total_metal_radicals,
                int max_assignments_per_state);

            PartialAssignmentFrontier EnumeratePartialAssignmentFrontier(
                const std::vector<MetalStateChoiceGroup> &available_state_search_groups,
                const std::optional<int> &max_mixed_valence_spread,
                const std::optional<int> &max_total_metal_radicals,
                int max_assignments_per_state);

            std::optional<ValenceBoundsKey> MergeValenceBounds(
                const ValenceBoundsKey &left_bounds,
                const ValenceBoundsKey &right_bounds,
                const std::optional<int> &max_mixed_valence_spread);

            RadicalBucketIndex BucketPartialAssignmentsByChargeRadicals(
                const PartialAssignmentFrontier &frontier);

            std::map<TargetBucket, std::vector<PartialMetalAssignment>> CombinePartialAssignmentFrontiers(
                const PartialAssignmentFrontier &left_frontier,
                const PartialAssignmentFrontier &right_frontier,
                int total_charge,
                int total_radical_electrons,
                const std::optional<int> &max_mixed_valence_spread,
                const std::optional<int> &max_total_metal_radicals,
                int max_assignments_per_target);

            TargetCandidateBuckets GroupCandidatesByTargetDp(
                const std::vector<std::string> &base_phase_history,
                const std::vector<std::vector<molgr::metal::MetalAtomPosition>> &available_valence_radical_states,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config);

            TargetCandidateBuckets GroupCandidatesByTargetDp(
                const std::vector<std::string> &base_phase_history,
                const std::vector<MetalStateChoiceGroup> &available_state_search_groups,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config);

            std::vector<TargetBucketTask> BuildTargetBucketTasks(TargetCandidateBuckets &&grouped_candidates);
        }
    }
}
