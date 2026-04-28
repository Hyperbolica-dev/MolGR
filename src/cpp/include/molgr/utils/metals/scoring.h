#pragma once

#include "molgr/config.h"
#include "molgr/state.h"

#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace molgr
{
    namespace metal
    {
        namespace scoring
        {
            int CandidateCombinationIndex(const molgr::state::MetalCandidateState &candidate);

            int MetadataInt(
                const molgr::state::MetalCandidateState &candidate,
                const std::string &key,
                int fallback = 0);

            double MetadataDouble(
                const molgr::state::MetalCandidateState &candidate,
                const std::string &key,
                double fallback = 0.0);

            int OrganicScoreBucketIndex(
                double score_value,
                double best_score,
                const molgr::config::MolGRConfig &config);

            std::optional<molgr::state::MetalCandidateState> SelectBestCandidate(
                const std::vector<molgr::state::MetalCandidateState> &candidates,
                const molgr::config::MolGRConfig &config);

            molgr::state::MetalCandidateState PrepareCandidateWithNoMetalState(
                const molgr::state::MetalCandidateState &candidate,
                const std::shared_ptr<molgr::state::ReconstructionState> &no_metal_state,
                const molgr::config::MolGRConfig &config);

            molgr::state::MetalCandidateState ScoreCandidateWithNoMetalState(
                const molgr::state::MetalCandidateState &candidate,
                const std::shared_ptr<molgr::state::ReconstructionState> &no_metal_state,
                const molgr::config::MolGRConfig &config);
        }
    }
}
