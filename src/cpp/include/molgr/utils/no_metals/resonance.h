#pragma once

#include "molgr/config.h"
#include "molgr/state.h"
#include "molgr/utils/perf.h"
#include "molgr/utils/resonance.h"

#include <cstddef>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace molgr
{
    namespace no_metals
    {
        namespace resonance
        {
            using RawStateKey =
                std::pair<molgr::resonance::ResonanceStateKey, int>;
            using ProcessedStateKey =
                std::pair<molgr::resonance::ProcessedResonanceKey, int>;

            struct ResonanceSearchSession
            {
                std::set<RawStateKey> seen_raw_states;
                std::map<RawStateKey, std::vector<std::pair<int, int>>> labels_by_state;
                std::set<ProcessedStateKey> seen_processed_states;
                std::size_t next_raw_index = 0;
            };

            std::vector<molgr::state::ReconstructionState> BuildResonanceSeedPool(
                std::vector<molgr::state::ReconstructionState> neighbor_seeds);

            std::vector<molgr::state::ReconstructionState> SearchResonanceCandidates(
                const std::vector<molgr::state::ReconstructionState> &states,
                const molgr::config::MolGRConfig &config,
                molgr::pipeline::perf::RunTimingReducer *timing_reducer = nullptr,
                ResonanceSearchSession *session = nullptr);
        }
    }
}
