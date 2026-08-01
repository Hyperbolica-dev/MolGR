#pragma once

#include "molgr/config.h"
#include "molgr/state.h"
#include "molgr/utils/organic_topology.h"

#include <tuple>
#include <vector>

namespace molgr
{
    namespace no_metals
    {
        namespace selection
        {
            using NoMetalTopologySelectionKey =
                std::tuple<int, int, int, double, double, double, double>;
            using NoMetalGraphTieBreakKey = std::vector<int>;

            molgr::organic_topology::OrganicTopologyMetrics AnnotateNoMetalCandidateTopology(
                molgr::state::ReconstructionState &candidate,
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());

            double ScoreReconstructionCandidate(
                molgr::state::ReconstructionState &candidate,
                const molgr::config::MolGRConfig &config);

            NoMetalTopologySelectionKey NoMetalCandidateTopologySelectionKey(
                molgr::state::ReconstructionState &candidate,
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());

            std::tuple<int, int, int, double, double, double, double, int, int, double>
            NoMetalCandidateSelectionKey(
                molgr::state::ReconstructionState &candidate,
                const molgr::config::MolGRConfig &config);

            NoMetalGraphTieBreakKey NoMetalCandidateGraphTieBreakKey(
                const molgr::state::ReconstructionState &candidate);
        }
    }
}
