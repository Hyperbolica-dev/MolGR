#pragma once

#include "molgr/config.h"
#include "molgr/state.h"
#include "molgr/utils/organic_topology.h"

#include <tuple>

namespace molgr
{
    namespace no_metals
    {
        namespace selection
        {
            using NoMetalTopologySelectionKey = std::tuple<double, int, double, double, double>;

            molgr::organic_topology::OrganicTopologyMetrics AnnotateNoMetalCandidateTopology(
                molgr::state::ReconstructionState &candidate,
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());

            double ScoreReconstructionCandidate(
                molgr::state::ReconstructionState &candidate,
                const molgr::config::MolGRConfig &config);

            NoMetalTopologySelectionKey NoMetalCandidateTopologySelectionKey(
                molgr::state::ReconstructionState &candidate,
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());

            std::tuple<double, int, double, double, double, int, double>
            NoMetalCandidateSelectionKey(
                molgr::state::ReconstructionState &candidate,
                const molgr::config::MolGRConfig &config);
        }
    }
}
