#pragma once

#include "molgr/state.h"

#include <vector>

namespace molgr
{
    namespace no_metals
    {
        namespace recovery
        {
            std::vector<molgr::state::ReconstructionState> EnumerateDeformedPiRecoverySeeds(
                const std::vector<molgr::state::ReconstructionState> &states);

            std::vector<molgr::state::ReconstructionState> EnumerateBondBreakRecoverySeeds(
                const std::vector<molgr::state::ReconstructionState> &states);
        }
    }
}
