#pragma once

#include "molgr/state.h"

#include <string>

namespace molgr
{
    namespace no_metals
    {
        namespace preparation
        {
            molgr::state::ReconstructionState SeedState(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons);

            molgr::state::ReconstructionState RunLinearPipeline(
                const molgr::state::ReconstructionState &state);
        }
    }
}
