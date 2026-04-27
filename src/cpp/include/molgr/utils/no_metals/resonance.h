#pragma once

#include "molgr/config.h"
#include "molgr/state.h"

#include <vector>

namespace molgr
{
    namespace no_metals
    {
        namespace resonance
        {
            std::vector<molgr::state::ReconstructionState> RecoverResonanceCandidates(
                const molgr::state::ReconstructionState &state,
                const molgr::config::MolGRConfig &config);
        }
    }
}
