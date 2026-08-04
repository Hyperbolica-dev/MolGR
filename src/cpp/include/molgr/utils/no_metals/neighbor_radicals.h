#pragma once

#include "molgr/state.h"

#include <openbabel/mol.h>

#include <optional>
#include <utility>
#include <vector>

namespace molgr
{
    namespace no_metals
    {
        namespace neighbor_radicals
        {
            std::vector<std::pair<int, int>> NeighborRadicalBondPairs(
                const OpenBabel::OBMol &mol);

            std::vector<molgr::state::ReconstructionState> EnumerateNeighborRadicalSeeds(
                const molgr::state::ReconstructionState &state,
                std::optional<int> exact_discrepancy = std::nullopt);
        }
    }
}
