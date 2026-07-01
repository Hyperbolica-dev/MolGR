#pragma once

#include "molgr/state.h"

#include <openbabel/mol.h>

#include <memory>
#include <string>

namespace molgr
{
    namespace no_metals
    {
        namespace preparation
        {
            std::shared_ptr<OpenBabel::OBMol> SeedOmolFromXyzBlock(
                const std::string &xyz_block);

            std::shared_ptr<OpenBabel::OBMol> NormalizeSeedOmolCopy(
                const OpenBabel::OBMol &seed_omol);

            molgr::state::ReconstructionState BuildSeedState(
                std::shared_ptr<OpenBabel::OBMol> omol,
                int total_charge,
                int total_radical_electrons);

            molgr::state::ReconstructionState SeedStateFromOmol(
                const OpenBabel::OBMol &seed_omol,
                int total_charge,
                int total_radical_electrons);

            molgr::state::ReconstructionState SeedState(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons);

            molgr::state::ReconstructionState RunLinearPipeline(
                const molgr::state::ReconstructionState &state);
        }
    }
}
