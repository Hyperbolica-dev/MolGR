#pragma once

#include "molgr/config.h"
#include "molgr/utils/utils.h"

#include <memory>
#include <string>

namespace molgr
{
    namespace pipeline
    {
        namespace reconstruct_with_metals
        {
            std::unique_ptr<molgr::utils::MoleculeData> Xyz2OmolMolData(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());
        }
    }
}
