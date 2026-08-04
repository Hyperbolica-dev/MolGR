#pragma once

#include "molgr/config.h"

#include <openbabel/elements.h>

namespace molgr::coordination
{
    inline double CoordinationDistanceCutoff(
        int metal_atomic_num,
        int ligand_atomic_num,
        const molgr::config::MetalScoringConfig &config)
    {
        const double metal_radius =
            OpenBabel::OBElements::GetCovalentRad(metal_atomic_num);
        const double ligand_radius =
            OpenBabel::OBElements::GetCovalentRad(ligand_atomic_num);
        return config.metal_access_radius_scale * (metal_radius + ligand_radius) +
               config.metal_coordination_extra_tolerance_angstrom;
    }
}
