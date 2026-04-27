#pragma once

#include "molgr/config.h"

#include <openbabel/atom.h>

#include <string>
#include <vector>

namespace molgr
{
    namespace metal_radical_inference
    {
        struct MetalRadicalInferenceResult
        {
            std::vector<int> radical_counts;
            int effective_d_electrons = 0;
            int residual_sp_electrons = 0;
            int remaining_f_electrons = 0;
            int coordination_number = 0;
            std::string geometry = "free_ion";
            double field_score = 0.0;
            std::string field_strength = "weak";
        };

        MetalRadicalInferenceResult InferMetalRadicalState(
            OpenBabel::OBAtom &metal_atom,
            int valence,
            const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());

        std::vector<int> InferMetalRadicalCounts(
            OpenBabel::OBAtom &metal_atom,
            int valence,
            const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());
    }
}
