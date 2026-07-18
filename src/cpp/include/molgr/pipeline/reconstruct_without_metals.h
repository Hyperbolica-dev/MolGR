#pragma once

#include "molgr/config.h"
#include "molgr/state.h"
#include "molgr/utils/perf.h"
#include "molgr/utils/utils.h"

#include <openbabel/mol.h>

#include <optional>
#include <string>
#include <vector>

namespace molgr
{
    namespace pipeline
    {
        namespace reconstruct_without_metals
        {
            struct DebugNoMetalCandidateSummary
            {
                std::string smiles;
                int resonance_index = -1;
                double score = 0.0;
                double aromatic_stability_score = 0.0;
                int aromatic_atom_count = 0;
                int max_conjugated_component_size = 0;
                int conjugated_atom_count = 0;
                int conjugated_bond_count = 0;
                int formal_charge_absolute_sum = 0;
                double conjugation_charge_penalty = 0.0;
                double adjusted_max_conjugated_component_size = 0.0;
                double adjusted_conjugated_atom_count = 0.0;
                double adjusted_conjugated_bond_count = 0.0;
                int excess_radical_labels = 0;
            };

            std::optional<molgr::state::ReconstructionState> XyzToOmolNoMetalState(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig(),
                perf::RunTimingReducer *timing_reducer = nullptr,
                bool preheat_score_bundle = true);

            std::optional<molgr::state::ReconstructionState> SeedOmolToOmolNoMetalState(
                const OpenBabel::OBMol &seed_omol,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig(),
                perf::RunTimingReducer *timing_reducer = nullptr,
                bool preheat_score_bundle = true);

            std::optional<molgr::state::ReconstructionState> SeedOmolCopyToOmolNoMetalState(
                std::shared_ptr<OpenBabel::OBMol> seed_omol,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig(),
                perf::RunTimingReducer *timing_reducer = nullptr,
                bool preheat_score_bundle = true);

            std::unique_ptr<molgr::utils::MoleculeData> XyzToMolDataNoMetal(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());

            std::vector<DebugNoMetalCandidateSummary> DebugNoMetalResonanceCandidateSummaries(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());
        }
    }
}
