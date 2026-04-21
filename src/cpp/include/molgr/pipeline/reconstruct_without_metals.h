#pragma once

#include "molgr/pipeline/resonance.h"
#include "molgr/state.h"
#include "molgr/stages/break_bond.h"
#include "molgr/stages/clean.h"
#include "molgr/stages/eliminate.h"
#include "molgr/stages/fresh.h"
#include "molgr/stages/preprocess.h"
#include "molgr/utils/utils.h"

#include <mutex>
#include <optional>
#include <string>

namespace molgr
{
    namespace pipeline
    {
        namespace perf
        {
            struct RunTimingBreakdown
            {
                double no_metal_pipeline_ms = 0.0;
                double resonance_handling_enumeration_ms = 0.0;
                double metal_enumeration_combination_ms = 0.0;
            };

            class RunTimingReducer
            {
            public:
                void AddNoMetalPipelineMs(double delta_ms);
                void AddResonanceHandlingEnumerationMs(double delta_ms);
                void AddMetalEnumerationCombinationMs(double delta_ms);
                RunTimingBreakdown Snapshot() const;

            private:
                mutable std::mutex mutex_;
                RunTimingBreakdown timing_;
            };

            class RunTimingScope
            {
            public:
                RunTimingReducer &Reducer();
                const RunTimingReducer &Reducer() const;
                ~RunTimingScope();

            private:
                RunTimingReducer reducer_;
            };

            RunTimingBreakdown GetRunTimingBreakdown();
            void SetRunTimingBreakdown(const RunTimingBreakdown &timing);
        }

        namespace reconstruct_without_metals
        {
            std::optional<molgr::state::ReconstructionState> XyzToOmolNoMetalState(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons,
                perf::RunTimingReducer *timing_reducer = nullptr,
                bool preheat_score_bundle = true);

            std::unique_ptr<molgr::utils::MoleculeData> XyzToMolDataNoMetal(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons);
        }
    }
}
