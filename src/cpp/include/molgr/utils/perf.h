#pragma once

#include <mutex>

namespace molgr
{
    namespace pipeline
    {
        namespace perf
        {
            struct RunTimingBreakdown
            {
                double no_metal_pipeline_ms = 0.0;
                double no_metal_linear_pipeline_ms = 0.0;
                double no_metal_validate_ms = 0.0;
                double resonance_handling_enumeration_ms = 0.0;
                double resonance_walk_ms = 0.0;
                double resonance_prepare_ms = 0.0;
                double resonance_dedup_score_ms = 0.0;
                double resonance_score_ms = 0.0;
                double resonance_topology_ms = 0.0;
                double resonance_raw_candidates = 0.0;
                double resonance_pruned_expansions = 0.0;
                double resonance_prepared_candidates = 0.0;
                double resonance_valid_candidates = 0.0;
                double resonance_dedup_candidates = 0.0;
                double metal_enumeration_combination_ms = 0.0;
                double force_field_total_ms = 0.0;
                double force_field_cache_key_ms = 0.0;
                double force_field_prepare_ms = 0.0;
                double force_field_setup_key_ms = 0.0;
                double force_field_setup_ms = 0.0;
                double force_field_energy_ms = 0.0;
                double force_field_calls = 0.0;
            };

            class RunTimingReducer
            {
            public:
                void AddNoMetalPipelineMs(double delta_ms);
                void AddNoMetalLinearPipelineMs(double delta_ms);
                void AddNoMetalValidateMs(double delta_ms);
                void AddResonanceHandlingEnumerationMs(double delta_ms);
                void AddResonanceWalkMs(double delta_ms);
                void AddResonancePrepareMs(double delta_ms);
                void AddResonanceDedupScoreMs(double delta_ms);
                void AddResonanceScoreMs(double delta_ms);
                void AddResonanceTopologyMs(double delta_ms);
                void AddResonanceRawCandidates(double delta_count);
                void AddResonancePrunedExpansions(double delta_count);
                void AddResonancePreparedCandidates(double delta_count);
                void AddResonanceValidCandidates(double delta_count);
                void AddResonanceDedupCandidates(double delta_count);
                void AddMetalEnumerationCombinationMs(double delta_ms);
                void AddForceFieldTotalMs(double delta_ms);
                void AddForceFieldCacheKeyMs(double delta_ms);
                void AddForceFieldPrepareMs(double delta_ms);
                void AddForceFieldSetupKeyMs(double delta_ms);
                void AddForceFieldSetupMs(double delta_ms);
                void AddForceFieldEnergyMs(double delta_ms);
                void AddForceFieldCalls(double delta_count);
                RunTimingBreakdown Snapshot() const;

            private:
                mutable std::mutex mutex_;
                RunTimingBreakdown timing_;
            };

            class RunTimingScope
            {
            public:
                RunTimingScope();
                RunTimingReducer &Reducer();
                const RunTimingReducer &Reducer() const;
                ~RunTimingScope();

            private:
                RunTimingReducer reducer_;
                RunTimingReducer *previous_ = nullptr;
            };

            class ActiveRunTimingReducerScope
            {
            public:
                explicit ActiveRunTimingReducerScope(RunTimingReducer *reducer);
                ~ActiveRunTimingReducerScope();

            private:
                RunTimingReducer *previous_ = nullptr;
            };

            RunTimingReducer *GetActiveRunTimingReducer();
            RunTimingBreakdown GetRunTimingBreakdown();
            void SetRunTimingBreakdown(const RunTimingBreakdown &timing);
        }
    }
}
