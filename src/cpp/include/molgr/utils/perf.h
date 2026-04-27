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
                double resonance_handling_enumeration_ms = 0.0;
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
                void AddResonanceHandlingEnumerationMs(double delta_ms);
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

            RunTimingReducer *GetActiveRunTimingReducer();
            RunTimingBreakdown GetRunTimingBreakdown();
            void SetRunTimingBreakdown(const RunTimingBreakdown &timing);
        }
    }
}
