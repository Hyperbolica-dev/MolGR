#include "molgr/utils/perf.h"

#include <mutex>

namespace molgr
{
    namespace pipeline
    {
        namespace perf
        {
            namespace
            {
                std::mutex t_last_run_timing_breakdown_mutex;
                molgr::pipeline::perf::RunTimingBreakdown t_last_run_timing_breakdown;
                thread_local molgr::pipeline::perf::RunTimingReducer *t_active_run_timing_reducer = nullptr;
            }

            void RunTimingReducer::AddNoMetalPipelineMs(double delta_ms)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                timing_.no_metal_pipeline_ms += delta_ms;
            }

            void RunTimingReducer::AddResonanceHandlingEnumerationMs(double delta_ms)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                timing_.resonance_handling_enumeration_ms += delta_ms;
            }

            void RunTimingReducer::AddMetalEnumerationCombinationMs(double delta_ms)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                timing_.metal_enumeration_combination_ms += delta_ms;
            }

            void RunTimingReducer::AddForceFieldTotalMs(double delta_ms)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                timing_.force_field_total_ms += delta_ms;
            }

            void RunTimingReducer::AddForceFieldCacheKeyMs(double delta_ms)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                timing_.force_field_cache_key_ms += delta_ms;
            }

            void RunTimingReducer::AddForceFieldPrepareMs(double delta_ms)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                timing_.force_field_prepare_ms += delta_ms;
            }

            void RunTimingReducer::AddForceFieldSetupKeyMs(double delta_ms)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                timing_.force_field_setup_key_ms += delta_ms;
            }

            void RunTimingReducer::AddForceFieldSetupMs(double delta_ms)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                timing_.force_field_setup_ms += delta_ms;
            }

            void RunTimingReducer::AddForceFieldEnergyMs(double delta_ms)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                timing_.force_field_energy_ms += delta_ms;
            }

            void RunTimingReducer::AddForceFieldCalls(double delta_count)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                timing_.force_field_calls += delta_count;
            }

            molgr::pipeline::perf::RunTimingBreakdown RunTimingReducer::Snapshot() const
            {
                std::lock_guard<std::mutex> lock(mutex_);
                return timing_;
            }

            RunTimingReducer &RunTimingScope::Reducer()
            {
                return reducer_;
            }

            const RunTimingReducer &RunTimingScope::Reducer() const
            {
                return reducer_;
            }

            RunTimingScope::RunTimingScope()
                : previous_(t_active_run_timing_reducer)
            {
                t_active_run_timing_reducer = &reducer_;
            }

            RunTimingScope::~RunTimingScope()
            {
                t_active_run_timing_reducer = previous_;
                SetRunTimingBreakdown(reducer_.Snapshot());
            }

            RunTimingReducer *GetActiveRunTimingReducer()
            {
                return t_active_run_timing_reducer;
            }

            molgr::pipeline::perf::RunTimingBreakdown GetRunTimingBreakdown()
            {
                std::lock_guard<std::mutex> lock(t_last_run_timing_breakdown_mutex);
                return t_last_run_timing_breakdown;
            }

            void SetRunTimingBreakdown(const RunTimingBreakdown &timing)
            {
                std::lock_guard<std::mutex> lock(t_last_run_timing_breakdown_mutex);
                t_last_run_timing_breakdown = timing;
            }
        }
    }
}
