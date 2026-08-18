#include "molgr/pipeline/reconstruct_batch.h"

#include "molgr/pipeline/reconstruct_with_metals.h"
#include "molgr/utils/parallel.h"
#include "molgr/utils/xyz.h"

#include <algorithm>
#include <condition_variable>
#include <deque>
#include <exception>
#include <limits>
#include <map>
#include <mutex>
#include <utility>
#include <vector>

namespace molgr::pipeline::reconstruct_batch
{
    namespace
    {
        bool ValidateElectronicTarget(
            const ReconstructionBatchRequest &request,
            molgr::diagnostics::ReconstructionDiagnostics *diagnostics)
        {
            if (request.total_radical_electrons < 0)
            {
                diagnostics->Fail(
                    "INVALID_ELECTRONIC_TARGET",
                    "input.electronic_state",
                    "The requested radical-electron count is invalid.");
                return false;
            }

            std::vector<int> atomic_numbers;
            if (!molgr::utils::ParseXyzAtomicNumbers(request.xyz_block, &atomic_numbers))
            {
                diagnostics->Fail(
                    "INVALID_XYZ",
                    "input.parse",
                    "The XYZ block could not be parsed.");
                return false;
            }

            int total_electrons = -request.total_charge;
            for (const int atomic_num : atomic_numbers)
            {
                total_electrons += atomic_num;
            }
            if (total_electrons < 0)
            {
                diagnostics->Fail(
                    "INVALID_ELECTRONIC_TARGET",
                    "input.electronic_state",
                    "The requested charge leaves a negative total electron count.");
                return false;
            }
            if (request.total_radical_electrons > total_electrons)
            {
                diagnostics->Fail(
                    "INVALID_ELECTRONIC_TARGET",
                    "input.electronic_state",
                    "The requested radical-electron count exceeds the total electron count.");
                return false;
            }
            if ((total_electrons % 2) != (request.total_radical_electrons % 2))
            {
                diagnostics->Fail(
                    "INVALID_ELECTRONIC_TARGET",
                    "input.electronic_state",
                    "The requested charge and radical-electron count have incompatible parity.");
                return false;
            }
            return true;
        }

        molgr::config::MolGRConfig BatchWorkerConfig(
            const molgr::config::MolGRConfig &config,
            std::size_t worker_count)
        {
            auto worker_config = config;
            if (worker_count > 1)
            {
                // The batch scheduler owns the outer worker budget. Nested
                // target-bucket/candidate pools would oversubscribe the process
                // and can deadlock the shared native executor.
                worker_config.cpp_backend.max_threads = 1;
                worker_config.cpp_backend.enable_target_bucket_parallelism = false;
                worker_config.cpp_backend.enable_candidate_scoring_parallelism = false;
                worker_config.cpp_backend.target_bucket_parallel_max_threads = 1;
            }
            return worker_config;
        }
    }

    struct ReconstructionBatchIterator::Impl
    {
        Impl(
            std::vector<ReconstructionBatchRequest> requests_in,
            const molgr::config::MolGRConfig &config,
            std::size_t max_workers,
            std::size_t queue_size_in,
            bool ordered_in)
            : requests(std::move(requests_in)),
              ordered(ordered_in),
              queue_size(std::max<std::size_t>(1, queue_size_in))
        {
            auto &executor = molgr::utils::parallel::detail::ParallelExecutor::Instance();
            const std::size_t requested_workers = max_workers == 0
                ? molgr::utils::parallel::ConfiguredParallelism(config, requests.size())
                : max_workers;
            const std::size_t bounded_workers = std::max<std::size_t>(1, requested_workers);
            worker_count = requests.empty()
                ? 0
                : std::min(
                      {bounded_workers, requests.size(), executor.HelperCapacity()});
            worker_config = BatchWorkerConfig(config, worker_count);
            if (queue_size > std::numeric_limits<std::size_t>::max() - worker_count)
            {
                buffered_capacity = std::numeric_limits<std::size_t>::max();
            }
            else
            {
                // Results held by running workers already consume memory even
                // before they enter the ready queue. Bound both populations
                // together without blocking a shared executor worker.
                buffered_capacity = queue_size + worker_count;
            }
            std::lock_guard<std::mutex> lock(mutex);
            ScheduleAvailableLocked();
        }

        ~Impl()
        {
            Close();
        }

        ReconstructionBatchResult Process(std::size_t index)
        {
            ReconstructionBatchResult result;
            result.index = index;
            const auto &request = requests[index];
            result.diagnostics.details["total_charge"] = std::to_string(request.total_charge);
            result.diagnostics.details["total_radical_electrons"] =
                std::to_string(request.total_radical_electrons);

            try
            {
                if (!ValidateElectronicTarget(request, &result.diagnostics))
                {
                    return result;
                }
                result.molecule = molgr::pipeline::reconstruct_with_metals::Xyz2OmolMolData(
                    request.xyz_block,
                    request.total_charge,
                    request.total_radical_electrons,
                    worker_config,
                    &result.diagnostics);
                if (!result.molecule && result.diagnostics.code.empty())
                {
                    result.diagnostics.Fail(
                        "NO_VALID_RECONSTRUCTION",
                        "reconstruction",
                        "The C++ reconstruction pipeline produced no molecule.");
                }
            }
            catch (const std::exception &exception)
            {
                result.diagnostics.details["cause_type"] = "std::exception";
                result.diagnostics.details["cause_message"] = exception.what();
                result.diagnostics.Fail(
                    "BACKEND_EXCEPTION",
                    "reconstruction",
                    "The C++ reconstruction pipeline raised an exception.");
            }
            catch (...)
            {
                result.diagnostics.Fail(
                    "BACKEND_EXCEPTION",
                    "reconstruction",
                    "The C++ reconstruction pipeline raised an unknown exception.");
            }
            return result;
        }

        void SetFatalExceptionLocked(std::exception_ptr exception)
        {
            if (!fatal_exception)
            {
                fatal_exception = std::move(exception);
            }
            stopped = true;
            result_cv.notify_all();
        }

        void ScheduleAvailableLocked()
        {
            auto &executor = molgr::utils::parallel::detail::ParallelExecutor::Instance();
            while (!stopped && next_input < requests.size() &&
                   active_tasks < worker_count &&
                   ready_count + active_tasks < buffered_capacity)
            {
                const std::size_t index = next_input++;
                ++active_tasks;
                try
                {
                    executor.Submit(
                        [this, index]()
                        {
                            try
                            {
                                Complete(Process(index));
                            }
                            catch (...)
                            {
                                CompleteFatal(std::current_exception());
                            }
                        });
                }
                catch (...)
                {
                    --next_input;
                    --active_tasks;
                    SetFatalExceptionLocked(std::current_exception());
                    workers_done_cv.notify_all();
                    return;
                }
            }
        }

        void Complete(ReconstructionBatchResult result)
        {
            std::lock_guard<std::mutex> lock(mutex);
            if (stopped)
            {
                --active_tasks;
                workers_done_cv.notify_all();
                return;
            }
            const auto index = result.index;
            if (ordered)
            {
                ordered_results.emplace(index, std::move(result));
            }
            else
            {
                unordered_results.push_back(std::move(result));
            }
            ++ready_count;
            ++completed_count;
            --active_tasks;
            ScheduleAvailableLocked();
            result_cv.notify_all();
            workers_done_cv.notify_all();
        }

        void CompleteFatal(std::exception_ptr exception)
        {
            std::lock_guard<std::mutex> lock(mutex);
            if (active_tasks > 0)
            {
                --active_tasks;
            }
            SetFatalExceptionLocked(std::move(exception));
            workers_done_cv.notify_all();
        }

        std::optional<ReconstructionBatchResult> Next()
        {
            std::unique_lock<std::mutex> lock(mutex);
            while (true)
            {
                if (fatal_exception)
                {
                    std::rethrow_exception(fatal_exception);
                }
                if (stopped)
                {
                    return std::nullopt;
                }
                if (ordered)
                {
                    auto result = ordered_results.find(next_output_index);
                    if (result != ordered_results.end())
                    {
                        auto output = std::move(result->second);
                        ordered_results.erase(result);
                        ++next_output_index;
                        --ready_count;
                        ScheduleAvailableLocked();
                        return output;
                    }
                }
                else if (!unordered_results.empty())
                {
                    auto output = std::move(unordered_results.front());
                    unordered_results.pop_front();
                    --ready_count;
                    ScheduleAvailableLocked();
                    return output;
                }
                if (ready_count == 0 && completed_count == requests.size())
                {
                    return std::nullopt;
                }
                result_cv.wait(lock);
            }
        }

        void Close()
        {
            std::unique_lock<std::mutex> lock(mutex);
            stopped = true;
            result_cv.notify_all();
            workers_done_cv.wait(lock, [this]() { return active_tasks == 0; });
        }

        std::vector<ReconstructionBatchRequest> requests;
        molgr::config::MolGRConfig worker_config;
        std::size_t worker_count = 0;
        bool ordered = false;
        std::size_t queue_size = 1;
        std::size_t buffered_capacity = 1;
        std::size_t next_input = 0;
        std::size_t next_output_index = 0;
        std::size_t ready_count = 0;
        std::size_t completed_count = 0;
        std::size_t active_tasks = 0;
        bool stopped = false;
        std::exception_ptr fatal_exception;
        std::mutex mutex;
        std::condition_variable result_cv;
        std::condition_variable workers_done_cv;
        std::deque<ReconstructionBatchResult> unordered_results;
        std::map<std::size_t, ReconstructionBatchResult> ordered_results;
    };

    ReconstructionBatchIterator::ReconstructionBatchIterator(
        std::vector<ReconstructionBatchRequest> requests,
        const molgr::config::MolGRConfig &config,
        std::size_t max_workers,
        std::size_t queue_size,
        bool ordered)
        : impl_(std::make_unique<Impl>(
              std::move(requests), config, max_workers, queue_size, ordered))
    {
    }

    ReconstructionBatchIterator::~ReconstructionBatchIterator() = default;

    std::optional<ReconstructionBatchResult> ReconstructionBatchIterator::Next()
    {
        return impl_->Next();
    }

    void ReconstructionBatchIterator::Close()
    {
        impl_->Close();
    }
}
