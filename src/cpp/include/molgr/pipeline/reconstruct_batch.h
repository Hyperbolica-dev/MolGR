#pragma once

#include "molgr/config.h"
#include "molgr/diagnostics.h"
#include "molgr/utils/utils.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace molgr::pipeline::reconstruct_batch
{
    struct ReconstructionBatchRequest
    {
        std::string xyz_block;
        int total_charge = 0;
        int total_radical_electrons = 0;
    };

    struct ReconstructionBatchResult
    {
        std::size_t index = 0;
        std::unique_ptr<molgr::utils::MoleculeData> molecule;
        molgr::diagnostics::ReconstructionDiagnostics diagnostics;
    };

    class ReconstructionBatchIterator
    {
    public:
        ReconstructionBatchIterator(
            std::vector<ReconstructionBatchRequest> requests,
            const molgr::config::MolGRConfig &config,
            std::size_t max_workers,
            std::size_t queue_size,
            bool ordered);
        ~ReconstructionBatchIterator();

        ReconstructionBatchIterator(const ReconstructionBatchIterator &) = delete;
        ReconstructionBatchIterator &operator=(const ReconstructionBatchIterator &) = delete;

        std::optional<ReconstructionBatchResult> Next();
        void Close();

    private:
        struct Impl;
        std::unique_ptr<Impl> impl_;
    };
}
