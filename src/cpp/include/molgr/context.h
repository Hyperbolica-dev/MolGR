#pragma once

#include "molgr/config.h"
#include "molgr/state.h"

#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace molgr
{
    namespace context
    {
        struct ReconstructionContext
        {
            std::string xyz_block;
            int total_charge = 0;
            int total_radical_electrons = 0;
            const molgr::config::MolGRConfig *config = &molgr::config::GetDefaultConfig();
            std::vector<std::string> phase_history;
            molgr::state::MetadataMap metadata;

            const molgr::config::MolGRConfig &Config() const
            {
                return config == nullptr ? molgr::config::GetDefaultConfig() : *config;
            }

            bool HasValidRadicalTarget() const
            {
                return total_radical_electrons >= 0;
            }
        };

        struct TargetBucket
        {
            int no_metal_charge = 0;
            int no_metal_radicals = 0;

            bool operator<(const TargetBucket &other) const
            {
                return std::tie(no_metal_charge, no_metal_radicals) <
                       std::tie(other.no_metal_charge, other.no_metal_radicals);
            }

            bool operator==(const TargetBucket &other) const
            {
                return no_metal_charge == other.no_metal_charge &&
                       no_metal_radicals == other.no_metal_radicals;
            }

            std::pair<int, int> AsPair() const
            {
                return {no_metal_charge, no_metal_radicals};
            }
        };
    }
}
