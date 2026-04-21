#pragma once

#include <string>
#include <tuple>
#include <vector>

namespace molgr
{
    struct MetalAtomPosition
    {
        int idx = 0;
        std::string symbol;
        int element_idx = 0;
        int valence = 0;
        int radical_num = 0;
        double position_x = 0.0;
        double position_y = 0.0;
        double position_z = 0.0;
    };

    using ChargedAtomSnapshot = std::tuple<int, double, double, double>;
    using ChargedAtomSnapshotList = std::vector<ChargedAtomSnapshot>;

    namespace metal
    {
        using MetalAtomPosition = ::molgr::MetalAtomPosition;
    }
}
