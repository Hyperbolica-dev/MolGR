#pragma once

#include <openbabel/mol.h>

#include <string>
#include <vector>

namespace molgr
{
    namespace utils
    {
        bool ReadXyzBlockToMol(const std::string &xyz_block, OpenBabel::OBMol *mol);

        // Parse the XYZ atom records without touching Open Babel format
        // plugins. This is used by batch-side electronic-state validation.
        bool ParseXyzAtomicNumbers(
            const std::string &xyz_block,
            std::vector<int> *atomic_numbers);
    }
}
