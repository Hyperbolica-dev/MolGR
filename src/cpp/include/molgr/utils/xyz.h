#pragma once

#include <openbabel/mol.h>

#include <string>

namespace molgr
{
    namespace utils
    {
        bool ReadXyzBlockToMol(const std::string &xyz_block, OpenBabel::OBMol *mol);
    }
}
