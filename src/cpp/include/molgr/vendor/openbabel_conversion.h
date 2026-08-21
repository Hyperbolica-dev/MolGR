#pragma once

#include <openbabel/mol.h>

#include <string>

namespace molgr
{
    namespace vendor
    {
        namespace openbabel_conversion
        {
            // Open Babel conversion mutates process-global option/plugin registries.
            // Keep the remaining SMILES compatibility surface in one adapter.
            bool ReadSmiles(const std::string &smiles, OpenBabel::OBMol *mol);
            std::string WriteSmilesFirstToken(const OpenBabel::OBMol &mol);
        }
    }
}
