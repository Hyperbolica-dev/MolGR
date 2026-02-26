#pragma once

#include <openbabel/mol.h>

#include <memory>
#include <string>

namespace molgr
{
    namespace reconstruct
    {
        std::unique_ptr<OpenBabel::OBMol> ReconstructFromXYZNoMetal(
            const std::string &xyz_block,
            int total_charge,
            int total_radical);

        void MakeConnections(OpenBabel::OBMol &mol, double factor = 1.4);
        void PreClean(OpenBabel::OBMol &mol);
        bool ValidateOmol(
            OpenBabel::OBMol &mol,
            int total_charge,
            int total_radical,
            bool emit_warnings = false);
    }
}
