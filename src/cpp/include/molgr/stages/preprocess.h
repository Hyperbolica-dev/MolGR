#pragma once

#include <openbabel/mol.h>

namespace molgr
{
    namespace reconstruct
    {
        bool MakeConnections(OpenBabel::OBMol &mol, double factor = 1.4);
        bool PreClean(OpenBabel::OBMol &mol);
        bool ValidateOmol(
            OpenBabel::OBMol &mol,
            int total_charge,
            int total_radical,
            bool emit_warnings = false);
    }
}
