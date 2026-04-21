#pragma once

#include <openbabel/mol.h>

namespace molgr
{
    namespace reconstruct
    {
        bool BreakDeformedEne(
            OpenBabel::OBMol &mol,
            int allowed_charge,
            int allowed_radical,
            double tolerance = 5.0);
        bool BreakOneBond(
            OpenBabel::OBMol &mol,
            int &current_charge_deficit,
            int allowed_radical);
    }
}
