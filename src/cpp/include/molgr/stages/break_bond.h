#pragma once

#include <openbabel/mol.h>

namespace molgr
{
    namespace reconstruct
    {
        void BreakDeformedEne(OpenBabel::OBMol &mol, int allowed_charge, int allowed_radical, double tolerance = 5.0);
        void BreakOneBond(OpenBabel::OBMol &mol, int &current_charge_deficit, int allowed_radical);
    }
}
