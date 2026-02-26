#pragma once

#include <openbabel/mol.h>

namespace molgr
{
    namespace reconstruct
    {
        void CleanCarbeneNeighborUnsaturated(OpenBabel::OBMol &mol);
        void CleanNeighborRadicals(OpenBabel::OBMol &mol);
        void CleanResonances(OpenBabel::OBMol &mol);
    }
}
