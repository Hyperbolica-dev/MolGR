#pragma once

#include <openbabel/mol.h>

namespace molgr
{
    namespace reconstruct
    {
        bool CleanCarbeneNeighborUnsaturated(OpenBabel::OBMol &mol);
        bool CleanNeighborRadicals(OpenBabel::OBMol &mol);
        bool CleanResonances(OpenBabel::OBMol &mol);
    }
}
