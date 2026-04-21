#pragma once

#include <openbabel/atom.h>
#include <openbabel/mol.h>

namespace molgr
{
    namespace reconstruct
    {
        int AssignRadicalDots(OpenBabel::OBAtom &atom);
        bool AssignChargeRadicalForAtom(OpenBabel::OBAtom &atom);
        bool FreshOmolChargeRadical(OpenBabel::OBMol &mol);
    }
}
