#pragma once

#include <openbabel/atom.h>
#include <openbabel/mol.h>

namespace molgr
{
    namespace reconstruct
    {
        int AssignRadicalDots(OpenBabel::OBAtom &atom);
        void AssignChargeRadicalForAtom(OpenBabel::OBAtom &atom);
        void FreshOmolChargeRadical(OpenBabel::OBMol &mol);
    }
}
