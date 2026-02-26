#pragma once

#include <openbabel/mol.h>

namespace molgr
{
    namespace scoring
    {
        double OmolScore(OpenBabel::OBMol &mol);
        double CalcSymmetryPenalty(OpenBabel::OBMol &mol);
        double CalculatePhysChemPenalty(OpenBabel::OBMol &mol);
        double CalculateMetalPenalty(OpenBabel::OBMol &mol);
        double GetDeviationScore(OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom);

    }
}
