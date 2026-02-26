#pragma once

#include "molgr/stages/clean.h"
#include "molgr/stages/eliminate.h"

#include <openbabel/mol.h>

#include <string>
#include <utility>
#include <vector>

namespace molgr
{
    namespace reconstruct
    {
        std::vector<OpenBabel::OBMol> GetOneStepResonance(const OpenBabel::OBMol &mol);
        std::vector<OpenBabel::OBMol> GetRadicalResonances(const OpenBabel::OBMol &mol);
        std::pair<OpenBabel::OBMol, int> ProcessResonance(const OpenBabel::OBMol &mol, int charge);
        std::string SmilesFirstToken(const OpenBabel::OBMol &mol);
    }
}
