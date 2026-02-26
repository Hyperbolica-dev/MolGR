#pragma once

#include "molgr/utils/utils.h"

namespace molgr
{
    namespace utils
    {
        OpenBabel::OBMol MolFromMoleculeData(const MoleculeData &data);
        MoleculeData MoleculeDataFromOBMol(const OpenBabel::OBMol &mol);
    }
}
