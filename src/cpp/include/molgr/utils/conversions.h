#pragma once

#include "molgr/utils/utils.h"

namespace molgr
{
    namespace utils
    {
        OpenBabel::OBMol CloneMolTopologyOnly(const OpenBabel::OBMol &mol);
        OpenBabel::OBMol MolFromMoleculeData(const MoleculeData &data);
        MoleculeData MoleculeDataFromOBMol(const OpenBabel::OBMol &mol);
    }
}
