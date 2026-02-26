#pragma once

#include "molgr/utils/consts.h"

#include <openbabel/mol.h>

#include <string>
#include <vector>

namespace molgr
{
    namespace reconstruct
    {
        const ElementInfo *GetElementInfo(int atomic_num);
        std::vector<int> GetFlatAtomList(OpenBabel::OBMol &mol, const std::string &smarts);
        bool ContainsAtomIdx(const std::vector<int> &atom_indices, int atom_idx);
        int PythonModulo(int value, int modulus);
    }
}
