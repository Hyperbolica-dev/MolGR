#pragma once

#include "molgr/utils/consts.h"
#include "molgr/utils/smarts.h"

#include <openbabel/mol.h>

#include <vector>

namespace molgr
{
    namespace reconstruct
    {
        const ElementInfo *GetElementInfo(int atomic_num);
        std::vector<int> GetFlatAtomList(OpenBabel::OBMol &mol, molgr::smarts::PatternId pattern_id);
        bool ContainsAtomIdx(const std::vector<int> &atom_indices, int atom_idx);
        int PythonModulo(int value, int modulus);
    }
}
