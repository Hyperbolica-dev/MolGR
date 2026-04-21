#pragma once

#include "molgr/types.h"
#include "molgr/utils/utils.h"

#include <openbabel/atom.h>
#include <openbabel/mol.h>

#include <memory>
#include <set>
#include <string>
#include <vector>

namespace molgr
{
    namespace metal
    {
        void ReinsertMetalStates(
            OpenBabel::OBMol &mol,
            const std::vector<MetalAtomPosition> &metals);
    }

    namespace pipeline
    {
        namespace reconstruct_with_metals
        {
            std::set<int> get_possible_metal_radicals(const std::string &metal_symbol, int valence);

            std::vector<molgr::metal::MetalAtomPosition> build_metal_states(const OpenBabel::OBAtom &obatom);

            void combine_metal_with_omol(
                OpenBabel::OBMol &mol,
                const std::vector<molgr::metal::MetalAtomPosition> &metals);

            std::unique_ptr<molgr::utils::MoleculeData> Xyz2OmolMolData(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons);
        }
    }
}
