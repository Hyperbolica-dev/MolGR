#pragma once

#include "molgr/config.h"
#include "molgr/state.h"
#include "molgr/types.h"

#include <openbabel/atom.h>
#include <openbabel/mol.h>

#include <set>
#include <string>
#include <vector>

namespace molgr
{
    namespace metal
    {
        namespace preparation
        {
            std::set<int> GetPossibleMetalRadicals(const std::string &metal_symbol, int valence);

            std::vector<molgr::metal::MetalAtomPosition> BuildMetalStates(
                OpenBabel::OBAtom &obatom,
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());

            void ReinsertMetalStates(
                OpenBabel::OBMol &mol,
                const std::vector<molgr::metal::MetalAtomPosition> &metals);

            void CombineMetalWithOmol(
                OpenBabel::OBMol &mol,
                const std::vector<molgr::metal::MetalAtomPosition> &metals);

            molgr::state::MetalPreparationState PrepareMetalState(
                const std::string &xyz_block,
                int total_charge,
                int total_radical_electrons,
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());
        }
    }
}
