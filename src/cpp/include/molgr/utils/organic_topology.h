#pragma once

#include <openbabel/bond.h>
#include <openbabel/mol.h>

#include <vector>

namespace molgr
{
    namespace organic_topology
    {
        struct OrganicTopologyMetrics
        {
            int aromatic_atom_count = 0;
            int aromatic_ring_count = 0;
            int conjugated_atom_count = 0;
            int conjugated_bond_count = 0;
            int max_conjugated_component_size = 0;
            std::vector<int> conjugated_atom_indices;
        };

        bool IsConjugatedBond(const OpenBabel::OBBond &bond);
        OrganicTopologyMetrics ComputeOrganicTopologyMetrics(const OpenBabel::OBMol &mol);
    }
}
