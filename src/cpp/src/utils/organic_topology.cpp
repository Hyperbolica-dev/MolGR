#include "molgr/utils/organic_topology.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/ring.h>
#include <openbabel/typer.h>

#include <algorithm>
#include <map>
#include <set>
#include <vector>

namespace
{
    bool AtomHasOddSpin(const OpenBabel::OBAtom &atom)
    {
        return atom.GetSpinMultiplicity() % 2 == 1;
    }

    bool AtomHasPiFeature(const OpenBabel::OBAtom &atom, int excluded_bond_idx = -1)
    {
        if (atom.IsAromatic() || atom.GetFormalCharge() != 0 || AtomHasOddSpin(atom))
        {
            return true;
        }
        FOR_BONDS_OF_ATOM(bond_iter, const_cast<OpenBabel::OBAtom *>(&atom))
        {
            OpenBabel::OBBond &bond = *bond_iter;
            if (excluded_bond_idx >= 0 && static_cast<int>(bond.GetIdx()) == excluded_bond_idx)
            {
                continue;
            }
            const OpenBabel::OBAtom *other = bond.GetNbrAtom(const_cast<OpenBabel::OBAtom *>(&atom));
            if (other == nullptr || other->GetAtomicNum() == 1)
            {
                continue;
            }
            if (bond.IsAromatic() || bond.GetBondOrder() >= 2)
            {
                return true;
            }
        }
        return false;
    }

    std::set<int> ValidatedConjugatedBondIndices(OpenBabel::OBMol &mol)
    {
        std::set<int> conjugated_bond_indices;
        FOR_BONDS_OF_MOL(bond_iter, mol)
        {
            OpenBabel::OBBond &bond = *bond_iter;
            OpenBabel::OBAtom *begin_atom = bond.GetBeginAtom();
            OpenBabel::OBAtom *end_atom = bond.GetEndAtom();
            if (begin_atom == nullptr || end_atom == nullptr ||
                begin_atom->GetAtomicNum() == 1 || end_atom->GetAtomicNum() == 1)
            {
                continue;
            }

            const int bond_idx = static_cast<int>(bond.GetIdx());
            const int bond_order = bond.GetBondOrder();
            if (bond.IsAromatic() || bond_order >= 2)
            {
                conjugated_bond_indices.insert(bond_idx);
                continue;
            }
            if (bond_order != 1)
            {
                continue;
            }
            if (!AtomHasPiFeature(*begin_atom, bond_idx) || !AtomHasPiFeature(*end_atom, bond_idx))
            {
                continue;
            }
            conjugated_bond_indices.insert(bond_idx);
        }
        return conjugated_bond_indices;
    }
}

namespace molgr
{
    namespace organic_topology
    {
        bool IsConjugatedBond(const OpenBabel::OBBond &bond)
        {
            if (bond.IsAromatic() || bond.GetBondOrder() >= 2)
            {
                return true;
            }
            if (bond.GetBondOrder() != 1)
            {
                return false;
            }
            const OpenBabel::OBAtom *begin_atom = bond.GetBeginAtom();
            const OpenBabel::OBAtom *end_atom = bond.GetEndAtom();
            if (begin_atom == nullptr || end_atom == nullptr)
            {
                return false;
            }
            return AtomHasPiFeature(*begin_atom, static_cast<int>(bond.GetIdx())) &&
                   AtomHasPiFeature(*end_atom, static_cast<int>(bond.GetIdx()));
        }

        OrganicTopologyMetrics ComputeOrganicTopologyMetrics(const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBMol working(mol);
            working.FindRingAtomsAndBonds();
            working.SetAromaticPerceived(false);
            OpenBabel::OBAromaticTyper().AssignAromaticFlags(working);

            OrganicTopologyMetrics metrics;
            FOR_ATOMS_OF_MOL(atom_iter, working)
            {
                OpenBabel::OBAtom &atom = *atom_iter;
                if (atom.GetAtomicNum() != 1 && atom.IsAromatic())
                {
                    ++metrics.aromatic_atom_count;
                }
            }

            for (OpenBabel::OBRing *ring : working.GetSSSR())
            {
                if (ring != nullptr && ring->IsAromatic())
                {
                    ++metrics.aromatic_ring_count;
                }
            }

            const std::set<int> conjugated_bond_indices = ValidatedConjugatedBondIndices(working);
            metrics.conjugated_bond_count = static_cast<int>(conjugated_bond_indices.size());

            std::map<int, std::set<int>> conjugated_neighbors;
            std::set<int> conjugated_atom_indices;
            FOR_BONDS_OF_MOL(bond_iter, working)
            {
                OpenBabel::OBBond &bond = *bond_iter;
                if (conjugated_bond_indices.find(static_cast<int>(bond.GetIdx())) ==
                    conjugated_bond_indices.end())
                {
                    continue;
                }
                OpenBabel::OBAtom *begin_atom = bond.GetBeginAtom();
                OpenBabel::OBAtom *end_atom = bond.GetEndAtom();
                if (begin_atom == nullptr || end_atom == nullptr ||
                    begin_atom->GetAtomicNum() == 1 || end_atom->GetAtomicNum() == 1)
                {
                    continue;
                }
                const int begin_idx = static_cast<int>(begin_atom->GetIdx()) - 1;
                const int end_idx = static_cast<int>(end_atom->GetIdx()) - 1;
                conjugated_atom_indices.insert(begin_idx);
                conjugated_atom_indices.insert(end_idx);
                conjugated_neighbors[begin_idx].insert(end_idx);
                conjugated_neighbors[end_idx].insert(begin_idx);
            }
            metrics.conjugated_atom_count = static_cast<int>(conjugated_atom_indices.size());
            metrics.conjugated_atom_indices.assign(
                conjugated_atom_indices.begin(),
                conjugated_atom_indices.end());

            std::set<int> visited;
            for (int atom_idx : conjugated_atom_indices)
            {
                if (visited.find(atom_idx) != visited.end())
                {
                    continue;
                }
                std::vector<int> stack{atom_idx};
                int component_size = 0;
                while (!stack.empty())
                {
                    const int current_idx = stack.back();
                    stack.pop_back();
                    if (!visited.insert(current_idx).second)
                    {
                        continue;
                    }
                    ++component_size;
                    const auto neighbors_it = conjugated_neighbors.find(current_idx);
                    if (neighbors_it == conjugated_neighbors.end())
                    {
                        continue;
                    }
                    for (int neighbor_idx : neighbors_it->second)
                    {
                        if (visited.find(neighbor_idx) == visited.end())
                        {
                            stack.push_back(neighbor_idx);
                        }
                    }
                }
                metrics.max_conjugated_component_size =
                    std::max(metrics.max_conjugated_component_size, component_size);
            }
            return metrics;
        }
    }
}
