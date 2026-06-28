#include "molgr/utils/organic_topology.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include "molgr/compat/openbabel_iter.h"
#include <openbabel/ring.h>
#include <openbabel/typer.h>

#include <algorithm>
#include <cstdlib>
#include <map>
#include <set>
#include <vector>

namespace
{
    constexpr int kAromaticRingFormalChargeAbsRejectionThreshold = 4;

    int RingFormalChargeSum(OpenBabel::OBMol &mol, const OpenBabel::OBRing &ring)
    {
        int charge_sum = 0;
        for (int atom_idx : ring._path)
        {
            const OpenBabel::OBAtom *atom = mol.GetAtom(atom_idx);
            if (atom == nullptr)
            {
                continue;
            }
            charge_sum += atom->GetFormalCharge();
        }
        return charge_sum;
    }

    bool IsChargeAcceptedAromaticRing(OpenBabel::OBMol &mol, OpenBabel::OBRing &ring)
    {
        if (!ring.IsAromatic())
        {
            return false;
        }
        return std::abs(RingFormalChargeSum(mol, ring)) <
               kAromaticRingFormalChargeAbsRejectionThreshold;
    }

    bool AtomHasOddSpin(const OpenBabel::OBAtom &atom)
    {
        return atom.GetSpinMultiplicity() % 2 == 1;
    }

    double AromaticRingStabilityWeight(
        OpenBabel::OBMol &mol,
        const OpenBabel::OBRing &ring,
        const molgr::config::OrganicTopologyConfig &config)
    {
        std::vector<const OpenBabel::OBAtom *> heavy_atoms;
        for (int atom_idx : ring._path)
        {
            const OpenBabel::OBAtom *atom = mol.GetAtom(atom_idx);
            if (atom == nullptr || atom->GetAtomicNum() == 1)
            {
                continue;
            }
            heavy_atoms.push_back(atom);
        }
        if (heavy_atoms.empty())
        {
            return 0.0;
        }

        const int ring_size = static_cast<int>(heavy_atoms.size());
        int hetero_count = 0;
        int charge_count = 0;
        int radical_count = 0;
        for (const OpenBabel::OBAtom *atom : heavy_atoms)
        {
            if (atom->GetAtomicNum() != 6)
            {
                ++hetero_count;
            }
            if (atom->GetFormalCharge() != 0)
            {
                ++charge_count;
            }
            if (AtomHasOddSpin(*atom))
            {
                ++radical_count;
            }
        }

        if (ring_size == 6 && hetero_count == 0 && charge_count == 0 && radical_count == 0)
        {
            return config.aromatic_stability_benzene_score;
        }

        const double size_factor =
            ring_size == 6
                ? config.aromatic_stability_ring_size_6_factor
                : (ring_size == 5
                       ? config.aromatic_stability_ring_size_5_factor
                       : config.aromatic_stability_other_ring_size_factor);
        const double hetero_factor = std::max(
            config.aromatic_stability_min_hetero_factor,
            1.0 - config.aromatic_stability_hetero_atom_penalty *
                      static_cast<double>(hetero_count));
        const double charge_factor = std::max(
            config.aromatic_stability_min_charge_factor,
            1.0 - config.aromatic_stability_formal_charge_penalty *
                      static_cast<double>(charge_count));
        const double radical_factor = std::max(
            config.aromatic_stability_min_radical_factor,
            1.0 - config.aromatic_stability_radical_penalty *
                      static_cast<double>(radical_count));
        return std::min(
            config.aromatic_stability_other_ring_max_score,
            size_factor * hetero_factor * charge_factor * radical_factor);
    }

    bool IsMultipleLikeBond(const OpenBabel::OBBond &bond)
    {
        return bond.IsAromatic() || bond.GetBondOrder() >= 2;
    }

    bool IsHeavyAtomBond(const OpenBabel::OBBond &bond)
    {
        const OpenBabel::OBAtom *begin_atom = bond.GetBeginAtom();
        const OpenBabel::OBAtom *end_atom = bond.GetEndAtom();
        return begin_atom != nullptr && end_atom != nullptr && begin_atom->GetAtomicNum() != 1 &&
               end_atom->GetAtomicNum() != 1;
    }

    std::set<int> ValidatedConjugatedBondIndices(OpenBabel::OBMol &mol)
    {
        std::vector<OpenBabel::OBBond *> heavy_bonds;
        std::set<int> atom_has_adjacent_multiple_like_bond;
        std::set<int> atom_has_adjacent_alternating_single_bond;
        std::set<int> aromatic_bond_indices;
        std::set<int> multiple_like_bond_indices;

        FOR_BONDS_OF_MOL(bond_iter, mol)
        {
            OpenBabel::OBBond &bond = *bond_iter;
            if (!IsHeavyAtomBond(bond))
            {
                continue;
            }
            heavy_bonds.push_back(&bond);
            if (!IsMultipleLikeBond(bond))
            {
                continue;
            }
            OpenBabel::OBAtom *begin_atom = bond.GetBeginAtom();
            OpenBabel::OBAtom *end_atom = bond.GetEndAtom();
            atom_has_adjacent_multiple_like_bond.insert(static_cast<int>(begin_atom->GetIdx()));
            atom_has_adjacent_multiple_like_bond.insert(static_cast<int>(end_atom->GetIdx()));
            if (bond.IsAromatic())
            {
                aromatic_bond_indices.insert(static_cast<int>(bond.GetIdx()));
            }
            else
            {
                multiple_like_bond_indices.insert(static_cast<int>(bond.GetIdx()));
            }
        }

        std::set<int> conjugated_bond_indices;
        std::set<int> alternating_single_bond_indices;
        for (OpenBabel::OBBond *bond : heavy_bonds)
        {
            if (bond == nullptr || bond->IsAromatic() || bond->GetBondOrder() != 1)
            {
                continue;
            }
            OpenBabel::OBAtom *begin_atom = bond->GetBeginAtom();
            OpenBabel::OBAtom *end_atom = bond->GetEndAtom();
            if (begin_atom == nullptr || end_atom == nullptr)
            {
                continue;
            }
            if (atom_has_adjacent_multiple_like_bond.find(static_cast<int>(begin_atom->GetIdx())) ==
                    atom_has_adjacent_multiple_like_bond.end() ||
                atom_has_adjacent_multiple_like_bond.find(static_cast<int>(end_atom->GetIdx())) ==
                    atom_has_adjacent_multiple_like_bond.end())
            {
                continue;
            }
            alternating_single_bond_indices.insert(static_cast<int>(bond->GetIdx()));
            atom_has_adjacent_alternating_single_bond.insert(static_cast<int>(begin_atom->GetIdx()));
            atom_has_adjacent_alternating_single_bond.insert(static_cast<int>(end_atom->GetIdx()));
        }

        conjugated_bond_indices.insert(aromatic_bond_indices.begin(), aromatic_bond_indices.end());
        conjugated_bond_indices.insert(
            alternating_single_bond_indices.begin(),
            alternating_single_bond_indices.end());
        for (OpenBabel::OBBond *bond : heavy_bonds)
        {
            if (bond == nullptr ||
                multiple_like_bond_indices.find(static_cast<int>(bond->GetIdx())) ==
                    multiple_like_bond_indices.end())
            {
                continue;
            }
            OpenBabel::OBAtom *begin_atom = bond->GetBeginAtom();
            OpenBabel::OBAtom *end_atom = bond->GetEndAtom();
            if (begin_atom == nullptr || end_atom == nullptr)
            {
                continue;
            }
            if (atom_has_adjacent_alternating_single_bond.find(static_cast<int>(begin_atom->GetIdx())) !=
                    atom_has_adjacent_alternating_single_bond.end() ||
                atom_has_adjacent_alternating_single_bond.find(static_cast<int>(end_atom->GetIdx())) !=
                    atom_has_adjacent_alternating_single_bond.end())
            {
                conjugated_bond_indices.insert(static_cast<int>(bond->GetIdx()));
            }
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
            const OpenBabel::OBMol *parent = const_cast<OpenBabel::OBBond &>(bond).GetParent();
            if (parent == nullptr)
            {
                return false;
            }
            OpenBabel::OBMol working(*parent);
            working.FindRingAtomsAndBonds();
            working.SetAromaticPerceived(false);
            OpenBabel::OBAromaticTyper().AssignAromaticFlags(working);
            const std::set<int> conjugated_bond_indices = ValidatedConjugatedBondIndices(working);
            return conjugated_bond_indices.find(static_cast<int>(bond.GetIdx())) !=
                   conjugated_bond_indices.end();
        }

        OrganicTopologyMetrics ComputeOrganicTopologyMetrics(
            const OpenBabel::OBMol &mol,
            const molgr::config::OrganicTopologyConfig &config)
        {
            OpenBabel::OBMol working(mol);
            working.FindRingAtomsAndBonds();
            working.SetAromaticPerceived(false);
            OpenBabel::OBAromaticTyper().AssignAromaticFlags(working);

            OrganicTopologyMetrics metrics;
            std::set<int> aromatic_atom_indices;
            for (OpenBabel::OBRing *ring : working.GetSSSR())
            {
                if (ring == nullptr || !IsChargeAcceptedAromaticRing(working, *ring))
                {
                    continue;
                }
                ++metrics.aromatic_ring_count;
                metrics.aromatic_stability_score += AromaticRingStabilityWeight(
                    working,
                    *ring,
                    config);
                for (int atom_idx : ring->_path)
                {
                    const OpenBabel::OBAtom *atom = working.GetAtom(atom_idx);
                    if (atom == nullptr || atom->GetAtomicNum() == 1)
                    {
                        continue;
                    }
                    aromatic_atom_indices.insert(atom_idx);
                }
            }
            metrics.aromatic_atom_count = static_cast<int>(aromatic_atom_indices.size());

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
