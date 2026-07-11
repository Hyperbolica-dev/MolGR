#include "molgr/utils/organic_topology.h"

#include "molgr/vendor/openbabel_threading.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include "molgr/compat/openbabel_iter.h"
#include <openbabel/ring.h>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <vector>

namespace
{
    constexpr int kAromaticRingFormalChargeAbsRejectionThreshold = 4;

    struct IntFlagSet
    {
        std::vector<unsigned char> present;
        std::vector<int> values;

        explicit IntFlagSet(std::size_t initial_size = 0)
            : present(initial_size, 0)
        {
        }

        bool Contains(int index) const
        {
            return index >= 0 &&
                   static_cast<std::size_t>(index) < present.size() &&
                   present[static_cast<std::size_t>(index)] != 0;
        }

        bool Add(int index)
        {
            if (index < 0)
            {
                return false;
            }
            const auto offset = static_cast<std::size_t>(index);
            if (offset >= present.size())
            {
                present.resize(offset + 1, 0);
            }
            if (present[offset] != 0)
            {
                return false;
            }
            present[offset] = 1;
            values.push_back(index);
            return true;
        }

        std::size_t Size() const
        {
            return values.size();
        }
    };

    struct HeavyBondRef
    {
        OpenBabel::OBBond *bond = nullptr;
        OpenBabel::OBAtom *begin_atom = nullptr;
        OpenBabel::OBAtom *end_atom = nullptr;
        int bond_idx = -1;
    };

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
        if (!molgr::vendor::openbabel_threading::RingIsAromatic(ring))
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

    int AdditionalAtomPiElectrons(
        const OpenBabel::OBAtom &atom,
        bool incident_to_ring_multiple_bond)
    {
        if (incident_to_ring_multiple_bond || atom.GetAtomicNum() == 1)
        {
            return 0;
        }
        if (atom.GetFormalCharge() < 0)
        {
            return 2;
        }
        if (atom.GetAtomicNum() != 6)
        {
            return 2;
        }
        if (AtomHasOddSpin(atom))
        {
            return 1;
        }
        return 0;
    }

    int RingPiElectronCount(OpenBabel::OBMol &mol, const OpenBabel::OBRing &ring)
    {
        if (ring._path.size() < 3)
        {
            return -1;
        }

        int pi_electron_count = 0;
        IntFlagSet atoms_incident_to_ring_multiple_bond(
            static_cast<std::size_t>(mol.NumAtoms()) + 1);
        for (std::size_t offset = 0; offset < ring._path.size(); ++offset)
        {
            const int begin_idx = ring._path[offset];
            const int end_idx = ring._path[(offset + 1) % ring._path.size()];
            OpenBabel::OBBond *bond = mol.GetBond(begin_idx, end_idx);
            if (bond == nullptr)
            {
                return -1;
            }
            if (bond->GetBondOrder() >= 2)
            {
                pi_electron_count += 2;
                atoms_incident_to_ring_multiple_bond.Add(begin_idx);
                atoms_incident_to_ring_multiple_bond.Add(end_idx);
            }
        }

        for (int atom_idx : ring._path)
        {
            const OpenBabel::OBAtom *atom = mol.GetAtom(atom_idx);
            if (atom == nullptr)
            {
                continue;
            }
            pi_electron_count += AdditionalAtomPiElectrons(
                *atom,
                atoms_incident_to_ring_multiple_bond.Contains(atom_idx));
        }
        return pi_electron_count;
    }

    bool HasHuckelPiElectronCount(int pi_electron_count)
    {
        return pi_electron_count >= 2 && (pi_electron_count - 2) % 4 == 0;
    }

    bool IsHuckelAcceptedAromaticRing(OpenBabel::OBMol &mol, OpenBabel::OBRing &ring)
    {
        return IsChargeAcceptedAromaticRing(mol, ring) &&
               HasHuckelPiElectronCount(RingPiElectronCount(mol, ring));
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
        OpenBabel::OBBond &mutable_bond = const_cast<OpenBabel::OBBond &>(bond);
        return molgr::vendor::openbabel_threading::BondIsAromatic(mutable_bond) ||
               bond.GetBondOrder() >= 2;
    }

    bool IsHeavyAtomBond(const OpenBabel::OBBond &bond)
    {
        const OpenBabel::OBAtom *begin_atom = bond.GetBeginAtom();
        const OpenBabel::OBAtom *end_atom = bond.GetEndAtom();
        return begin_atom != nullptr && end_atom != nullptr && begin_atom->GetAtomicNum() != 1 &&
               end_atom->GetAtomicNum() != 1;
    }

    IntFlagSet ValidatedConjugatedBondIndices(OpenBabel::OBMol &mol)
    {
        std::vector<HeavyBondRef> heavy_bonds;
        heavy_bonds.reserve(static_cast<std::size_t>(mol.NumBonds()));
        IntFlagSet atom_has_adjacent_multiple_like_bond(
            static_cast<std::size_t>(mol.NumAtoms()) + 1);
        IntFlagSet atom_has_adjacent_alternating_single_bond(
            static_cast<std::size_t>(mol.NumAtoms()) + 1);
        IntFlagSet aromatic_bond_indices(
            static_cast<std::size_t>(mol.NumBonds()) + 1);
        IntFlagSet multiple_like_bond_indices(
            static_cast<std::size_t>(mol.NumBonds()) + 1);

        FOR_BONDS_OF_MOL(bond_iter, mol)
        {
            OpenBabel::OBBond &bond = *bond_iter;
            if (!IsHeavyAtomBond(bond))
            {
                continue;
            }
            OpenBabel::OBAtom *begin_atom = bond.GetBeginAtom();
            OpenBabel::OBAtom *end_atom = bond.GetEndAtom();
            const int bond_idx = static_cast<int>(bond.GetIdx());
            heavy_bonds.push_back(HeavyBondRef{
                &bond,
                begin_atom,
                end_atom,
                bond_idx,
            });
            if (!IsMultipleLikeBond(bond))
            {
                continue;
            }
            atom_has_adjacent_multiple_like_bond.Add(static_cast<int>(begin_atom->GetIdx()));
            atom_has_adjacent_multiple_like_bond.Add(static_cast<int>(end_atom->GetIdx()));
            if (molgr::vendor::openbabel_threading::BondIsAromatic(bond))
            {
                aromatic_bond_indices.Add(bond_idx);
            }
            else
            {
                multiple_like_bond_indices.Add(bond_idx);
            }
        }

        IntFlagSet conjugated_bond_indices(
            static_cast<std::size_t>(mol.NumBonds()) + 1);
        IntFlagSet alternating_single_bond_indices(
            static_cast<std::size_t>(mol.NumBonds()) + 1);
        for (const HeavyBondRef &heavy_bond : heavy_bonds)
        {
            OpenBabel::OBBond *bond = heavy_bond.bond;
            OpenBabel::OBAtom *begin_atom = heavy_bond.begin_atom;
            OpenBabel::OBAtom *end_atom = heavy_bond.end_atom;
            if (bond == nullptr || begin_atom == nullptr || end_atom == nullptr ||
                molgr::vendor::openbabel_threading::BondIsAromatic(*bond) ||
                bond->GetBondOrder() != 1)
            {
                continue;
            }
            if (!atom_has_adjacent_multiple_like_bond.Contains(static_cast<int>(begin_atom->GetIdx())) ||
                !atom_has_adjacent_multiple_like_bond.Contains(static_cast<int>(end_atom->GetIdx())))
            {
                continue;
            }
            alternating_single_bond_indices.Add(heavy_bond.bond_idx);
            atom_has_adjacent_alternating_single_bond.Add(static_cast<int>(begin_atom->GetIdx()));
            atom_has_adjacent_alternating_single_bond.Add(static_cast<int>(end_atom->GetIdx()));
        }

        for (int bond_idx : aromatic_bond_indices.values)
        {
            conjugated_bond_indices.Add(bond_idx);
        }
        for (int bond_idx : alternating_single_bond_indices.values)
        {
            conjugated_bond_indices.Add(bond_idx);
        }
        for (const HeavyBondRef &heavy_bond : heavy_bonds)
        {
            OpenBabel::OBBond *bond = heavy_bond.bond;
            OpenBabel::OBAtom *begin_atom = heavy_bond.begin_atom;
            OpenBabel::OBAtom *end_atom = heavy_bond.end_atom;
            if (bond == nullptr || begin_atom == nullptr || end_atom == nullptr ||
                !multiple_like_bond_indices.Contains(heavy_bond.bond_idx))
            {
                continue;
            }
            if (atom_has_adjacent_alternating_single_bond.Contains(static_cast<int>(begin_atom->GetIdx())) ||
                atom_has_adjacent_alternating_single_bond.Contains(static_cast<int>(end_atom->GetIdx())))
            {
                conjugated_bond_indices.Add(heavy_bond.bond_idx);
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
            molgr::vendor::openbabel_threading::ResetAndAssignAromaticFlags(working);
            const IntFlagSet conjugated_bond_indices = ValidatedConjugatedBondIndices(working);
            return conjugated_bond_indices.Contains(static_cast<int>(bond.GetIdx()));
        }

        OrganicTopologyMetrics ComputeOrganicTopologyMetrics(
            const OpenBabel::OBMol &mol,
            const molgr::config::OrganicTopologyConfig &config)
        {
            OpenBabel::OBMol working(mol);
            molgr::vendor::openbabel_threading::ResetAndAssignAromaticFlags(working);

            OrganicTopologyMetrics metrics;
            IntFlagSet aromatic_atom_indices(
                static_cast<std::size_t>(working.NumAtoms()) + 1);
            FOR_RINGS_OF_MOL(ring_iter, working)
            {
                OpenBabel::OBRing *ring = &(*ring_iter);
                if (ring == nullptr || !IsHuckelAcceptedAromaticRing(working, *ring))
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
                    aromatic_atom_indices.Add(atom_idx);
                }
            }
            metrics.aromatic_atom_count = static_cast<int>(aromatic_atom_indices.Size());

            const IntFlagSet conjugated_bond_indices = ValidatedConjugatedBondIndices(working);
            metrics.conjugated_bond_count = static_cast<int>(conjugated_bond_indices.Size());

            std::vector<std::vector<int>> conjugated_neighbors(
                static_cast<std::size_t>(working.NumAtoms()));
            IntFlagSet conjugated_atom_indices(
                static_cast<std::size_t>(working.NumAtoms()));
            FOR_BONDS_OF_MOL(bond_iter, working)
            {
                OpenBabel::OBBond &bond = *bond_iter;
                if (!conjugated_bond_indices.Contains(static_cast<int>(bond.GetIdx())))
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
                if (begin_idx < 0 || end_idx < 0)
                {
                    continue;
                }
                conjugated_atom_indices.Add(begin_idx);
                conjugated_atom_indices.Add(end_idx);
                const auto begin_offset = static_cast<std::size_t>(begin_idx);
                const auto end_offset = static_cast<std::size_t>(end_idx);
                if (begin_offset >= conjugated_neighbors.size() ||
                    end_offset >= conjugated_neighbors.size())
                {
                    continue;
                }
                conjugated_neighbors[begin_offset].push_back(end_idx);
                conjugated_neighbors[end_offset].push_back(begin_idx);
            }
            metrics.conjugated_atom_count = static_cast<int>(conjugated_atom_indices.Size());
            metrics.conjugated_atom_indices = conjugated_atom_indices.values;
            std::sort(metrics.conjugated_atom_indices.begin(), metrics.conjugated_atom_indices.end());

            IntFlagSet visited(static_cast<std::size_t>(working.NumAtoms()));
            for (int atom_idx : metrics.conjugated_atom_indices)
            {
                if (visited.Contains(atom_idx))
                {
                    continue;
                }
                std::vector<int> stack{atom_idx};
                int component_size = 0;
                while (!stack.empty())
                {
                    const int current_idx = stack.back();
                    stack.pop_back();
                    if (!visited.Add(current_idx))
                    {
                        continue;
                    }
                    ++component_size;
                    if (current_idx < 0 ||
                        static_cast<std::size_t>(current_idx) >= conjugated_neighbors.size())
                    {
                        continue;
                    }
                    for (int neighbor_idx : conjugated_neighbors[static_cast<std::size_t>(current_idx)])
                    {
                        if (!visited.Contains(neighbor_idx))
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
