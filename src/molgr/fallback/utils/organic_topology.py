from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Set, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.config import OrganicTopologyConfig, resolve_config


_AROMATIC_RING_FORMAL_CHARGE_ABS_REJECTION_THRESHOLD = 4


@dataclass(frozen=True)
class OrganicTopologyMetrics:
    aromatic_atom_count: int
    aromatic_ring_count: int
    aromatic_stability_score: float
    conjugated_atom_count: int
    conjugated_bond_count: int
    max_conjugated_component_size: int
    conjugated_atom_indices: Tuple[int, ...]


def _prepare_topology_working_molecule(omol: pybel.Molecule) -> pybel.Molecule:
    working = omol.clone
    working_obmol = cast(ob.OBMol, working.OBMol)
    working_obmol.FindRingAtomsAndBonds()
    working_obmol.SetAromaticPerceived(False)
    ob.OBAromaticTyper().AssignAromaticFlags(working_obmol)
    return working


def _ring_atom_indices(ring: ob.OBRing) -> Tuple[int, ...]:
    return tuple(int(atom_idx) for atom_idx in getattr(ring, "_path", ()))


def _ring_formal_charge_sum(obmol: ob.OBMol, ring: ob.OBRing) -> int:
    charge_sum = 0
    for atom_idx in _ring_atom_indices(ring):
        atom = obmol.GetAtom(atom_idx)
        if atom is None:
            continue
        charge_sum += cast(int, atom.GetFormalCharge())
    return charge_sum


def _is_charge_accepted_aromatic_ring(obmol: ob.OBMol, ring: ob.OBRing) -> bool:
    if not ring.IsAromatic():
        return False
    return (
        abs(_ring_formal_charge_sum(obmol, ring))
        < _AROMATIC_RING_FORMAL_CHARGE_ABS_REJECTION_THRESHOLD
    )


def _aromatic_ring_stability_weight(
    obmol: ob.OBMol,
    ring_atom_indices: Iterable[int],
    config: OrganicTopologyConfig,
) -> float:
    atoms = [obmol.GetAtom(atom_idx) for atom_idx in ring_atom_indices]
    heavy_atoms = [atom for atom in atoms if atom is not None and atom.GetAtomicNum() != 1]
    if not heavy_atoms:
        return 0.0

    ring_size = len(heavy_atoms)
    hetero_count = sum(1 for atom in heavy_atoms if int(atom.GetAtomicNum()) != 6)
    charge_count = sum(1 for atom in heavy_atoms if int(atom.GetFormalCharge()) != 0)
    radical_count = sum(1 for atom in heavy_atoms if _atom_has_odd_spin(atom))

    if ring_size == 6 and hetero_count == 0 and charge_count == 0 and radical_count == 0:
        return config.aromatic_stability_benzene_score

    size_factor = (
        config.aromatic_stability_ring_size_6_factor
        if ring_size == 6
        else (
            config.aromatic_stability_ring_size_5_factor
            if ring_size == 5
            else config.aromatic_stability_other_ring_size_factor
        )
    )
    hetero_factor = max(
        config.aromatic_stability_min_hetero_factor,
        1.0 - config.aromatic_stability_hetero_atom_penalty * hetero_count,
    )
    charge_factor = max(
        config.aromatic_stability_min_charge_factor,
        1.0 - config.aromatic_stability_formal_charge_penalty * charge_count,
    )
    radical_factor = max(
        config.aromatic_stability_min_radical_factor,
        1.0 - config.aromatic_stability_radical_penalty * radical_count,
    )
    return min(
        config.aromatic_stability_other_ring_max_score,
        size_factor * hetero_factor * charge_factor * radical_factor,
    )


def _atom_has_odd_spin(atom: ob.OBAtom) -> bool:
    return cast(int, atom.GetSpinMultiplicity()) % 2 == 1


def _atom_has_pi_feature(atom: ob.OBAtom, *, exclude_bond_idx: int | None = None) -> bool:
    if atom.IsAromatic() or cast(int, atom.GetFormalCharge()) != 0 or _atom_has_odd_spin(atom):
        return True
    for bond_iter in ob.OBAtomBondIter(atom):
        bond = cast(ob.OBBond, bond_iter)
        if exclude_bond_idx is not None and cast(int, bond.GetIdx()) == exclude_bond_idx:
            continue
        other_atom = cast(ob.OBAtom, bond.GetNbrAtom(atom))
        if other_atom.GetAtomicNum() == 1:
            continue
        if bond.IsAromatic() or cast(int, bond.GetBondOrder()) >= 2:
            return True
    return False


def _validated_conjugated_bond_indices(omol: pybel.Molecule) -> Set[int]:
    obmol = cast(ob.OBMol, omol.OBMol)
    conjugated_bond_indices: Set[int] = set()
    for bond_iter in ob.OBMolBondIter(obmol):
        bond = cast(ob.OBBond, bond_iter)
        begin_atom = cast(ob.OBAtom, bond.GetBeginAtom())
        end_atom = cast(ob.OBAtom, bond.GetEndAtom())
        if begin_atom.GetAtomicNum() == 1 or end_atom.GetAtomicNum() == 1:
            continue
        bond_idx = cast(int, bond.GetIdx())
        bond_order = cast(int, bond.GetBondOrder())
        if bond.IsAromatic() or bond_order >= 2:
            conjugated_bond_indices.add(bond_idx)
            continue
        if bond_order != 1:
            continue
        if not _atom_has_pi_feature(begin_atom, exclude_bond_idx=bond_idx):
            continue
        if not _atom_has_pi_feature(end_atom, exclude_bond_idx=bond_idx):
            continue
        conjugated_bond_indices.add(bond_idx)
    return conjugated_bond_indices


def is_conjugated_bond(bond: ob.OBBond) -> bool:
    parent = bond.GetParent()
    if parent is None:
        return False
    parent_mol = _prepare_topology_working_molecule(pybel.Molecule(parent))
    return cast(int, bond.GetIdx()) in _validated_conjugated_bond_indices(parent_mol)


def compute_organic_topology_metrics(
    omol: pybel.Molecule,
    config: OrganicTopologyConfig | None = None,
) -> OrganicTopologyMetrics:
    try:
        topology_config = resolve_config().organic_topology if config is None else config
        working_omol = _prepare_topology_working_molecule(omol)
        obmol = cast(ob.OBMol, working_omol.OBMol)

        aromatic_ring_count = 0
        aromatic_stability_score = 0.0
        aromatic_atom_indices: Set[int] = set()
        for ring_iter in ob.OBMolRingIter(obmol):
            ring = cast(ob.OBRing, ring_iter)
            if not _is_charge_accepted_aromatic_ring(obmol, ring):
                continue
            aromatic_ring_count += 1
            ring_atom_indices = _ring_atom_indices(ring)
            aromatic_stability_score += _aromatic_ring_stability_weight(
                obmol,
                ring_atom_indices,
                topology_config,
            )
            for atom_idx in ring_atom_indices:
                atom = obmol.GetAtom(atom_idx)
                if atom is None or atom.GetAtomicNum() == 1:
                    continue
                aromatic_atom_indices.add(atom_idx)
        aromatic_atom_count = len(aromatic_atom_indices)

        conjugated_neighbors: Dict[int, Set[int]] = {}
        conjugated_atom_indices: Set[int] = set()
        conjugated_bond_indices = _validated_conjugated_bond_indices(working_omol)
        for bond_iter in ob.OBMolBondIter(obmol):
            bond = cast(ob.OBBond, bond_iter)
            begin_atom = cast(ob.OBAtom, bond.GetBeginAtom())
            end_atom = cast(ob.OBAtom, bond.GetEndAtom())
            if begin_atom.GetAtomicNum() == 1 or end_atom.GetAtomicNum() == 1:
                continue
            if cast(int, bond.GetIdx()) not in conjugated_bond_indices:
                continue
            begin_idx = begin_atom.GetIdx() - 1
            end_idx = end_atom.GetIdx() - 1
            conjugated_atom_indices.add(begin_idx)
            conjugated_atom_indices.add(end_idx)
            conjugated_neighbors.setdefault(begin_idx, set()).add(end_idx)
            conjugated_neighbors.setdefault(end_idx, set()).add(begin_idx)

        max_conjugated_component_size = 0
        visited: Set[int] = set()
        for atom_idx in conjugated_atom_indices:
            if atom_idx in visited:
                continue
            stack = [atom_idx]
            component_size = 0
            while stack:
                current_idx = stack.pop()
                if current_idx in visited:
                    continue
                visited.add(current_idx)
                component_size += 1
                stack.extend(
                    neighbor_idx
                    for neighbor_idx in conjugated_neighbors.get(current_idx, set())
                    if neighbor_idx not in visited
                )
            max_conjugated_component_size = max(max_conjugated_component_size, component_size)

        return OrganicTopologyMetrics(
            aromatic_atom_count=aromatic_atom_count,
            aromatic_ring_count=aromatic_ring_count,
            aromatic_stability_score=aromatic_stability_score,
            conjugated_atom_count=len(conjugated_atom_indices),
            conjugated_bond_count=len(conjugated_bond_indices),
            max_conjugated_component_size=max_conjugated_component_size,
            conjugated_atom_indices=tuple(sorted(conjugated_atom_indices)),
        )
    except Exception:  # noqa: BLE001
        return OrganicTopologyMetrics(
            aromatic_atom_count=0,
            aromatic_ring_count=0,
            aromatic_stability_score=0.0,
            conjugated_atom_count=0,
            conjugated_bond_count=0,
            max_conjugated_component_size=0,
            conjugated_atom_indices=(),
        )


__all__ = ["OrganicTopologyMetrics", "compute_organic_topology_metrics", "is_conjugated_bond"]
