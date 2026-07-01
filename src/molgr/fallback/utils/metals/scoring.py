"""Scoring and selection utilities for metal-aware fallback reconstruction."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import astuple, dataclass
from typing import List, Optional, Set, Tuple, cast

import numpy as np
from openbabel import openbabel as ob
from openbabel import pybel

from molgr.config import CONFIG, MetalScoringConfig, MolGRConfig
from molgr.fallback.state import (
    MetalCandidateState,
    MetalCandidateStateMachine,
    ReconstructionState,
)
from molgr.fallback.utils import consts, dataclasses
from molgr.fallback.utils.organic_topology import compute_organic_topology_metrics
from molgr.utils.coordination_visibility import (
    CoordinationBlockerArrays,
    Point3D,
    coordination_visibility_mask,
    empty_coordination_blocker_arrays,
)


_INNER_VISIBLE_DIRADICAL_EXEMPT_ATOMIC_NUMS = frozenset({15, 16, 17, 35, 53})
_NEGATIVE_METAL_DISCORDANCE_PENALTY = 0.5


@dataclass(frozen=True)
class _OrganicElectronicStateMetrics:
    aromatic_atom_count: int
    aromatic_ring_count: int
    aromatic_stability_score: float
    conjugated_atom_count: int
    conjugated_bond_count: int
    max_conjugated_component_size: int
    radical_localization_penalty: float
    charge_localization_penalty: float


def _distance_to_metal(
    atom: ob.OBAtom,
    metal_state: dataclasses.MetalAtomPosition,
) -> float:
    dx = float(atom.GetX()) - metal_state.position_x
    dy = float(atom.GetY()) - metal_state.position_y
    dz = float(atom.GetZ()) - metal_state.position_z
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def _charge_sign(charge: int) -> int:
    if charge > 0:
        return 1
    if charge < 0:
        return -1
    return 0


def _is_inner_visible_diradical_discordance_atom(atom: ob.OBAtom) -> bool:
    return (
        int(atom.GetSpinMultiplicity()) >= 2
        and int(atom.GetAtomicNum()) not in _INNER_VISIBLE_DIRADICAL_EXEMPT_ATOMIC_NUMS
    )


def _nearest_nonzero_metal_charge_sign_to_bond(
    begin_atom: ob.OBAtom,
    end_atom: ob.OBAtom,
    metal_states: Sequence[dataclasses.MetalAtomPosition],
) -> int:
    midpoint_x = (float(begin_atom.GetX()) + float(end_atom.GetX())) * 0.5
    midpoint_y = (float(begin_atom.GetY()) + float(end_atom.GetY())) * 0.5
    midpoint_z = (float(begin_atom.GetZ()) + float(end_atom.GetZ())) * 0.5
    best_distance_sq = float("inf")
    best_charge_sign = 0
    for metal_state in metal_states:
        metal_charge_sign = _charge_sign(int(metal_state.valence))
        if metal_charge_sign == 0:
            continue
        dx = midpoint_x - metal_state.position_x
        dy = midpoint_y - metal_state.position_y
        dz = midpoint_z - metal_state.position_z
        distance_sq = dx * dx + dy * dy + dz * dz
        if distance_sq < best_distance_sq:
            best_distance_sq = distance_sq
            best_charge_sign = metal_charge_sign
    return best_charge_sign


def _atom_coordinates(atom: ob.OBAtom) -> Point3D:
    return float(atom.GetX()), float(atom.GetY()), float(atom.GetZ())


def _build_coordination_blocker_arrays(
    obmol: ob.OBMol,
    *,
    metal_scoring_config: MetalScoringConfig,
) -> CoordinationBlockerArrays:
    blocker_indices: list[int] = []
    blocker_coordinates: list[Point3D] = []
    blocker_radii: list[float] = []
    for blocker_iter in ob.OBMolAtomIter(obmol):
        blocker = cast(ob.OBAtom, blocker_iter)
        blocker_radius = (
            metal_scoring_config.metal_access_radius_scale
            * float(ob.GetCovalentRad(int(blocker.GetAtomicNum())))
            + metal_scoring_config.metal_access_clearance_angstrom
        )
        if blocker_radius <= 0.0:
            continue
        blocker_indices.append(int(blocker.GetIdx()))
        blocker_coordinates.append(_atom_coordinates(blocker))
        blocker_radii.append(blocker_radius)
    if not blocker_indices:
        return empty_coordination_blocker_arrays()
    return (
        np.asarray(blocker_indices, dtype=np.int64),
        np.asarray(blocker_coordinates, dtype=np.float64),
        np.asarray(blocker_radii, dtype=np.float64),
    )


def _coordination_radius_cutoff(
    atom: ob.OBAtom,
    metal_state: dataclasses.MetalAtomPosition,
    *,
    metal_scoring_config: MetalScoringConfig,
) -> float:
    atom_radius = float(ob.GetCovalentRad(int(atom.GetAtomicNum())))
    metal_radius = float(ob.GetCovalentRad(int(metal_state.element_idx)))
    return (
        metal_scoring_config.metal_access_radius_scale * (atom_radius + metal_radius)
        + metal_scoring_config.metal_coordination_extra_tolerance_angstrom
    )


def _is_inner_sphere_atom(
    atom: ob.OBAtom,
    metal_state: dataclasses.MetalAtomPosition,
    *,
    metal_scoring_config: MetalScoringConfig,
) -> bool:
    distance = _distance_to_metal(atom, metal_state)
    return (
        0.0
        < distance
        <= _coordination_radius_cutoff(
            atom,
            metal_state,
            metal_scoring_config=metal_scoring_config,
        )
    )


def _non_metal_atom_entries(obmol: ob.OBMol) -> tuple[ob.OBAtom, ...]:
    atoms: list[ob.OBAtom] = []
    for atom_iter in ob.OBMolAtomIter(obmol):
        atom = cast(ob.OBAtom, atom_iter)
        if atom.IsMetal():
            continue
        atoms.append(atom)
    return tuple(atoms)


def _inner_visible_atoms_to_metal(
    atoms: Sequence[ob.OBAtom],
    metal_state: dataclasses.MetalAtomPosition,
    *,
    metal_scoring_config: MetalScoringConfig,
    blocker_arrays: CoordinationBlockerArrays,
) -> tuple[ob.OBAtom, ...]:
    if not atoms:
        return ()

    atom_indices = np.asarray([int(atom.GetIdx()) for atom in atoms], dtype=np.int64)
    atom_coordinates = np.asarray([_atom_coordinates(atom) for atom in atoms], dtype=np.float64)
    metal_coordinates = np.asarray(
        (
            float(metal_state.position_x),
            float(metal_state.position_y),
            float(metal_state.position_z),
        ),
        dtype=np.float64,
    )
    atom_radii = np.asarray(
        [float(ob.GetCovalentRad(int(atom.GetAtomicNum()))) for atom in atoms],
        dtype=np.float64,
    )
    metal_radius = float(ob.GetCovalentRad(int(metal_state.element_idx)))
    cutoff = float(metal_scoring_config.metal_access_radius_scale) * (
        atom_radii + metal_radius
    ) + float(metal_scoring_config.metal_coordination_extra_tolerance_angstrom)
    delta = atom_coordinates - metal_coordinates
    distance_sq = np.einsum("ij,ij->i", delta, delta)
    inner_mask = (distance_sq > 0.0) & (distance_sq <= cutoff * cutoff)
    if not np.any(inner_mask):
        return ()

    inner_indices = atom_indices[inner_mask]
    inner_coordinates = atom_coordinates[inner_mask]
    visible_mask = coordination_visibility_mask(
        inner_indices,
        inner_coordinates,
        metal_coordinates,
        blocker_arrays,
    )
    inner_atoms = [atom for atom, is_inner in zip(atoms, inner_mask) if is_inner]
    return tuple(atom for atom, is_visible in zip(inner_atoms, visible_mask) if is_visible)


def _bond_order(bond: ob.OBBond) -> int:
    if bond.IsAromatic():
        return 2
    return cast(int, bond.GetBondOrder())


def _other_bond_atom(bond: ob.OBBond, atom: ob.OBAtom) -> ob.OBAtom:
    return (
        cast(ob.OBAtom, bond.GetEndAtom())
        if cast(ob.OBAtom, bond.GetBeginAtom()).GetIdx() == atom.GetIdx()
        else cast(ob.OBAtom, bond.GetBeginAtom())
    )


def _bond_between_atoms(lhs: ob.OBAtom, rhs: ob.OBAtom) -> Optional[ob.OBBond]:
    rhs_idx = rhs.GetIdx()
    for bond_iter in ob.OBAtomBondIter(lhs):
        bond = cast(ob.OBBond, bond_iter)
        if _other_bond_atom(bond, lhs).GetIdx() == rhs_idx:
            return bond
    return None


def _has_conjugated_bridge_between_charged_carbons(
    begin_atom: ob.OBAtom,
    end_atom: ob.OBAtom,
) -> bool:
    for begin_bond_iter in ob.OBAtomBondIter(begin_atom):
        begin_bond = cast(ob.OBBond, begin_bond_iter)
        if _bond_order(begin_bond) != 1:
            continue
        begin_neighbor = _other_bond_atom(begin_bond, begin_atom)
        if begin_neighbor.GetAtomicNum() == 1:
            continue
        for middle_bond_iter in ob.OBAtomBondIter(begin_neighbor):
            middle_bond = cast(ob.OBBond, middle_bond_iter)
            if (
                int(middle_bond.GetIdx()) == int(begin_bond.GetIdx())
                or _bond_order(middle_bond) != 2
            ):
                continue
            middle_neighbor = _other_bond_atom(middle_bond, begin_neighbor)
            if middle_neighbor.GetIdx() == begin_atom.GetIdx():
                continue
            bridge_end_bond = _bond_between_atoms(middle_neighbor, end_atom)
            if bridge_end_bond is None:
                continue
            if _bond_order(bridge_end_bond) == 1:
                return True
    return False


def _inner_visible_conjugated_charged_carbon_pair_count(
    obmol: ob.OBMol,
    visible_inner_atom_indices: Set[int],
) -> int:
    charged_carbons: List[ob.OBAtom] = []
    for atom_iter in ob.OBMolAtomIter(obmol):
        atom = cast(ob.OBAtom, atom_iter)
        if int(atom.GetIdx()) not in visible_inner_atom_indices:
            continue
        if atom.GetAtomicNum() != 6 or int(atom.GetFormalCharge()) == 0:
            continue
        charged_carbons.append(atom)

    count = 0
    for index, begin_atom in enumerate(charged_carbons):
        for end_atom in charged_carbons[index + 1 :]:
            if _charge_sign(int(begin_atom.GetFormalCharge())) != _charge_sign(
                int(end_atom.GetFormalCharge())
            ):
                continue
            if _bond_between_atoms(begin_atom, end_atom) is not None:
                continue
            if _has_conjugated_bridge_between_charged_carbons(begin_atom, end_atom):
                count += 1
    return count


def _has_outer_sphere_proton(
    obmol: ob.OBMol,
    metal_states: Sequence[dataclasses.MetalAtomPosition],
    *,
    metal_scoring_config: MetalScoringConfig,
) -> bool:
    for atom_iter in ob.OBMolAtomIter(obmol):
        atom = cast(ob.OBAtom, atom_iter)
        if atom.IsMetal() or int(atom.GetAtomicNum()) != 1 or int(atom.GetFormalCharge()) <= 0:
            continue
        if any(
            _is_inner_sphere_atom(
                atom,
                metal_state,
                metal_scoring_config=metal_scoring_config,
            )
            for metal_state in metal_states
        ):
            continue
        return True
    return False


def _negative_metal_discordance_count(
    obmol: ob.OBMol,
    metal_states: Sequence[dataclasses.MetalAtomPosition],
    *,
    metal_scoring_config: MetalScoringConfig,
) -> Tuple[int, bool, bool]:
    negative_metal_count = sum(
        abs(int(metal_state.valence))
        for metal_state in metal_states
        if int(metal_state.valence) < 0
    )
    if negative_metal_count == 0:
        return 0, False, False

    has_outer_sphere_proton = _has_outer_sphere_proton(
        obmol,
        metal_states,
        metal_scoring_config=metal_scoring_config,
    )
    has_positive_metal_counterion = any(
        int(metal_state.valence) > 0 for metal_state in metal_states
    )
    if has_outer_sphere_proton or has_positive_metal_counterion:
        return 0, has_outer_sphere_proton, has_positive_metal_counterion
    return negative_metal_count, False, False


def _zero_valent_metals_with_organic_cation_count(
    obmol: ob.OBMol,
    metal_states: Sequence[dataclasses.MetalAtomPosition],
    *,
    total_charge: int,
) -> int:
    if total_charge > 0:
        return 0
    if not metal_states or any(int(metal_state.valence) != 0 for metal_state in metal_states):
        return 0
    for atom_iter in ob.OBMolAtomIter(obmol):
        atom = cast(ob.OBAtom, atom_iter)
        if atom.IsMetal() or int(atom.GetFormalCharge()) <= 0:
            continue
        if int(atom.GetAtomicNum()) > 0 and not _is_locally_zwitterionic_organic_cation(atom):
            return 1
    return 0


def _is_unsaturated_organic_cation(atom: ob.OBAtom) -> bool:
    if atom.IsMetal() or int(atom.GetFormalCharge()) <= 0:
        return False

    total_degree = int(atom.GetTotalDegree())
    total_valence = int(atom.GetTotalValence())
    typical_valence = int(
        ob.GetTypicalValence(
            int(atom.GetAtomicNum()),
            total_valence,
            int(atom.GetFormalCharge()),
        )
    )
    return total_degree < total_valence or total_valence < typical_valence


def _has_adjacent_formal_charge_cancellation(atom: ob.OBAtom) -> bool:
    formal_charge = int(atom.GetFormalCharge())
    if formal_charge == 0:
        return False
    adjacent_charge = sum(
        int(cast(ob.OBAtom, neighbor_iter).GetFormalCharge())
        for neighbor_iter in ob.OBAtomAtomIter(atom)
        if not cast(ob.OBAtom, neighbor_iter).IsMetal()
    )
    return formal_charge + adjacent_charge == 0


def _is_locally_zwitterionic_organic_cation(atom: ob.OBAtom) -> bool:
    return int(atom.GetFormalCharge()) > 0 and _has_adjacent_formal_charge_cancellation(atom)


def _nonnegative_metal_unsaturated_organic_cation_count(
    obmol: ob.OBMol,
    metal_states: Sequence[dataclasses.MetalAtomPosition],
) -> int:
    if not any(int(metal_state.valence) >= 0 for metal_state in metal_states):
        return 0

    for atom_iter in ob.OBMolAtomIter(obmol):
        atom = cast(ob.OBAtom, atom_iter)
        if _is_unsaturated_organic_cation(atom) and not _is_locally_zwitterionic_organic_cation(
            atom
        ):
            return 1
    return 0


def _charge_localization_penalty_for_atom(
    atom: ob.OBAtom,
    *,
    is_conjugated: bool,
) -> float:
    formal_charge = cast(int, atom.GetFormalCharge())
    if formal_charge == 0:
        return 0.0

    magnitude = float(abs(formal_charge))
    atomic_num = cast(int, atom.GetAtomicNum())
    is_aromatic = bool(atom.IsAromatic())
    radical_electrons = cast(int, atom.GetSpinMultiplicity()) % 2

    def _normalize(value: float | None, *, lower: float, upper: float, fallback: float) -> float:
        if value is None or upper <= lower:
            return fallback
        clamped = min(max(value, lower), upper)
        return (clamped - lower) / (upper - lower)

    def _approx_local_electron_count(element_info: dataclasses.ElementInfo) -> int:
        return int(element_info.num_outer_electrons) - formal_charge + int(atom.GetTotalValence())

    def _generic_penalty() -> float:
        if atomic_num == 1:
            penalty = 4.0 * magnitude
        elif formal_charge < 0:
            if atomic_num in {9, 17, 35, 53}:
                penalty = 0.2 * magnitude
            elif atomic_num in {8, 16, 34, 52}:
                penalty = 0.3 * magnitude
            elif atomic_num in {7, 15, 33, 51}:
                penalty = 0.6 * magnitude
            elif atomic_num == 6:
                penalty = (1.5 if is_conjugated or is_aromatic else 4.0) * magnitude
            else:
                penalty = (1.0 if is_conjugated or is_aromatic else 2.0) * magnitude
        else:
            if atomic_num in {7, 15, 33, 51}:
                penalty = 0.4 * magnitude
            elif atomic_num in {8, 16, 34, 52}:
                penalty = 1.0 * magnitude
            elif atomic_num == 6:
                penalty = (1.2 if is_conjugated or is_aromatic else 3.0) * magnitude
            else:
                penalty = (0.8 if is_conjugated or is_aromatic else 1.8) * magnitude
        if not is_conjugated and not is_aromatic:
            penalty += 0.5 * magnitude
        return penalty

    def _main_group_charge_penalty() -> float | None:
        element_info = consts.NON_METAL_DICT.get(atomic_num)
        if element_info is None:
            return None

        local_electron_count = _approx_local_electron_count(element_info)
        total_valence = int(atom.GetTotalValence())
        default_valence = int(element_info.default_valence)
        shell_target = 2 if atomic_num == 1 else 8
        coordination_excess = max(total_valence - default_valence, 0)
        shell_deficiency = max(shell_target - local_electron_count, 0)
        shell_surplus = max(local_electron_count - shell_target, 0)
        ionization_norm = _normalize(
            consts.NON_METAL_FIRST_IONIZATION_ENERGY_EV.get(atomic_num),
            lower=7.5,
            upper=18.0,
            fallback=0.5,
        )
        electronegativity_norm = _normalize(
            consts.NON_METAL_PAULING_ELECTRONEGATIVITY.get(atomic_num),
            lower=1.9,
            upper=4.0,
            fallback=0.5,
        )
        local_environment_penalty = 0.05 * magnitude if not (is_conjugated or is_aromatic) else 0.0

        if formal_charge > 0 and shell_deficiency > 0:
            return magnitude * (
                0.45
                + 0.55 * float(shell_deficiency)
                + 0.9 * ionization_norm
                + local_environment_penalty
            )

        if local_electron_count >= shell_target and coordination_excess > 0:
            return magnitude * (
                0.18
                + 0.12 * float(coordination_excess)
                + 0.04 * float(shell_surplus)
                + local_environment_penalty
            )

        if formal_charge < 0 and local_electron_count >= shell_target:
            return magnitude * (
                0.10 + 0.60 * (1.0 - electronegativity_norm) + local_environment_penalty
            )

        return None

    penalty = _main_group_charge_penalty()
    if penalty is None:
        penalty = _generic_penalty()
    if radical_electrons > 0 and atomic_num == 6 and not (is_conjugated or is_aromatic):
        penalty += 2.0 * float(radical_electrons)
    return penalty


def _radical_localization_penalty_for_atom(
    atom: ob.OBAtom,
    *,
    is_conjugated: bool,
) -> float:
    radical_electrons = cast(int, atom.GetSpinMultiplicity()) % 2
    if radical_electrons <= 0:
        return 0.0

    magnitude = float(radical_electrons)
    atomic_num = cast(int, atom.GetAtomicNum())
    is_aromatic = bool(atom.IsAromatic())

    if atomic_num == 1:
        return 4.0 * magnitude
    if atomic_num == 6:
        return (0.6 if is_conjugated or is_aromatic else 2.5) * magnitude
    if atomic_num in {7, 15, 33, 51}:
        return (1.0 if is_conjugated or is_aromatic else 2.5) * magnitude
    if atomic_num in {8, 16, 34, 52}:
        return (1.5 if is_conjugated or is_aromatic else 3.0) * magnitude
    return (1.2 if is_conjugated or is_aromatic else 2.5) * magnitude


def _compute_organic_electronic_state_metrics(
    omol: pybel.Molecule,
    config: MolGRConfig | None = None,
) -> _OrganicElectronicStateMetrics:
    try:
        topology_metrics = compute_organic_topology_metrics(
            omol,
            (CONFIG if config is None else config).organic_topology,
        )
        obmol = cast(ob.OBMol, omol.OBMol)
        conjugated_atom_indices = set(topology_metrics.conjugated_atom_indices)
        radical_localization_penalty = 0.0
        charge_localization_penalty = 0.0
        for atom_iter in ob.OBMolAtomIter(obmol):
            atom = cast(ob.OBAtom, atom_iter)
            atom_idx = atom.GetIdx() - 1
            radical_localization_penalty += _radical_localization_penalty_for_atom(
                atom,
                is_conjugated=(atom_idx in conjugated_atom_indices),
            )
            charge_localization_penalty += _charge_localization_penalty_for_atom(
                atom,
                is_conjugated=(atom_idx in conjugated_atom_indices),
            )

        return _OrganicElectronicStateMetrics(
            aromatic_atom_count=topology_metrics.aromatic_atom_count,
            aromatic_ring_count=topology_metrics.aromatic_ring_count,
            aromatic_stability_score=topology_metrics.aromatic_stability_score,
            conjugated_atom_count=topology_metrics.conjugated_atom_count,
            conjugated_bond_count=topology_metrics.conjugated_bond_count,
            max_conjugated_component_size=topology_metrics.max_conjugated_component_size,
            radical_localization_penalty=radical_localization_penalty,
            charge_localization_penalty=charge_localization_penalty,
        )
    except Exception:  # noqa: BLE001
        return _OrganicElectronicStateMetrics(
            aromatic_atom_count=0,
            aromatic_ring_count=0,
            aromatic_stability_score=0.0,
            conjugated_atom_count=0,
            conjugated_bond_count=0,
            max_conjugated_component_size=0,
            radical_localization_penalty=float("inf"),
            charge_localization_penalty=float("inf"),
        )


def _annotate_candidate_discordance_features(
    candidate: MetalCandidateState,
    *,
    config: MolGRConfig | None = None,
) -> float:
    no_metal_state = candidate.no_metal_state
    if no_metal_state is None:
        raise ValueError("MetalCandidateState requires no_metal_state before discordance scoring")
    resolved_config = CONFIG if config is None else config
    metal_scoring_config = resolved_config.metal_scoring
    obmol = cast(ob.OBMol, no_metal_state.omol.OBMol)
    blocker_arrays = _build_coordination_blocker_arrays(
        obmol,
        metal_scoring_config=metal_scoring_config,
    )
    non_metal_atoms = _non_metal_atom_entries(obmol)

    inner_visible_diradical_count = 0
    inner_visible_same_sign_charge_count = 0
    visible_inner_atom_indices: set[int] = set()

    for metal_state in candidate.metal_states:
        metal_charge_sign = _charge_sign(int(metal_state.valence))
        for atom in _inner_visible_atoms_to_metal(
            non_metal_atoms,
            metal_state,
            metal_scoring_config=metal_scoring_config,
            blocker_arrays=blocker_arrays,
        ):
            atom_idx = int(atom.GetIdx())
            visible_inner_atom_indices.add(atom_idx)
            if _is_inner_visible_diradical_discordance_atom(atom):
                inner_visible_diradical_count += 1

            formal_charge = int(atom.GetFormalCharge())
            atom_charge_sign = _charge_sign(formal_charge)
            has_inner_same_sign_charge = (metal_charge_sign == 0 and formal_charge > 0) or (
                metal_charge_sign != 0 and atom_charge_sign == metal_charge_sign
            )
            if has_inner_same_sign_charge and not _has_adjacent_formal_charge_cancellation(atom):
                inner_visible_same_sign_charge_count += 1

    outer_or_invisible_adjacent_double_charge_count = 0
    outer_or_invisible_adjacent_same_sign_double_charge_count = 0
    outer_or_invisible_adjacent_opposite_sign_double_charge_count = 0
    outer_or_invisible_adjacent_unknown_metal_sign_double_charge_count = 0
    inner_visible_adjacent_carbanion_pair_count = 0
    inner_visible_conjugated_carbanion_pair_count = (
        _inner_visible_conjugated_charged_carbon_pair_count(obmol, visible_inner_atom_indices)
    )
    (
        negative_metal_count,
        negative_metal_has_outer_sphere_cation_exception,
        negative_metal_has_positive_metal_counterion_exception,
    ) = _negative_metal_discordance_count(
        obmol,
        candidate.metal_states,
        metal_scoring_config=metal_scoring_config,
    )
    zero_valent_metals_with_organic_cation_count = _zero_valent_metals_with_organic_cation_count(
        obmol,
        candidate.metal_states,
        total_charge=int(no_metal_state.total_charge)
        + sum(int(metal_state.valence) for metal_state in candidate.metal_states),
    )
    nonnegative_metal_unsaturated_organic_cation_count = (
        _nonnegative_metal_unsaturated_organic_cation_count(obmol, candidate.metal_states)
    )
    for bond_iter in ob.OBMolBondIter(obmol):
        bond = cast(ob.OBBond, bond_iter)
        begin_atom = cast(ob.OBAtom, bond.GetBeginAtom())
        end_atom = cast(ob.OBAtom, bond.GetEndAtom())
        begin_charge = int(begin_atom.GetFormalCharge())
        end_charge = int(end_atom.GetFormalCharge())
        begin_charge_sign = _charge_sign(begin_charge)
        end_charge_sign = _charge_sign(end_charge)
        if begin_charge_sign == 0 or end_charge_sign == 0:
            continue
        if begin_charge_sign != end_charge_sign:
            continue

        begin_idx = int(begin_atom.GetIdx())
        end_idx = int(end_atom.GetIdx())
        pair_both_inner_visible = (
            begin_idx in visible_inner_atom_indices and end_idx in visible_inner_atom_indices
        )
        is_inner_visible_adjacent_carbanion_pair = (
            pair_both_inner_visible
            and begin_atom.GetAtomicNum() == 6
            and end_atom.GetAtomicNum() == 6
            and begin_charge_sign == end_charge_sign
        )
        if is_inner_visible_adjacent_carbanion_pair:
            inner_visible_adjacent_carbanion_pair_count += 1
        elif not pair_both_inner_visible:
            outer_or_invisible_adjacent_double_charge_count += 1
            metal_charge_sign = _nearest_nonzero_metal_charge_sign_to_bond(
                begin_atom,
                end_atom,
                candidate.metal_states,
            )
            if metal_charge_sign == 0:
                outer_or_invisible_adjacent_unknown_metal_sign_double_charge_count += 1
            elif begin_charge_sign == metal_charge_sign:
                outer_or_invisible_adjacent_same_sign_double_charge_count += 1
            else:
                outer_or_invisible_adjacent_opposite_sign_double_charge_count += 1

    negative_metal_penalty = _NEGATIVE_METAL_DISCORDANCE_PENALTY * negative_metal_count
    discordance_count = (
        inner_visible_diradical_count
        + outer_or_invisible_adjacent_double_charge_count
        + inner_visible_adjacent_carbanion_pair_count
        + inner_visible_conjugated_carbanion_pair_count
        + inner_visible_same_sign_charge_count
        + negative_metal_penalty
        + zero_valent_metals_with_organic_cation_count
        + nonnegative_metal_unsaturated_organic_cation_count
    )
    candidate.metadata["metal_discordance_structural_count"] = discordance_count
    candidate.metadata["metal_discordance_aromatic_ring_deficit_count"] = 0
    candidate.metadata["metal_discordance_count"] = discordance_count
    candidate.metadata["metal_discordance_inner_visible_diradical_count"] = (
        inner_visible_diradical_count
    )
    candidate.metadata["metal_discordance_outer_or_invisible_adjacent_double_charge_count"] = (
        outer_or_invisible_adjacent_double_charge_count
    )
    candidate.metadata[
        "metal_discordance_outer_or_invisible_adjacent_same_sign_double_charge_count"
    ] = outer_or_invisible_adjacent_same_sign_double_charge_count
    candidate.metadata[
        "metal_discordance_outer_or_invisible_adjacent_opposite_sign_double_charge_count"
    ] = outer_or_invisible_adjacent_opposite_sign_double_charge_count
    candidate.metadata[
        "metal_discordance_outer_or_invisible_adjacent_unknown_metal_sign_double_charge_count"
    ] = outer_or_invisible_adjacent_unknown_metal_sign_double_charge_count
    candidate.metadata["metal_discordance_inner_visible_adjacent_carbanion_pair_count"] = (
        inner_visible_adjacent_carbanion_pair_count
    )
    candidate.metadata["metal_discordance_inner_visible_conjugated_carbanion_pair_count"] = (
        inner_visible_conjugated_carbanion_pair_count
    )
    candidate.metadata["metal_discordance_inner_visible_same_sign_charge_count"] = (
        inner_visible_same_sign_charge_count
    )
    candidate.metadata["metal_discordance_negative_metal_count"] = negative_metal_count
    candidate.metadata["metal_discordance_negative_metal_penalty"] = negative_metal_penalty
    candidate.metadata["metal_discordance_negative_metal_outer_sphere_cation_exception"] = (
        negative_metal_has_outer_sphere_cation_exception
    )
    candidate.metadata["metal_discordance_negative_metal_positive_metal_counterion_exception"] = (
        negative_metal_has_positive_metal_counterion_exception
    )
    candidate.metadata["metal_discordance_zero_valent_metals_with_organic_cation_count"] = (
        zero_valent_metals_with_organic_cation_count
    )
    candidate.metadata["metal_discordance_nonnegative_metal_unsaturated_organic_cation_count"] = (
        nonnegative_metal_unsaturated_organic_cation_count
    )
    return discordance_count


def _annotate_organic_electronic_state_consistency(
    candidate: MetalCandidateState,
    *,
    config: MolGRConfig | None = None,
) -> None:
    no_metal_state = candidate.no_metal_state
    if no_metal_state is None:
        raise ValueError("MetalCandidateState requires no_metal_state before organic-state scoring")

    resolved_config = CONFIG if config is None else config
    topology_config = resolved_config.organic_topology
    cache_key = f"organic_electronic_state_metrics:{astuple(topology_config)!r}"
    cached_metrics = cast(
        Optional[_OrganicElectronicStateMetrics],
        no_metal_state.metadata.get(cache_key),
    )
    if cached_metrics is None:
        cached_metrics = no_metal_state.get_cached_omol_value(
            cache_key,
            lambda omol: _compute_organic_electronic_state_metrics(omol, config=config),
        )
        no_metal_state.metadata[cache_key] = cached_metrics
        no_metal_state.metadata["organic_electronic_state_metrics"] = cached_metrics

    candidate.metadata["organic_aromatic_atom_count"] = cached_metrics.aromatic_atom_count
    candidate.metadata["organic_aromatic_ring_count"] = cached_metrics.aromatic_ring_count
    candidate.metadata["organic_aromatic_stability_score"] = cached_metrics.aromatic_stability_score
    candidate.metadata["organic_conjugated_atom_count"] = cached_metrics.conjugated_atom_count
    candidate.metadata["organic_conjugated_bond_count"] = cached_metrics.conjugated_bond_count
    candidate.metadata["organic_max_conjugated_component_size"] = (
        cached_metrics.max_conjugated_component_size
    )
    candidate.metadata["organic_radical_localization_penalty"] = (
        cached_metrics.radical_localization_penalty
    )
    candidate.metadata["organic_charge_localization_penalty"] = (
        cached_metrics.charge_localization_penalty
    )


def _prepare_candidate_with_no_metal_state(
    candidate: MetalCandidateState,
    no_metal_state: ReconstructionState,
    *,
    config: MolGRConfig | None = None,
) -> MetalCandidateState:
    """Attach the shared no-metal reconstruction and count discordance features."""

    candidate_machine = MetalCandidateStateMachine.from_candidate_state(candidate)
    candidate_machine.set_no_metal_state("reconstruct_no_metal", no_metal_state)
    candidate_machine.annotate("score_candidate")
    scored_candidate = candidate_machine.freeze()
    if config is None:
        scored_candidate.combined_score()
    else:
        scored_candidate.combined_score(config=config)
    _annotate_candidate_discordance_features(scored_candidate, config=config)
    return scored_candidate


def _annotate_selected_candidate_metrics(
    candidate: MetalCandidateState,
    *,
    config: MolGRConfig | None = None,
) -> None:
    if candidate.score is None:
        if config is None:
            candidate.combined_score()
        else:
            candidate.combined_score(config=config)
    if "organic_aromatic_atom_count" not in candidate.metadata:
        _annotate_organic_electronic_state_consistency(candidate, config=config)


def _ensure_candidate_organic_metrics(
    candidate: MetalCandidateState,
    *,
    config: MolGRConfig | None = None,
) -> None:
    if "organic_aromatic_ring_count" not in candidate.metadata:
        _annotate_organic_electronic_state_consistency(candidate, config=config)
    if "organic_aromatic_stability_score" not in candidate.metadata:
        candidate.metadata["organic_aromatic_stability_score"] = float(
            candidate.metadata.get("organic_aromatic_ring_count", 0)
        )


def _annotate_candidate_set_discordance_features(
    candidates: Sequence[MetalCandidateState],
    *,
    config: MolGRConfig | None = None,
) -> None:
    if not candidates:
        return

    max_aromatic_ring_count = 0
    max_aromatic_stability_score = 0.0
    for candidate in candidates:
        _ensure_candidate_organic_metrics(candidate, config=config)
        max_aromatic_ring_count = max(
            max_aromatic_ring_count,
            int(candidate.metadata.get("organic_aromatic_ring_count", 0)),
        )
        max_aromatic_stability_score = max(
            max_aromatic_stability_score,
            float(candidate.metadata.get("organic_aromatic_stability_score", 0.0)),
        )

    for candidate in candidates:
        structural_discordance_count = float(
            candidate.metadata.get(
                "metal_discordance_structural_count",
                candidate.metadata.get("metal_discordance_count", 0),
            )
        )
        aromatic_ring_deficit_count = max(
            0,
            max_aromatic_ring_count - int(candidate.metadata.get("organic_aromatic_ring_count", 0)),
        )
        aromatic_stability_deficit = max(
            0.0,
            max_aromatic_stability_score
            - float(candidate.metadata.get("organic_aromatic_stability_score", 0.0)),
        )
        candidate.metadata["metal_discordance_structural_count"] = structural_discordance_count
        candidate.metadata["metal_discordance_max_aromatic_ring_count"] = max_aromatic_ring_count
        candidate.metadata["metal_discordance_max_aromatic_stability_score"] = (
            max_aromatic_stability_score
        )
        candidate.metadata["metal_discordance_aromatic_ring_deficit_count"] = (
            aromatic_ring_deficit_count
        )
        candidate.metadata["metal_discordance_aromatic_stability_deficit"] = (
            aromatic_stability_deficit
        )
        candidate.metadata["metal_discordance_count"] = (
            structural_discordance_count + aromatic_stability_deficit
        )


def select_best_candidate(
    scored_candidates: Sequence[MetalCandidateState],
    *,
    config: MolGRConfig | None = None,
) -> Optional[MetalCandidateState]:
    if not scored_candidates:
        return None

    _annotate_candidate_set_discordance_features(scored_candidates, config=config)
    min_discordance_count = min(
        float(candidate.metadata.get("metal_discordance_count", 0))
        for candidate in scored_candidates
    )
    discordance_filtered_candidates: list[MetalCandidateState] = []
    for scored_candidate in scored_candidates:
        discordance_count = float(scored_candidate.metadata.get("metal_discordance_count", 0))
        passes_discordance_filter = discordance_count == min_discordance_count
        scored_candidate.metadata["passes_metal_discordance_filter"] = passes_discordance_filter
        if passes_discordance_filter:
            discordance_filtered_candidates.append(scored_candidate)

    for scored_candidate in discordance_filtered_candidates:
        _annotate_selected_candidate_metrics(scored_candidate, config=config)

    best_candidate: Optional[MetalCandidateState] = None
    best_selection_key: Optional[Tuple[float, float, int]] = None
    for scored_candidate in discordance_filtered_candidates:
        score_value = (
            cast(float, scored_candidate.score)
            if scored_candidate.score is not None
            else scored_candidate.combined_score(config=config)
        )
        selection_key = (
            float(min_discordance_count),
            float(score_value),
            int(scored_candidate.metadata.get("combination_index", 0)),
        )
        scored_candidate.metadata["selection_key"] = selection_key
        if best_selection_key is not None and selection_key >= best_selection_key:
            continue
        best_selection_key = selection_key
        best_candidate = scored_candidate
    return best_candidate


__all__ = [
    "_charge_localization_penalty_for_atom",
    "_annotate_candidate_set_discordance_features",
    "_radical_localization_penalty_for_atom",
    "_prepare_candidate_with_no_metal_state",
    "select_best_candidate",
]
