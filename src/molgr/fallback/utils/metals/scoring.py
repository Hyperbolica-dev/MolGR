"""Scoring and selection utilities for metal-aware fallback reconstruction."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.config import MetalScoringConfig, MolGRConfig, resolve_config
from molgr.fallback.state import (
    MetalCandidateState,
    MetalCandidateStateMachine,
    ReconstructionState,
)
from molgr.fallback.utils import consts, dataclasses
from molgr.fallback.utils.organic_topology import compute_organic_topology_metrics


_MetalSiteAnchorKey = Tuple[int, int, int, int, int]
_CoordinationBlocker = Tuple[int, Tuple[float, float, float], float]
_METAL_SITE_CACHE_COORDINATE_SCALE = 1_000_000
_WEIGHTED_SELECTION_FIELD_NAMES = (
    "organic_aromatic_ring_loss",
    "organic_max_conjugated_component_loss",
    "organic_charge_localization_penalty",
    "organic_aromatic_atom_loss",
    "organic_conjugated_atom_loss",
    "organic_radical_localization_penalty",
    "metal_coordination_access_penalty",
    "metal_same_element_valence_spread_penalty",
    "metal_electrostatic_penalty",
    "metal_donor_penalty",
    "metal_prior_penalty",
)


@dataclass(frozen=True)
class _MetalSiteHeuristicScore:
    electrostatic_support: float
    anionic_donor_support: float
    neutral_donor_support: float
    coordination_access_penalty: float
    visible_coordination_reward: float
    negative_metal_visible_coordination_penalty: float
    obstructed_opposite_charge_penalty: float
    electrostatic_penalty: float
    donor_penalty: float


@dataclass(frozen=True)
class _MetalSiteEnvironmentProfile:
    electrostatic_support: float
    visible_anionic_donor_support: float
    visible_neutral_donor_support: float
    visible_effective_donor_support: float
    obstructed_negative_effective_donor_support: float


@dataclass(frozen=True)
class _OrganicElectronicStateMetrics:
    aromatic_atom_count: int
    aromatic_ring_count: int
    conjugated_atom_count: int
    conjugated_bond_count: int
    max_conjugated_component_size: int
    radical_localization_penalty: float
    charge_localization_penalty: float


@dataclass(frozen=True)
class _OrganicElectronicStateSelectionContext:
    max_aromatic_atom_count: int
    max_aromatic_ring_count: int
    max_conjugated_atom_count: int
    max_conjugated_component_size: int


@dataclass(frozen=True)
class _WeightedSelectionContext:
    field_names: Tuple[str, ...]
    best_values: Tuple[float, ...]
    weights: Tuple[float, ...]
    scales: Tuple[float, ...]


def _distance_to_metal(
    atom: ob.OBAtom,
    metal_state: dataclasses.MetalAtomPosition,
) -> float:
    dx = float(atom.GetX()) - metal_state.position_x
    dy = float(atom.GetY()) - metal_state.position_y
    dz = float(atom.GetZ()) - metal_state.position_z
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def _distance_weight(distance: float, cutoff: float, min_distance_angstrom: float) -> float:
    if distance <= 0.0 or distance >= cutoff:
        return 0.0
    scaled = distance / cutoff
    attenuation = max(0.0, 1.0 - scaled * scaled)
    return (attenuation * attenuation) / max(distance, min_distance_angstrom)


def _atom_coordinates(atom: ob.OBAtom) -> Tuple[float, float, float]:
    return float(atom.GetX()), float(atom.GetY()), float(atom.GetZ())


def _quantized_metal_site_coordinate(value: float) -> int:
    return int(round(float(value) * _METAL_SITE_CACHE_COORDINATE_SCALE))


def _metal_site_anchor_key(
    metal_state: dataclasses.MetalAtomPosition,
) -> _MetalSiteAnchorKey:
    return (
        int(metal_state.idx),
        int(metal_state.element_idx),
        _quantized_metal_site_coordinate(metal_state.position_x),
        _quantized_metal_site_coordinate(metal_state.position_y),
        _quantized_metal_site_coordinate(metal_state.position_z),
    )


def _distance_point_to_segment(
    point: Tuple[float, float, float],
    segment_start: Tuple[float, float, float],
    segment_end: Tuple[float, float, float],
) -> float:
    vx = segment_end[0] - segment_start[0]
    vy = segment_end[1] - segment_start[1]
    vz = segment_end[2] - segment_start[2]
    seg_len_sq = vx * vx + vy * vy + vz * vz
    if seg_len_sq <= 1.0e-12:
        dx = point[0] - segment_start[0]
        dy = point[1] - segment_start[1]
        dz = point[2] - segment_start[2]
        return math.sqrt(dx * dx + dy * dy + dz * dz)

    wx = point[0] - segment_start[0]
    wy = point[1] - segment_start[1]
    wz = point[2] - segment_start[2]
    projection = max(0.0, min(1.0, (wx * vx + wy * vy + wz * vz) / seg_len_sq))
    closest = (
        segment_start[0] + projection * vx,
        segment_start[1] + projection * vy,
        segment_start[2] + projection * vz,
    )
    dx = point[0] - closest[0]
    dy = point[1] - closest[1]
    dz = point[2] - closest[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _build_coordination_blockers(
    obmol: ob.OBMol,
    *,
    metal_scoring_config: MetalScoringConfig,
) -> Tuple[_CoordinationBlocker, ...]:
    blockers: List[_CoordinationBlocker] = []
    for blocker_iter in ob.OBMolAtomIter(obmol):
        blocker = cast(ob.OBAtom, blocker_iter)
        blocker_radius = (
            metal_scoring_config.metal_access_radius_scale
            * float(ob.GetCovalentRad(int(blocker.GetAtomicNum())))
            + metal_scoring_config.metal_access_clearance_angstrom
        )
        if blocker_radius <= 0.0:
            continue
        blockers.append((int(blocker.GetIdx()), _atom_coordinates(blocker), blocker_radius))
    return tuple(blockers)


def _has_unobstructed_coordination_path_from_blockers(
    atom_idx: int,
    atom_coordinates: Tuple[float, float, float],
    segment_start: Tuple[float, float, float],
    blockers: Sequence[_CoordinationBlocker],
) -> bool:
    for blocker_idx, blocker_coordinates, blocker_radius in blockers:
        if blocker_idx == atom_idx:
            continue
        if _distance_point_to_segment(blocker_coordinates, segment_start, atom_coordinates) < blocker_radius:
            return False
    return True


def _bond_order(bond: ob.OBBond) -> int:
    if bond.IsAromatic():
        return 2
    return cast(int, bond.GetBondOrder())


def _has_multiple_bond_to_atomic_num(atom: ob.OBAtom, atomic_nums: Set[int]) -> bool:
    for bond in ob.OBAtomBondIter(atom):
        obbond = cast(ob.OBBond, bond)
        if _bond_order(obbond) < 2:
            continue
        neighbor = (
            cast(ob.OBAtom, obbond.GetEndAtom())
            if cast(ob.OBAtom, obbond.GetBeginAtom()).GetIdx() == atom.GetIdx()
            else cast(ob.OBAtom, obbond.GetBeginAtom())
        )
        if cast(int, neighbor.GetAtomicNum()) in atomic_nums:
            return True
    return False


def _is_carbonyl_like_oxygen(atom: ob.OBAtom) -> bool:
    if cast(int, atom.GetAtomicNum()) not in {8, 16, 34, 52}:
        return False
    return _has_multiple_bond_to_atomic_num(atom, {6})


def _is_amide_like_nitrogen(atom: ob.OBAtom) -> bool:
    if cast(int, atom.GetAtomicNum()) not in {7, 15, 33, 51}:
        return False
    for bond in ob.OBAtomBondIter(atom):
        obbond = cast(ob.OBBond, bond)
        if _bond_order(obbond) != 1:
            continue
        neighbor = (
            cast(ob.OBAtom, obbond.GetEndAtom())
            if cast(ob.OBAtom, obbond.GetBeginAtom()).GetIdx() == atom.GetIdx()
            else cast(ob.OBAtom, obbond.GetBeginAtom())
        )
        if cast(int, neighbor.GetAtomicNum()) != 6:
            continue
        if _has_multiple_bond_to_atomic_num(neighbor, {8, 16, 34, 52}):
            return True
    return False


def _is_nitrile_like_nitrogen(atom: ob.OBAtom) -> bool:
    if cast(int, atom.GetAtomicNum()) not in {7, 15, 33, 51}:
        return False
    return any(_bond_order(cast(ob.OBBond, bond)) >= 3 for bond in ob.OBAtomBondIter(atom))


def _atom_donor_support(atom: ob.OBAtom) -> Tuple[float, float]:
    atomic_num = cast(int, atom.GetAtomicNum())
    formal_charge = cast(int, atom.GetFormalCharge())
    if atomic_num == 1 or atom.IsMetal() or formal_charge > 0:
        return 0.0, 0.0

    if formal_charge < 0:
        magnitude = float(abs(formal_charge))
        if atomic_num in {9, 17, 35, 53}:
            return 2.2 * magnitude, 0.0
        if atomic_num in {8, 16, 34, 52}:
            return 2.0 * magnitude, 0.0
        if atomic_num in {7, 15, 33, 51}:
            return 1.6 * magnitude, 0.0
        if atomic_num == 6:
            return 1.0 * magnitude, 0.0
        return 1.2 * magnitude, 0.0

    if atomic_num in {9, 17, 35, 53}:
        return 0.0, 0.0
    if atomic_num in {8, 16, 34, 52}:
        return 0.0, 0.7 if _is_carbonyl_like_oxygen(atom) else 1.0
    if atomic_num in {7, 15, 33, 51}:
        if _is_amide_like_nitrogen(atom):
            return 0.0, 0.2
        if _is_nitrile_like_nitrogen(atom):
            return 0.0, 0.4
        if atom.IsAromatic():
            return 0.0, 0.7
        return 0.0, 0.9
    return 0.0, 0.0


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
                0.10
                + 0.60 * (1.0 - electronegativity_norm)
                + local_environment_penalty
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


def _compute_organic_electronic_state_metrics(omol: pybel.Molecule) -> _OrganicElectronicStateMetrics:
    try:
        topology_metrics = compute_organic_topology_metrics(omol)
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
            conjugated_atom_count=0,
            conjugated_bond_count=0,
            max_conjugated_component_size=0,
            radical_localization_penalty=float("inf"),
            charge_localization_penalty=float("inf"),
        )


def _build_metal_site_environment_profile(
    omol: pybel.Molecule,
    metal_state: dataclasses.MetalAtomPosition,
    *,
    metal_scoring_config: MetalScoringConfig,
) -> _MetalSiteEnvironmentProfile:
    obmol = cast(ob.OBMol, omol.OBMol)
    blockers = _build_coordination_blockers(obmol, metal_scoring_config=metal_scoring_config)
    segment_start = (
        float(metal_state.position_x),
        float(metal_state.position_y),
        float(metal_state.position_z),
    )

    electrostatic_support = 0.0
    visible_anionic_donor_support = 0.0
    visible_neutral_donor_support = 0.0
    visible_effective_donor_support = 0.0
    obstructed_negative_effective_donor_support = 0.0

    for atom_iter in ob.OBMolAtomIter(obmol):
        atom = cast(ob.OBAtom, atom_iter)
        if atom.IsMetal():
            continue

        atom_idx = int(atom.GetIdx())
        atom_coordinates = _atom_coordinates(atom)
        distance = _distance_to_metal(atom, metal_state)
        electrostatic_weight = _distance_weight(
            distance,
            metal_scoring_config.metal_local_potential_cutoff_angstrom,
            metal_scoring_config.min_distance_angstrom,
        )
        formal_charge = float(cast(int, atom.GetFormalCharge()))
        if electrostatic_weight > 0.0:
            electrostatic_support += -formal_charge * electrostatic_weight

        coordination_weight = _distance_weight(
            distance,
            metal_scoring_config.metal_donor_cutoff_angstrom,
            metal_scoring_config.min_distance_angstrom,
        )
        if coordination_weight <= 0.0:
            continue

        atom_anionic_support, atom_neutral_support = _atom_donor_support(atom)
        atom_effective_donor_support = (
            atom_anionic_support
            + metal_scoring_config.local_neutral_donor_weight * atom_neutral_support
        )
        if atom_effective_donor_support <= 0.0:
            continue

        weighted_effective_support = atom_effective_donor_support * coordination_weight
        if _has_unobstructed_coordination_path_from_blockers(
            atom_idx,
            atom_coordinates,
            segment_start,
            blockers,
        ):
            visible_anionic_donor_support += atom_anionic_support * coordination_weight
            visible_neutral_donor_support += atom_neutral_support * coordination_weight
            visible_effective_donor_support += weighted_effective_support
            continue

        if formal_charge < 0.0:
            obstructed_negative_effective_donor_support += weighted_effective_support

    return _MetalSiteEnvironmentProfile(
        electrostatic_support=electrostatic_support,
        visible_anionic_donor_support=visible_anionic_donor_support,
        visible_neutral_donor_support=visible_neutral_donor_support,
        visible_effective_donor_support=visible_effective_donor_support,
        obstructed_negative_effective_donor_support=obstructed_negative_effective_donor_support,
    )


def _score_metal_site_environment_from_profile(
    profile: _MetalSiteEnvironmentProfile,
    metal_state: dataclasses.MetalAtomPosition,
    *,
    metal_scoring_config: MetalScoringConfig,
) -> _MetalSiteHeuristicScore:
    anionic_donor_support = profile.visible_anionic_donor_support
    neutral_donor_support = profile.visible_neutral_donor_support
    coordination_access_penalty = 0.0
    visible_coordination_reward = 0.0
    negative_metal_visible_coordination_penalty = 0.0
    obstructed_opposite_charge_penalty = 0.0
    metal_valence = float(metal_state.valence)

    if metal_valence >= 0.0 and profile.visible_effective_donor_support > 0.0:
        visible_coordination_reward = (
            metal_scoring_config.visible_coordination_reward_weight
            * profile.visible_effective_donor_support
        )
        coordination_access_penalty -= visible_coordination_reward

    if metal_valence < 0.0 and profile.visible_effective_donor_support > 0.0:
        negative_metal_visible_coordination_penalty = (
            metal_scoring_config.negative_metal_visible_coordination_penalty_weight
            * max(abs(metal_valence), 1.0)
            * profile.visible_effective_donor_support
        )
        coordination_access_penalty += negative_metal_visible_coordination_penalty

    if metal_valence > 0.0 and profile.obstructed_negative_effective_donor_support > 0.0:
        obstructed_opposite_charge_penalty = (
            metal_scoring_config.obstructed_opposite_charge_penalty_weight
            * max(abs(metal_valence), 1.0)
            * profile.obstructed_negative_effective_donor_support
        )
        coordination_access_penalty += obstructed_opposite_charge_penalty

    target_valence = float(max(metal_state.valence, 0))
    electrostatic_target = metal_scoring_config.local_potential_target_per_valence * target_valence
    electrostatic_under = max(electrostatic_target - profile.electrostatic_support, 0.0)
    electrostatic_over = max(profile.electrostatic_support - electrostatic_target, 0.0)
    electrostatic_penalty = (
        electrostatic_under
        + metal_scoring_config.local_potential_oversupport_weight * electrostatic_over
    )

    effective_donor_support = (
        anionic_donor_support
        + metal_scoring_config.local_neutral_donor_weight * neutral_donor_support
    )
    donor_target = metal_scoring_config.local_donor_target_per_valence * target_valence
    donor_under = max(donor_target - effective_donor_support, 0.0)
    donor_over = max(effective_donor_support - donor_target, 0.0)
    donor_penalty = donor_under + metal_scoring_config.local_donor_oversupport_weight * donor_over

    return _MetalSiteHeuristicScore(
        electrostatic_support=profile.electrostatic_support,
        anionic_donor_support=anionic_donor_support,
        neutral_donor_support=neutral_donor_support,
        coordination_access_penalty=coordination_access_penalty,
        visible_coordination_reward=visible_coordination_reward,
        negative_metal_visible_coordination_penalty=negative_metal_visible_coordination_penalty,
        obstructed_opposite_charge_penalty=obstructed_opposite_charge_penalty,
        electrostatic_penalty=electrostatic_penalty,
        donor_penalty=donor_penalty,
    )


def _same_element_valence_spread_penalty(
    metal_states: Sequence[dataclasses.MetalAtomPosition],
    *,
    metal_scoring_config: MetalScoringConfig,
) -> float:
    grouped_valences: Dict[str, List[int]] = defaultdict(list)
    for metal_state in metal_states:
        grouped_valences[metal_state.symbol].append(metal_state.valence)

    penalty = 0.0
    for valences in grouped_valences.values():
        if len(valences) >= 2:
            penalty += metal_scoring_config.same_element_valence_spread_weight * float(
                max(valences) - min(valences)
            )
    return penalty


def _get_cached_metal_site_environment_profile(
    no_metal_state: ReconstructionState,
    metal_state: dataclasses.MetalAtomPosition,
    *,
    metal_scoring_config: MetalScoringConfig,
) -> _MetalSiteEnvironmentProfile:
    profiles = no_metal_state.get_cached_revision_value(
        f"metal_site_environment_profiles::{metal_scoring_config!r}",
        dict,
    )
    profile_key = _metal_site_anchor_key(metal_state)
    cached_profile = cast(Optional[_MetalSiteEnvironmentProfile], profiles.get(profile_key))
    if cached_profile is not None:
        return cached_profile

    profile = _build_metal_site_environment_profile(
        no_metal_state.omol,
        metal_state,
        metal_scoring_config=metal_scoring_config,
    )
    profiles[profile_key] = profile
    return profile


def _metal_state_assignment_penalty(metal_state: dataclasses.MetalAtomPosition) -> float:
    penalty = 0.0
    if metal_state.valence <= 0:
        penalty += 10.0 * max(abs(metal_state.valence), 1)

    prior_list = consts.METAL_VALENCE_AVAILABLE_PRIOR.get(metal_state.symbol, [])
    minor_list = consts.METAL_VALENCE_AVAILABLE_MINOR.get(metal_state.symbol, [])
    if metal_state.valence not in prior_list:
        if metal_state.valence in minor_list:
            penalty += 10.0
        else:
            penalty += 20.0
    return penalty


def _annotate_metal_environment_consistency(
    candidate: MetalCandidateState,
    *,
    config: MolGRConfig | None = None,
) -> None:
    no_metal_state = candidate.no_metal_state
    if no_metal_state is None:
        raise ValueError("MetalCandidateState requires no_metal_state before metal scoring")
    metal_scoring_config = resolve_config(config).metal_scoring

    total_prior_penalty = 0.0
    total_coordination_access_penalty = 0.0
    total_electrostatic_penalty = 0.0
    total_donor_penalty = 0.0
    site_breakdown: List[Dict[str, float | int | str]] = []

    for metal_state in candidate.metal_states:
        prior_penalty = _metal_state_assignment_penalty(metal_state)
        site_score = _score_metal_site_environment_from_profile(
            _get_cached_metal_site_environment_profile(
                no_metal_state,
                metal_state,
                metal_scoring_config=metal_scoring_config,
            ),
            metal_state,
            metal_scoring_config=metal_scoring_config,
        )
        total_prior_penalty += prior_penalty
        total_coordination_access_penalty += site_score.coordination_access_penalty
        total_electrostatic_penalty += site_score.electrostatic_penalty
        total_donor_penalty += site_score.donor_penalty
        site_breakdown.append(
            {
                "idx": metal_state.idx,
                "symbol": metal_state.symbol,
                "valence": metal_state.valence,
                "radical_num": metal_state.radical_num,
                "prior_penalty": prior_penalty,
                "electrostatic_support": site_score.electrostatic_support,
                "anionic_donor_support": site_score.anionic_donor_support,
                "neutral_donor_support": site_score.neutral_donor_support,
                "coordination_access_penalty": site_score.coordination_access_penalty,
                "visible_coordination_reward": site_score.visible_coordination_reward,
                "negative_metal_visible_coordination_penalty": site_score.negative_metal_visible_coordination_penalty,
                "obstructed_opposite_charge_penalty": site_score.obstructed_opposite_charge_penalty,
                "electrostatic_penalty": site_score.electrostatic_penalty,
                "donor_penalty": site_score.donor_penalty,
            }
        )

    candidate.metadata["metal_prior_penalty"] = total_prior_penalty
    candidate.metadata["metal_coordination_access_penalty"] = total_coordination_access_penalty
    candidate.metadata["metal_same_element_valence_spread_penalty"] = _same_element_valence_spread_penalty(
        candidate.metal_states,
        metal_scoring_config=metal_scoring_config,
    )
    candidate.metadata["metal_electrostatic_penalty"] = total_electrostatic_penalty
    candidate.metadata["metal_donor_penalty"] = total_donor_penalty
    candidate.metadata["metal_environment_breakdown"] = site_breakdown


def _annotate_organic_electronic_state_consistency(candidate: MetalCandidateState) -> None:
    no_metal_state = candidate.no_metal_state
    if no_metal_state is None:
        raise ValueError("MetalCandidateState requires no_metal_state before organic-state scoring")

    cached_metrics = cast(
        Optional[_OrganicElectronicStateMetrics],
        no_metal_state.metadata.get("organic_electronic_state_metrics"),
    )
    if cached_metrics is None:
        cached_metrics = no_metal_state.get_cached_omol_value(
            "organic_electronic_state_metrics",
            _compute_organic_electronic_state_metrics,
        )
        no_metal_state.metadata["organic_electronic_state_metrics"] = cached_metrics

    candidate.metadata["organic_aromatic_atom_count"] = cached_metrics.aromatic_atom_count
    candidate.metadata["organic_aromatic_ring_count"] = cached_metrics.aromatic_ring_count
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


def _score_candidate_with_no_metal_state(
    candidate: MetalCandidateState,
    no_metal_state: ReconstructionState,
    *,
    config: MolGRConfig | None = None,
) -> MetalCandidateState:
    """Attach the shared no-metal reconstruction and score the metal candidate."""

    candidate_machine = MetalCandidateStateMachine.from_candidate_state(candidate)
    candidate_machine.set_no_metal_state("reconstruct_no_metal", no_metal_state)
    candidate_machine.annotate("score_candidate")
    scored_candidate = candidate_machine.freeze()
    if config is None:
        scored_candidate.combined_score()
    else:
        scored_candidate.combined_score(config=config)
    _annotate_organic_electronic_state_consistency(scored_candidate)
    _annotate_metal_environment_consistency(scored_candidate, config=config)
    return scored_candidate


def _build_organic_electronic_state_selection_context(
    candidates: Sequence[MetalCandidateState],
) -> _OrganicElectronicStateSelectionContext:
    return _OrganicElectronicStateSelectionContext(
        max_aromatic_atom_count=max(int(candidate.metadata.get("organic_aromatic_atom_count", 0)) for candidate in candidates),
        max_aromatic_ring_count=max(int(candidate.metadata.get("organic_aromatic_ring_count", 0)) for candidate in candidates),
        max_conjugated_atom_count=max(int(candidate.metadata.get("organic_conjugated_atom_count", 0)) for candidate in candidates),
        max_conjugated_component_size=max(
            int(candidate.metadata.get("organic_max_conjugated_component_size", 0))
            for candidate in candidates
        ),
    )


def _organic_electronic_state_key(
    candidate: MetalCandidateState,
    selection_context: _OrganicElectronicStateSelectionContext,
) -> Tuple[int, int, float, int, int, float]:
    aromatic_atom_loss = selection_context.max_aromatic_atom_count - int(
        candidate.metadata.get("organic_aromatic_atom_count", 0)
    )
    aromatic_ring_loss = selection_context.max_aromatic_ring_count - int(
        candidate.metadata.get("organic_aromatic_ring_count", 0)
    )
    conjugated_atom_loss = selection_context.max_conjugated_atom_count - int(
        candidate.metadata.get("organic_conjugated_atom_count", 0)
    )
    max_conjugated_component_loss = selection_context.max_conjugated_component_size - int(
        candidate.metadata.get("organic_max_conjugated_component_size", 0)
    )
    radical_localization_penalty = float(
        candidate.metadata.get("organic_radical_localization_penalty", float("inf"))
    )
    charge_localization_penalty = float(
        candidate.metadata.get("organic_charge_localization_penalty", float("inf"))
    )

    candidate.metadata["organic_aromatic_atom_loss"] = aromatic_atom_loss
    candidate.metadata["organic_aromatic_ring_loss"] = aromatic_ring_loss
    candidate.metadata["organic_conjugated_atom_loss"] = conjugated_atom_loss
    candidate.metadata["organic_max_conjugated_component_loss"] = max_conjugated_component_loss
    candidate.metadata["organic_electronic_state_key"] = (
        aromatic_ring_loss,
        max_conjugated_component_loss,
        charge_localization_penalty,
        aromatic_atom_loss,
        conjugated_atom_loss,
        radical_localization_penalty,
    )
    return cast(Tuple[int, int, float, int, int, float], candidate.metadata["organic_electronic_state_key"])


def _organic_score_bucket_index(
    score_value: float,
    best_force_field_score: float,
    *,
    config: MolGRConfig | None = None,
) -> int:
    if score_value <= best_force_field_score:
        return 0
    metal_scoring_config = resolve_config(config).metal_scoring
    baseline_scale = max(abs(best_force_field_score), 1.0)
    relative_excess = (score_value - best_force_field_score) / baseline_scale
    return int(relative_excess // metal_scoring_config.organic_score_bucket_relative_ratio)


def _passes_organic_force_field_guard(
    score_value: float,
    best_force_field_score: float,
    *,
    config: MolGRConfig | None = None,
) -> bool:
    metal_scoring_config = resolve_config(config).metal_scoring
    hard_max_ratio = float(metal_scoring_config.organic_force_field_hard_max_ratio)
    if hard_max_ratio <= 0.0 or best_force_field_score <= 0.0:
        return True
    return score_value <= best_force_field_score * hard_max_ratio


def _raw_candidate_selection_key(
    candidate: MetalCandidateState,
    organic_selection_context: _OrganicElectronicStateSelectionContext,
    best_force_field_score: float,
    *,
    config: MolGRConfig | None = None,
) -> Tuple[int, int, float, int, int, float, float, float, float, float, float, int, float, int]:
    score_value = cast(float, candidate.score)
    organic_bucket = _organic_score_bucket_index(score_value, best_force_field_score, config=config)
    candidate.metadata["organic_score_bucket"] = organic_bucket
    organic_state_key = _organic_electronic_state_key(candidate, organic_selection_context)
    return (
        *organic_state_key,
        float(candidate.metadata.get("metal_coordination_access_penalty", 0.0)),
        float(candidate.metadata.get("metal_same_element_valence_spread_penalty", 0.0)),
        float(candidate.metadata.get("metal_electrostatic_penalty", 0.0)),
        float(candidate.metadata.get("metal_donor_penalty", 0.0)),
        float(candidate.metadata.get("metal_prior_penalty", candidate.metadata.get("metal_assignment_rank", 0.0))),
        organic_bucket,
        score_value,
        int(candidate.metadata.get("combination_index", 0)),
    )


def _normalize_weighted_selection_values(
    raw_values: Sequence[float],
    *,
    default_value: float,
) -> Tuple[float, ...]:
    if not raw_values:
        return (default_value,) * len(_WEIGHTED_SELECTION_FIELD_NAMES)
    normalized = [max(float(value), 0.0) for value in raw_values]
    if len(normalized) < len(_WEIGHTED_SELECTION_FIELD_NAMES):
        normalized.extend([normalized[-1]] * (len(_WEIGHTED_SELECTION_FIELD_NAMES) - len(normalized)))
    return tuple(normalized[: len(_WEIGHTED_SELECTION_FIELD_NAMES)])


def _build_weighted_selection_context(
    candidates: Sequence[MetalCandidateState],
    organic_selection_context: _OrganicElectronicStateSelectionContext,
    best_force_field_score: float,
    *,
    config: MolGRConfig | None = None,
) -> _WeightedSelectionContext:
    for candidate in candidates:
        raw_selection_key = _raw_candidate_selection_key(
            candidate,
            organic_selection_context,
            best_force_field_score,
            config=config,
        )
        candidate.metadata["raw_selection_key"] = raw_selection_key
        candidate.metadata["weighted_selection_metric_values"] = tuple(
            float(raw_selection_key[idx]) for idx in range(len(_WEIGHTED_SELECTION_FIELD_NAMES))
        )
    best_values = tuple(
        min(
            float(cast(Tuple[float, ...], candidate.metadata["weighted_selection_metric_values"])[idx])
            for candidate in candidates
        )
        for idx in range(len(_WEIGHTED_SELECTION_FIELD_NAMES))
    )
    metal_scoring_config = resolve_config(config).metal_scoring
    return _WeightedSelectionContext(
        field_names=_WEIGHTED_SELECTION_FIELD_NAMES,
        best_values=best_values,
        weights=_normalize_weighted_selection_values(
            metal_scoring_config.selection_weight_values,
            default_value=1.0,
        ),
        scales=_normalize_weighted_selection_values(
            metal_scoring_config.selection_scale_values,
            default_value=1.0,
        ),
    )


def select_best_candidate(
    scored_candidates: Sequence[MetalCandidateState],
    *,
    config: MolGRConfig | None = None,
) -> Optional[MetalCandidateState]:
    if not scored_candidates:
        return None

    best_force_field_score = min(cast(float, candidate.score) for candidate in scored_candidates)
    eligible_candidates: list[MetalCandidateState] = []
    for scored_candidate in scored_candidates:
        score_value = cast(float, scored_candidate.score)
        passes_force_field_guard = _passes_organic_force_field_guard(
            score_value,
            best_force_field_score,
            config=config,
        )
        scored_candidate.metadata["passes_organic_force_field_guard"] = passes_force_field_guard
        if passes_force_field_guard:
            eligible_candidates.append(scored_candidate)
    if not eligible_candidates:
        return None

    organic_selection_context = _build_organic_electronic_state_selection_context(eligible_candidates)
    weighted_selection_context = _build_weighted_selection_context(
        eligible_candidates,
        organic_selection_context,
        best_force_field_score,
        config=config,
    )
    best_candidate: Optional[MetalCandidateState] = None
    best_selection_key: Optional[Tuple[float, ...]] = None
    for scored_candidate in eligible_candidates:
        metric_values = cast(
            Tuple[float, ...], scored_candidate.metadata["weighted_selection_metric_values"]
        )
        regrets = tuple(
            max(0.0, (value - best_value) / max(scale, 1e-12))
            for value, best_value, scale in zip(
                metric_values,
                weighted_selection_context.best_values,
                weighted_selection_context.scales,
            )
        )
        weighted_components = tuple(
            regret * weight for regret, weight in zip(regrets, weighted_selection_context.weights)
        )
        weighted_score = sum(weighted_components)
        scored_candidate.metadata["weighted_selection_fields"] = weighted_selection_context.field_names
        scored_candidate.metadata["weighted_selection_best_values"] = weighted_selection_context.best_values
        scored_candidate.metadata["weighted_selection_weights"] = weighted_selection_context.weights
        scored_candidate.metadata["weighted_selection_scales"] = weighted_selection_context.scales
        scored_candidate.metadata["weighted_selection_regrets"] = regrets
        scored_candidate.metadata["weighted_selection_components"] = weighted_components
        scored_candidate.metadata["weighted_selection_score"] = weighted_score
        raw_selection_key = cast(Tuple[float, ...], scored_candidate.metadata["raw_selection_key"])
        selection_key = (weighted_score, *raw_selection_key)
        scored_candidate.metadata["selection_key"] = selection_key
        if best_selection_key is not None and selection_key >= best_selection_key:
            continue
        best_selection_key = selection_key
        best_candidate = scored_candidate
    return best_candidate


__all__ = [
    "_build_metal_site_environment_profile",
    "_build_organic_electronic_state_selection_context",
    "_charge_localization_penalty_for_atom",
    "_organic_score_bucket_index",
    "_passes_organic_force_field_guard",
    "_radical_localization_penalty_for_atom",
    "_score_candidate_with_no_metal_state",
    "select_best_candidate",
]
