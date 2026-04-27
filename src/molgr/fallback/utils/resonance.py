"""Resonance helpers shared by fallback pipeline orchestration."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Protocol, Sequence, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.config import MolGRConfig, resolve_config
from molgr.fallback.stages.clean import clean_neighbor_radicals, clean_resonances
from molgr.fallback.stages.eliminate import (
    eliminate_1_3_dipole,
    eliminate_negative_charges,
    eliminate_positive_charges,
)
from molgr.fallback.state import OmolStateMachine
from molgr.fallback.utils import consts
from molgr.fallback.utils.force_field import OmolForceFieldContext, selection_force_field_energy
from molgr.fallback.utils.tools import typed_lru_cache

from . import smarts


ResonanceAtomKey = Tuple[int, int, int, bool]
ResonanceBondKey = Tuple[int, int, int, bool]
ResonanceStateKey = Tuple[Tuple[ResonanceAtomKey, ...], Tuple[ResonanceBondKey, ...]]
ProcessedResonanceKey = str
ResonanceBondIndexMap = Dict[Tuple[int, int], int]
_DEFAULT_RESONANCE_MOVE_SCORE_CACHE_MAXSIZE = 4096
_DIRECT_GAIN_CONJUGATION_SCORE_WEIGHT = 4.0
_DIRECT_GAIN_DEVIATION_SCORE_WEIGHT = 10.0

_BondOrderOverrides = Dict[Tuple[int, int], int]
_DirectGainMetrics = Tuple[float, float, float, float]


@dataclass(frozen=True)
class ResonanceTraversalContext:
    """Read-only state exposed to traversal policies for the current expansion step."""

    root_omol: pybel.Molecule
    current_omol: pybel.Molecule
    current_state_key: ResonanceStateKey
    depth: int
    max_depth: int


@dataclass(frozen=True)
class ResonanceTraversalMove:
    """A one-step radical resonance move and the incremental key of its output state."""

    path: Tuple[int, int, int]
    next_state_key: ResonanceStateKey


@dataclass(frozen=True)
class ResonanceSearchNode:
    """A resonance state yielded during traversal."""

    omol: pybel.Molecule
    state_key: ResonanceStateKey
    depth: int


@dataclass(frozen=True)
class _IndexedResonanceTraversalMove:
    idxs: Tuple[int, int, int]
    next_state_key: ResonanceStateKey


class ResonanceTraversalPolicy(Protocol):
    def __call__(
        self,
        context: ResonanceTraversalContext,
        moves: Sequence[ResonanceTraversalMove],
    ) -> Optional[Sequence[ResonanceTraversalMove]]: ...


@dataclass(frozen=True)
class _LimitedDiscrepancyForceFieldTraversalPolicy:
    max_discrepancy: int
    config: MolGRConfig | None = None

    def __call__(
        self,
        context: ResonanceTraversalContext,
        moves: Sequence[ResonanceTraversalMove],
    ) -> Optional[Sequence[ResonanceTraversalMove]]:
        ordered_moves = _order_force_field_moves(
            context.current_omol,
            moves,
            config=self.config,
        )
        return [move for move, _score in ordered_moves]


@dataclass(frozen=True)
class _LimitedDiscrepancyDirectGainTraversalPolicy:
    max_discrepancy: int

    def __call__(
        self,
        context: ResonanceTraversalContext,
        moves: Sequence[ResonanceTraversalMove],
    ) -> Optional[Sequence[ResonanceTraversalMove]]:
        ordered_moves = _order_direct_gain_moves(context.current_omol, moves)
        return [move for move, _score in ordered_moves]


@dataclass(frozen=True)
class _LimitedDiscrepancyInputOrderTraversalPolicy:
    max_discrepancy: int

    def __call__(
        self,
        context: ResonanceTraversalContext,
        moves: Sequence[ResonanceTraversalMove],
    ) -> Optional[Sequence[ResonanceTraversalMove]]:
        return list(moves)


_LIMITED_DISCREPANCY_POLICY_TYPES = (
    _LimitedDiscrepancyForceFieldTraversalPolicy,
    _LimitedDiscrepancyDirectGainTraversalPolicy,
    _LimitedDiscrepancyInputOrderTraversalPolicy,
)


def process_resonance(
    resonance: pybel.Molecule,
    charge: int,
) -> tuple[pybel.Molecule, int, bool]:
    """Normalize one resonance candidate before validation and scoring."""

    machine = OmolStateMachine(resonance, charge)
    hit = machine.run_omol_charge_stage(None, eliminate_1_3_dipole)
    hit = machine.run_omol_charge_stage(None, eliminate_positive_charges) or hit
    hit = machine.run_omol_charge_stage(None, eliminate_negative_charges) or hit
    hit = machine.run_omol_stage(None, clean_neighbor_radicals) or hit
    hit = machine.run_omol_stage(None, clean_resonances) or hit
    return machine.omol, machine.given_charge, hit


def build_resonance_state_key(omol: pybel.Molecule) -> ResonanceStateKey:
    """Build the structural key used to deduplicate raw resonance states."""

    obmol = cast(ob.OBMol, omol.OBMol)
    atom_keys: List[ResonanceAtomKey] = []
    for atom in ob.OBMolAtomIter(obmol):
        atom = cast(ob.OBAtom, atom)
        atom_keys.append(
            (
                atom.GetAtomicNum(),
                atom.GetFormalCharge(),
                atom.GetSpinMultiplicity(),
                bool(atom.IsAromatic()),
            )
        )

    bond_keys: List[ResonanceBondKey] = []
    for bond in ob.OBMolBondIter(obmol):
        begin_idx = cast(ob.OBAtom, bond.GetBeginAtom()).GetIdx()
        end_idx = cast(ob.OBAtom, bond.GetEndAtom()).GetIdx()
        if begin_idx > end_idx:
            begin_idx, end_idx = end_idx, begin_idx
        bond_keys.append((begin_idx, end_idx, bond.GetBondOrder(), bool(bond.IsAromatic())))
    return tuple(atom_keys), tuple(bond_keys)


def build_processed_resonance_key(omol: pybel.Molecule) -> ProcessedResonanceKey:
    """Build the dedup key after `process_resonance` has normalized the candidate."""

    return omol.write("molreport") or ""


def _build_resonance_search_context(
    omol: pybel.Molecule,
) -> Tuple[ResonanceStateKey, ResonanceBondIndexMap]:
    obmol = cast(ob.OBMol, omol.OBMol)
    atom_keys: List[ResonanceAtomKey] = []
    for atom in ob.OBMolAtomIter(obmol):
        atom = cast(ob.OBAtom, atom)
        atom_keys.append(
            (
                atom.GetAtomicNum(),
                atom.GetFormalCharge(),
                atom.GetSpinMultiplicity(),
                bool(atom.IsAromatic()),
            )
        )

    bond_index_map: ResonanceBondIndexMap = {}
    bond_keys: List[ResonanceBondKey] = []
    for bond in ob.OBMolBondIter(obmol):
        begin_idx = cast(ob.OBAtom, bond.GetBeginAtom()).GetIdx()
        end_idx = cast(ob.OBAtom, bond.GetEndAtom()).GetIdx()
        if begin_idx > end_idx:
            begin_idx, end_idx = end_idx, begin_idx
        bond_index_map[(begin_idx, end_idx)] = len(bond_keys)
        bond_keys.append((begin_idx, end_idx, bond.GetBondOrder(), bool(bond.IsAromatic())))
    return (tuple(atom_keys), tuple(bond_keys)), bond_index_map


def _increment_resonance_state_key(
    state_key: ResonanceStateKey,
    bond_index_map: ResonanceBondIndexMap,
    idxs: Tuple[int, int, int],
) -> ResonanceStateKey:
    atom_keys, bond_keys = state_key
    atom_key_list = list(atom_keys)
    bond_key_list = list(bond_keys)

    atom1_idx = idxs[0] - 1
    atom3_idx = idxs[2] - 1

    atom1 = atom_key_list[atom1_idx]
    atom3 = atom_key_list[atom3_idx]
    atom_key_list[atom1_idx] = (atom1[0], atom1[1], atom1[2] - 1, atom1[3])
    atom_key_list[atom3_idx] = (atom3[0], atom3[1], atom3[2] + 1, atom3[3])

    bond1_pair = (min(idxs[0], idxs[1]), max(idxs[0], idxs[1]))
    bond2_pair = (min(idxs[1], idxs[2]), max(idxs[1], idxs[2]))
    bond1_idx = bond_index_map[bond1_pair]
    bond2_idx = bond_index_map[bond2_pair]
    bond1 = bond_key_list[bond1_idx]
    bond2 = bond_key_list[bond2_idx]
    bond_key_list[bond1_idx] = (bond1[0], bond1[1], bond1[2] + 1, bond1[3])
    bond_key_list[bond2_idx] = (bond2[0], bond2[1], bond2[2] - 1, bond2[3])

    return tuple(atom_key_list), tuple(bond_key_list)


def _enumerate_one_step_resonance_moves(
    omol: pybel.Molecule,
    state_key: ResonanceStateKey,
    bond_index_map: ResonanceBondIndexMap,
) -> List[_IndexedResonanceTraversalMove]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int]] = list(smarts.RESONANCE_ONE_STEP.findall(omol))
    result: List[_IndexedResonanceTraversalMove] = []
    for idxs in res:
        atom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        bond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        bond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if (
            atom1.GetSpinMultiplicity() == 1
            and atom3.GetSpinMultiplicity() == 0
            and bond1.GetBondOrder() <= 2
            and bond2.GetBondOrder() >= 2
        ):
            result.append(
                _IndexedResonanceTraversalMove(
                    idxs=idxs,
                    next_state_key=_increment_resonance_state_key(state_key, bond_index_map, idxs),
                )
            )
    return result


def _materialize_one_step_resonance(
    omol: pybel.Molecule,
    idxs: Tuple[int, int, int],
) -> pybel.Molecule:
    new_omol = omol.clone
    new_obmol = cast(ob.OBMol, new_omol.OBMol)
    atom1_clone = cast(ob.OBAtom, new_obmol.GetAtom(idxs[0]))
    atom3_clone = cast(ob.OBAtom, new_obmol.GetAtom(idxs[2]))
    bond1_clone = cast(ob.OBBond, new_obmol.GetBond(idxs[0], idxs[1]))
    bond2_clone = cast(ob.OBBond, new_obmol.GetBond(idxs[1], idxs[2]))
    bond1_clone.SetBondOrder(bond1_clone.GetBondOrder() + 1)
    bond2_clone.SetBondOrder(bond2_clone.GetBondOrder() - 1)
    atom1_clone.SetSpinMultiplicity(atom1_clone.GetSpinMultiplicity() - 1)
    atom3_clone.SetSpinMultiplicity(atom3_clone.GetSpinMultiplicity() + 1)
    return new_omol


@typed_lru_cache(maxsize=_DEFAULT_RESONANCE_MOVE_SCORE_CACHE_MAXSIZE, typed=True)
def _score_one_step_resonance_with_force_field_cached(
    context: OmolForceFieldContext,
    move_path: Tuple[int, int, int],
    config: MolGRConfig,
) -> float:
    moved_omol = _materialize_one_step_resonance(context.omol, move_path)
    try:
        return selection_force_field_energy(moved_omol, config=config)
    except ValueError:
        return float("inf")


def _score_one_step_resonance_with_force_field(
    omol: pybel.Molecule,
    move_path: Tuple[int, int, int],
    *,
    config: MolGRConfig | None = None,
) -> float:
    return _score_one_step_resonance_with_force_field_cached(
        OmolForceFieldContext(omol),
        move_path,
        resolve_config(config),
    )


def _bond_pair(atom_idx_a: int, atom_idx_b: int) -> Tuple[int, int]:
    if atom_idx_a <= atom_idx_b:
        return atom_idx_a, atom_idx_b
    return atom_idx_b, atom_idx_a


def _get_bond_order_with_overrides(
    bond: ob.OBBond,
    bond_order_overrides: _BondOrderOverrides,
) -> int:
    override = bond_order_overrides.get(
        _bond_pair(
            cast(int, cast(ob.OBAtom, bond.GetBeginAtom()).GetIdx()),
            cast(int, cast(ob.OBAtom, bond.GetEndAtom()).GetIdx()),
        )
    )
    if override is not None:
        return override
    return cast(int, bond.GetBondOrder())


def _conjugated_bond_kind(
    bond: ob.OBBond,
    bond_order_overrides: _BondOrderOverrides,
) -> int:
    if bond.IsAromatic():
        return 2
    bond_order = _get_bond_order_with_overrides(bond, bond_order_overrides)
    if bond_order <= 0:
        return 0
    if bond_order == 1:
        return 1
    return 2


def _calculate_radical_penalty(atom: ob.OBAtom) -> float:
    radical_num = cast(int, atom.GetSpinMultiplicity())
    if radical_num == 0:
        return 0.0
    if cast(int, atom.GetAtomicNum()) in consts.HETEROATOM:
        return float(radical_num) * 10.0
    return (3.0 - float(atom.GetHvyDegree())) * 1.5


def _estimate_radical_conjugation_size(
    obmol: ob.OBMol,
    radical_atom_idx: int,
    bond_order_overrides: Optional[_BondOrderOverrides] = None,
) -> int:
    overrides = bond_order_overrides or {}
    visited_states = {(radical_atom_idx, 0)}
    visited_atoms = {radical_atom_idx}
    frontier = [(radical_atom_idx, 0)]

    while frontier:
        atom_idx, previous_bond_kind = frontier.pop()
        atom = cast(ob.OBAtom, obmol.GetAtom(atom_idx))
        if atom is None:
            continue
        for bond_iter in ob.OBAtomBondIter(atom):
            bond = cast(ob.OBBond, bond_iter)
            bond_kind = _conjugated_bond_kind(bond, overrides)
            if bond_kind == 0:
                continue
            if previous_bond_kind != 0 and bond_kind == previous_bond_kind:
                continue
            neighbor = cast(ob.OBAtom, bond.GetNbrAtom(atom))
            neighbor_idx = cast(int, neighbor.GetIdx())
            state = (neighbor_idx, bond_kind)
            if state in visited_states:
                continue
            visited_states.add(state)
            visited_atoms.add(neighbor_idx)
            frontier.append(state)

    return len(visited_atoms)


def _calculate_all_double_bond_carbon_bonus(
    obmol: ob.OBMol,
    atom_idx: int,
    bond_order_overrides: Optional[_BondOrderOverrides] = None,
) -> float:
    atom = cast(ob.OBAtom, obmol.GetAtom(atom_idx))
    if atom is None or cast(int, atom.GetAtomicNum()) != 6:
        return 0.0

    overrides = bond_order_overrides or {}
    saw_bond = False
    for bond_iter in ob.OBAtomBondIter(atom):
        saw_bond = True
        bond = cast(ob.OBBond, bond_iter)
        if _get_bond_order_with_overrides(bond, overrides) != 2:
            return 0.0
    return 5.0 if saw_bond else 0.0


def _atom_coords(atom: ob.OBAtom) -> Tuple[float, float, float]:
    return float(atom.GetX()), float(atom.GetY()), float(atom.GetZ())


def _vector_sub(
    left: Tuple[float, float, float],
    right: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    return left[0] - right[0], left[1] - right[1], left[2] - right[2]


def _vector_dot(
    left: Tuple[float, float, float],
    right: Tuple[float, float, float],
) -> float:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def _vector_cross(
    left: Tuple[float, float, float],
    right: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _vector_length_sq(vector: Tuple[float, float, float]) -> float:
    return _vector_dot(vector, vector)


def _calculate_tetrahedron_volume(
    p1: Tuple[float, float, float],
    p2: Tuple[float, float, float],
    p3: Tuple[float, float, float],
    p4: Tuple[float, float, float],
) -> float:
    v1 = _vector_sub(p1, p4)
    v2 = _vector_sub(p2, p4)
    v3 = _vector_sub(p3, p4)
    return abs(_vector_dot(v1, _vector_cross(v2, v3))) / 6.0


def _calculate_shape_quality(
    p1: Tuple[float, float, float],
    p2: Tuple[float, float, float],
    p3: Tuple[float, float, float],
    p4: Tuple[float, float, float],
) -> float:
    volume = _calculate_tetrahedron_volume(p1, p2, p3, p4)
    if volume < 1e-9:
        return 0.0

    edges_sq_sum = 0.0
    edges_sq_sum += _vector_length_sq(_vector_sub(p1, p2))
    edges_sq_sum += _vector_length_sq(_vector_sub(p1, p3))
    edges_sq_sum += _vector_length_sq(_vector_sub(p1, p4))
    edges_sq_sum += _vector_length_sq(_vector_sub(p2, p3))
    edges_sq_sum += _vector_length_sq(_vector_sub(p2, p4))
    edges_sq_sum += _vector_length_sq(_vector_sub(p3, p4))

    l_rms_squared = edges_sq_sum / 6.0
    l_rms_cubed = math.pow(l_rms_squared, 1.5)
    if l_rms_cubed < 1e-9:
        return 0.0

    quality = 6.0 * math.sqrt(2.0) * (volume / l_rms_cubed)
    return max(0.0, min(1.0, quality))


def _get_deviation_score(obmol: ob.OBMol, atom: Optional[ob.OBAtom]) -> float:
    if atom is None:
        return 0.0

    neighbors = [cast(ob.OBAtom, neighbor) for neighbor in ob.OBAtomAtomIter(atom)]
    if len(neighbors) == 2:
        angle = float(obmol.GetAngle(neighbors[0], atom, neighbors[1]))
        return abs(angle - 108.0) / 108.0
    if len(neighbors) == 3:
        return 1.0 - _calculate_shape_quality(
            _atom_coords(neighbors[0]),
            _atom_coords(neighbors[1]),
            _atom_coords(neighbors[2]),
            _atom_coords(atom),
        )
    return 0.0


def _compute_direct_gain_resonance_metrics(
    omol: pybel.Molecule,
    move_path: Tuple[int, int, int],
) -> _DirectGainMetrics:
    obmol = cast(ob.OBMol, omol.OBMol)
    old_radical_idx, center_idx, new_radical_idx = move_path

    bond_old_center = cast(ob.OBBond, obmol.GetBond(old_radical_idx, center_idx))
    bond_center_new = cast(ob.OBBond, obmol.GetBond(center_idx, new_radical_idx))
    old_atom = cast(ob.OBAtom, obmol.GetAtom(old_radical_idx))
    new_atom = cast(ob.OBAtom, obmol.GetAtom(new_radical_idx))
    if bond_old_center is None or bond_center_new is None or old_atom is None or new_atom is None:
        return 0.0, 0.0, 0.0, 0.0

    bond_order_overrides = {
        _bond_pair(old_radical_idx, center_idx): cast(int, bond_old_center.GetBondOrder()) + 1,
        _bond_pair(center_idx, new_radical_idx): cast(int, bond_center_new.GetBondOrder()) - 1,
    }
    conjugation_gain = float(
        _estimate_radical_conjugation_size(obmol, new_radical_idx, bond_order_overrides)
        - _estimate_radical_conjugation_size(obmol, old_radical_idx)
    )
    deviation_gain = _get_deviation_score(obmol, old_atom) - _get_deviation_score(obmol, new_atom)
    radical_penalty_gain = _calculate_radical_penalty(old_atom) - _calculate_radical_penalty(new_atom)
    double_bond_bonus_gain = (
        _calculate_all_double_bond_carbon_bonus(obmol, new_radical_idx, bond_order_overrides)
        - _calculate_all_double_bond_carbon_bonus(obmol, old_radical_idx)
    )
    return (
        conjugation_gain,
        deviation_gain,
        radical_penalty_gain,
        double_bond_bonus_gain,
    )


def _direct_gain_move_score(metrics: _DirectGainMetrics) -> float:
    return -(
        metrics[0] * _DIRECT_GAIN_CONJUGATION_SCORE_WEIGHT
        + metrics[1] * _DIRECT_GAIN_DEVIATION_SCORE_WEIGHT
        + metrics[2]
        + metrics[3]
    )


def resonance_move_score_cache_info() -> Tuple[int, int, int]:
    info = _score_one_step_resonance_with_force_field_cached.cache_info()
    return info.hits, info.misses, info.currsize


def resonance_move_score_cache_clear() -> None:
    _score_one_step_resonance_with_force_field_cached.cache_clear()


def _order_force_field_moves(
    omol: pybel.Molecule,
    moves: Sequence[ResonanceTraversalMove],
    *,
    config: MolGRConfig | None = None,
) -> List[Tuple[ResonanceTraversalMove, float]]:
    ordered_moves: List[Tuple[ResonanceTraversalMove, float]] = []
    for move in moves:
        ordered_moves.append(
            (
                move,
                _score_one_step_resonance_with_force_field(
                    omol,
                    move.path,
                    config=config,
                ),
            )
        )
    ordered_moves.sort(key=lambda item: (item[1], item[0].path))
    return ordered_moves


def _order_direct_gain_moves(
    omol: pybel.Molecule,
    moves: Sequence[ResonanceTraversalMove],
) -> List[Tuple[ResonanceTraversalMove, float]]:
    ordered_moves: List[Tuple[ResonanceTraversalMove, float]] = []
    for move in moves:
        ordered_moves.append((move, _direct_gain_move_score(_compute_direct_gain_resonance_metrics(omol, move.path))))
    ordered_moves.sort(key=lambda item: (item[1], item[0].path))
    return ordered_moves


def make_limited_discrepancy_force_field_traversal_policy(
    *,
    max_discrepancy: int = 2,
    config: MolGRConfig | None = None,
) -> ResonanceTraversalPolicy:
    """Prefer lower one-step force-field-energy moves, with bounded discrepancy search."""

    return _LimitedDiscrepancyForceFieldTraversalPolicy(
        max_discrepancy=max(0, max_discrepancy),
        config=config,
    )


def make_limited_discrepancy_direct_gain_traversal_policy(
    *,
    max_discrepancy: int = 2,
) -> ResonanceTraversalPolicy:
    """Prefer larger direct-gain improvements, with bounded discrepancy search."""

    return _LimitedDiscrepancyDirectGainTraversalPolicy(
        max_discrepancy=max(0, max_discrepancy),
    )


def make_limited_discrepancy_input_order_traversal_policy(
    *,
    max_discrepancy: int = 2,
) -> ResonanceTraversalPolicy:
    """Use raw resonance move order, with bounded discrepancy search."""

    return _LimitedDiscrepancyInputOrderTraversalPolicy(
        max_discrepancy=max(0, max_discrepancy),
    )


def _apply_resonance_traversal_policy(
    context: ResonanceTraversalContext,
    moves: Sequence[_IndexedResonanceTraversalMove],
    traversal_policy: Optional[ResonanceTraversalPolicy],
) -> Sequence[_IndexedResonanceTraversalMove]:
    if traversal_policy is None or not moves:
        return moves

    public_moves = tuple(
        ResonanceTraversalMove(
            path=move.idxs,
            next_state_key=move.next_state_key,
        )
        for move in moves
    )
    selected_moves = traversal_policy(context, public_moves)
    if selected_moves is None:
        return moves

    move_by_key = {(move.idxs, move.next_state_key): move for move in moves}
    selected_indexed_moves: List[_IndexedResonanceTraversalMove] = []
    seen_selected_keys = set()
    for move in selected_moves:
        move_key = (move.path, move.next_state_key)
        if move_key in seen_selected_keys:
            continue
        indexed_move = move_by_key.get(move_key)
        if indexed_move is None:
            raise ValueError("resonance traversal policy returned a move outside the candidate frontier")
        seen_selected_keys.add(move_key)
        selected_indexed_moves.append(indexed_move)
    return selected_indexed_moves


__all__ = [
    "ProcessedResonanceKey",
    "ResonanceBondIndexMap",
    "ResonanceSearchNode",
    "ResonanceStateKey",
    "ResonanceTraversalContext",
    "ResonanceTraversalMove",
    "ResonanceTraversalPolicy",
    "build_processed_resonance_key",
    "build_resonance_state_key",
    "make_limited_discrepancy_direct_gain_traversal_policy",
    "make_limited_discrepancy_force_field_traversal_policy",
    "make_limited_discrepancy_input_order_traversal_policy",
    "process_resonance",
    "resonance_move_score_cache_clear",
    "resonance_move_score_cache_info",
]
