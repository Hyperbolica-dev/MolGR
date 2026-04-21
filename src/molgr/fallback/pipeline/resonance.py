"""Radical resonance search and post-processing for the fallback pipeline.

The v2 mainline only keeps two layers of behavior here:
1. Enumerate resonance states with a bounded search policy.
2. Normalize each resonance state with charge/radical cleanup before scoring.

Older experimental traversal policies and duplicate collection wrappers have been
removed so the module matches the current production path more directly.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from heapq import heappop, heappush
from typing import Callable, Deque, Dict, List, Optional, Protocol, Sequence, Set, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.stages.clean import clean_neighbor_radicals, clean_resonances
from molgr.fallback.stages.eliminate import (
    eliminate_1_3_dipole,
    eliminate_negative_charges,
    eliminate_positive_charges,
)
from molgr.fallback.state import OmolStateMachine
from molgr.fallback.utils import smarts
from molgr.fallback.utils.scoring import calculate_radical_penalty, get_deviation_score


ResonanceAtomKey = Tuple[int, int, int, bool]
ResonanceBondKey = Tuple[int, int, int, bool]
ResonanceStateKey = Tuple[Tuple[ResonanceAtomKey, ...], Tuple[ResonanceBondKey, ...]]
ProcessedResonanceKey = str
ResonanceBondIndexMap = Dict[Tuple[int, int], int]
DirectGainMetrics = Tuple[float, float, float, float]
_ZERO_DIRECT_GAIN_METRICS: DirectGainMetrics = (0.0, 0.0, 0.0, 0.0)
_DIRECT_GAIN_CONJUGATION_SCORE_WEIGHT = 4.0
_DIRECT_GAIN_DEVIATION_SCORE_WEIGHT = 10.0
_DIRECT_GAIN_BRANCH_BOUND_STEP_SLACK = 10.0


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
class _LimitedDiscrepancyDirectGainTraversalPolicy:
    max_discrepancy: int
    fallback_to_full_frontier: bool

    def __call__(
        self,
        context: ResonanceTraversalContext,
        moves: Sequence[ResonanceTraversalMove],
    ) -> Optional[Sequence[ResonanceTraversalMove]]:
        positive_moves = _order_direct_gain_moves(
            context.current_omol,
            moves,
            positive_only=True,
        )
        if positive_moves:
            return [move for move, _metrics in positive_moves]

        if self.fallback_to_full_frontier:
            return [
                move
                for move, _metrics in _order_direct_gain_moves(
                    context.current_omol,
                    moves,
                    positive_only=False,
                )
            ]
        return ()


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
        bond_keys.append(
            (begin_idx, end_idx, bond.GetBondOrder(), bool(bond.IsAromatic()))
        )
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


def _build_bond_index_map_from_state_key(state_key: ResonanceStateKey) -> ResonanceBondIndexMap:
    return {
        (begin_idx, end_idx): idx
        for idx, (begin_idx, end_idx, _bond_order, _is_aromatic) in enumerate(state_key[1])
    }


def _bond_pair(atom_idx_a: int, atom_idx_b: int) -> Tuple[int, int]:
    if atom_idx_a <= atom_idx_b:
        return atom_idx_a, atom_idx_b
    return atom_idx_b, atom_idx_a


def _get_bond_order_with_overrides(
    bond: ob.OBBond,
    bond_order_overrides: Dict[Tuple[int, int], int],
) -> int:
    begin_idx = cast(ob.OBAtom, bond.GetBeginAtom()).GetIdx()
    end_idx = cast(ob.OBAtom, bond.GetEndAtom()).GetIdx()
    return bond_order_overrides.get(_bond_pair(begin_idx, end_idx), cast(int, bond.GetBondOrder()))


def _conjugated_bond_kind(
    bond: ob.OBBond,
    bond_order_overrides: Dict[Tuple[int, int], int],
) -> int:
    if bond.IsAromatic():
        return 2
    bond_order = _get_bond_order_with_overrides(bond, bond_order_overrides)
    if bond_order <= 0:
        return 0
    if bond_order == 1:
        return 1
    return 2


def _estimate_radical_conjugation_size(
    omol: pybel.Molecule,
    radical_atom_idx: int,
    bond_order_overrides: Optional[Dict[Tuple[int, int], int]] = None,
) -> int:
    obmol = cast(ob.OBMol, omol.OBMol)
    overrides = {} if bond_order_overrides is None else bond_order_overrides
    visited_states: Set[Tuple[int, int]] = {(radical_atom_idx, 0)}
    visited_atoms = {radical_atom_idx}
    frontier: List[Tuple[int, int]] = [(radical_atom_idx, 0)]

    while frontier:
        atom_idx, previous_bond_kind = frontier.pop()
        atom = cast(ob.OBAtom, obmol.GetAtom(atom_idx))
        for bond in ob.OBAtomBondIter(atom):
            bond_kind = _conjugated_bond_kind(cast(ob.OBBond, bond), overrides)
            if bond_kind == 0:
                continue
            if previous_bond_kind != 0 and bond_kind == previous_bond_kind:
                continue
            begin_idx = cast(ob.OBAtom, bond.GetBeginAtom()).GetIdx()
            end_idx = cast(ob.OBAtom, bond.GetEndAtom()).GetIdx()
            neighbor_idx = end_idx if begin_idx == atom_idx else begin_idx
            state = (neighbor_idx, bond_kind)
            if state in visited_states:
                continue
            visited_states.add(state)
            visited_atoms.add(neighbor_idx)
            frontier.append(state)

    return len(visited_atoms)


def _calculate_all_double_bond_carbon_bonus(
    omol: pybel.Molecule,
    atom_idx: int,
    bond_order_overrides: Optional[Dict[Tuple[int, int], int]] = None,
) -> float:
    obmol = cast(ob.OBMol, omol.OBMol)
    atom = cast(ob.OBAtom, obmol.GetAtom(atom_idx))
    if atom.GetAtomicNum() != 6:
        return 0.0

    saw_bond = False
    overrides = {} if bond_order_overrides is None else bond_order_overrides
    for bond in ob.OBAtomBondIter(atom):
        saw_bond = True
        if _get_bond_order_with_overrides(cast(ob.OBBond, bond), overrides) != 2:
            return 0.0
    return 5.0 if saw_bond else 0.0


def _compute_direct_gain_resonance_metrics(
    omol: pybel.Molecule,
    move_path: Tuple[int, int, int],
) -> DirectGainMetrics:
    old_radical_idx, center_idx, new_radical_idx = move_path
    obmol = cast(ob.OBMol, omol.OBMol)
    bond_old_center = cast(ob.OBBond, obmol.GetBond(old_radical_idx, center_idx))
    bond_center_new = cast(ob.OBBond, obmol.GetBond(center_idx, new_radical_idx))
    bond_order_overrides = {
        _bond_pair(old_radical_idx, center_idx): cast(int, bond_old_center.GetBondOrder()) + 1,
        _bond_pair(center_idx, new_radical_idx): cast(int, bond_center_new.GetBondOrder()) - 1,
    }

    old_atom = cast(ob.OBAtom, obmol.GetAtom(old_radical_idx))
    new_atom = cast(ob.OBAtom, obmol.GetAtom(new_radical_idx))

    conjugation_gain = _estimate_radical_conjugation_size(
        omol,
        new_radical_idx,
        bond_order_overrides,
    ) - _estimate_radical_conjugation_size(omol, old_radical_idx)
    deviation_gain = get_deviation_score(omol, old_radical_idx) - get_deviation_score(omol, new_radical_idx)
    radical_penalty_gain = calculate_radical_penalty(old_atom) - calculate_radical_penalty(new_atom)
    double_bond_bonus_gain = _calculate_all_double_bond_carbon_bonus(
        omol,
        new_radical_idx,
        bond_order_overrides,
    ) - _calculate_all_double_bond_carbon_bonus(omol, old_radical_idx)
    return conjugation_gain, deviation_gain, radical_penalty_gain, double_bond_bonus_gain


def _has_positive_direct_gain(metrics: DirectGainMetrics) -> bool:
    return any(metric > 0 for metric in metrics)


def _direct_gain_move_sort_key(
    move_path: Tuple[int, int, int],
    metrics: DirectGainMetrics,
) -> Tuple[int, float, float, float, float, Tuple[int, int, int]]:
    return (
        0 if _has_positive_direct_gain(metrics) else 1,
        -metrics[0],
        -metrics[1],
        -metrics[2],
        -metrics[3],
        move_path,
    )


def _order_direct_gain_moves(
    omol: pybel.Molecule,
    moves: Sequence[ResonanceTraversalMove],
    *,
    positive_only: bool,
) -> List[Tuple[ResonanceTraversalMove, DirectGainMetrics]]:
    ordered_moves: List[Tuple[ResonanceTraversalMove, DirectGainMetrics]] = []
    for move in moves:
        metrics = _compute_direct_gain_resonance_metrics(omol, move.path)
        if positive_only and not _has_positive_direct_gain(metrics):
            continue
        ordered_moves.append((move, metrics))
    ordered_moves.sort(key=lambda item: _direct_gain_move_sort_key(item[0].path, item[1]))
    return ordered_moves


def _add_direct_gain_metrics(
    left: DirectGainMetrics,
    right: DirectGainMetrics,
) -> DirectGainMetrics:
    return tuple(left[idx] + right[idx] for idx in range(len(left)))  # type: ignore[return-value]


def _componentwise_max_direct_gain_metrics(
    left: DirectGainMetrics,
    right: DirectGainMetrics,
) -> DirectGainMetrics:
    return tuple(max(left[idx], right[idx]) for idx in range(len(left)))  # type: ignore[return-value]


def _compute_n_step_direct_gain_upper_bound(
    omol: pybel.Molecule,
    state_key: ResonanceStateKey,
    remaining_steps: int,
    cache: Dict[Tuple[ResonanceStateKey, int], DirectGainMetrics],
) -> DirectGainMetrics:
    if remaining_steps <= 0:
        return _ZERO_DIRECT_GAIN_METRICS

    cache_key = (state_key, remaining_steps)
    cached = cache.get(cache_key)
    if cached is not None:
        return cached

    upper_bound = _ZERO_DIRECT_GAIN_METRICS
    bond_index_map = _build_bond_index_map_from_state_key(state_key)
    for move in _enumerate_one_step_resonance_moves(omol, state_key, bond_index_map):
        direct_metrics = _compute_direct_gain_resonance_metrics(omol, move.idxs)
        if not _has_positive_direct_gain(direct_metrics):
            continue
        total_metrics = direct_metrics
        if remaining_steps > 1:
            next_omol = _materialize_one_step_resonance(omol, move.idxs)
            future_upper_bound = _compute_n_step_direct_gain_upper_bound(
                next_omol,
                move.next_state_key,
                remaining_steps - 1,
                cache,
            )
            total_metrics = _add_direct_gain_metrics(direct_metrics, future_upper_bound)
        upper_bound = _componentwise_max_direct_gain_metrics(upper_bound, total_metrics)

    cache[cache_key] = upper_bound
    return upper_bound


def _estimate_direct_gain_score_improvement_upper_bound(
    metrics: DirectGainMetrics,
    *,
    remaining_steps: int,
) -> float:
    # The direct-gain metrics only cover the dominant local score terms. Add a per-step
    # slack so branch pruning stays conservative with respect to process_resonance side
    # effects and any smaller score terms that are not modeled here.
    return max(
        0.0,
        metrics[0] * _DIRECT_GAIN_CONJUGATION_SCORE_WEIGHT
        + metrics[1] * _DIRECT_GAIN_DEVIATION_SCORE_WEIGHT
        + metrics[2]
        + metrics[3]
        + max(remaining_steps, 0) * _DIRECT_GAIN_BRANCH_BOUND_STEP_SLACK,
    )


def estimate_remaining_resonance_score_improvement_upper_bound(
    omol: pybel.Molecule,
    state_key: ResonanceStateKey,
    remaining_steps: int,
    cache: Optional[Dict[Tuple[ResonanceStateKey, int], DirectGainMetrics]] = None,
) -> float:
    """Estimate a conservative upper bound for the score gain still reachable from a node."""

    if remaining_steps <= 0:
        return 0.0
    optimistic_metrics = _compute_n_step_direct_gain_upper_bound(
        omol,
        state_key,
        remaining_steps,
        {} if cache is None else cache,
    )
    return _estimate_direct_gain_score_improvement_upper_bound(
        optimistic_metrics,
        remaining_steps=remaining_steps,
    )


def make_limited_discrepancy_direct_gain_traversal_policy(
    *,
    max_discrepancy: int = 2,
    fallback_to_full_frontier: bool = True,
) -> ResonanceTraversalPolicy:
    """Prefer direct-gain moves first, but allow a bounded number of lower-ranked detours."""

    return _LimitedDiscrepancyDirectGainTraversalPolicy(
        max_discrepancy=max(0, max_discrepancy),
        fallback_to_full_frontier=fallback_to_full_frontier,
    )


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


def get_radical_resonances(
    omol: pybel.Molecule,
    max_depth: int = 2,
    *,
    traversal_policy: Optional[ResonanceTraversalPolicy] = None,
) -> List[pybel.Molecule]:
    """Collect raw resonance states reachable within the configured search depth."""

    resonances: List[pybel.Molecule] = []

    def _collect(node: ResonanceSearchNode) -> bool:
        resonances.append(node.omol)
        return True

    walk_radical_resonances(
        omol,
        max_depth=max_depth,
        traversal_policy=traversal_policy,
        visit=_collect,
    )
    return resonances


def walk_radical_resonances(
    omol: pybel.Molecule,
    max_depth: int = 2,
    *,
    traversal_policy: Optional[ResonanceTraversalPolicy] = None,
    visit: Optional[Callable[[ResonanceSearchNode], bool]] = None,
) -> None:
    """Traverse resonance states and let the caller decide whether each node should expand."""

    if isinstance(traversal_policy, _LimitedDiscrepancyDirectGainTraversalPolicy):
        _walk_radical_resonances_limited_discrepancy(
            omol,
            max_depth=max_depth,
            traversal_policy=traversal_policy,
            visit=visit,
        )
        return

    _walk_radical_resonances_bfs(
        omol,
        max_depth=max_depth,
        traversal_policy=traversal_policy,
        visit=visit,
    )


def _walk_radical_resonances_bfs(
    omol: pybel.Molecule,
    *,
    max_depth: int,
    traversal_policy: Optional[ResonanceTraversalPolicy],
    visit: Optional[Callable[[ResonanceSearchNode], bool]],
) -> None:
    visitor = visit if visit is not None else (lambda _node: True)
    root_key, bond_index_map = _build_resonance_search_context(omol)
    seen = {root_key}
    frontier: Deque[Tuple[pybel.Molecule, ResonanceStateKey, int]] = deque([(omol, root_key, 0)])

    while frontier:
        current, current_key, depth = frontier.popleft()
        should_expand = visitor(
            ResonanceSearchNode(
                omol=current,
                state_key=current_key,
                depth=depth,
            )
        )
        if depth >= max_depth or not should_expand:
            continue

        indexed_moves = _enumerate_one_step_resonance_moves(current, current_key, bond_index_map)
        selected_moves = _apply_resonance_traversal_policy(
            ResonanceTraversalContext(
                root_omol=omol,
                current_omol=current,
                current_state_key=current_key,
                depth=depth,
                max_depth=max_depth,
            ),
            indexed_moves,
            traversal_policy,
        )
        for move in selected_moves:
            state_key = move.next_state_key
            if state_key in seen:
                continue
            seen.add(state_key)
            new_resonance = _materialize_one_step_resonance(current, move.idxs)
            frontier.append((new_resonance, state_key, depth + 1))


def _walk_radical_resonances_limited_discrepancy(
    omol: pybel.Molecule,
    *,
    max_depth: int,
    traversal_policy: _LimitedDiscrepancyDirectGainTraversalPolicy,
    visit: Optional[Callable[[ResonanceSearchNode], bool]],
) -> None:
    visitor = visit if visit is not None else (lambda _node: True)
    root_key, bond_index_map = _build_resonance_search_context(omol)
    best_discrepancy_by_state = {root_key: 0}
    emitted_states = set()
    frontier: List[Tuple[int, int, int, pybel.Molecule, ResonanceStateKey]] = [
        (0, 0, 0, omol, root_key)
    ]
    push_order = 0

    while frontier:
        discrepancy, depth, _order, current, current_key = heappop(frontier)
        if discrepancy != best_discrepancy_by_state.get(current_key):
            continue
        if current_key in emitted_states:
            should_expand = True
        else:
            emitted_states.add(current_key)
            should_expand = visitor(
                ResonanceSearchNode(
                    omol=current,
                    state_key=current_key,
                    depth=depth,
                )
            )
        if depth >= max_depth or not should_expand:
            continue

        indexed_moves = _enumerate_one_step_resonance_moves(current, current_key, bond_index_map)
        selected_moves = _apply_resonance_traversal_policy(
            ResonanceTraversalContext(
                root_omol=omol,
                current_omol=current,
                current_state_key=current_key,
                depth=depth,
                max_depth=max_depth,
            ),
            indexed_moves,
            traversal_policy,
        )
        for move_rank, move in enumerate(selected_moves):
            next_discrepancy = discrepancy + move_rank
            if next_discrepancy > traversal_policy.max_discrepancy:
                break
            state_key = move.next_state_key
            best_known_discrepancy = best_discrepancy_by_state.get(state_key)
            if (
                best_known_discrepancy is not None
                and best_known_discrepancy <= next_discrepancy
            ):
                continue
            best_discrepancy_by_state[state_key] = next_discrepancy
            next_omol = _materialize_one_step_resonance(current, move.idxs)
            push_order += 1
            heappush(
                frontier,
                (next_discrepancy, depth + 1, push_order, next_omol, state_key),
            )


__all__ = [
    "ResonanceSearchNode",
    "ResonanceTraversalContext",
    "ResonanceTraversalMove",
    "ResonanceTraversalPolicy",
    "build_processed_resonance_key",
    "build_resonance_state_key",
    "estimate_remaining_resonance_score_improvement_upper_bound",
    "get_radical_resonances",
    "make_limited_discrepancy_direct_gain_traversal_policy",
    "process_resonance",
    "walk_radical_resonances",
]
