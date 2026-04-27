"""Radical resonance traversal flow for the fallback pipeline."""

from __future__ import annotations

from collections import deque
from heapq import heappop, heappush
from typing import Callable, Deque, List, Optional, Tuple

from openbabel import pybel

from molgr.fallback.utils import resonance as resonance_utils
from molgr.fallback.utils.resonance import (
    ProcessedResonanceKey,
    ResonanceBondIndexMap,
    ResonanceSearchNode,
    ResonanceStateKey,
    ResonanceTraversalContext,
    ResonanceTraversalMove,
    ResonanceTraversalPolicy,
    build_processed_resonance_key,
    build_resonance_state_key,
    make_limited_discrepancy_direct_gain_traversal_policy,
    make_limited_discrepancy_force_field_traversal_policy,
    make_limited_discrepancy_input_order_traversal_policy,
    process_resonance,
    resonance_move_score_cache_clear,
    resonance_move_score_cache_info,
)


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

    if isinstance(traversal_policy, resonance_utils._LIMITED_DISCREPANCY_POLICY_TYPES):
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
    root_key, bond_index_map = resonance_utils._build_resonance_search_context(omol)
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

        indexed_moves = resonance_utils._enumerate_one_step_resonance_moves(
            current,
            current_key,
            bond_index_map,
        )
        selected_moves = resonance_utils._apply_resonance_traversal_policy(
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
            if move.next_state_key in seen:
                continue
            seen.add(move.next_state_key)
            frontier.append(
                (
                    resonance_utils._materialize_one_step_resonance(current, move.idxs),
                    move.next_state_key,
                    depth + 1,
                )
            )


def _walk_radical_resonances_limited_discrepancy(
    omol: pybel.Molecule,
    *,
    max_depth: int,
    traversal_policy: ResonanceTraversalPolicy,
    visit: Optional[Callable[[ResonanceSearchNode], bool]],
) -> None:
    visitor = visit if visit is not None else (lambda _node: True)
    root_key, bond_index_map = resonance_utils._build_resonance_search_context(omol)
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

        indexed_moves = resonance_utils._enumerate_one_step_resonance_moves(
            current,
            current_key,
            bond_index_map,
        )
        selected_moves = resonance_utils._apply_resonance_traversal_policy(
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
            best_known_discrepancy = best_discrepancy_by_state.get(move.next_state_key)
            if best_known_discrepancy is not None and best_known_discrepancy <= next_discrepancy:
                continue
            best_discrepancy_by_state[move.next_state_key] = next_discrepancy
            push_order += 1
            heappush(
                frontier,
                (
                    next_discrepancy,
                    depth + 1,
                    push_order,
                    resonance_utils._materialize_one_step_resonance(current, move.idxs),
                    move.next_state_key,
                ),
            )


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
    "get_radical_resonances",
    "make_limited_discrepancy_direct_gain_traversal_policy",
    "make_limited_discrepancy_force_field_traversal_policy",
    "make_limited_discrepancy_input_order_traversal_policy",
    "process_resonance",
    "resonance_move_score_cache_clear",
    "resonance_move_score_cache_info",
    "walk_radical_resonances",
]
