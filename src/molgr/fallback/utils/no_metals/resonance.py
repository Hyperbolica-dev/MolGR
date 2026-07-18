"""Resonance recovery helpers for no-metal fallback reconstruction."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import Callable, Optional, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.config import CONFIG, MolGRConfig
from molgr.fallback.pipeline.resonance import ResonanceTraversalPolicy, walk_radical_resonances
from molgr.fallback.stages.clean import (
    clean_carbene_neighbor_unsaturated,
    clean_resonances,
)
from molgr.fallback.stages.eliminate import (
    assign_negative_charges_from_radicals,
)
from molgr.fallback.stages.fresh import fresh_omol_charge_radical
from molgr.fallback.stages.preprocess import validate_omol
from molgr.fallback.state import OmolStateMachine, ReconstructionState
from molgr.fallback.utils import resonance as resonance_utils

from .selection import _annotate_no_metal_candidate_topology, _score_reconstruction_candidate


_OmolStage = Callable[[pybel.Molecule], Tuple[pybel.Molecule, bool]]
_OmolChargeStage = Callable[
    [pybel.Molecule, int],
    Tuple[pybel.Molecule, int, bool],
]


def _default_resonance_traversal_policy(
    config: MolGRConfig | None = None,
) -> ResonanceTraversalPolicy:
    resolved_config = CONFIG if config is None else config
    max_discrepancy = max(0, int(resolved_config.resonance.limited_discrepancy_max_discrepancy))
    traversal_score = resolved_config.resonance.traversal_score
    if traversal_score == "uff_lite_gain":
        return resonance_utils.make_limited_discrepancy_uff_lite_gain_traversal_policy(
            max_discrepancy=max_discrepancy,
        )
    if traversal_score == "input_order":
        return resonance_utils.make_limited_discrepancy_input_order_traversal_policy(
            max_discrepancy=max_discrepancy,
        )
    raise ValueError(
        "Unsupported resonance traversal_score. Expected 'uff_lite_gain' or 'input_order'."
    )


def _resonance_max_depth(config: MolGRConfig | None = None) -> int:
    resolved_config = CONFIG if config is None else config
    return max(0, int(resolved_config.resonance.max_depth))


def _clone_omol(omol: pybel.Molecule) -> pybel.Molecule:
    return pybel.Molecule(ob.OBMol(cast(ob.OBMol, omol.OBMol)))


@dataclass
class _RawResonanceCandidate:
    seed: ReconstructionState
    seed_index: int
    resonance_index: int
    raw_index: int
    omol: pybel.Molecule


@dataclass
class _ResonanceSearchSession:
    seen_raw_states: set[tuple[object, int]] = field(default_factory=set)
    labels_by_state: dict[tuple[object, int], list[tuple[int, int]]] = field(default_factory=dict)
    seen_processed_states: set[tuple[str, int]] = field(default_factory=set)
    next_raw_index: int = 0


def _candidate_key(state: ReconstructionState) -> tuple[object, int]:
    return resonance_utils.build_resonance_state_key(state.omol), state.given_charge


def _deduplicate_states(states: Iterable[ReconstructionState]) -> list[ReconstructionState]:
    unique: dict[tuple[object, int], ReconstructionState] = {}
    for state in states:
        unique.setdefault(_candidate_key(state), state)
    return list(unique.values())


def build_resonance_seed_pool(
    neighbor_seeds: Sequence[ReconstructionState],
) -> list[ReconstructionState]:
    """Add optional electronic-state normalizations to the neighbor-radical seeds."""

    pool = _deduplicate_states(neighbor_seeds)
    stages = (
        ("relocate_carbene_radical_for_resonance", clean_carbene_neighbor_unsaturated, False),
        (
            "assign_negative_charges_from_radicals",
            assign_negative_charges_from_radicals,
            True,
        ),
        ("refresh_electronic_labels_for_resonance", fresh_omol_charge_radical, False),
    )
    for phase, stage, uses_charge_budget in stages:
        additions: list[ReconstructionState] = []
        for state in pool:
            machine = OmolStateMachine.from_reconstruction_state(state).branch(
                None,
                omol=_clone_omol(state.omol),
            )
            if uses_charge_budget:
                hit = machine.run_omol_charge_stage(
                    phase,
                    cast(_OmolChargeStage, stage),
                )
            else:
                hit = machine.run_omol_stage(phase, cast(_OmolStage, stage))
            if hit:
                additions.append(machine.freeze_like(state))
        pool = _deduplicate_states((*pool, *additions))
    return pool


def _register_traversal_label(
    labels_by_state: dict[tuple[object, int], list[tuple[int, int]]],
    state_key: tuple[object, int],
    label: tuple[int, int],
) -> bool:
    labels = labels_by_state.setdefault(state_key, [])
    if any(depth <= label[0] and discrepancy <= label[1] for depth, discrepancy in labels):
        return False
    labels[:] = [
        (depth, discrepancy)
        for depth, discrepancy in labels
        if not (label[0] <= depth and label[1] <= discrepancy)
    ]
    labels.append(label)
    return True


def search_resonance_candidates(
    states: Sequence[ReconstructionState],
    *,
    resonance_traversal_policy: Optional[ResonanceTraversalPolicy] = None,
    config: MolGRConfig | None = None,
    session: _ResonanceSearchSession | None = None,
) -> list[ReconstructionState]:
    """Search new states and normalize them once across an optional shared session."""

    search_session = _ResonanceSearchSession() if session is None else session
    raw_candidates: list[_RawResonanceCandidate] = []
    use_dominance_pruning = resonance_traversal_policy is None or isinstance(
        resonance_traversal_policy,
        resonance_utils._LIMITED_DISCREPANCY_POLICY_TYPES,
    )
    resonance_max_depth = _resonance_max_depth(config)
    for seed_index, state in enumerate(states):
        resonance_index = 0

        def collect(
            node: resonance_utils.ResonanceSearchNode,
            state: ReconstructionState = state,
            seed_index: int = seed_index,
        ) -> bool:
            nonlocal resonance_index
            raw_key = (node.state_key, state.given_charge)
            should_expand = True
            if use_dominance_pruning and node.depth < resonance_max_depth:
                should_expand = _register_traversal_label(
                    search_session.labels_by_state,
                    raw_key,
                    (node.depth, node.discrepancy),
                )
            if raw_key not in search_session.seen_raw_states:
                search_session.seen_raw_states.add(raw_key)
                raw_candidates.append(
                    _RawResonanceCandidate(
                        seed=state,
                        seed_index=seed_index,
                        resonance_index=resonance_index,
                        raw_index=search_session.next_raw_index,
                        omol=node.omol,
                    )
                )
                search_session.next_raw_index += 1
            resonance_index += 1
            return should_expand

        walk_radical_resonances(
            state.omol,
            max_depth=resonance_max_depth,
            traversal_policy=resonance_traversal_policy,
            visit=collect,
        )

    candidates: list[ReconstructionState] = []
    for raw in raw_candidates:
        base_machine = OmolStateMachine.from_reconstruction_state(raw.seed)
        variant_specs = (
            ("resonance_rule_normalization", clean_resonances, False),
            ("full_resonance_normalization", resonance_utils.process_resonance, True),
        )
        for normalization, stage, uses_charge_budget in variant_specs:
            machine = base_machine.branch(
                "branch_resonance_candidate",
                omol=_clone_omol(raw.omol),
            )
            if uses_charge_budget:
                machine.run_omol_charge_stage(
                    normalization,
                    cast(_OmolChargeStage, stage),
                )
            else:
                machine.run_omol_stage(
                    normalization,
                    cast(_OmolStage, stage),
                )
            processed_key = (
                resonance_utils.build_processed_resonance_key(machine.omol),
                machine.given_charge,
            )
            if processed_key in search_session.seen_processed_states:
                continue
            search_session.seen_processed_states.add(processed_key)
            if not validate_omol(
                machine.omol,
                raw.seed.total_charge,
                raw.seed.total_radical_electrons,
            ):
                continue
            machine.annotate(
                "validate_no_metal_candidate",
                resonance_seed_index=raw.seed_index,
                resonance_index=raw.resonance_index,
                resonance_raw_index=raw.raw_index,
                resonance_normalization=normalization,
            )
            candidate = machine.freeze_like(raw.seed)
            try:
                _score_reconstruction_candidate(candidate, config=config)
            except ValueError:
                continue
            _annotate_no_metal_candidate_topology(candidate, config=config)
            candidates.append(candidate)

    return candidates


__all__ = [
    "ResonanceTraversalPolicy",
    "_default_resonance_traversal_policy",
    "build_resonance_seed_pool",
    "search_resonance_candidates",
]
