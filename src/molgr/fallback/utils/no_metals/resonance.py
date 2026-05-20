"""Resonance recovery helpers for no-metal fallback reconstruction."""

from __future__ import annotations

from typing import List, Optional

from molgr.config import MolGRConfig, resolve_config
from molgr.fallback.pipeline.resonance import ResonanceTraversalPolicy, get_radical_resonances
from molgr.fallback.stages.preprocess import validate_omol
from molgr.fallback.state import OmolStateMachine, ReconstructionState
from molgr.fallback.utils import resonance as resonance_utils

from .selection import _annotate_no_metal_candidate_topology, _score_reconstruction_candidate


def _default_resonance_traversal_policy(
    config: MolGRConfig | None = None,
) -> ResonanceTraversalPolicy:
    resolved_config = resolve_config(config)
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
    return max(0, int(resolve_config(config).resonance.max_depth))


def _recover_resonance_candidates(
    state: ReconstructionState,
    *,
    resonance_traversal_policy: Optional[ResonanceTraversalPolicy] = None,
    config: MolGRConfig | None = None,
) -> List[ReconstructionState]:
    """Enumerate, normalize, validate, and score resonance candidates."""

    candidates: List[ReconstructionState] = []
    seen_processed_states = set()
    base_machine = OmolStateMachine.from_reconstruction_state(state)
    resonance_max_depth = _resonance_max_depth(config)
    if resonance_traversal_policy is None:
        if config is None and resonance_max_depth == 2:
            resonance_iterable = get_radical_resonances(state.omol)
        else:
            resonance_iterable = get_radical_resonances(
                state.omol,
                max_depth=resonance_max_depth,
            )
    else:
        if config is None and resonance_max_depth == 2:
            resonance_iterable = get_radical_resonances(
                state.omol,
                traversal_policy=resonance_traversal_policy,
            )
        else:
            resonance_iterable = get_radical_resonances(
                state.omol,
                max_depth=resonance_max_depth,
                traversal_policy=resonance_traversal_policy,
            )

    for resonance_index, resonance in enumerate(resonance_iterable):
        candidate_machine = base_machine.branch("branch_resonance_candidate", omol=resonance)
        candidate_machine.run_omol_charge_stage(
            "process_resonance", resonance_utils.process_resonance
        )
        processed_state_key = candidate_machine.get_cached_omol_value(
            "resonance_state_key",
            resonance_utils.build_processed_resonance_key,
        )
        if processed_state_key in seen_processed_states:
            continue
        seen_processed_states.add(processed_state_key)
        if not validate_omol(
            candidate_machine.omol,
            state.total_charge,
            state.total_radical_electrons,
        ):
            continue
        candidate_machine.annotate("validate_resonance_candidate", resonance_index=resonance_index)
        candidate = candidate_machine.freeze_like(state)
        if config is None:
            candidate.score("organic_core")
        else:
            candidate.score("organic_core", config=config)
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
    "_recover_resonance_candidates",
]
