"""No-metal reconstruction pipeline for fallback.

The flow is intentionally linear until validation fails:
1. Apply deterministic preprocess / eliminate / clean / bond-breaking stages.
2. If the structure is already valid, score it directly.
3. Otherwise enumerate resonance candidates, normalize them, and choose the winner.
"""

from __future__ import annotations

from typing import List, Optional, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.pipeline import resonance as resonance_module
from molgr.fallback.pipeline.resonance import (
    ResonanceTraversalPolicy,
    build_processed_resonance_key,
    get_radical_resonances,
    make_limited_discrepancy_direct_gain_traversal_policy,
    process_resonance,
)
from molgr.fallback.stages.break_bond import break_deformed_ene, break_one_bond
from molgr.fallback.stages.clean import (
    clean_carbene_neighbor_unsaturated,
    clean_neighbor_radicals,
    clean_resonances,
)
from molgr.fallback.stages.eliminate import (
    eliminate_carbene_neighbor_heteroatom,
    eliminate_carboxyl,
    eliminate_charge_spliting,
    eliminate_CN_in_doubt,
    eliminate_high_positive_charge_atoms,
    eliminate_NNN,
)
from molgr.fallback.stages.fresh import fresh_omol_charge_radical
from molgr.fallback.stages.preprocess import make_connections, pre_clean, validate_omol
from molgr.fallback.state import OmolStateMachine, ReconstructionState
from molgr.fallback.utils.tools import typed_lru_cache


_DEFAULT_RESONANCE_TRAVERSAL_POLICY = make_limited_discrepancy_direct_gain_traversal_policy(
    max_discrepancy=1,
    fallback_to_full_frontier=True,
)
_ENABLE_RESONANCE_INCUMBENT_BOUND = True
_RESONANCE_SEARCH_MAX_DEPTH = 2
_RESONANCE_INCUMBENT_PRUNE_MARGIN = 5.0


def _seed_state(
    xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
) -> ReconstructionState:
    """Create the initial reconstruction state from the input XYZ block."""

    return ReconstructionState(
        omol=pybel.readstring("xyz", xyz_block),
        given_charge=0,
        total_charge=total_charge,
        total_radical_electrons=total_radical_electrons,
        phase_history=("read_xyz",),
        metadata={"source": "xyz_to_omol_no_metal_state"},
    )


def _run_linear_pipeline(state: ReconstructionState) -> ReconstructionState:
    """Run the deterministic no-metal stage sequence before resonance recovery."""

    machine = OmolStateMachine.from_reconstruction_state(state)

    machine.run_omol_stage("make_connections", make_connections)
    machine.run_omol_stage("pre_clean", pre_clean)
    machine.run_omol_stage("fresh_omol_charge_radical_initial", fresh_omol_charge_radical)

    machine.set_given_charge(
        "initialize_charge_budget",
        state.total_charge
        - sum(
            cast(ob.OBAtom, atom.OBAtom).GetFormalCharge() for atom in machine.omol.atoms
        ),
    )

    machine.run_omol_charge_stage("eliminate_NNN_negative", eliminate_NNN, False)
    machine.run_omol_charge_stage(
        "eliminate_high_positive_charge_atoms",
        eliminate_high_positive_charge_atoms,
    )
    machine.run_omol_charge_stage("eliminate_CN_in_doubt", eliminate_CN_in_doubt)
    machine.run_omol_charge_stage("eliminate_NNN_positive", eliminate_NNN, True)
    machine.run_omol_charge_stage("eliminate_carboxyl", eliminate_carboxyl)
    machine.run_omol_stage(
        "clean_carbene_neighbor_unsaturated_first",
        clean_carbene_neighbor_unsaturated,
    )
    machine.run_omol_charge_stage(
        "eliminate_carbene_neighbor_heteroatom",
        eliminate_carbene_neighbor_heteroatom,
    )
    machine.run_omol_stage("clean_neighbor_radicals", clean_neighbor_radicals)
    machine.run_omol_stage(
        "clean_carbene_neighbor_unsaturated_second",
        clean_carbene_neighbor_unsaturated,
    )
    machine.run_omol_charge_stage("eliminate_charge_spliting", eliminate_charge_spliting)
    machine.run_omol_stage(
        "break_deformed_ene",
        break_deformed_ene,
        machine.given_charge,
        state.total_radical_electrons,
    )
    machine.run_omol_charge_stage(
        "break_one_bond",
        break_one_bond,
        state.total_radical_electrons,
    )
    machine.run_omol_stage("fresh_omol_charge_radical_final", fresh_omol_charge_radical)

    return machine.freeze_like(state)


def _recover_resonance_candidates(
    state: ReconstructionState,
    *,
    resonance_traversal_policy: Optional[ResonanceTraversalPolicy] = None,
) -> List[ReconstructionState]:
    """Enumerate, normalize, validate, and score resonance candidates."""

    if _ENABLE_RESONANCE_INCUMBENT_BOUND and isinstance(
        resonance_traversal_policy,
        resonance_module._LimitedDiscrepancyDirectGainTraversalPolicy,
    ):
        return _recover_resonance_candidates_with_incumbent_bound(
            state,
            resonance_traversal_policy=resonance_traversal_policy,
        )

    candidates: List[ReconstructionState] = []
    seen_processed_states = set()
    base_machine = OmolStateMachine.from_reconstruction_state(state)
    if resonance_traversal_policy is None:
        resonance_iterable = get_radical_resonances(state.omol)
    else:
        resonance_iterable = get_radical_resonances(
            state.omol,
            traversal_policy=resonance_traversal_policy,
        )

    for resonance_index, resonance in enumerate(resonance_iterable):
        candidate_machine = base_machine.branch("branch_resonance_candidate", omol=resonance)
        candidate_machine.run_omol_charge_stage("process_resonance", process_resonance)
        processed_state_key = candidate_machine.get_cached_omol_value(
            "resonance_state_key",
            build_processed_resonance_key,
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
        candidate.score("organic_core")
        candidate.post_reinsertion_base_components()
        candidate.full_score()
        candidates.append(candidate)

    return candidates


def _recover_resonance_candidates_with_incumbent_bound(
    state: ReconstructionState,
    *,
    resonance_traversal_policy: resonance_module._LimitedDiscrepancyDirectGainTraversalPolicy,
) -> List[ReconstructionState]:
    """Traverse resonance states with incumbent-based branch pruning enabled."""

    # Score validated processed states once, then prune branches whose optimistic
    # remaining improvement still cannot beat the current incumbent.
    candidates: List[ReconstructionState] = []
    processed_state_scores = {}
    base_machine = OmolStateMachine.from_reconstruction_state(state)
    best_score = float("inf")
    bound_cache = {}
    resonance_index = 0

    def _visit(node: resonance_module.ResonanceSearchNode) -> bool:
        nonlocal best_score, resonance_index

        current_resonance_index = resonance_index
        resonance_index += 1

        candidate_machine = base_machine.branch(
            "branch_resonance_candidate",
            omol=node.omol.clone,
        )
        candidate_machine.run_omol_charge_stage("process_resonance", process_resonance)
        processed_state_key = candidate_machine.get_cached_omol_value(
            "resonance_state_key",
            build_processed_resonance_key,
        )
        cached_score = processed_state_scores.get(processed_state_key)
        if cached_score is None and processed_state_key not in processed_state_scores:
            processed_state_scores[processed_state_key] = None
            if validate_omol(
                candidate_machine.omol,
                state.total_charge,
                state.total_radical_electrons,
            ):
                candidate_machine.annotate(
                    "validate_resonance_candidate",
                    resonance_index=current_resonance_index,
                )
                candidate = candidate_machine.freeze_like(state)
                candidate.score("organic_core")
                candidate.post_reinsertion_base_components()
                cached_score = candidate.full_score()
                processed_state_scores[processed_state_key] = cached_score
                candidates.append(candidate)
                if cached_score < best_score:
                    best_score = cached_score
            else:
                cached_score = None

        remaining_steps = _RESONANCE_SEARCH_MAX_DEPTH - node.depth
        if remaining_steps <= 0 or cached_score is None:
            return node.depth < _RESONANCE_SEARCH_MAX_DEPTH
        if cached_score < best_score + _RESONANCE_INCUMBENT_PRUNE_MARGIN:
            return True

        optimistic_improvement = (
            resonance_module.estimate_remaining_resonance_score_improvement_upper_bound(
                node.omol,
                node.state_key,
                remaining_steps,
                bound_cache,
            )
        )
        return cached_score - optimistic_improvement < best_score

    resonance_module.walk_radical_resonances(
        state.omol,
        max_depth=_RESONANCE_SEARCH_MAX_DEPTH,
        traversal_policy=resonance_traversal_policy,
        visit=_visit,
    )
    return candidates


@typed_lru_cache(maxsize=1024, typed=True)
def xyz_to_omol_no_metal_state(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
) -> Optional[ReconstructionState]:
    """Return the best no-metal reconstruction state for the requested charge/radicals."""

    if total_radical_electrons < 0:
        return None

    state = _seed_state(xyz_block, total_charge, total_radical_electrons)
    state = _run_linear_pipeline(state)

    if validate_omol(state.omol, total_charge, total_radical_electrons):
        result_machine = OmolStateMachine.from_reconstruction_state(state)
        result_machine.annotate("validate_direct_candidate")
        result_machine.run_omol_stage("clean_resonances", clean_resonances)
        result = result_machine.freeze_like(state)
        result.score("organic_core")
        result.post_reinsertion_base_components()
        return result

    resonance_candidates = _recover_resonance_candidates(
        state,
        resonance_traversal_policy=_DEFAULT_RESONANCE_TRAVERSAL_POLICY,
    )
    if not resonance_candidates:
        return None

    best_candidate: Optional[ReconstructionState] = None
    best_score = float("inf")
    for candidate in resonance_candidates:
        score = candidate.full_score()
        if score >= best_score:
            continue
        best_score = score
        best_candidate = candidate

    if best_candidate is None:
        return None

    result_machine = OmolStateMachine.from_reconstruction_state(best_candidate)
    result_machine.annotate("select_best_resonance_candidate")
    return result_machine.freeze_like(best_candidate)


def xyz_to_omol_no_metal(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
) -> Optional[pybel.Molecule]:
    """Materialize the winning no-metal reconstruction."""

    state = xyz_to_omol_no_metal_state(xyz_block, total_charge, total_radical_electrons)
    if state is None:
        return None
    return state.omol


__all__ = ["xyz_to_omol_no_metal", "xyz_to_omol_no_metal_state"]
