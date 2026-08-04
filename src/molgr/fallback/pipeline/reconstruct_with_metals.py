"""Metal-aware reconstruction orchestration for fallback.

The production path is:
1. Strip metals and enumerate a small set of metal valence/radical assignments.
2. Group assignments by the no-metal target they induce via meet-in-the-middle DP.
3. Reconstruct each no-metal target once, reuse that state across the bucket, and
   score the organic core with fixed UFF scoring.
4. Across no-metal charge states, prefer candidates that preserve aromaticity,
   conjugation, and charge localization on chemically plausible sites.
5. Within that organic electronic-state preference, prefer metal assignments whose
   oxidation states are better supported by the fixed local electrostatics and donor field.
"""

from __future__ import annotations

from typing import List, Optional, cast

from openbabel import pybel

from molgr.config import MolGRConfig
from molgr.fallback.state import MetalCandidateState, MetalCandidateStateMachine
from molgr.fallback.utils.metals import preparation, scoring, search
from molgr.fallback.utils.no_metals import preparation as no_metal_preparation

from . import reconstruct_without_metals


def _candidate_matches_global_electronic_state(
    candidate: MetalCandidateState,
    total_charge: int,
    total_radical_electrons: int,
) -> bool:
    return (
        candidate.no_metal_charge_target
        + sum(int(state.valence) for state in candidate.metal_states)
        == total_charge
        and candidate.no_metal_radical_target
        + sum(int(state.radical_num) for state in candidate.metal_states)
        == total_radical_electrons
    )


def xyz2omol_state(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
    *,
    config: MolGRConfig | None = None,
) -> Optional[MetalCandidateState]:
    """Return the best scored metal candidate state for the input XYZ block."""

    base_state = preparation.prepare_metal_state(
        xyz_block,
        total_charge,
        total_radical_electrons,
        config=config,
    )
    state_search_groups = search._build_metal_state_search_groups(
        base_state.available_valence_radical_states,
        config=config,
    )
    layered_state_search_groups = search._build_layered_metal_state_search_groups(
        state_search_groups,
        total_radical_electrons,
        config=config,
    )
    no_metal_seed_omol: Optional[pybel.Molecule] = None
    scored_candidates: List[MetalCandidateState] = []
    winning_layer_index = 0
    for layer_index, available_valence_radical_states in enumerate(layered_state_search_groups):
        grouped_candidates = search._group_candidates_by_target_dp(
            base_state.phase_history,
            available_valence_radical_states,
            total_charge,
            total_radical_electrons,
            config=config,
        )
        if not grouped_candidates:
            continue

        if no_metal_seed_omol is None:
            try:
                no_metal_seed_omol = no_metal_preparation._seed_omol_from_xyz(
                    base_state.no_metal_xyz_block
                )
            except (OSError, ValueError):
                return None

        current_layer_scored_candidates: list[MetalCandidateState] = []
        for candidates in grouped_candidates.values():
            if not candidates:
                continue
            prototype = candidates[0]
            try:
                no_metal_state = reconstruct_without_metals._seed_omol_to_omol_no_metal_state(
                    no_metal_seed_omol,
                    prototype.no_metal_charge_target,
                    prototype.no_metal_radical_target,
                    config=config,
                )
            except (OSError, ValueError):
                continue
            if no_metal_state is None:
                continue

            for candidate in candidates:
                try:
                    scored_candidate = scoring._prepare_candidate_with_no_metal_state(
                        candidate,
                        no_metal_state,
                        config=config,
                    )
                except ValueError:
                    continue
                if cast(Optional[float], scored_candidate.score) is None:
                    continue
                current_layer_scored_candidates.append(scored_candidate)

        if not current_layer_scored_candidates:
            continue
        scored_candidates = current_layer_scored_candidates
        winning_layer_index = layer_index
        break

    if not scored_candidates:
        return None

    for scored_candidate in scored_candidates:
        scored_candidate.metadata["search_layer_index"] = winning_layer_index

    best_candidate = scoring.select_best_candidate(scored_candidates, config=config)
    if best_candidate is None:
        return None
    if not _candidate_matches_global_electronic_state(
        best_candidate,
        total_charge,
        total_radical_electrons,
    ):
        return None
    if best_candidate.combined_omol is None:
        best_candidate.materialize_combined_omol(preparation.combine_metal_with_omol)
        winner_machine = MetalCandidateStateMachine.from_candidate_state(best_candidate)
        winner_machine.annotate("combine_metal_with_omol")
        best_candidate = winner_machine.freeze()

    winner_machine = MetalCandidateStateMachine.from_candidate_state(best_candidate)
    winner_machine.annotate("select_best_candidate")
    return winner_machine.freeze()


def xyz2omol(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
    *,
    config: MolGRConfig | None = None,
) -> Optional[pybel.Molecule]:
    """Materialize the winning metal-aware reconstruction."""

    candidate = xyz2omol_state(
        xyz_block,
        total_charge,
        total_radical_electrons,
        config=config,
    )
    if candidate is None:
        return None
    return candidate.combined_omol


__all__ = [
    "xyz2omol",
    "xyz2omol_state",
]
