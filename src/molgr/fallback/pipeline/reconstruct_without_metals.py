"""No-metal reconstruction pipeline for fallback.

The flow is intentionally linear until validation fails:
1. Apply deterministic preprocess / eliminate / clean / bond-breaking stages.
2. If the structure is already valid, score it directly.
3. Otherwise enumerate resonance candidates, normalize them, and choose the winner.
"""

from __future__ import annotations

from dataclasses import astuple, dataclass, field
from typing import Optional

from openbabel import pybel

from molgr.config import CONFIG, MolGRConfig
from molgr.fallback.state import OmolStateMachine, ReconstructionState
from molgr.fallback.utils.no_metals import (
    preparation,
    selection,
)
from molgr.fallback.utils.no_metals import (
    resonance as no_metal_resonance,
)
from molgr.fallback.utils.tools import typed_lru_cache


_DEFAULT_NO_METAL_STATE_CACHE_MAXSIZE = 1024


@dataclass(frozen=True)
class _ConfigCacheToken:
    key: tuple[object, ...]
    config: MolGRConfig = field(compare=False, hash=False)


def _config_cache_token(config: MolGRConfig) -> _ConfigCacheToken:
    return _ConfigCacheToken(astuple(config), config)


def _run_no_metal_pipeline_from_state(
    seed_state: ReconstructionState,
    *,
    config: MolGRConfig | None = None,
) -> Optional[ReconstructionState]:
    if seed_state.total_radical_electrons < 0:
        return None

    state = preparation._run_linear_pipeline(seed_state)

    if preparation.validate_omol(
        state.omol,
        seed_state.total_charge,
        seed_state.total_radical_electrons,
    ):
        result_machine = OmolStateMachine.from_reconstruction_state(state)
        result_machine.annotate("validate_direct_candidate")
        result_machine.run_omol_stage("clean_resonances", preparation.clean_resonances)
        result = result_machine.freeze_like(state)
        try:
            selection._score_reconstruction_candidate(result, config=config)
        except ValueError:
            return None
        selection._annotate_no_metal_candidate_topology(result, config=config)
        return result

    resonance_candidates = no_metal_resonance._recover_resonance_candidates(
        state,
        resonance_traversal_policy=no_metal_resonance._default_resonance_traversal_policy(config),
        config=config,
    )
    if not resonance_candidates:
        return None

    best_candidate: Optional[ReconstructionState] = None
    best_selection_key: Optional[tuple[float, int, float, float, float, float]] = None
    for candidate in resonance_candidates:
        selection._score_reconstruction_candidate(candidate, config=config)
        selection_key = selection._no_metal_candidate_selection_key(candidate, config=config)
        if best_selection_key is not None and selection_key >= best_selection_key:
            continue
        best_selection_key = selection_key
        best_candidate = candidate

    if best_candidate is None:
        return None

    result_machine = OmolStateMachine.from_reconstruction_state(best_candidate)
    result_machine.annotate("select_best_resonance_candidate")
    return result_machine.freeze_like(best_candidate)


@typed_lru_cache(maxsize=_DEFAULT_NO_METAL_STATE_CACHE_MAXSIZE, typed=True)
def _run_no_metal_pipeline_cached(
    seed_state: ReconstructionState,
    config_token: _ConfigCacheToken,
) -> Optional[ReconstructionState]:
    return _run_no_metal_pipeline_from_state(seed_state, config=config_token.config)


def xyz_to_omol_no_metal_state(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
    *,
    config: MolGRConfig | None = None,
) -> Optional[ReconstructionState]:
    """Return the best no-metal reconstruction state for the requested charge/radicals."""

    seed_state = preparation._seed_state(
        xyz_block,
        total_charge,
        total_radical_electrons,
    )
    resolved_config = CONFIG if config is None else config
    return _run_no_metal_pipeline_cached(seed_state, _config_cache_token(resolved_config))


def _seed_omol_to_omol_no_metal_state(
    seed_omol: pybel.Molecule,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
    *,
    config: MolGRConfig | None = None,
) -> Optional[ReconstructionState]:
    """Return the best no-metal state from a shared parsed seed molecule."""

    seed_state = preparation._seed_state_from_omol(
        seed_omol,
        total_charge,
        total_radical_electrons,
    )
    resolved_config = CONFIG if config is None else config
    return _run_no_metal_pipeline_cached(seed_state, _config_cache_token(resolved_config))


def xyz_to_omol_no_metal(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
    *,
    config: MolGRConfig | None = None,
) -> Optional[pybel.Molecule]:
    """Materialize the winning no-metal reconstruction."""

    state = xyz_to_omol_no_metal_state(
        xyz_block,
        total_charge,
        total_radical_electrons,
        config=config,
    )
    if state is None:
        return None
    return state.omol


__all__ = [
    "xyz_to_omol_no_metal",
    "xyz_to_omol_no_metal_state",
]
