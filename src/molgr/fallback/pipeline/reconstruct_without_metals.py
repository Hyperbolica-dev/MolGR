"""No-metal reconstruction through unified resonance seeds and recovery tiers."""

from __future__ import annotations

from dataclasses import astuple, dataclass, field
from typing import Optional

from openbabel import pybel

from molgr.config import CONFIG, MolGRConfig
from molgr.fallback.state import OmolStateMachine, ReconstructionState
from molgr.fallback.utils.no_metals import neighbor_radicals, preparation, recovery, selection
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

    prepared_seed = preparation.prepare_no_metal_seed(seed_state)
    traversal_policy = no_metal_resonance._default_resonance_traversal_policy(config)
    resolved_config = CONFIG if config is None else config
    max_discrepancy = max(
        0,
        int(resolved_config.resonance.limited_discrepancy_max_discrepancy),
    )
    resonance_seeds: list[ReconstructionState] = []
    candidate_pool: list[ReconstructionState] = []
    search_session = no_metal_resonance._ResonanceSearchSession()
    for discrepancy in range(max_discrepancy + 1):
        neighbor_seeds = neighbor_radicals.enumerate_neighbor_radical_seeds(
            prepared_seed,
            exact_discrepancy=discrepancy,
        )
        layer_seeds = no_metal_resonance.build_resonance_seed_pool(neighbor_seeds)
        resonance_seeds.extend(layer_seeds)
        candidate_pool = no_metal_resonance.search_resonance_candidates(
            layer_seeds,
            resonance_traversal_policy=traversal_policy,
            config=config,
            session=search_session,
        )
        if candidate_pool:
            break

    deformed_pi_seeds: list[ReconstructionState] = []
    if not candidate_pool:
        deformed_pi_seeds = recovery.enumerate_deformed_pi_recovery_seeds(resonance_seeds)
        if deformed_pi_seeds:
            candidate_pool = no_metal_resonance.search_resonance_candidates(
                no_metal_resonance.build_resonance_seed_pool(deformed_pi_seeds),
                resonance_traversal_policy=traversal_policy,
                config=config,
            )

    if not candidate_pool:
        bond_break_seeds = recovery.enumerate_bond_break_recovery_seeds(
            (*resonance_seeds, *deformed_pi_seeds)
        )
        if bond_break_seeds:
            candidate_pool = no_metal_resonance.search_resonance_candidates(
                no_metal_resonance.build_resonance_seed_pool(bond_break_seeds),
                resonance_traversal_policy=traversal_policy,
                config=config,
            )

    if not candidate_pool:
        return None

    best_candidate: Optional[ReconstructionState] = None
    best_selection_key: Optional[tuple[float, int, float, float, float, int, float]] = None
    for candidate in candidate_pool:
        selection_key = selection._no_metal_candidate_selection_key(candidate, config=config)
        if best_selection_key is not None and selection_key >= best_selection_key:
            continue
        best_selection_key = selection_key
        best_candidate = candidate

    if best_candidate is None:
        return None

    result_machine = OmolStateMachine.from_reconstruction_state(best_candidate)
    result_machine.annotate("select_best_no_metal_candidate")
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
