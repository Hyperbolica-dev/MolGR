from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import Optional, Tuple


ForceFieldChoiceConfig = Tuple[str, ...]


@dataclass(frozen=True)
class CacheConfig:
    pass


@dataclass(frozen=True)
class ForceFieldConfig:
    auto_force_fields_metal_free: ForceFieldChoiceConfig = ("uff",)
    auto_force_fields_with_metals: ForceFieldChoiceConfig = ("uff",)
    organic_force_field: str = "auto"
    selection_force_field: str = "auto"
    combined_force_field: str = "uff"


@dataclass(frozen=True)
class ResonanceConfig:
    max_depth: int = 2
    limited_discrepancy_max_discrepancy: int = 1
    traversal_score: str = "force_field"


@dataclass(frozen=True)
class CppBackendConfig:
    max_threads: Optional[int] = None
    enable_target_bucket_parallelism: bool = True
    enable_candidate_scoring_parallelism: bool = False
    enable_resonance_candidate_parallelism: bool = True
    enable_uff_atom_typing_cache: bool = True
    resonance_candidate_parallel_threshold: int = 8
    candidate_score_parallel_threshold: int = 32


@dataclass(frozen=True)
class MetalScoringConfig:
    organic_score_bucket_relative_ratio: float = 0.20
    organic_force_field_hard_max_ratio: float = 2.5
    open_shell_multimetal_state_penalty_window: float = 10.0
    open_shell_multimetal_min_state_options: int = 6
    same_element_multimetal_unify_threshold: int = 3
    max_mixed_valence_spread: Optional[int] = 3
    max_assignments_per_target: int = 64
    selection_weight_values: Tuple[float, ...] = (
        8.0,
        6.0,
        6.0,
        1.5,
        1.5,
        2.0,
        2.0,
        0.5,
        1.2,
        1.2,
        2.0,
    )
    selection_scale_values: Tuple[float, ...] = (
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        2.0,
        1.0,
        1.0,
        10.0,
    )
    metal_local_potential_cutoff_angstrom: float = 6.0
    metal_donor_cutoff_angstrom: float = 3.4
    min_distance_angstrom: float = 1.2
    metal_access_radius_scale: float = 1.0
    metal_access_clearance_angstrom: float = 0.0
    local_potential_target_per_valence: float = 0.20
    local_potential_oversupport_weight: float = 0.25
    local_donor_target_per_valence: float = 0.80
    local_donor_oversupport_weight: float = 0.35
    local_neutral_donor_weight: float = 0.35
    visible_coordination_reward_weight: float = 1.5
    negative_metal_visible_coordination_penalty_weight: float = 2.0
    obstructed_opposite_charge_penalty_weight: float = 12.0
    same_element_valence_spread_weight: float = 2.0


@dataclass(frozen=True)
class MetalRadicalInferenceConfig:
    coordination_cutoff_angstrom: float = 3.2
    max_considered_donors: int = 6
    square_planar_planarity_tolerance_angstrom: float = 0.45
    trigonal_planar_planarity_tolerance_angstrom: float = 0.35
    linear_angle_min_degrees: float = 150.0
    strong_field_threshold: float = 1.10
    weak_field_threshold: float = 0.75


@dataclass(frozen=True)
class MolGRConfig:
    cache: CacheConfig = field(default_factory=CacheConfig)
    force_field: ForceFieldConfig = field(default_factory=ForceFieldConfig)
    resonance: ResonanceConfig = field(default_factory=ResonanceConfig)
    cpp_backend: CppBackendConfig = field(default_factory=CppBackendConfig)
    metal_scoring: MetalScoringConfig = field(default_factory=MetalScoringConfig)
    metal_radical_inference: MetalRadicalInferenceConfig = field(
        default_factory=MetalRadicalInferenceConfig
    )


def make_default_config() -> MolGRConfig:
    return MolGRConfig()


DEFAULT_MOLGR_CONFIG = make_default_config()
_ACTIVE_MOLGR_CONFIG = make_default_config()


def _sync_cpp_backend_default_config(config: MolGRConfig) -> None:
    core = sys.modules.get("molgr._core")
    if core is None:
        return
    set_default_config = getattr(core, "set_default_config", None)
    if set_default_config is None:
        return
    set_default_config(config)


def get_config() -> MolGRConfig:
    return _ACTIVE_MOLGR_CONFIG


def set_config(config: MolGRConfig) -> None:
    global _ACTIVE_MOLGR_CONFIG
    _ACTIVE_MOLGR_CONFIG = config
    _sync_cpp_backend_default_config(config)


def reset_config() -> None:
    global _ACTIVE_MOLGR_CONFIG
    _ACTIVE_MOLGR_CONFIG = make_default_config()
    _sync_cpp_backend_default_config(_ACTIVE_MOLGR_CONFIG)


def resolve_config(config: Optional[MolGRConfig] = None) -> MolGRConfig:
    return get_config() if config is None else config


def sync_cpp_backend_default_config(config: Optional[MolGRConfig] = None) -> None:
    _sync_cpp_backend_default_config(resolve_config(config))


def is_default_config(config: Optional[MolGRConfig] = None) -> bool:
    return resolve_config(config) == make_default_config()


def force_field_config_cache_key(config: Optional[MolGRConfig] = None) -> str:
    return repr(resolve_config(config).force_field)


__all__ = [
    "DEFAULT_MOLGR_CONFIG",
    "CacheConfig",
    "CppBackendConfig",
    "ForceFieldChoiceConfig",
    "ForceFieldConfig",
    "MetalRadicalInferenceConfig",
    "MetalScoringConfig",
    "MolGRConfig",
    "ResonanceConfig",
    "force_field_config_cache_key",
    "get_config",
    "is_default_config",
    "make_default_config",
    "reset_config",
    "resolve_config",
    "set_config",
    "sync_cpp_backend_default_config",
]
