from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import Literal, Optional


ResonanceTraversalScore = Literal[
    "uff_lite_gain",
    "input_order",
]
ReconstructionFailurePolicy = Literal[
    "raise",
    "return_suspicious",
]


@dataclass(frozen=True)
class ResonanceConfig:
    max_depth: int = 2
    limited_discrepancy_max_discrepancy: int = 1
    traversal_score: ResonanceTraversalScore = "uff_lite_gain"


@dataclass(frozen=True)
class CppBackendConfig:
    max_threads: Optional[int] = None
    enable_target_bucket_parallelism: bool = True
    enable_candidate_scoring_parallelism: bool = False
    enable_uff_atom_typing_cache: bool = True
    candidate_score_parallel_threshold: int = 32


@dataclass(frozen=True)
class OrganicTopologyConfig:
    aromatic_stability_benzene_score: float = 1.0
    aromatic_stability_other_ring_max_score: float = 0.99
    aromatic_stability_ring_size_6_factor: float = 0.92
    aromatic_stability_ring_size_5_factor: float = 0.84
    aromatic_stability_other_ring_size_factor: float = 0.76
    aromatic_stability_hetero_atom_penalty: float = 0.10
    aromatic_stability_min_hetero_factor: float = 0.62
    aromatic_stability_formal_charge_penalty: float = 0.16
    aromatic_stability_min_charge_factor: float = 0.50
    aromatic_stability_radical_penalty: float = 0.20
    aromatic_stability_min_radical_factor: float = 0.50


@dataclass(frozen=True)
class MetalScoringConfig:
    open_shell_multimetal_state_penalty_window: float = 10.0
    open_shell_multimetal_min_state_options: int = 6
    same_element_multimetal_unify_threshold: int = 3
    max_mixed_valence_spread: Optional[int] = 3
    max_assignments_per_target: int = 64
    metal_coordination_extra_tolerance_angstrom: float = 0.35
    pi_dative_distance_difference_tolerance_angstrom: float = 0.10
    metal_access_radius_scale: float = 1.0
    metal_access_clearance_angstrom: float = 0.0


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
class PythonInterfaceConfig:
    reconstruction_failure_policy: ReconstructionFailurePolicy = "raise"


@dataclass(frozen=True)
class MolGRConfig:
    resonance: ResonanceConfig = field(default_factory=ResonanceConfig)
    cpp_backend: CppBackendConfig = field(default_factory=CppBackendConfig)
    organic_topology: OrganicTopologyConfig = field(default_factory=OrganicTopologyConfig)
    metal_scoring: MetalScoringConfig = field(default_factory=MetalScoringConfig)
    metal_radical_inference: MetalRadicalInferenceConfig = field(
        default_factory=MetalRadicalInferenceConfig
    )
    interface: PythonInterfaceConfig = field(default_factory=PythonInterfaceConfig)


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


__all__ = [
    "DEFAULT_MOLGR_CONFIG",
    "CppBackendConfig",
    "MetalRadicalInferenceConfig",
    "MetalScoringConfig",
    "MolGRConfig",
    "OrganicTopologyConfig",
    "PythonInterfaceConfig",
    "ReconstructionFailurePolicy",
    "ResonanceConfig",
    "ResonanceTraversalScore",
    "get_config",
    "is_default_config",
    "make_default_config",
    "reset_config",
    "resolve_config",
    "set_config",
    "sync_cpp_backend_default_config",
]
