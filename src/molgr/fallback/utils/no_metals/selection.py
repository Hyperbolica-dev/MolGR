"""Selection helpers for no-metal fallback reconstruction."""

from __future__ import annotations

from molgr.config import MolGRConfig
from molgr.fallback.state import ReconstructionState
from molgr.fallback.utils.organic_topology import (
    OrganicTopologyMetrics,
    compute_organic_topology_metrics,
)


def _score_reconstruction_candidate(
    candidate: ReconstructionState,
    *,
    config: MolGRConfig | None = None,
) -> float:
    if config is None:
        return candidate.full_score()
    return candidate.full_score(config=config)


def _annotate_no_metal_candidate_topology(
    candidate: ReconstructionState,
) -> OrganicTopologyMetrics:
    metrics = candidate.get_cached_omol_value(
        "organic_topology_metrics",
        compute_organic_topology_metrics,
    )
    candidate.metadata["organic_aromatic_atom_count"] = metrics.aromatic_atom_count
    candidate.metadata["organic_aromatic_ring_count"] = metrics.aromatic_ring_count
    candidate.metadata["organic_conjugated_atom_count"] = metrics.conjugated_atom_count
    candidate.metadata["organic_conjugated_bond_count"] = metrics.conjugated_bond_count
    candidate.metadata["organic_max_conjugated_component_size"] = (
        metrics.max_conjugated_component_size
    )
    return metrics


def _no_metal_candidate_selection_key(
    candidate: ReconstructionState,
) -> tuple[int, int, int, int, float]:
    metrics = _annotate_no_metal_candidate_topology(candidate)
    score = float(candidate.metadata.get("score", float("inf")))
    selection_key = (
        -metrics.aromatic_atom_count,
        -metrics.max_conjugated_component_size,
        -metrics.conjugated_atom_count,
        -metrics.conjugated_bond_count,
        score,
    )
    candidate.metadata["organic_topology_selection_key"] = selection_key
    return selection_key


__all__ = [
    "_annotate_no_metal_candidate_topology",
    "_no_metal_candidate_selection_key",
    "_score_reconstruction_candidate",
]
