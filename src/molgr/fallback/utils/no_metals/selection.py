"""Selection helpers for no-metal fallback reconstruction."""

from __future__ import annotations

from dataclasses import astuple

from molgr.config import MolGRConfig, resolve_config
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
    *,
    config: MolGRConfig | None = None,
) -> OrganicTopologyMetrics:
    topology_config = resolve_config(config).organic_topology
    metrics = candidate.get_cached_omol_value(
        f"organic_topology_metrics:{astuple(topology_config)!r}",
        lambda omol: compute_organic_topology_metrics(
            omol,
            topology_config,
        ),
    )
    candidate.metadata["organic_aromatic_atom_count"] = metrics.aromatic_atom_count
    candidate.metadata["organic_aromatic_ring_count"] = metrics.aromatic_ring_count
    candidate.metadata["organic_aromatic_stability_score"] = metrics.aromatic_stability_score
    candidate.metadata["organic_conjugated_atom_count"] = metrics.conjugated_atom_count
    candidate.metadata["organic_conjugated_bond_count"] = metrics.conjugated_bond_count
    candidate.metadata["organic_max_conjugated_component_size"] = (
        metrics.max_conjugated_component_size
    )
    return metrics


def _no_metal_candidate_selection_key(
    candidate: ReconstructionState,
    *,
    config: MolGRConfig | None = None,
) -> tuple[float, int, int, int, int, float]:
    metrics = _annotate_no_metal_candidate_topology(candidate, config=config)
    score = float(candidate.metadata.get("score", float("inf")))
    selection_key = (
        -metrics.aromatic_stability_score,
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
