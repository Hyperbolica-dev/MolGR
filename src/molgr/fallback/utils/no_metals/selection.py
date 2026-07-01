"""Selection helpers for no-metal fallback reconstruction."""

from __future__ import annotations

from dataclasses import astuple
from typing import cast

from openbabel import openbabel as ob

from molgr.config import CONFIG, MolGRConfig
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
    resolved_config = CONFIG if config is None else config
    topology_config = resolved_config.organic_topology
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


def _formal_charge_absolute_sum(candidate: ReconstructionState) -> int:
    return int(
        candidate.get_cached_omol_value(
            "organic_formal_charge_absolute_sum",
            lambda omol: sum(
                abs(int(cast(ob.OBAtom, atom_iter).GetFormalCharge()))
                for atom_iter in ob.OBMolAtomIter(cast(ob.OBMol, omol.OBMol))
            ),
        )
    )


def _no_metal_candidate_selection_key(
    candidate: ReconstructionState,
    *,
    config: MolGRConfig | None = None,
) -> tuple[float, int, float, float, float, float]:
    metrics = _annotate_no_metal_candidate_topology(candidate, config=config)
    score = float(candidate.metadata.get("score", float("inf")))
    formal_charge_absolute_sum = _formal_charge_absolute_sum(candidate)
    conjugation_charge_penalty = formal_charge_absolute_sum / 2.0
    adjusted_max_conjugated_component_size = (
        metrics.max_conjugated_component_size - conjugation_charge_penalty
    )
    adjusted_conjugated_atom_count = metrics.conjugated_atom_count - conjugation_charge_penalty
    adjusted_conjugated_bond_count = metrics.conjugated_bond_count - conjugation_charge_penalty
    candidate.metadata["organic_formal_charge_absolute_sum"] = formal_charge_absolute_sum
    candidate.metadata["organic_conjugation_charge_penalty"] = conjugation_charge_penalty
    candidate.metadata["organic_adjusted_max_conjugated_component_size"] = (
        adjusted_max_conjugated_component_size
    )
    candidate.metadata["organic_adjusted_conjugated_atom_count"] = adjusted_conjugated_atom_count
    candidate.metadata["organic_adjusted_conjugated_bond_count"] = adjusted_conjugated_bond_count
    selection_key = (
        -metrics.aromatic_stability_score,
        -metrics.aromatic_atom_count,
        -adjusted_max_conjugated_component_size,
        -adjusted_conjugated_atom_count,
        -adjusted_conjugated_bond_count,
        score,
    )
    candidate.metadata["organic_topology_selection_key"] = selection_key
    return selection_key


__all__ = [
    "_annotate_no_metal_candidate_topology",
    "_no_metal_candidate_selection_key",
    "_score_reconstruction_candidate",
]
