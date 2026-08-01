"""Selection helpers for no-metal fallback reconstruction."""

from __future__ import annotations

from dataclasses import astuple
from typing import Optional, Sequence, cast

from openbabel import openbabel as ob

from molgr.config import CONFIG, MolGRConfig
from molgr.fallback.state import ReconstructionState
from molgr.fallback.utils.electrons import (
    get_lone_pair_count,
    get_unpaired_electron_count,
    has_unresolved_two_electron_center,
)
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
    candidate.metadata["organic_hyperconjugative_donor_count"] = (
        metrics.hyperconjugative_donor_count
    )
    candidate.metadata["organic_hyperconjugation_score"] = metrics.hyperconjugation_score
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


def _excess_radical_labels(candidate: ReconstructionState) -> int:
    radical_sum = int(
        candidate.get_cached_omol_value(
            "organic_radical_label_sum",
            lambda omol: sum(
                get_unpaired_electron_count(cast(ob.OBAtom, atom_iter))
                for atom_iter in ob.OBMolAtomIter(cast(ob.OBMol, omol.OBMol))
            ),
        )
    )
    excess = max(0, radical_sum - candidate.total_radical_electrons)
    candidate.metadata["organic_radical_label_sum"] = radical_sum
    candidate.metadata["organic_excess_radical_labels"] = excess
    return excess


def _no_metal_candidate_selection_key(
    candidate: ReconstructionState,
    *,
    config: MolGRConfig | None = None,
) -> tuple[int, int, int, float, float, float, float, int, int, float]:
    metrics = _annotate_no_metal_candidate_topology(candidate, config=config)
    score = float(candidate.metadata.get("score", float("inf")))
    formal_charge_absolute_sum = _formal_charge_absolute_sum(candidate)
    conjugation_charge_penalty = formal_charge_absolute_sum / 2.0
    adjusted_max_conjugated_component_size = (
        metrics.max_conjugated_component_size - conjugation_charge_penalty
    )
    adjusted_conjugated_atom_count = metrics.conjugated_atom_count - conjugation_charge_penalty
    adjusted_conjugated_bond_count = metrics.conjugated_bond_count - conjugation_charge_penalty
    excess_radical_labels = _excess_radical_labels(candidate)
    candidate.metadata["organic_formal_charge_absolute_sum"] = formal_charge_absolute_sum
    candidate.metadata["organic_conjugation_charge_penalty"] = conjugation_charge_penalty
    candidate.metadata["organic_adjusted_max_conjugated_component_size"] = (
        adjusted_max_conjugated_component_size
    )
    candidate.metadata["organic_adjusted_conjugated_atom_count"] = adjusted_conjugated_atom_count
    candidate.metadata["organic_adjusted_conjugated_bond_count"] = adjusted_conjugated_bond_count
    selection_key = (
        formal_charge_absolute_sum,
        -metrics.aromatic_atom_count,
        -metrics.aromatic_ring_count,
        -metrics.aromatic_stability_score,
        -adjusted_max_conjugated_component_size,
        -adjusted_conjugated_atom_count,
        -adjusted_conjugated_bond_count,
        excess_radical_labels,
        -metrics.hyperconjugation_score,
        score,
    )
    candidate.metadata["organic_topology_selection_key"] = selection_key
    return selection_key


def _no_metal_candidate_graph_tie_break_key(
    candidate: ReconstructionState,
) -> tuple[int, ...]:
    """Return a backend-independent key for exact no-metal score ties.

    The key deliberately uses atom indices, explicit electron labels, and the
    sorted bond table only. It does not use Open Babel aromaticity flags,
    iterator order, coordinates, or resonance indices.
    """

    atoms = tuple(
        value
        for atom in candidate.omol
        for value in (
            int(atom.idx),
            int(atom.OBAtom.GetAtomicNum()),
            int(atom.OBAtom.GetFormalCharge()),
            int(get_unpaired_electron_count(atom.OBAtom)),
            int(get_lone_pair_count(atom.OBAtom)),
            int(has_unresolved_two_electron_center(atom.OBAtom)),
        )
    )
    bonds = tuple(
        value
        for begin_idx, end_idx, order in sorted(
            (
                min(int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())),
                max(int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())),
                int(bond.GetBondOrder()),
            )
            for bond in ob.OBMolBondIter(candidate.omol.OBMol)
        )
        for value in (begin_idx, end_idx, order)
    )
    return (len(candidate.omol.atoms), *atoms, len(bonds) // 3, *bonds)


def select_best_no_metal_candidate(
    candidates: Sequence[ReconstructionState],
    *,
    config: MolGRConfig | None = None,
) -> Optional[ReconstructionState]:
    """Return the lowest-ranked candidate using chemistry first and graph ties last."""

    if not candidates:
        return None
    return min(
        candidates,
        key=lambda candidate: (
            _no_metal_candidate_selection_key(candidate, config=config),
            _no_metal_candidate_graph_tie_break_key(candidate),
        ),
    )


__all__ = [
    "_annotate_no_metal_candidate_topology",
    "_no_metal_candidate_graph_tie_break_key",
    "_no_metal_candidate_selection_key",
    "_score_reconstruction_candidate",
    "select_best_no_metal_candidate",
]
