"""Preparation helpers for no-metal fallback reconstruction."""

from __future__ import annotations

from typing import cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.state import ReconstructionState
from molgr.fallback.utils.electrons import (
    set_lone_pair_count,
    set_unpaired_electron_count,
    set_unresolved_two_electron_center,
)


def _normalize_seed_electronic_labels(omol: pybel.Molecule) -> pybel.Molecule:
    """Clear all inferred seed charges and explicit MolGR electron classifications.

    XYZ/Open Babel perception is not trusted as an electronic-state source. The
    no-metal pipeline therefore resets formal charge, real unpaired electrons,
    active lone pairs, and unresolved-center markers before rebuilding them from
    topology and the supplied global budgets.
    """

    for atom in omol.atoms:
        atom.OBAtom.SetFormalCharge(0)
        set_unpaired_electron_count(atom.OBAtom, 0)
        set_lone_pair_count(atom.OBAtom, 0)
        set_unresolved_two_electron_center(atom.OBAtom, False)
    return omol


def _seed_state(
    xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
) -> ReconstructionState:
    """Create the initial reconstruction state from the input XYZ block."""

    return _seed_state_from_omol(
        _seed_omol_from_xyz(xyz_block),
        total_charge,
        total_radical_electrons,
    )


def _seed_omol_from_xyz(xyz_block: str) -> pybel.Molecule:
    """Parse and normalize the shared initial no-metal molecule from XYZ."""

    return _normalize_seed_electronic_labels(pybel.readstring("xyz", xyz_block))


def _seed_state_from_omol(
    seed_omol: pybel.Molecule,
    total_charge: int,
    total_radical_electrons: int,
) -> ReconstructionState:
    """Create an independent initial state by cloning a normalized seed molecule."""

    omol = _normalize_seed_electronic_labels(
        pybel.Molecule(ob.OBMol(cast(ob.OBMol, seed_omol.OBMol)))
    )
    return ReconstructionState(
        omol=omol,
        given_charge=0,
        total_charge=total_charge,
        total_radical_electrons=total_radical_electrons,
        phase_history=("read_xyz", "normalize_seed_electronic_labels"),
        metadata={"source": "xyz_to_omol_no_metal_state"},
    )


def prepare_no_metal_seed(state: ReconstructionState) -> ReconstructionState:
    """Compatibility entry point for the pipeline-owned linear reconstruction."""

    from molgr.fallback.pipeline.reconstruct_without_metals import _run_linear_preparation

    return _run_linear_preparation(state)


__all__ = [
    "_seed_omol_from_xyz",
    "_seed_state",
    "_seed_state_from_omol",
    "prepare_no_metal_seed",
]
