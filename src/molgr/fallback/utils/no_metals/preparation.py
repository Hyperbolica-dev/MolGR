"""Preparation helpers for no-metal fallback reconstruction."""

from __future__ import annotations

from typing import cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.stages.clean import clean_carbene_neighbor_unsaturated
from molgr.fallback.stages.eliminate import (
    eliminate_carbene_neighbor_heteroatom,
    eliminate_carboxyl,
    eliminate_CN_in_doubt,
    eliminate_high_positive_charge_atoms,
    eliminate_NNN,
)
from molgr.fallback.stages.fresh import fresh_omol_charge_radical
from molgr.fallback.stages.preprocess import make_connections, pre_clean
from molgr.fallback.state import OmolStateMachine, ReconstructionState


def _normalize_seed_electronic_labels(omol: pybel.Molecule) -> pybel.Molecule:
    """Clear Open Babel's seed-time charge/radical guesses before reconstruction."""

    for atom in omol.atoms:
        atom.OBAtom.SetFormalCharge(0)
        atom.OBAtom.SetSpinMultiplicity(0)
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
    """Run deterministic connectivity and charge cleanup before seed enumeration."""

    machine = OmolStateMachine.from_reconstruction_state(state)

    machine.run_omol_stage("make_connections", make_connections)
    machine.run_omol_stage("pre_clean", pre_clean)
    machine.run_omol_stage("fresh_omol_charge_radical_initial", fresh_omol_charge_radical)

    machine.set_given_charge(
        "initialize_charge_budget",
        state.total_charge
        - sum(cast(ob.OBAtom, atom.OBAtom).GetFormalCharge() for atom in machine.omol.atoms),
    )

    machine.run_omol_charge_stage("eliminate_NNN_negative", eliminate_NNN, False)
    machine.run_omol_charge_stage(
        "eliminate_high_positive_charge_atoms",
        eliminate_high_positive_charge_atoms,
    )
    machine.run_omol_charge_stage("eliminate_CN_in_doubt", eliminate_CN_in_doubt)
    machine.run_omol_charge_stage("eliminate_NNN_positive", eliminate_NNN, True)
    machine.run_omol_charge_stage("eliminate_carboxyl", eliminate_carboxyl)
    machine.run_omol_stage(
        "clean_carbene_neighbor_unsaturated_first",
        clean_carbene_neighbor_unsaturated,
    )
    machine.run_omol_charge_stage(
        "eliminate_carbene_neighbor_heteroatom",
        eliminate_carbene_neighbor_heteroatom,
    )

    machine.annotate("prepare_no_metal_seed")
    return machine.freeze_like(state)


__all__ = [
    "_seed_omol_from_xyz",
    "_seed_state",
    "_seed_state_from_omol",
    "prepare_no_metal_seed",
]
