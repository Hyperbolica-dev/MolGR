"""Preparation helpers for no-metal fallback reconstruction."""

from __future__ import annotations

from typing import cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.stages.break_bond import break_deformed_ene, break_one_bond
from molgr.fallback.stages.clean import (
    clean_carbene_neighbor_unsaturated,
    clean_neighbor_radicals,
    clean_resonances,
)
from molgr.fallback.stages.eliminate import (
    eliminate_carbene_neighbor_heteroatom,
    eliminate_carboxyl,
    eliminate_charge_spliting,
    eliminate_CN_in_doubt,
    eliminate_high_positive_charge_atoms,
    eliminate_NNN,
)
from molgr.fallback.stages.fresh import fresh_omol_charge_radical
from molgr.fallback.stages.preprocess import make_connections, pre_clean, validate_omol
from molgr.fallback.state import OmolStateMachine, ReconstructionState


NEIGHBOR_RADICAL_RESOLUTION_STRATEGY_KEY = "neighbor_radical_resolution_strategy"


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


def _clone_omol(omol: pybel.Molecule) -> pybel.Molecule:
    return pybel.Molecule(ob.OBMol(cast(ob.OBMol, omol.OBMol)))


def _neighbor_radical_bond_pairs(omol: pybel.Molecule) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    for bond in ob.OBMolBondIter(cast(ob.OBMol, omol.OBMol)):
        begin_atom = cast(ob.OBAtom, bond.GetBeginAtom())
        end_atom = cast(ob.OBAtom, bond.GetEndAtom())
        if begin_atom.GetSpinMultiplicity() > 0 and end_atom.GetSpinMultiplicity() > 0:
            pairs.append((begin_atom.GetIdx(), end_atom.GetIdx()))
    return pairs


def _clean_neighbor_radicals_charge_split(
    omol: pybel.Molecule,
    begin_charge_sign: int,
) -> tuple[pybel.Molecule, bool]:
    """Convert adjacent radical pairs into one charge-separated orientation."""

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False
    for begin_idx, end_idx in _neighbor_radical_bond_pairs(omol):
        begin_atom = cast(ob.OBAtom, obmol.GetAtom(begin_idx))
        end_atom = cast(ob.OBAtom, obmol.GetAtom(end_idx))
        if begin_atom is None or end_atom is None:
            continue
        spin_to_consume = min(
            begin_atom.GetSpinMultiplicity(),
            end_atom.GetSpinMultiplicity(),
        )
        if spin_to_consume <= 0:
            continue
        begin_atom.SetSpinMultiplicity(begin_atom.GetSpinMultiplicity() - spin_to_consume)
        end_atom.SetSpinMultiplicity(end_atom.GetSpinMultiplicity() - spin_to_consume)
        begin_atom.SetFormalCharge(
            begin_atom.GetFormalCharge() + begin_charge_sign * spin_to_consume
        )
        end_atom.SetFormalCharge(
            end_atom.GetFormalCharge() - begin_charge_sign * spin_to_consume
        )
        hit = True
    return omol, hit


def _run_deterministic_pre_resolution_stages(state: ReconstructionState) -> OmolStateMachine:
    """Run deterministic no-metal stages before neighbor-radical resolution."""

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

    return machine


def _run_direct_neighbor_radical_resolution(machine: OmolStateMachine) -> OmolStateMachine:
    """Consume adjacent radicals by increasing bond order on the current path."""

    machine.metadata[NEIGHBOR_RADICAL_RESOLUTION_STRATEGY_KEY] = "direct"
    machine.run_omol_stage("clean_neighbor_radicals", clean_neighbor_radicals)
    return machine


def _neighbor_radical_charge_split_resolution(
    machine: OmolStateMachine,
    *,
    strategy: str,
    phase: str,
    begin_charge_sign: int,
) -> OmolStateMachine:
    machine.metadata[NEIGHBOR_RADICAL_RESOLUTION_STRATEGY_KEY] = strategy
    machine.run_omol_stage(
        phase,
        _clean_neighbor_radicals_charge_split,
        begin_charge_sign,
    )
    return machine


def _enumerate_neighbor_radical_resolution_machines(
    machine: OmolStateMachine,
) -> list[OmolStateMachine]:
    """Enumerate explicit neighbor-radical resolution strategies for candidates."""

    if not _neighbor_radical_bond_pairs(machine.omol):
        return [_run_direct_neighbor_radical_resolution(machine)]

    direct_machine = machine.branch(None, omol=_clone_omol(machine.omol))
    direct_machine = _run_direct_neighbor_radical_resolution(direct_machine)

    begin_positive_machine = machine.branch(None, omol=_clone_omol(machine.omol))
    begin_positive_machine = _neighbor_radical_charge_split_resolution(
        begin_positive_machine,
        strategy="charge_begin_positive",
        phase="clean_neighbor_radicals_charge_begin_positive",
        begin_charge_sign=1,
    )

    begin_negative_machine = machine.branch(None, omol=_clone_omol(machine.omol))
    begin_negative_machine = _neighbor_radical_charge_split_resolution(
        begin_negative_machine,
        strategy="charge_begin_negative",
        phase="clean_neighbor_radicals_charge_begin_negative",
        begin_charge_sign=-1,
    )

    return [direct_machine, begin_positive_machine, begin_negative_machine]


def _run_deterministic_post_resolution_stages(
    machine: OmolStateMachine,
    state: ReconstructionState,
) -> ReconstructionState:
    """Run deterministic no-metal stages after one neighbor-radical strategy."""

    machine.run_omol_stage(
        "clean_carbene_neighbor_unsaturated_second",
        clean_carbene_neighbor_unsaturated,
    )
    machine.run_omol_charge_stage("eliminate_charge_spliting", eliminate_charge_spliting)
    machine.run_omol_stage(
        "break_deformed_ene",
        break_deformed_ene,
        machine.given_charge,
        state.total_radical_electrons,
        5.0,
    )
    machine.run_omol_charge_stage(
        "break_one_bond",
        break_one_bond,
        state.total_radical_electrons,
    )
    machine.run_omol_stage("fresh_omol_charge_radical_final", fresh_omol_charge_radical)

    return machine.freeze_like(state)


def _enumerate_no_metal_candidate_states(state: ReconstructionState) -> list[ReconstructionState]:
    """Run deterministic stages around explicit neighbor-radical strategies."""

    machine = _run_deterministic_pre_resolution_stages(state)
    candidates: list[ReconstructionState] = []
    for strategy_machine in _enumerate_neighbor_radical_resolution_machines(machine):
        candidates.append(_run_deterministic_post_resolution_stages(strategy_machine, state))
    return candidates


def _run_linear_pipeline(state: ReconstructionState) -> ReconstructionState:
    """Run the deterministic direct no-metal stage sequence."""

    machine = _run_deterministic_pre_resolution_stages(state)
    machine = _run_direct_neighbor_radical_resolution(machine)
    return _run_deterministic_post_resolution_stages(machine, state)


__all__ = [
    "NEIGHBOR_RADICAL_RESOLUTION_STRATEGY_KEY",
    "_enumerate_neighbor_radical_resolution_machines",
    "_enumerate_no_metal_candidate_states",
    "_run_linear_pipeline",
    "_seed_omol_from_xyz",
    "_seed_state",
    "_seed_state_from_omol",
    "clean_resonances",
    "validate_omol",
]
