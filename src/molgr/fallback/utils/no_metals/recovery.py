"""Explicit fallback recovery tiers for no-metal reconstruction."""

from __future__ import annotations

from typing import Callable, Iterable, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.stages.break_bond import break_deformed_ene, break_one_bond
from molgr.fallback.state import OmolStateMachine, ReconstructionState
from molgr.fallback.utils import resonance as resonance_utils


def _clone_omol(omol: pybel.Molecule) -> pybel.Molecule:
    return pybel.Molecule(ob.OBMol(cast(ob.OBMol, omol.OBMol)))


def _deduplicate(states: Iterable[ReconstructionState]) -> list[ReconstructionState]:
    unique: dict[tuple[object, int], ReconstructionState] = {}
    for state in states:
        key = resonance_utils.build_resonance_state_key(state.omol), state.given_charge
        unique.setdefault(key, state)
    return list(unique.values())


def enumerate_deformed_pi_recovery_seeds(
    states: Iterable[ReconstructionState],
    *,
    break_stage: Callable[..., tuple[pybel.Molecule, bool]] | None = None,
) -> list[ReconstructionState]:
    """Lower deformed bonds whose stage explicitly updates endpoint electrons."""

    break_stage = break_deformed_ene if break_stage is None else break_stage
    recovered: list[ReconstructionState] = []
    for state in states:
        machine = OmolStateMachine.from_reconstruction_state(state).branch(
            None,
            omol=_clone_omol(state.omol),
        )
        hit = machine.run_omol_stage(
            "recover_deformed_pi_bonds",
            break_stage,
            machine.given_charge,
            state.total_radical_electrons,
            5.0,
        )
        if not hit:
            continue
        machine.annotate(None, recovery_tier=1, recovery_strategy="deformed_pi_bonds")
        recovered.append(machine.freeze_like(state))
    return _deduplicate(recovered)


def enumerate_bond_break_recovery_seeds(
    states: Iterable[ReconstructionState],
    *,
    break_stage: Callable[..., tuple[pybel.Molecule, int, bool]] | None = None,
) -> list[ReconstructionState]:
    """Apply bond-breaking rules whose stages explicitly update endpoint electrons."""

    break_stage = break_one_bond if break_stage is None else break_stage
    recovered: list[ReconstructionState] = []
    for state in states:
        machine = OmolStateMachine.from_reconstruction_state(state).branch(
            None,
            omol=_clone_omol(state.omol),
        )
        hit = machine.run_omol_charge_stage(
            "recover_by_breaking_bonds",
            break_stage,
            state.total_radical_electrons,
        )
        if not hit:
            continue
        machine.annotate(None, recovery_tier=2, recovery_strategy="bond_break")
        recovered.append(machine.freeze_like(state))
    return _deduplicate(recovered)


__all__ = [
    "enumerate_bond_break_recovery_seeds",
    "enumerate_deformed_pi_recovery_seeds",
]
