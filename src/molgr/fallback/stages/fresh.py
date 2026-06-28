"""Charge and radical refresh rules for fallback."""

from __future__ import annotations

from typing import cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.utils import consts


def assign_radical_dots(atom: ob.OBAtom) -> int:
    """Estimate how many radical electrons are needed to satisfy the atom valence."""
    if atom.GetSpinMultiplicity() != 0:
        return atom.GetSpinMultiplicity()

    return max(
        0,
        cast(
            int,
            ob.GetTypicalValence(
                atom.GetAtomicNum(), atom.GetTotalValence(), atom.GetFormalCharge()
            ),
        )
        - cast(int, atom.GetTotalValence()),
    )


def assign_charge_radical_for_atom(atom: ob.OBAtom) -> bool:
    """Refresh one atom's formal charge / radical count from its local valence state."""

    old_charge = atom.GetFormalCharge()
    old_spin = atom.GetSpinMultiplicity()
    if assign_radical_dots(atom):
        atom.SetSpinMultiplicity(assign_radical_dots(atom))
    else:
        if (
            consts.NON_METAL_DICT[atom.GetAtomicNum()].num_outer_electrons == 3
            and atom.GetTotalValence() == 4
        ):
            atom.SetFormalCharge(-1)
        else:
            low_valence_total_elec = (
                consts.NON_METAL_DICT[atom.GetAtomicNum()].num_outer_electrons
                + atom.GetTotalValence()
                + atom.GetSpinMultiplicity()
                - atom.GetFormalCharge()
            ) % 8
            high_valence_total_elec = (
                consts.NON_METAL_DICT[atom.GetAtomicNum()].num_outer_electrons
                - atom.GetTotalValence()
                + atom.GetSpinMultiplicity()
                - atom.GetFormalCharge()
            ) % 2
            if low_valence_total_elec == 0:
                return False
            if low_valence_total_elec <= high_valence_total_elec:
                atom.SetFormalCharge(low_valence_total_elec)
            else:
                atom.SetSpinMultiplicity(
                    consts.NON_METAL_DICT[atom.GetAtomicNum()].num_outer_electrons
                    - atom.GetTotalValence()
                    + atom.GetSpinMultiplicity()
                    - atom.GetFormalCharge()
                )
    return old_charge != atom.GetFormalCharge() or old_spin != atom.GetSpinMultiplicity()


def fresh_omol_charge_radical(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    """Refresh formal charges and radical counts for the whole molecule."""

    hit = False
    for atom in omol.atoms:
        hit = assign_charge_radical_for_atom(atom.OBAtom) or hit
    return omol, hit


__all__ = [
    "assign_charge_radical_for_atom",
    "assign_radical_dots",
    "fresh_omol_charge_radical",
]
