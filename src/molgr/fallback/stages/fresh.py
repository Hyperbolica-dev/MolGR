from typing import cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.utils import consts


def assign_radical_dots(atom: ob.OBAtom) -> int:
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


def assign_charge_radical_for_atom(atom: ob.OBAtom):
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
                return
            if low_valence_total_elec <= high_valence_total_elec:
                atom.SetFormalCharge(low_valence_total_elec)
            else:
                atom.SetSpinMultiplicity(
                    consts.NON_METAL_DICT[atom.GetAtomicNum()].num_outer_electrons
                    - atom.GetTotalValence()
                    + atom.GetSpinMultiplicity()
                    - atom.GetFormalCharge()
                )


def fresh_omol_charge_radical(omol: pybel.Molecule) -> pybel.Molecule:
    for atom in omol.atoms:
        assign_charge_radical_for_atom(atom.OBAtom)
    return omol
