"""Bond-breaking heuristics used when the deterministic cleanup path still stalls."""

from __future__ import annotations

from typing import List, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.utils import smarts


def break_deformed_ene(
    omol: pybel.Molecule,
    given_charge: int = 0,
    given_radical: int = 0,
    tolerance: float = 10,
) -> tuple[pybel.Molecule, bool]:
    """Break strongly deformed ene-like bonds to seed new radical sites."""

    possible_ene_pairs: List[Tuple[int, int, float]] = []
    obmol = cast(ob.OBMol, omol.OBMol)
    matches: List[Tuple[int, int, int, int]] = list(smarts.BREAK_DEFORMED_ENE_A.findall(omol))
    while len(matches):
        idxs = matches.pop(0)
        bond2_1 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if bond2_1.IsRotor() or bond2_1.GetBondOrder() == 1:
            continue
        torsion_angle: float = abs(obmol.GetTorsion(*idxs))
        torsion_angle = min(torsion_angle, 180 - torsion_angle)
        if torsion_angle > tolerance:
            possible_ene_pairs.append((idxs[1], idxs[2], torsion_angle))

    matches = list(smarts.BREAK_DEFORMED_ENE_B.findall(omol))
    while len(matches):
        idxs = matches.pop(0)
        bond2_2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if bond2_2.IsRotor() or bond2_2.GetBondOrder() == 1:
            continue
        torsion_angle = abs(obmol.GetTorsion(*idxs))
        torsion_angle = min(torsion_angle, 180 - torsion_angle)
        if torsion_angle > tolerance:
            possible_ene_pairs.append((idxs[1], idxs[2], torsion_angle))

    possible_ene_pairs.sort(key=lambda x: x[2], reverse=True)
    hit = False
    for idx1, idx2, _ in possible_ene_pairs:
        if (
            sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(obmol))
            >= abs(given_charge) + given_radical
        ):
            return omol, hit
        bond = cast(ob.OBBond, obmol.GetBond(idx1, idx2))
        if bond.IsRotor() or bond.GetBondOrder() == 1:
            continue
        bond.SetBondOrder(bond.GetBondOrder() - 1)
        begin_atom = cast(ob.OBAtom, bond.GetBeginAtom())
        end_atom = cast(ob.OBAtom, bond.GetEndAtom())
        begin_atom.SetSpinMultiplicity(begin_atom.GetSpinMultiplicity() + 1)
        end_atom.SetSpinMultiplicity(end_atom.GetSpinMultiplicity() + 1)
        hit = True
    return omol, hit


def break_one_bond(
    omol: pybel.Molecule,
    given_charge: int = 0,
    given_radical: int = 0,
) -> tuple[pybel.Molecule, int, bool]:
    """Apply the last-resort bond-breaking rules that create radical candidates."""

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False
    while res := smarts.BREAK_ONE_BOND_MULTIPLE.findall(omol):
        if (
            sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(obmol))
            >= abs(given_charge) + given_radical
        ):
            return omol, given_charge, hit
        idxs = res.pop(0)
        bond = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        bond.SetBondOrder(bond.GetBondOrder() - 1)
        begin_atom = cast(ob.OBAtom, bond.GetBeginAtom())
        end_atom = cast(ob.OBAtom, bond.GetEndAtom())
        begin_atom.SetSpinMultiplicity(begin_atom.GetSpinMultiplicity() + 1)
        end_atom.SetSpinMultiplicity(end_atom.GetSpinMultiplicity() + 1)
        hit = True

    while res := smarts.BREAK_ONE_BOND_CATION.findall(omol):
        if (
            sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(obmol))
            >= abs(given_charge) + given_radical
        ):
            return omol, given_charge, hit
        idxs = res.pop(0)
        bond2 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        bond2.SetBondOrder(bond2.GetBondOrder() - 1)
        begin_atom2 = cast(ob.OBAtom, bond2.GetBeginAtom())
        end_atom2 = cast(ob.OBAtom, bond2.GetEndAtom())
        end_atom2.SetSpinMultiplicity(end_atom2.GetSpinMultiplicity() + 1)
        begin_atom2.SetFormalCharge(int(begin_atom2.GetFormalCharge() - 1))
        given_charge += 1
        hit = True

    res = list(smarts.BREAK_ONE_BOND_AROMATIC.findall(omol))
    if res:
        if (
            sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(obmol))
            >= abs(given_charge) + given_radical
        ):
            return omol, given_charge, hit
        for idxs in res:
            bond3 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
            if bond3.GetBondOrder() == 1:
                continue
            bond3.SetBondOrder(bond3.GetBondOrder() - 1)
            begin_atom3 = cast(ob.OBAtom, bond3.GetBeginAtom())
            end_atom3 = cast(ob.OBAtom, bond3.GetEndAtom())
            begin_atom3.SetSpinMultiplicity(begin_atom3.GetSpinMultiplicity() + 1)
            end_atom3.SetSpinMultiplicity(end_atom3.GetSpinMultiplicity() + 1)
            hit = True
            break

    if all(cast(ob.OBBond, bond).GetBondOrder() == 1 for bond in ob.OBMolBondIter(obmol)):
        single_bonds = list(ob.OBMolBondIter(obmol))
        if single_bonds:
            if (
                sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(obmol))
                >= abs(given_charge) + given_radical
            ):
                return omol, given_charge, hit
            single_bond = cast(ob.OBBond, single_bonds[0])
            single_begin_atom = cast(ob.OBAtom, single_bond.GetBeginAtom())
            single_end_atom = cast(ob.OBAtom, single_bond.GetEndAtom())
            single_begin_atom.SetSpinMultiplicity(single_begin_atom.GetSpinMultiplicity() + 1)
            single_end_atom.SetSpinMultiplicity(single_end_atom.GetSpinMultiplicity() + 1)
            obmol.DeleteBond(single_bond)
            hit = True
    return omol, given_charge, hit


__all__ = ["break_deformed_ene", "break_one_bond"]
