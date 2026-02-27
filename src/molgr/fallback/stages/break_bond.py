"""
Author: TMJ
Date: 2026-02-21 22:53:05
LastEditors: TMJ
LastEditTime: 2026-02-22 18:28:35
Description: 请填写简介
"""

from __future__ import annotations

from typing import List, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel


def break_deformed_ene(
    omol: pybel.Molecule, given_charge: int = 0, given_radical: int = 0, tolerance: float = 5
):
    """
    假设：在一个烯烃形成的平面上，取代基偏离平面可能意味着这个双键实际上是一个旋转键
    情况1：异侧检验
    情况2：同侧检验
    """
    possible_ene_pairs: List[Tuple[int, int, float]] = []
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[*]~[*+0]=,:[*+0]~[*]")
    matches: List[Tuple[int, int, int, int]] = list(smarts.findall(omol))
    while len(matches):
        idxs = matches.pop(0)
        bond2_1 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if bond2_1.IsRotor() or bond2_1.GetBondOrder() == 1:
            continue
        torsion_angle: float = abs(obmol.GetTorsion(*idxs))
        torsion_angle = min(torsion_angle, 180 - torsion_angle)
        if torsion_angle > tolerance:
            possible_ene_pairs.append((idxs[1], idxs[2], torsion_angle))

    smarts = pybel.Smarts("[*]~[*+0](=,:[*+0])~[*]")
    matches = list(smarts.findall(omol))
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
    # 优先断开偏离程度较大的键
    for idx1, idx2, _ in possible_ene_pairs:
        if (
            sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(obmol))
            >= abs(given_charge) + given_radical
        ):
            return omol
        bond = cast(ob.OBBond, obmol.GetBond(idx1, idx2))
        if bond.IsRotor() or bond.GetBondOrder() == 1:
            continue
        bond.SetBondOrder(bond.GetBondOrder() - 1)
        begin_atom = cast(ob.OBAtom, bond.GetBeginAtom())
        end_atom = cast(ob.OBAtom, bond.GetEndAtom())
        begin_atom.SetSpinMultiplicity(begin_atom.GetSpinMultiplicity() + 1)
        end_atom.SetSpinMultiplicity(end_atom.GetSpinMultiplicity() + 1)
    return omol


def break_one_bond(omol: pybel.Molecule, given_charge: int = 0, given_radical: int = 0):
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[*+0]#,=[*+0]")
    while res := smarts.findall(omol):
        if (
            sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(obmol))
            >= abs(given_charge) + given_radical
        ):
            return omol, given_charge
        idxs = res.pop(0)
        bond = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        bond.SetBondOrder(bond.GetBondOrder() - 1)
        begin_atom = cast(ob.OBAtom, bond.GetBeginAtom())
        end_atom = cast(ob.OBAtom, bond.GetEndAtom())
        begin_atom.SetSpinMultiplicity(begin_atom.GetSpinMultiplicity() + 1)
        end_atom.SetSpinMultiplicity(end_atom.GetSpinMultiplicity() + 1)

    smarts = pybel.Smarts("[#7+1,#15+1]=[*+0]")
    while res := smarts.findall(omol):
        if (
            sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(obmol))
            >= abs(given_charge) + given_radical
        ):
            return omol, given_charge
        idxs = res.pop(0)
        bond2 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        bond2.SetBondOrder(bond2.GetBondOrder() - 1)
        begin_atom2 = cast(ob.OBAtom, bond2.GetBeginAtom())
        end_atom2 = cast(ob.OBAtom, bond2.GetEndAtom())
        end_atom2.SetSpinMultiplicity(end_atom2.GetSpinMultiplicity() + 1)
        begin_atom2.SetFormalCharge(int(begin_atom2.GetFormalCharge() - 1))
        given_charge += 1

    smarts = pybel.Smarts("[*+0]:[*+0]")
    while res := smarts.findall(omol):
        if (
            sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(obmol))
            >= abs(given_charge) + given_radical
        ):
            return omol, given_charge
        idxs = res.pop(0)
        bond3 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        bond3.SetBondOrder(bond3.GetBondOrder() - 1)
        begin_atom3 = cast(ob.OBAtom, bond3.GetBeginAtom())
        end_atom3 = cast(ob.OBAtom, bond3.GetEndAtom())
        begin_atom3.SetSpinMultiplicity(begin_atom3.GetSpinMultiplicity() + 1)
        end_atom3.SetSpinMultiplicity(end_atom3.GetSpinMultiplicity() + 1)

    if all(cast(ob.OBBond, bond).GetBondOrder() == 1 for bond in ob.OBMolBondIter(obmol)):
        for single_bond in list(ob.OBMolBondIter(obmol)):
            if (
                sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(obmol))
                >= abs(given_charge) + given_radical
            ):
                return omol, given_charge
            single_begin_atom = cast(ob.OBAtom, single_bond.GetBeginAtom())
            single_end_atom = cast(ob.OBAtom, single_bond.GetEndAtom())
            single_begin_atom.SetSpinMultiplicity(single_begin_atom.GetSpinMultiplicity() + 1)
            single_end_atom.SetSpinMultiplicity(single_end_atom.GetSpinMultiplicity() + 1)
            obmol.DeleteBond(single_bond)
    return omol, given_charge


__all__ = [
    "break_deformed_ene",
    "break_one_bond",
]
