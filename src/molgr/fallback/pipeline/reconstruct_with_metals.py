"""
Author: TMJ
Date: 2026-02-21 22:30:55
LastEditors: TMJ
LastEditTime: 2026-02-22 13:29:26
Description: 请填写简介
"""

from __future__ import annotations

import itertools
from collections.abc import Sequence
from typing import List, Optional, Set, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.pipeline.reconstruct_without_metals import xyz_to_omol_no_metal
from molgr.fallback.utils import consts, dataclasses, scoring


def get_possible_metal_radicals(metal: str, valence: int) -> Set[int]:
    f_d_s_p = consts.METAL_F_D_S_P_ELECTRONS[metal]
    f, d, s, p = f_d_s_p.f, f_d_s_p.d, f_d_s_p.s, f_d_s_p.p
    if valence <= s + p:
        return {(f + s + p - valence) % 2 + dd for dd in consts.D_ELECTRONS_SPIN[d]}
    if valence <= s + p + d:
        return {f % 2 + dd for dd in consts.D_ELECTRONS_SPIN[d - valence + s + p]}
    if valence <= s + p + d + f:
        return {f % 2}
    raise ValueError("Valence is too high for this metal")


def combine_metal_with_omol(
    omol: pybel.Molecule, metal_list: Sequence[dataclasses.MetalAtomPosition]
) -> pybel.Molecule:
    obmol = cast(ob.OBMol, omol.clone.OBMol)
    for metal in metal_list:
        atom = cast(ob.OBAtom, obmol.NewAtom())
        atom.SetAtomicNum(metal.element_idx)
        atom.SetFormalCharge(metal.valence)
        atom.SetSpinMultiplicity(metal.radical_num)
        atom.SetVector(metal.position_x, metal.position_y, metal.position_z)
        obmol.RenumberAtoms(
            list(range(1, metal.idx)) + [atom.GetIdx()] + list(range(metal.idx + 1, atom.GetIdx()))
        )
    return pybel.Molecule(obmol)


def xyz2omol(
    xyz_block: str, total_charge: int = 0, total_radical_electrons: int = 0
) -> Optional[pybel.Molecule]:
    omol = pybel.readstring("xyz", xyz_block)
    removable_metal_atoms: List[ob.OBAtom] = [
        cast(ob.OBAtom, atom.OBAtom) for atom in omol.atoms if atom.OBAtom.IsMetal()
    ]
    available_valence_radical_states = [
        [
            dataclasses.MetalAtomPosition(
                idx=obatom.GetIdx(),
                symbol=ob.GetSymbol(obatom.GetAtomicNum()),
                element_idx=obatom.GetAtomicNum(),
                valence=valence,
                radical_num=radical_num,
                position_x=obatom.GetX(),
                position_y=obatom.GetY(),
                position_z=obatom.GetZ(),
            )
            for valence in consts.METAL_VALENCE_AVAILABLE_PRIOR[ob.GetSymbol(obatom.GetAtomicNum())]
            + consts.METAL_VALENCE_AVAILABLE_MINOR[ob.GetSymbol(obatom.GetAtomicNum())]
            for radical_num in get_possible_metal_radicals(
                ob.GetSymbol(obatom.GetAtomicNum()), valence
            )
        ]
        for obatom in removable_metal_atoms
    ]
    for obatom in removable_metal_atoms:
        omol.OBMol.DeleteAtom(obatom)
    no_metal_xyz = cast(str, omol.write("xyz"))
    possible_metal_valence_radical_product = (
        product
        for product in itertools.product(*available_valence_radical_states)
        if sum(metal_pos.radical_num for metal_pos in product) <= total_radical_electrons
    )
    possible_omols: List[pybel.Molecule] = []
    for metal_atom_product in possible_metal_valence_radical_product:
        total_metal_charge = sum(metal_pos.valence for metal_pos in metal_atom_product)
        total_metal_radical_electrons = sum(
            metal_pos.radical_num for metal_pos in metal_atom_product
        )
        possible_omol = xyz_to_omol_no_metal(
            no_metal_xyz,
            total_charge - total_metal_charge,
            total_radical_electrons - total_metal_radical_electrons,
        )
        if possible_omol is None:
            continue
        possible_omols.append(combine_metal_with_omol(possible_omol, metal_atom_product))
    if not possible_omols:
        return None
    scored_omols = [(scoring.omol_score(res), res) for res in possible_omols]
    scored_omols.sort(key=lambda x: x[0])
    return scored_omols[0][1]


__all__ = ["xyz2omol"]
