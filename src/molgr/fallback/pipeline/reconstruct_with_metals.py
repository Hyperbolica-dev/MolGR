"""
Author: TMJ
Date: 2026-02-21 22:30:55
LastEditors: TMJ
LastEditTime: 2026-03-16 22:31:07
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
    f_d_s_p = consts.METAL_F_D_S_P_ELECTRONS.get(metal)
    if f_d_s_p is None:
        return set()

    f, d, s, p = f_d_s_p.f, f_d_s_p.d, f_d_s_p.s, f_d_s_p.p

    if valence <= s + p:
        if 0 <= d < len(consts.D_ELECTRONS_SPIN):
            base = (f + s + p - valence) % 2
            return {base + dd for dd in consts.D_ELECTRONS_SPIN[d]}
        return set()

    if valence <= s + p + d:
        idx = d - valence + s + p
        if 0 <= idx < len(consts.D_ELECTRONS_SPIN):
            return {f % 2 + dd for dd in consts.D_ELECTRONS_SPIN[idx]}
        return set()

    if valence <= s + p + d + f:
        return {f % 2}

    return set()


def _build_metal_states(obatom: ob.OBAtom) -> List[dataclasses.MetalAtomPosition]:
    symbol = ob.GetSymbol(obatom.GetAtomicNum())

    def _default_state() -> dataclasses.MetalAtomPosition:
        return dataclasses.MetalAtomPosition(
            idx=obatom.GetIdx(),
            symbol=symbol,
            element_idx=obatom.GetAtomicNum(),
            valence=0,
            radical_num=0,
            position_x=obatom.GetX(),
            position_y=obatom.GetY(),
            position_z=obatom.GetZ(),
        )

    prior = consts.METAL_VALENCE_AVAILABLE_PRIOR.get(symbol, [])
    minor = consts.METAL_VALENCE_AVAILABLE_MINOR.get(symbol, [])
    seen_valences = set()
    valences: List[int] = []
    for valence in prior + minor:
        if valence in seen_valences:
            continue
        seen_valences.add(valence)
        valences.append(valence)
    if not valences:
        valences = [0]

    if symbol not in consts.METAL_F_D_S_P_ELECTRONS:
        return [_default_state()]

    states: List[dataclasses.MetalAtomPosition] = []
    for valence in valences:
        try:
            radicals = get_possible_metal_radicals(symbol, valence)
        except ValueError:
            continue
        for radical_num in sorted(radicals):
            states.append(
                dataclasses.MetalAtomPosition(
                    idx=obatom.GetIdx(),
                    symbol=symbol,
                    element_idx=obatom.GetAtomicNum(),
                    valence=valence,
                    radical_num=radical_num,
                    position_x=obatom.GetX(),
                    position_y=obatom.GetY(),
                    position_z=obatom.GetZ(),
                )
            )

    if not states:
        return [_default_state()]
    return states


def combine_metal_with_omol(
    omol: pybel.Molecule, metal_list: Sequence[dataclasses.MetalAtomPosition]
) -> pybel.Molecule:
    obmol = cast(ob.OBMol, omol.clone.OBMol)
    obmol.BeginModify()
    try:
        num_organic = obmol.NumAtoms()
        num_metals = len(metal_list)
        total_atoms = num_organic + num_metals

        for metal in metal_list:
            atom = cast(ob.OBAtom, obmol.NewAtom())
            atom.SetAtomicNum(metal.element_idx)
            atom.SetFormalCharge(metal.valence)
            atom.SetSpinMultiplicity(metal.radical_num)
            atom.SetVector(metal.position_x, metal.position_y, metal.position_z)

        new_order = [0] * total_atoms
        has_error = False
        for i, metal in enumerate(metal_list):
            current_idx = num_organic + 1 + i
            target_slot = metal.idx - 1
            if target_slot < 0 or target_slot >= total_atoms:
                has_error = True
                continue
            if new_order[target_slot] != 0:
                has_error = True
                continue
            new_order[target_slot] = current_idx

        if not has_error:
            current_organic_idx = 1
            for i in range(total_atoms):
                if new_order[i] != 0:
                    continue
                if current_organic_idx > num_organic:
                    has_error = True
                    break
                new_order[i] = current_organic_idx
                current_organic_idx += 1

        if not has_error and all(idx > 0 for idx in new_order):
            obmol.RenumberAtoms(new_order)
    finally:
        obmol.EndModify()
    return pybel.Molecule(obmol)


def xyz2omol(
    xyz_block: str, total_charge: int = 0, total_radical_electrons: int = 0
) -> Optional[pybel.Molecule]:
    omol = pybel.readstring("xyz", xyz_block)
    removable_metal_atoms: List[ob.OBAtom] = [
        cast(ob.OBAtom, atom.OBAtom) for atom in omol.atoms if atom.OBAtom.IsMetal()
    ]
    available_valence_radical_states = [
        _build_metal_states(obatom) for obatom in removable_metal_atoms
    ]
    for obatom in removable_metal_atoms:
        omol.OBMol.DeleteAtom(obatom)
    no_metal_xyz = cast(str, omol.write("xyz"))
    possible_omols: List[pybel.Molecule] = []
    possible_metal_valence_radical_product = (
        product for product in itertools.product(*available_valence_radical_states)
    )
    for metal_atom_product in possible_metal_valence_radical_product:
        total_metal_charge = sum(metal_pos.valence for metal_pos in metal_atom_product)
        total_metal_radical_electrons = sum(
            metal_pos.radical_num for metal_pos in metal_atom_product
        )
        try:
            possible_omol = xyz_to_omol_no_metal(
                no_metal_xyz,
                total_charge - total_metal_charge,
                total_radical_electrons - total_metal_radical_electrons,
            )
        except (OSError, ValueError):
            continue

        if possible_omol is None:
            try:
                possible_omol = xyz_to_omol_no_metal(
                    no_metal_xyz,
                    total_charge - total_metal_charge,
                    total_radical_electrons - total_metal_radical_electrons % 2,
                )
            except (OSError, ValueError):
                continue
        if possible_omol is None:
            continue
        possible_omols.append(combine_metal_with_omol(possible_omol, metal_atom_product))
    if not possible_omols:
        return None
    scored_omols = [(scoring.omol_score(res), res) for res in possible_omols]
    scored_omols.sort(key=lambda x: x[0])
    return scored_omols[0][1]


__all__ = ["xyz2omol"]
