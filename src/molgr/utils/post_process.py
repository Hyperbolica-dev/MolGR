"""
Author: TMJ
Date: 2026-02-27 23:35:43
LastEditors: TMJ
LastEditTime: 2026-02-28 00:43:00
Description: 请填写简介
"""

import numpy as np
from rdkit import Chem

from molgr.fallback.utils.consts import NON_METAL_DICT


pt = Chem.GetPeriodicTable()


def make_dative_bond(rdmol: Chem.Mol, extra_tolerance: float = 0.35) -> Chem.Mol:
    rwmol = Chem.RWMol(rdmol)
    metal_atom_ids = [
        atom_id
        for atom_id in range(rwmol.GetNumAtoms())
        if rwmol.GetAtomWithIdx(atom_id).GetAtomicNum() not in NON_METAL_DICT
    ]
    distance_matrix = Chem.Get3DDistanceMatrix(rwmol)

    for metal_atom_id in metal_atom_ids:
        metal_atom = rwmol.GetAtomWithIdx(metal_atom_id)
        distance_to_non_metal_atoms = distance_matrix[metal_atom_id]
        close_non_metal_atom_ids = np.argsort(distance_to_non_metal_atoms)

        for close_non_metal_atom_id in close_non_metal_atom_ids:
            non_metal_atom = rwmol.GetAtomWithIdx(int(close_non_metal_atom_id))
            if non_metal_atom.GetAtomicNum() not in NON_METAL_DICT:
                continue
            if distance_to_non_metal_atoms[
                close_non_metal_atom_id
            ] > extra_tolerance + pt.GetRcovalent(non_metal_atom.GetAtomicNum()) + pt.GetRcovalent(
                metal_atom.GetAtomicNum()
            ):
                break
            if (
                pt.GetNOuterElecs(non_metal_atom.GetAtomicNum()) == 3
                and non_metal_atom.GetFormalCharge() < 0
            ):
                continue
            if non_metal_atom.GetFormalCharge() < 0:
                rwmol.AddBond(
                    int(close_non_metal_atom_id), int(metal_atom_id), Chem.BondType.DATIVE
                )
            if (
                pt.GetNOuterElecs(non_metal_atom.GetAtomicNum()) in (5, 6, 7)
                and non_metal_atom.GetFormalCharge() == 0
            ):
                rwmol.AddBond(
                    int(close_non_metal_atom_id), int(metal_atom_id), Chem.BondType.DATIVE
                )
    return rwmol.GetMol()
