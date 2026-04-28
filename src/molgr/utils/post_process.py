"""
Author: TMJ
Date: 2026-02-27 23:35:43
LastEditors: TMJ
LastEditTime: 2026-04-28 23:23:34
Description: 请填写简介
"""

from __future__ import annotations

from typing import Optional

import numpy as np
from rdkit import Chem

from molgr.config import MolGRConfig, resolve_config
from molgr.fallback.utils.consts import NON_METAL_DICT


pt = Chem.GetPeriodicTable()


def make_dative_bond(
    rdmol: Chem.Mol,
    extra_tolerance: Optional[float] = None,
    *,
    config: MolGRConfig | None = None,
) -> Chem.Mol:
    if extra_tolerance is None:
        extra_tolerance = resolve_config(
            config
        ).metal_scoring.metal_coordination_extra_tolerance_angstrom
    metal_atom_ids = [
        atom_id
        for atom_id in range(rdmol.GetNumAtoms())
        if rdmol.GetAtomWithIdx(atom_id).GetAtomicNum() not in NON_METAL_DICT
    ]
    if not metal_atom_ids:
        return rdmol

    rwmol = Chem.RWMol(rdmol)
    distance_matrix = Chem.Get3DDistanceMatrix(rwmol)

    for metal_atom_id in metal_atom_ids:
        metal_atom = rwmol.GetAtomWithIdx(metal_atom_id)
        distance_to_non_metal_atoms = distance_matrix[metal_atom_id]
        close_non_metal_atom_ids = np.argsort(distance_to_non_metal_atoms)

        for close_non_metal_atom_id in close_non_metal_atom_ids:
            if (
                rwmol.GetBondBetweenAtoms(int(close_non_metal_atom_id), int(metal_atom_id))
                is not None
            ):
                continue
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
                if non_metal_atom.GetAtomicNum() == 1:
                    rwmol.AddBond(
                        int(close_non_metal_atom_id), int(metal_atom_id), Chem.BondType.SINGLE
                    )
                    non_metal_atom.SetFormalCharge(0)
                    metal_atom.SetFormalCharge(metal_atom.GetFormalCharge() - 1)
                else:
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


def make_stereochemistry(rdmol: Chem.Mol) -> Chem.Mol:
    rdmol.UpdatePropertyCache(strict=False)
    Chem.SetAromaticity(rdmol)
    Chem.DetectBondStereochemistry(rdmol)
    Chem.SetBondStereoFromDirections(rdmol)
    Chem.AssignAtomChiralTagsFromStructure(rdmol)
    Chem.AssignStereochemistryFrom3D(rdmol)
    Chem.AssignCIPLabels(rdmol)
    for bond_idx in range(rdmol.GetNumBonds()):
        rd_bond = rdmol.GetBondWithIdx(bond_idx)
        if rd_bond.GetStereo() == Chem.BondStereo.STEREONONE:
            rd_bond.SetBondDir(Chem.BondDir.NONE)
    return rdmol
