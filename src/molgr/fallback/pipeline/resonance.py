"""
Author: TMJ
Date: 2026-02-23 22:15:39
LastEditors: TMJ
LastEditTime: 2026-02-28 02:43:11
Description: 请填写简介
"""
# pyright: reportMissingImports=false, reportAttributeAccessIssue=false

from __future__ import annotations

from typing import List, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.stages.clean import clean_neighbor_radicals, clean_resonances
from molgr.fallback.stages.eliminate import (
    eliminate_1_3_dipole,
    eliminate_negative_charges,
    eliminate_positive_charges,
)
from molgr.fallback.utils.tools import typed_lru_cache


@typed_lru_cache(maxsize=1024, typed=True)
def process_resonance(resonance: pybel.Molecule, charge: int) -> Tuple[pybel.Molecule, int]:
    resonance, charge = eliminate_1_3_dipole(resonance, charge)
    resonance, charge = eliminate_positive_charges(resonance, charge)
    resonance, charge = eliminate_negative_charges(resonance, charge)
    resonance = clean_neighbor_radicals(resonance)
    resonance = clean_resonances(resonance)
    return resonance, charge


def get_one_step_resonance(omol: pybel.Molecule) -> List[pybel.Molecule]:
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[*]-,=,:[*]=,#,:[*]")
    res: List[Tuple[int, int, int]] = list(smarts.findall(omol))
    result: List[pybel.Molecule] = []
    for idxs in res:
        atom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        bond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        bond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if (
            atom1.GetSpinMultiplicity() == 1
            and atom3.GetSpinMultiplicity() == 0
            and bond1.GetBondOrder() <= 2
            and bond2.GetBondOrder() >= 2
        ):
            new_omol = omol.clone
            new_obmol = cast(ob.OBMol, new_omol.OBMol)
            atom1_clone = cast(ob.OBAtom, new_obmol.GetAtom(idxs[0]))
            atom3_clone = cast(ob.OBAtom, new_obmol.GetAtom(idxs[2]))
            bond1_clone = cast(ob.OBBond, new_obmol.GetBond(idxs[0], idxs[1]))
            bond2_clone = cast(ob.OBBond, new_obmol.GetBond(idxs[1], idxs[2]))
            bond1_clone.SetBondOrder(bond1_clone.GetBondOrder() + 1)
            bond2_clone.SetBondOrder(bond2_clone.GetBondOrder() - 1)
            atom1_clone.SetSpinMultiplicity(atom1_clone.GetSpinMultiplicity() - 1)
            atom3_clone.SetSpinMultiplicity(atom3_clone.GetSpinMultiplicity() + 1)
            result.append(new_omol)
    return result


def get_radical_resonances(omol: pybel.Molecule) -> List[pybel.Molecule]:
    resonances = {cast(str, omol.write("smi")): omol}
    new_resonances = get_one_step_resonance(omol)
    for new_resonance in new_resonances:
        resonances[cast(str, new_resonance.write("smi"))] = new_resonance
    for temp_omol in new_resonances:
        new_new_resonances = get_one_step_resonance(temp_omol)
        for new_resonance in new_new_resonances:
            resonances[cast(str, new_resonance.write("smi"))] = new_resonance
    return list(resonances.values())


__all__ = [
    "get_one_step_resonance",
    "get_radical_resonances",
    "process_resonance",
]
