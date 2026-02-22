"""
Author: TMJ
Date: 2026-02-21 22:30:55
LastEditors: TMJ
LastEditTime: 2026-02-22 16:54:57
Description: 请填写简介
"""

from __future__ import annotations

from typing import List, Optional, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.pipeline.resonance import get_radical_resonances, process_resonance
from molgr.fallback.stages.break_bond import break_deformed_ene, break_one_bond
from molgr.fallback.stages.clean import (
    clean_carbene_neighbor_unsaturated,
    clean_neighbor_radicals,
    clean_resonances,
)
from molgr.fallback.stages.eliminate import (
    eliminate_carbene_neighbor_heteroatom,
    eliminate_carboxyl,
    eliminate_charge_spliting,
    eliminate_CN_in_doubt,
    eliminate_high_positive_charge_atoms,
    eliminate_NNN,
)
from molgr.fallback.stages.fresh import fresh_omol_charge_radical
from molgr.fallback.stages.preprocess import (
    make_connections,
    pre_clean,
    validate_omol,
)
from molgr.fallback.utils.scoring import omol_score
from molgr.fallback.utils.tools import typed_lru_cache


@typed_lru_cache(maxsize=1024, typed=True)
def xyz_to_omol_no_metal(
    xyz_block: str, total_charge: int = 0, total_radical_electrons: int = 0
) -> Optional[pybel.Molecule]:
    if total_radical_electrons < 0:
        return None

    omol = pybel.readstring("xyz", xyz_block)
    omol = make_connections(omol)
    omol = pre_clean(omol)
    omol = fresh_omol_charge_radical(omol)

    given_charge = total_charge - sum(
        cast(ob.OBAtom, atom.OBAtom).GetFormalCharge() for atom in omol.atoms
    )

    omol, given_charge = eliminate_NNN(omol, given_charge)
    omol, given_charge = eliminate_high_positive_charge_atoms(omol, given_charge)
    omol, given_charge = eliminate_CN_in_doubt(omol, given_charge)
    omol, given_charge = eliminate_carboxyl(omol, given_charge)
    omol = clean_carbene_neighbor_unsaturated(omol)
    omol, given_charge = eliminate_carbene_neighbor_heteroatom(omol, given_charge)
    omol = clean_neighbor_radicals(omol)
    omol = clean_carbene_neighbor_unsaturated(omol)
    omol, given_charge = eliminate_charge_spliting(omol, given_charge)
    omol = break_deformed_ene(omol, given_charge, total_radical_electrons)
    omol, given_charge = break_one_bond(omol, given_charge, total_radical_electrons)
    omol = fresh_omol_charge_radical(omol)

    if validate_omol(omol, total_charge, total_radical_electrons):
        return clean_resonances(omol)
    possible_resonances = get_radical_resonances(omol)
    recovered_resonances: List[pybel.Molecule] = []
    for resonance in possible_resonances:
        charge = given_charge
        resonance, charge = process_resonance(resonance, charge)
        if validate_omol(resonance, total_charge, total_radical_electrons):
            recovered_resonances.append(resonance)
    if len(recovered_resonances) == 0:
        return None
    scored_resonances = [(omol_score(res), res) for res in recovered_resonances]
    scored_resonances.sort(key=lambda x: x[0])
    return scored_resonances[0][1]


__all__ = [
    "xyz_to_omol_no_metal",
]
