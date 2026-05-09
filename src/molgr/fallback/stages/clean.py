"""Radical/resonance cleanup stages for fallback."""

from __future__ import annotations

from typing import List, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.stages.fresh import fresh_omol_charge_radical
from molgr.fallback.utils import consts, smarts


def clean_carbene_neighbor_unsaturated(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    """Shift a carbene radical toward an adjacent unsaturated bond when favorable."""

    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int]] = list(smarts.CLEAN_CARBENE_NEIGHBOR_UNSAT.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        atom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        if atom1.GetSpinMultiplicity() == 2 and atom3.GetSpinMultiplicity() == 0:
            cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2])).SetBondOrder(
                int(cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2])).GetBondOrder() - 1)
            )
            cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1])).SetBondOrder(
                int(cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1])).GetBondOrder() + 1)
            )
            atom1.SetSpinMultiplicity(atom1.GetSpinMultiplicity() - 1)
            atom3.SetSpinMultiplicity(atom3.GetSpinMultiplicity() + 1)
            hit = True
    return omol, hit


def clean_neighbor_radicals(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    """Fuse adjacent radicals into higher bond order whenever both endpoints allow it."""

    hit = False
    for bond in list(ob.OBMolBondIter(omol.OBMol)):
        begin_atom = cast(ob.OBAtom, bond.GetBeginAtom())
        end_atom = cast(ob.OBAtom, bond.GetEndAtom())
        if begin_atom.GetSpinMultiplicity() and end_atom.GetSpinMultiplicity():
            bond_to_add = min(begin_atom.GetSpinMultiplicity(), end_atom.GetSpinMultiplicity())
            bond.SetBondOrder(bond.GetBondOrder() + bond_to_add)
            begin_atom.SetSpinMultiplicity(begin_atom.GetSpinMultiplicity() - bond_to_add)
            end_atom.SetSpinMultiplicity(end_atom.GetSpinMultiplicity() - bond_to_add)
            if bond_to_add > 0:
                hit = True
    return omol, hit


def clean_resonances_0(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int, int]] = list(smarts.CLEAN_RESONANCE_0.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom4 = cast(ob.OBAtom, obmol.GetAtom(idxs[-1]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        obbond3 = cast(ob.OBBond, obmol.GetBond(idxs[2], idxs[3]))
        if (
            consts.NON_METAL_DICT[obatom1.GetAtomicNum()].default_valence
            > obatom1.GetTotalValence()
            and consts.NON_METAL_DICT[obatom4.GetAtomicNum()].default_valence
            > obatom4.GetTotalValence()
        ):
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obbond3.SetBondOrder(obbond3.GetBondOrder() + 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() + 1)
            obatom4.SetFormalCharge(obatom4.GetFormalCharge() - 1)
            hit = True
    return omol, hit


def clean_resonances_1(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int]] = list(smarts.CLEAN_RESONANCE_1.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[-1]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
        obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
        obatom1.SetFormalCharge(obatom1.GetFormalCharge() + 1)
        obatom3.SetFormalCharge(obatom3.GetFormalCharge() - 1)
        hit = True
    return omol, hit


def clean_resonances_2(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int, int, int, int]] = list(smarts.CLEAN_RESONANCE_2.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom5 = cast(ob.OBAtom, obmol.GetAtom(idxs[-1]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[3]))
        obbond3 = cast(ob.OBBond, obmol.GetBond(idxs[3], idxs[4]))
        obbond4 = cast(ob.OBBond, obmol.GetBond(idxs[4], idxs[5]))
        if (
            obbond4.GetBondOrder() == 1
            and obbond3.GetBondOrder() == 2
            and obbond2.GetBondOrder() == 1
            and obbond1.GetBondOrder() == 2
            and obatom1.GetFormalCharge() == 0
            and obatom5.GetFormalCharge() == -1
        ):
            obbond4.SetBondOrder(obbond4.GetBondOrder() + 1)
            obbond3.SetBondOrder(obbond3.GetBondOrder() - 1)
            obbond2.SetBondOrder(obbond2.GetBondOrder() + 1)
            obbond1.SetBondOrder(obbond1.GetBondOrder() - 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() - 1)
            obatom5.SetFormalCharge(obatom5.GetFormalCharge() + 1)
            hit = True
    return omol, hit


def clean_resonances_3(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int, int, int]] = list(smarts.CLEAN_RESONANCE_3.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom5 = cast(ob.OBAtom, obmol.GetAtom(idxs[-1]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        obbond3 = cast(ob.OBBond, obmol.GetBond(idxs[2], idxs[3]))
        obbond4 = cast(ob.OBBond, obmol.GetBond(idxs[3], idxs[4]))
        obbond1.SetBondOrder(obbond1.GetBondOrder() - 1)
        obbond2.SetBondOrder(obbond2.GetBondOrder() + 1)
        obbond3.SetBondOrder(obbond3.GetBondOrder() - 1)
        obbond4.SetBondOrder(obbond4.GetBondOrder() + 1)
        obatom1.SetFormalCharge(obatom1.GetFormalCharge() - 1)
        obatom5.SetFormalCharge(obatom5.GetFormalCharge() + 1)
        omol, _ = fresh_omol_charge_radical(omol)
        hit = True
    return omol, hit


def clean_resonances_4(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int]] = list(smarts.CLEAN_RESONANCE_4.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[-1]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        obbond2.SetBondOrder(obbond2.GetBondOrder() + 1)
        obbond1.SetBondOrder(obbond1.GetBondOrder() - 1)
        obatom1.SetFormalCharge(obatom1.GetFormalCharge() - 1)
        obatom3.SetFormalCharge(obatom3.GetFormalCharge() + 1)
        hit = True
    return omol, hit


def clean_resonances_5(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int]] = list(smarts.CLEAN_RESONANCE_5.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[-1]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if (
            obbond1.GetBondOrder() == 2
            and obbond2.GetBondOrder() == 1
            and obatom3.GetFormalCharge() == -1
            and obatom1.GetFormalCharge() == 0
        ):
            obbond2.SetBondOrder(obbond2.GetBondOrder() + 1)
            obbond1.SetBondOrder(obbond1.GetBondOrder() - 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() - 1)
            obatom3.SetFormalCharge(obatom3.GetFormalCharge() + 1)
            hit = True
    return omol, hit


def clean_resonances_6(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int]] = list(smarts.CLEAN_RESONANCE_6.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[-1]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        obbond2.SetBondOrder(obbond2.GetBondOrder() + 1)
        obbond1.SetBondOrder(obbond1.GetBondOrder() - 1)
        obatom1.SetFormalCharge(obatom1.GetFormalCharge() - 1)
        obatom3.SetFormalCharge(obatom3.GetFormalCharge() + 1)
        hit = True
    return omol, hit


def clean_resonances_7(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int, int, int, int, int]] = list(
        smarts.CLEAN_RESONANCE_7.findall(omol)
    )
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if (
            obbond1.GetBondOrder() == 1
            and obbond2.GetBondOrder() == 2
            and obatom1.GetFormalCharge() == -1
            and obatom3.GetFormalCharge() == 0
        ):
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() + 1)
            obatom3.SetFormalCharge(obatom3.GetFormalCharge() - 1)
            hit = True
    return omol, hit


def clean_resonances_8(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    obmol.SetAromaticPerceived(False)
    res = smarts.CLEAN_RESONANCE_8.findall(omol)
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom5 = cast(ob.OBAtom, obmol.GetAtom(idxs[4]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        obbond3 = cast(ob.OBBond, obmol.GetBond(idxs[2], idxs[3]))
        obbond4 = cast(ob.OBBond, obmol.GetBond(idxs[3], idxs[4]))
        if (
            obbond1.GetBondOrder() == 1
            and obbond2.GetBondOrder() == 2
            and obbond3.GetBondOrder() == 1
            and obbond4.GetBondOrder() == 2
        ):
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obbond3.SetBondOrder(obbond3.GetBondOrder() + 1)
            obbond4.SetBondOrder(obbond4.GetBondOrder() - 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() + 1)
            obatom5.SetFormalCharge(obatom5.GetFormalCharge() - 1)
            hit = True
    return omol, hit


def clean_resonances_9(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int]] = list(smarts.CLEAN_RESONANCE_9.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        if (
            consts.NON_METAL_DICT[obatom1.GetAtomicNum()].default_valence
            - obatom1.GetTotalValence()
            >= 1
            and consts.NON_METAL_DICT[obatom2.GetAtomicNum()].default_valence
            - obatom2.GetTotalValence()
            >= 1
        ):
            bond_to_add = min(
                consts.NON_METAL_DICT[obatom1.GetAtomicNum()].default_valence
                - obatom1.GetTotalValence(),
                consts.NON_METAL_DICT[obatom2.GetAtomicNum()].default_valence
                - obatom2.GetTotalValence(),
            )
            obbond1.SetBondOrder(obbond1.GetBondOrder() + bond_to_add)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() - bond_to_add)
            obatom2.SetFormalCharge(obatom2.GetFormalCharge() + bond_to_add)
            hit = True
    return omol, hit


def clean_resonances_10(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int, int]] = list(smarts.CLEAN_RESONANCE_10.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[-1]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        obbond3 = cast(ob.OBBond, obmol.GetBond(idxs[2], idxs[3]))
        if obatom1.GetSpinMultiplicity() == 1 and obatom3.GetSpinMultiplicity() == 1:
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond3.SetBondOrder(obbond3.GetBondOrder() + 1)
            obatom1.SetSpinMultiplicity(obatom1.GetSpinMultiplicity() - 1)
            obatom3.SetSpinMultiplicity(obatom3.GetSpinMultiplicity() - 1)
            hit = True
    return omol, hit


def clean_resonances_11(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int]] = list(smarts.CLEAN_RESONANCE_11.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        if consts.NON_METAL_DICT[
            obatom2.GetAtomicNum()
        ].default_valence - obatom2.GetTotalValence() >= 1 and (
            obbond1.GetBondOrder() == 1 or obbond1.GetBondOrder() == 2
        ):
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() + 1)
            obatom2.SetFormalCharge(obatom2.GetFormalCharge() - 1)
            hit = True
    return omol, hit


def clean_resonances_12(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int, int]] = list(smarts.CLEAN_RESONANCE_12.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom4 = cast(ob.OBAtom, obmol.GetAtom(idxs[3]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        obbond3 = cast(ob.OBBond, obmol.GetBond(idxs[2], idxs[3]))
        if consts.NON_METAL_DICT[
            obatom4.GetAtomicNum()
        ].default_valence - obatom4.GetTotalValence() >= 1 and (
            obbond1.GetBondOrder() == 1
            and obbond2.GetBondOrder() == 2
            and obbond3.GetBondOrder() == 1
        ):
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obbond3.SetBondOrder(obbond3.GetBondOrder() + 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() + 1)
            obatom4.SetFormalCharge(obatom4.GetFormalCharge() - 1)
            hit = True
    return omol, hit


def clean_resonances_13(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int]] = list(smarts.CLEAN_RESONANCE_13.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if (
            consts.NON_METAL_DICT[obatom1.GetAtomicNum()].default_valence
            - obatom1.GetTotalValence()
            >= 1
            and obbond1.GetBondOrder() == 1
        ):
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() + 1)
            obatom3.SetFormalCharge(obatom3.GetFormalCharge() - 1)
            hit = True
    return omol, hit


def clean_resonances(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    """Run the ordered resonance normalization rule set after candidate generation."""

    hit = False
    omol, stage_hit = clean_resonances_11(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_0(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_1(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_2(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_3(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_4(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_9(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_5(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_6(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_7(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_8(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_9(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_10(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_12(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_13(omol)
    hit = stage_hit or hit
    return omol, hit


__all__ = ["clean_carbene_neighbor_unsaturated", "clean_neighbor_radicals", "clean_resonances"]
