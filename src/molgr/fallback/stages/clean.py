from __future__ import annotations

from typing import List, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.stages.fresh import fresh_omol_charge_radical
from molgr.fallback.utils import consts


def clean_carbene_neighbor_unsaturated(omol: pybel.Molecule) -> pybel.Molecule:
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[*]-[*]=[*]")
    res: List[Tuple[int, int, int]] = list(smarts.findall(omol))
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
    return omol


def clean_neighbor_radicals(omol: pybel.Molecule) -> pybel.Molecule:
    for bond in list(ob.OBMolBondIter(omol.OBMol)):
        begin_atom = cast(ob.OBAtom, bond.GetBeginAtom())
        end_atom = cast(ob.OBAtom, bond.GetEndAtom())
        if begin_atom.GetSpinMultiplicity() and end_atom.GetSpinMultiplicity():
            bond_to_add = min(begin_atom.GetSpinMultiplicity(), end_atom.GetSpinMultiplicity())
            bond.SetBondOrder(bond.GetBondOrder() + bond_to_add)
            begin_atom.SetSpinMultiplicity(begin_atom.GetSpinMultiplicity() - bond_to_add)
            end_atom.SetSpinMultiplicity(end_atom.GetSpinMultiplicity() - bond_to_add)
    return omol


def clean_resonances_0(omol: pybel.Molecule) -> pybel.Molecule:
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[*-]-[*]=[*]~[*+]")
    res: List[Tuple[int, int, int, int]] = list(smarts.findall(omol))
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
    return omol


def clean_resonances_1(omol: pybel.Molecule) -> pybel.Molecule:
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[*-]=[*+]=[*+0]")
    res: List[Tuple[int, int, int]] = list(smarts.findall(omol))
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
    return omol


def clean_resonances_2(omol: pybel.Molecule) -> pybel.Molecule:
    obmol = cast(ob.OBMol, omol.OBMol)

    smarts = pybel.Smarts("[#8]=[#6](-[!-])-[*]=[*]-[#7-,#6-]")
    res: List[Tuple[int, int, int, int, int, int]] = list(smarts.findall(omol))
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
    return omol


def clean_resonances_3(omol: pybel.Molecule) -> pybel.Molecule:
    """
    净结果产生一个氮宾
    """
    obmol = cast(ob.OBMol, omol.OBMol)

    smarts = pybel.Smarts("[#7v2+]=[*]-[*]=[*]-[#8-]")
    res: List[Tuple[int, int, int, int, int]] = list(smarts.findall(omol))
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
        omol = fresh_omol_charge_radical(omol)
    return omol


def clean_resonances_4(omol: pybel.Molecule) -> pybel.Molecule:
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[#7+,#8+]=[*]-[#6-,#7-,#8-]")
    res: List[Tuple[int, int, int]] = list(smarts.findall(omol))
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
    return omol


def clean_resonances_5(omol: pybel.Molecule) -> pybel.Molecule:
    """
    1,3负离子共振，净结果产生一个相对更稳定的阴离子
    """
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[#7+0,#8+0,#16+0]=[*+0]-[#6-,#7-]")
    res: List[Tuple[int, int, int]] = list(smarts.findall(omol))
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
    return omol


def clean_resonances_6(omol: pybel.Molecule) -> pybel.Molecule:
    """
    净结果产生一个炔丙基（氰基亚甲基）负离子
    """
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[#6]=[#6]=[#6-,#7-]")
    res: List[Tuple[int, int, int]] = list(smarts.findall(omol))
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
    return omol


def clean_resonances_7(omol: pybel.Molecule) -> pybel.Molecule:
    """
    酚基邻位类型芳构化，净结果产生一个酚基（或等电子体）负离子
    """
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[*-]1-,:[*](=,:[*])-,:[*]=,:[*]-,:[*]=,:[*]1")
    res: List[Tuple[int, int, int, int, int, int, int]] = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if (
            obbond1.GetBondOrder() == 1
            and obbond2.GetBondOrder() == 2
            and obatom3.GetFormalCharge() == -1
            and obatom1.GetFormalCharge() == 0
        ):
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() + 1)
            obatom3.SetFormalCharge(obatom3.GetFormalCharge() - 1)
    return omol


def clean_resonances_8(omol: pybel.Molecule) -> pybel.Molecule:
    """
    酚基对位类型芳构化，净结果产生一个酚基（或等电子体）负离子
    """
    obmol = cast(ob.OBMol, omol.OBMol)
    obmol.SetAromaticPerceived(False)
    smarts = pybel.Smarts("[*-]1-,:[*]=,:[*]-,:[*](=,:[*])-,:[*]=,:[*]1")
    res = smarts.findall(omol)
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
    return omol


def clean_resonances_9(omol: pybel.Molecule) -> pybel.Molecule:
    """
    净结果消除相邻相反电荷对
    """
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[*+,*+2,*+3]-,=[*-,*-2,*-3]")
    res: List[Tuple[int, int]] = list(smarts.findall(omol))
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
    return omol


def clean_resonances_10(omol: pybel.Molecule) -> pybel.Molecule:
    """
    净结果消除离域的两个自由基并形成共轭
    """
    obmol = cast(ob.OBMol, omol.OBMol)

    smarts = pybel.Smarts("[*]-[*]=,#[*]-[*]")
    res: List[Tuple[int, int, int, int]] = list(smarts.findall(omol))
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
    return omol


def clean_resonances_11(omol: pybel.Molecule) -> pybel.Molecule:
    """
    净结果形成更稳定的鎓离子
    """
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[#7v3+0,#8v2+0,#16v2+0]-,=,:[*+1]")
    res: List[Tuple[int, int]] = list(smarts.findall(omol))
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
    return omol


def clean_resonances_12(omol: pybel.Molecule) -> pybel.Molecule:
    """
    通过离域形成更稳定的鎓离子
    """
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[#7v3+0,#8v2+0,#16v2+0]-,:[*]=,:[*]-,:[*+1]")
    res: List[Tuple[int, int, int, int]] = list(smarts.findall(omol))
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
    return omol


def clean_resonances_13(omol: pybel.Molecule) -> pybel.Molecule:
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[*-]:[*]=[#7+0,#8+0]")
    res: List[Tuple[int, int, int]] = list(smarts.findall(omol))
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
    return omol


def clean_resonances(omol: pybel.Molecule) -> pybel.Molecule:
    omol = clean_resonances_0(omol)
    omol = clean_resonances_1(omol)
    omol = clean_resonances_2(omol)
    omol = clean_resonances_3(omol)
    omol = clean_resonances_4(omol)
    omol = clean_resonances_9(omol)
    omol = clean_resonances_5(omol)
    omol = clean_resonances_6(omol)
    omol = clean_resonances_7(omol)
    omol = clean_resonances_8(omol)
    omol = clean_resonances_9(omol)
    omol = clean_resonances_10(omol)
    omol = clean_resonances_11(omol)
    omol = clean_resonances_12(omol)
    omol = clean_resonances_13(omol)
    return omol


__all__ = [
    "clean_carbene_neighbor_unsaturated",
    "clean_neighbor_radicals",
    "clean_resonances",
    "clean_resonances_0",
    "clean_resonances_1",
    "clean_resonances_10",
    "clean_resonances_11",
    "clean_resonances_12",
    "clean_resonances_13",
    "clean_resonances_2",
    "clean_resonances_3",
    "clean_resonances_4",
    "clean_resonances_5",
    "clean_resonances_6",
    "clean_resonances_7",
    "clean_resonances_8",
    "clean_resonances_9",
]
