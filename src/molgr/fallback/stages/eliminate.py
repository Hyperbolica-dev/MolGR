from __future__ import annotations

from typing import List, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.utils import consts


def eliminate_high_positive_charge_atoms(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[*+1,*+2,*+3]-[Ov1+0,Nv2+0,Sv1+0]")
    while res := smarts.findall(omol):
        idxs = cast(List[Tuple[int, int]], res.pop(0))
        atom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        atom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
        if (
            -sum(cast(ob.OBAtom, atom).GetFormalCharge() for atom in ob.OBAtomAtomIter(atom1))
            >= atom1.GetFormalCharge()
        ):
            break
        atom2.SetSpinMultiplicity(atom2.GetSpinMultiplicity() - 1)
        atom2.SetFormalCharge(atom2.GetFormalCharge() - 1)
        given_charge += 1
    return omol, given_charge


def eliminate_CN_in_doubt(omol: pybel.Molecule, given_charge: int) -> Tuple[pybel.Molecule, int]:
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[#6v4+0]=,#[#7v4+1,#15v4+1]")
    doubt_pair: List[Tuple[int, int]] = smarts.findall(omol)
    cn_in_doubt = len(doubt_pair)
    if cn_in_doubt % 2 == 0 and cn_in_doubt > 0:
        for atom_1_idx, atom_2_idx in doubt_pair[: cn_in_doubt // 2]:
            atom_1 = cast(ob.OBAtom, obmol.GetAtom(atom_1_idx))
            atom_2 = cast(ob.OBAtom, obmol.GetAtom(atom_2_idx))
            bond = cast(ob.OBBond, obmol.GetBond(atom_1_idx, atom_2_idx))
            atom_1.SetFormalCharge(-1)
            bond.SetBondOrder(bond.GetBondOrder() - 1)
            atom_2.SetFormalCharge(0)
            given_charge += 2
    return omol, given_charge


def eliminate_carboxyl(omol: pybel.Molecule, given_charge: int) -> Tuple[pybel.Molecule, int]:
    obmol = cast(ob.OBMol, omol.OBMol)
    smarts = pybel.Smarts("[Ov1+0]-C=O")
    res: List[Tuple[int, int, int]] = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        atom_1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        atom_1.SetSpinMultiplicity(atom_1.GetSpinMultiplicity() - 1)
        atom_1.SetFormalCharge(atom_1.GetFormalCharge() - 1)
        given_charge += 1
    return omol, given_charge


def eliminate_carbene_neighbor_heteroatom(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    for atom in omol.atoms:
        obatom = cast(ob.OBAtom, atom.OBAtom)
        if obatom.GetSpinMultiplicity() == 2:
            for neighbor in ob.OBAtomAtomIter(obatom):
                if cast(ob.OBAtom, neighbor).GetSpinMultiplicity():
                    return omol, given_charge
            for neighbor in ob.OBAtomAtomIter(obatom):
                if (
                    cast(ob.OBAtom, neighbor).GetAtomicNum() in consts.HETEROATOM
                    and cast(ob.OBAtom, neighbor).GetFormalCharge() == 0
                    and cast(ob.OBAtom, neighbor).GetSpinMultiplicity() == 0
                ):
                    bond = cast(ob.OBBond, obatom.GetBond(neighbor))
                    bond.SetBondOrder(bond.GetBondOrder() + 1)
                    obatom.SetSpinMultiplicity(0)
                    obatom.SetFormalCharge(obatom.GetFormalCharge() - 1)
                    cast(ob.OBAtom, neighbor).SetFormalCharge(
                        cast(ob.OBAtom, neighbor).GetFormalCharge() + 1
                    )
                    break
    return omol, given_charge


def eliminate_NNN(omol: pybel.Molecule, given_charge: int) -> Tuple[pybel.Molecule, int]:
    obmol = cast(ob.OBMol, omol.OBMol)

    smarts = pybel.Smarts("[#7v1+0]-[#7v2+0]-[#7v1+0]")
    while res := smarts.findall(omol):
        idxs = cast(List[Tuple[int, int, int]], res.pop(0))
        atom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        atom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
        atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        bond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        bond1.SetBondOrder(bond1.GetBondOrder() + 1)
        bond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        bond2.SetBondOrder(bond2.GetBondOrder() + 1)
        atom1.SetSpinMultiplicity(atom1.GetSpinMultiplicity() - 1)
        atom1.SetFormalCharge(atom1.GetFormalCharge() - 1)
        atom2.SetSpinMultiplicity(atom2.GetSpinMultiplicity() - 1)
        atom2.SetFormalCharge(atom2.GetFormalCharge() + 1)
        atom3.SetSpinMultiplicity(atom3.GetSpinMultiplicity() - 1)
        atom3.SetFormalCharge(atom3.GetFormalCharge() - 1)
        given_charge += 1
    return omol, given_charge


def eliminate_charge_spliting(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    obmol = cast(ob.OBMol, omol.OBMol)

    if (
        all(cast(ob.OBAtom, atom).GetFormalCharge() == 0 for atom in ob.OBMolAtomIter(obmol))
        and sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(obmol))
        >= 2
    ):
        radical_atoms: List[ob.OBAtom] = [
            atom for atom in ob.OBMolAtomIter(obmol) if cast(ob.OBAtom, atom).GetSpinMultiplicity()
        ]
        while len(radical_atoms) > abs(given_charge) + 1:
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 8:
                    atom.SetSpinMultiplicity(atom.GetSpinMultiplicity() - 1)
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                    given_charge += 1
                    radical_atoms.remove(atom)
                    break
            else:
                break
        while len(radical_atoms) > abs(given_charge) + 1:
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 7:
                    atom.SetSpinMultiplicity(atom.GetSpinMultiplicity() - 1)
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                    given_charge += 1
                    radical_atoms.remove(atom)
                    break
            else:
                break
        while len(radical_atoms) > abs(given_charge) + 1:
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 6 and not any(
                    _atom
                    for _atom in ob.OBAtomAtomIter(atom)
                    if cast(ob.OBAtom, _atom).GetAtomicNum() in consts.HETEROATOM
                ):
                    atom.SetSpinMultiplicity(atom.GetSpinMultiplicity() - 1)
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                    given_charge += 1
                    radical_atoms.remove(atom)
                    break
            else:
                break
        while len(radical_atoms) > abs(given_charge) + 1:
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 6:
                    atom.SetSpinMultiplicity(atom.GetSpinMultiplicity() - 1)
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                    given_charge += 1
                    radical_atoms.remove(atom)
                    break
            else:
                break
    return omol, given_charge


def eliminate_1_3_dipole(omol: pybel.Molecule, given_charge: int) -> Tuple[pybel.Molecule, int]:
    obmol = cast(ob.OBMol, omol.OBMol)

    smarts = pybel.Smarts("[*-1]-,=[N+0,O+0]-,=[*]")
    res: List[Tuple[int, int, int]] = list(smarts.findall(omol))
    while len(res):
        idxs = res.pop(0)
        atom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
        atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        bond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if (
            atom3.GetSpinMultiplicity()
            and consts.NON_METAL_DICT[atom2.GetAtomicNum()].num_outer_electrons
            + atom2.GetTotalValence()
            == 8
        ):
            atom2.SetFormalCharge(atom2.GetFormalCharge() + 1)
            bond2.SetBondOrder(int(bond2.GetBondOrder() + 1))
            atom3.SetSpinMultiplicity(atom3.GetSpinMultiplicity() - 1)
            given_charge -= 1
    return omol, given_charge


def eliminate_positive_charges(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    obmol = cast(ob.OBMol, omol.OBMol)

    smarts = pybel.Smarts("[Nv3+0]=[Nv2+0]")
    while given_charge > 0 and (res := smarts.findall(omol)):
        idxs_1 = cast(Tuple[int, int], res.pop(0))
        abatom = cast(ob.OBAtom, obmol.GetAtom(idxs_1[1]))
        abatom.SetSpinMultiplicity(abatom.GetSpinMultiplicity() - 1)
        abatom.SetFormalCharge(1)
        given_charge -= 1

    smarts = pybel.Smarts("[#6v3+0,#6v2+0,#1v0+0]")
    while given_charge > 0 and (res := smarts.findall(omol)):
        idxs_2 = cast(Tuple[int, int, int], res.pop(0))
        abatom2 = cast(ob.OBAtom, obmol.GetAtom(idxs_2[0]))
        abatom2.SetSpinMultiplicity(abatom2.GetSpinMultiplicity() - 1)
        abatom2.SetFormalCharge(1)
        given_charge -= 1
    for atom in omol.atoms:
        obatom = cast(ob.OBAtom, atom.OBAtom)
        if given_charge <= 0:
            break
        if obatom.GetSpinMultiplicity() >= 1 and obatom.GetFormalCharge() == 0:
            to_add = min(obatom.GetSpinMultiplicity(), given_charge)
            obatom.SetSpinMultiplicity(obatom.GetSpinMultiplicity() - to_add)
            obatom.SetFormalCharge(to_add)
            given_charge -= to_add
    return omol, given_charge


def eliminate_negative_charges(
    omol: pybel.Molecule, given_charge: int
) -> Tuple[pybel.Molecule, int]:
    obmol = cast(ob.OBMol, omol.OBMol)
    possible_heteroatoms: List[Tuple[ob.OBAtom, int]] = []
    for atom in omol.atoms:
        obatom = cast(ob.OBAtom, atom.OBAtom)
        if (
            obatom.GetAtomicNum() in consts.HETEROATOM
            and obatom.GetFormalCharge() == 0
            and obatom.GetSpinMultiplicity() >= 1
        ):
            possible_heteroatoms.append((obatom, consts.HETEROATOM.index(obatom.GetAtomicNum())))
    possible_heteroatoms.sort(key=lambda x: x[1])
    for obatom, _ in possible_heteroatoms:
        if given_charge >= 0:
            break
        to_add = min(obatom.GetSpinMultiplicity(), abs(given_charge))
        obatom.SetSpinMultiplicity(obatom.GetSpinMultiplicity() - to_add)
        obatom.SetFormalCharge(-to_add)
        given_charge += to_add

    smarts = pybel.Smarts("[#6v3+0]")
    while given_charge < 0 and (res := smarts.findall(omol)):
        idxs = cast(Tuple[int], res.pop(0))
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        to_add = min(obatom1.GetSpinMultiplicity(), abs(given_charge))
        obatom1.SetSpinMultiplicity(obatom1.GetSpinMultiplicity() - to_add)
        obatom1.SetFormalCharge(-to_add)
        given_charge += to_add

    smarts = pybel.Smarts("[#1v0+0]")
    while given_charge < 0 and (res := smarts.findall(omol)):
        idxs = cast(Tuple[int], res.pop(0))
        obatom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        to_add = min(obatom2.GetSpinMultiplicity(), abs(given_charge))
        obatom2.SetSpinMultiplicity(obatom2.GetSpinMultiplicity() - to_add)
        obatom2.SetFormalCharge(-to_add)
        given_charge += to_add

    smarts = pybel.Smarts("[#6v2+0,#6v1+0,#6v0+0]")
    while given_charge < 0 and (res := smarts.findall(omol)):
        idxs = cast(Tuple[int], res.pop(0))
        obatom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        to_add = min(obatom3.GetSpinMultiplicity(), abs(given_charge))
        obatom3.SetSpinMultiplicity(obatom3.GetSpinMultiplicity() - to_add)
        obatom3.SetFormalCharge(-to_add)
        given_charge += to_add
    return omol, given_charge


__all__ = [
    "eliminate_1_3_dipole",
    "eliminate_CN_in_doubt",
    "eliminate_NNN",
    "eliminate_carboxyl",
    "eliminate_carbene_neighbor_heteroatom",
    "eliminate_charge_spliting",
    "eliminate_high_positive_charge_atoms",
    "eliminate_negative_charges",
    "eliminate_positive_charges",
]
