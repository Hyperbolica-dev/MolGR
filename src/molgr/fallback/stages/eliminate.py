"""Charge-elimination heuristics shared by the linear and resonance cleanup paths."""

from __future__ import annotations

from typing import List, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.utils import consts, smarts


def eliminate_high_positive_charge_atoms(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    """Neutralize unstable highly positive atoms by borrowing electrons from neighbors."""

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False
    while res := smarts.ELIM_HIGH_POSITIVE.findall(omol):
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
        hit = True
    return omol, given_charge, hit


def eliminate_CN_in_doubt(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    """Resolve ambiguous C/N charge assignments in paired motifs."""

    obmol = cast(ob.OBMol, omol.OBMol)
    doubt_pair: List[Tuple[int, int]] = smarts.ELIM_CN_IN_DOUBT.findall(omol)
    cn_in_doubt = len(doubt_pair)
    # confirm that all atoms in doubt_pair are unique
    if len({atom_id for pair in doubt_pair for atom_id in pair}) != cn_in_doubt * 2:
        return omol, given_charge, False
    hit = False
    if cn_in_doubt % 2 == 0 and cn_in_doubt > 0:
        for atom_1_idx, atom_2_idx in doubt_pair[: cn_in_doubt // 2]:
            atom_1 = cast(ob.OBAtom, obmol.GetAtom(atom_1_idx))
            atom_2 = cast(ob.OBAtom, obmol.GetAtom(atom_2_idx))
            bond = cast(ob.OBBond, obmol.GetBond(atom_1_idx, atom_2_idx))
            atom_1.SetFormalCharge(-1)
            bond.SetBondOrder(bond.GetBondOrder() - 1)
            atom_2.SetFormalCharge(0)
            given_charge += 2
            hit = True
    return omol, given_charge, hit


def eliminate_carboxyl(omol: pybel.Molecule, given_charge: int) -> tuple[pybel.Molecule, int, bool]:
    """Collapse carboxyl-like radical patterns into their charged form."""

    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int]] = list(smarts.ELIM_CARBOXYL.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        atom_1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        atom_1.SetSpinMultiplicity(atom_1.GetSpinMultiplicity() - 1)
        atom_1.SetFormalCharge(atom_1.GetFormalCharge() - 1)
        given_charge += 1
        hit = True
    return omol, given_charge, hit


def eliminate_carbene_neighbor_heteroatom(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    """Push carbene radical density onto a neighboring heteroatom when possible."""

    hit = False
    for atom in omol.atoms:
        obatom = cast(ob.OBAtom, atom.OBAtom)
        if obatom.GetSpinMultiplicity() == 2:
            for neighbor in ob.OBAtomAtomIter(obatom):
                if cast(ob.OBAtom, neighbor).GetSpinMultiplicity():
                    return omol, given_charge, hit
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
                    hit = True
                    break
    return omol, given_charge, hit


def eliminate_NNN(
    omol: pybel.Molecule, given_charge: int, positive: bool = False
) -> tuple[pybel.Molecule, int, bool]:
    """Resolve the two N-N-N motifs that are common charge/radical ambiguities."""

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False
    if not positive:
        while res := smarts.ELIM_NNN_NEGATIVE.findall(omol):
            idxs = cast(List[Tuple[int, int, int]], res.pop(0))
            atom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
            atom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
            atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
            bond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
            bond1.SetBondOrder(bond1.GetBondOrder() + 1)
            bond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
            bond2.SetBondOrder(bond2.GetBondOrder() + 1)
            atom1.SetSpinMultiplicity(atom1.GetSpinMultiplicity() - 2)
            atom1.SetFormalCharge(atom1.GetFormalCharge() - 1)
            atom2.SetSpinMultiplicity(atom2.GetSpinMultiplicity() - 1)
            atom2.SetFormalCharge(atom2.GetFormalCharge() + 1)
            atom3.SetSpinMultiplicity(atom3.GetSpinMultiplicity() - 2)
            atom3.SetFormalCharge(atom3.GetFormalCharge() - 1)
            given_charge += 1
            hit = True
    else:
        while res := smarts.ELIM_NNN_POSITIVE.findall(omol):
            idxs = cast(List[Tuple[int, int, int]], res.pop(0))
            atom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
            atom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
            atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
            bond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
            bond1.SetBondOrder(bond1.GetBondOrder() + 1)
            atom1.SetFormalCharge(atom1.GetFormalCharge() + 1)
            atom2.SetSpinMultiplicity(atom2.GetSpinMultiplicity() - 1)
            given_charge -= 1
            hit = True
    return omol, given_charge, hit


def eliminate_charge_spliting(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    """Reduce overly split radical charge patterns before resonance expansion."""

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False

    if (
        all(cast(ob.OBAtom, atom).GetFormalCharge() == 0 for atom in ob.OBMolAtomIter(obmol))
        and sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(obmol))
        >= 2
    ):
        radical_atoms: List[ob.OBAtom] = [
            atom for atom in ob.OBMolAtomIter(obmol) if cast(ob.OBAtom, atom).GetSpinMultiplicity()
        ]
        total_radicals = sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in radical_atoms)
        while total_radicals > abs(given_charge):
            for atom in radical_atoms:
                if atom.GetAtomicNum() in (8, 9, 17, 35, 53):
                    atom.SetFormalCharge(atom.GetFormalCharge() - atom.GetSpinMultiplicity())
                    given_charge += atom.GetSpinMultiplicity()
                    total_radicals -= atom.GetSpinMultiplicity()
                    atom.SetSpinMultiplicity(0)
                    radical_atoms.remove(atom)
                    hit = True
                    break
            else:
                break
        while total_radicals > abs(given_charge):
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 16:
                    atom.SetFormalCharge(atom.GetFormalCharge() - atom.GetSpinMultiplicity())
                    given_charge += atom.GetSpinMultiplicity()
                    total_radicals -= atom.GetSpinMultiplicity()
                    atom.SetSpinMultiplicity(0)
                    radical_atoms.remove(atom)
                    hit = True
                    break
            else:
                break
        while total_radicals > abs(given_charge):
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 7:
                    atom.SetFormalCharge(atom.GetFormalCharge() - atom.GetSpinMultiplicity())
                    given_charge += atom.GetSpinMultiplicity()
                    total_radicals -= atom.GetSpinMultiplicity()
                    atom.SetSpinMultiplicity(0)
                    radical_atoms.remove(atom)
                    hit = True
                    break
            else:
                break
        while total_radicals > abs(given_charge):
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 6 and not any(
                    _atom
                    for _atom in ob.OBAtomAtomIter(atom)
                    if cast(ob.OBAtom, _atom).GetAtomicNum() in consts.HETEROATOM
                ):
                    atom.SetFormalCharge(atom.GetFormalCharge() - atom.GetSpinMultiplicity())
                    given_charge += atom.GetSpinMultiplicity()
                    total_radicals -= atom.GetSpinMultiplicity()
                    atom.SetSpinMultiplicity(0)
                    radical_atoms.remove(atom)
                    hit = True
                    break
            else:
                break
        while total_radicals > abs(given_charge):
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 6:
                    atom.SetFormalCharge(atom.GetFormalCharge() - atom.GetSpinMultiplicity())
                    given_charge += atom.GetSpinMultiplicity()
                    total_radicals -= atom.GetSpinMultiplicity()
                    atom.SetSpinMultiplicity(0)
                    radical_atoms.remove(atom)
                    hit = True
                    break
            else:
                break
    return omol, given_charge, hit


def eliminate_1_3_dipole(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    """Collapse simple 1,3-dipole motifs during resonance post-processing."""

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False

    res: List[Tuple[int, int, int]] = list(smarts.ELIM_1_3_DIPOLE.findall(omol))
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
            hit = True
    return omol, given_charge, hit


def eliminate_positive_charges(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False

    while given_charge > 0 and (res := smarts.ELIM_POSITIVE_N.findall(omol)):
        idxs_1 = cast(Tuple[int, int], res.pop(0))
        abatom = cast(ob.OBAtom, obmol.GetAtom(idxs_1[1]))
        abatom.SetSpinMultiplicity(abatom.GetSpinMultiplicity() - 1)
        abatom.SetFormalCharge(1)
        given_charge -= 1
        hit = True

    while given_charge > 0 and (res := smarts.ELIM_POSITIVE_C_H.findall(omol)):
        idxs_2 = cast(Tuple[int, int, int], res.pop(0))
        abatom2 = cast(ob.OBAtom, obmol.GetAtom(idxs_2[0]))
        abatom2.SetSpinMultiplicity(abatom2.GetSpinMultiplicity() - 1)
        abatom2.SetFormalCharge(1)
        given_charge -= 1
        hit = True
    for atom in omol.atoms:
        obatom = cast(ob.OBAtom, atom.OBAtom)
        if given_charge <= 0:
            break
        if obatom.GetSpinMultiplicity() >= 1 and obatom.GetFormalCharge() == 0:
            to_add = min(obatom.GetSpinMultiplicity(), given_charge)
            obatom.SetSpinMultiplicity(obatom.GetSpinMultiplicity() - to_add)
            obatom.SetFormalCharge(to_add)
            given_charge -= to_add
            if to_add > 0:
                hit = True
    return omol, given_charge, hit


def eliminate_negative_charges(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False
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
        if to_add > 0:
            hit = True

    while given_charge < 0 and (res := smarts.ELIM_NEGATIVE_C_V3.findall(omol)):
        idxs = cast(Tuple[int], res.pop(0))
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        to_add = min(obatom1.GetSpinMultiplicity(), abs(given_charge))
        obatom1.SetSpinMultiplicity(obatom1.GetSpinMultiplicity() - to_add)
        obatom1.SetFormalCharge(-to_add)
        given_charge += to_add
        if to_add > 0:
            hit = True

    while given_charge < 0 and (res := smarts.ELIM_NEGATIVE_H.findall(omol)):
        idxs = cast(Tuple[int], res.pop(0))
        obatom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        to_add = min(obatom2.GetSpinMultiplicity(), abs(given_charge))
        obatom2.SetSpinMultiplicity(obatom2.GetSpinMultiplicity() - to_add)
        obatom2.SetFormalCharge(-to_add)
        given_charge += to_add
        if to_add > 0:
            hit = True

    while given_charge < 0 and (res := smarts.ELIM_NEGATIVE_C_LOW.findall(omol)):
        idxs = cast(Tuple[int], res.pop(0))
        obatom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        to_add = min(obatom3.GetSpinMultiplicity(), abs(given_charge))
        obatom3.SetSpinMultiplicity(obatom3.GetSpinMultiplicity() - to_add)
        obatom3.SetFormalCharge(-to_add)
        given_charge += to_add
        if to_add > 0:
            hit = True

    while given_charge < 0:
        for atom in ob.OBMolAtomIter(obmol):
            atom = cast(ob.OBAtom, atom)
            if atom.GetSpinMultiplicity() >= 1 and atom.GetFormalCharge() == 0:
                to_add = min(atom.GetSpinMultiplicity(), abs(given_charge))
                atom.SetSpinMultiplicity(atom.GetSpinMultiplicity() - to_add)
                atom.SetFormalCharge(-to_add)
                given_charge += to_add
                if to_add > 0:
                    hit = True
        else:
            break
    return omol, given_charge, hit


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
