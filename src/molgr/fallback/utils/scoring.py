from __future__ import annotations

import math
from collections.abc import Sequence
from typing import List, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.utils import consts


def calculate_tetrahedron_volume(
    p1: Sequence[float], p2: Sequence[float], p3: Sequence[float], p4: Sequence[float]
) -> float:
    a = [p1[0] - p4[0], p1[1] - p4[1], p1[2] - p4[2]]
    b = [p2[0] - p4[0], p2[1] - p4[1], p2[2] - p4[2]]
    c = [p3[0] - p4[0], p3[1] - p4[1], p3[2] - p4[2]]
    det = (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )
    return abs(det) / 6.0


def calculate_shape_quality(
    p1: Sequence[float], p2: Sequence[float], p3: Sequence[float], p4: Sequence[float]
) -> float:
    volume = calculate_tetrahedron_volume(p1, p2, p3, p4)
    if math.isclose(volume, 0.0):
        return 0.0

    def d2(a: Sequence[float], b: Sequence[float]) -> float:
        return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2

    edges_sq = [d2(p1, p2), d2(p1, p3), d2(p1, p4), d2(p2, p3), d2(p2, p4), d2(p3, p4)]
    sum_edges_sq = sum(edges_sq)
    l_rms_cubed = (sum_edges_sq / 6.0) ** 1.5
    if math.isclose(l_rms_cubed, 0.0):
        return 0.0
    quality = (6.0 * math.sqrt(2.0)) * (volume / l_rms_cubed)
    return max(0.0, min(1.0, quality))


def get_deviation_score(omol: pybel.Molecule, atom_idx: int) -> float:
    def _get_coords(a: ob.OBAtom) -> Sequence[float]:
        vec = cast(ob.vector3, a.GetVector())
        return [vec.GetX(), vec.GetY(), vec.GetZ()]

    obmol = cast(ob.OBMol, omol.OBMol)
    atom = cast(ob.OBAtom, obmol.GetAtom(atom_idx))
    neighbor_atoms: List[ob.OBAtom] = list(ob.OBAtomAtomIter(atom))
    if len(neighbor_atoms) == 2:
        angle = obmol.GetAngle(neighbor_atoms[0], atom, neighbor_atoms[1])
        return abs(angle - 108) / 108.0
    if len(neighbor_atoms) == 3:
        quality = calculate_shape_quality(
            _get_coords(neighbor_atoms[0]),
            _get_coords(neighbor_atoms[1]),
            _get_coords(neighbor_atoms[2]),
            _get_coords(atom),
        )
        return 1.0 - quality
    return 0.0


def calc_symmetry_penalty(omol: pybel.Molecule) -> float:
    obmol = cast(ob.OBMol, omol.OBMol)
    gs = ob.OBGraphSym(obmol)
    symmetry_ids_vec = ob.vectorUnsignedInt()
    gs.GetSymmetry(symmetry_ids_vec)
    return (len(set(symmetry_ids_vec)) - obmol.NumAtoms()) * 2.0


def calculate_charge_penalty(atom: ob.OBAtom) -> float:
    charge = cast(int, atom.GetFormalCharge())
    if charge == 0:
        return 0.0
    total_electrons = cast(
        int,
        consts.NON_METAL_DICT[atom.GetAtomicNum()].num_outer_electrons
        + atom.GetTotalValence()
        - atom.GetFormalCharge()
        + atom.GetSpinMultiplicity(),
    )
    if total_electrons == 8 or total_electrons == 2:
        return 0.0
    en = cast(float, ob.GetElectroNeg(atom.GetAtomicNum()))
    if charge > 0:
        return abs(charge) * max(0.0, en - 2) * 3.0
    return abs(charge) * max(0.0, 4 - en) * 3.0


def calculate_radical_penalty(atom: ob.OBAtom) -> float:
    radical_num = cast(int, atom.GetSpinMultiplicity())
    if radical_num == 0:
        return 0.0
    if cast(int, atom.GetAtomicNum()) in consts.HETEROATOM:
        return radical_num * 2.0
    return (3 - cast(int, atom.GetHvyDegree())) * 1.5


def calculate_coulombic_penalty(bond: ob.OBBond) -> float:
    q1 = cast(int, cast(ob.OBAtom, bond.GetBeginAtom()).GetFormalCharge())
    q2 = cast(int, cast(ob.OBAtom, bond.GetEndAtom()).GetFormalCharge())
    if q1 == 0 or q2 == 0:
        return 0.0
    if q1 * q2 > 0:
        return 15.0
    return -0.5


def calculate_physchem_penalty(omol: pybel.Molecule) -> float:
    total_penalty = 0.0
    obmol = cast(ob.OBMol, omol.OBMol)
    for atom in ob.OBMolAtomIter(obmol):
        if cast(ob.OBAtom, atom).IsMetal():
            continue
        total_penalty += calculate_charge_penalty(atom)
        total_penalty += calculate_radical_penalty(atom)
    for bond in ob.OBMolBondIter(obmol):
        total_penalty += calculate_coulombic_penalty(bond)
    return total_penalty


def get_metal_coordination_sphere(
    omol: pybel.Molecule, metal_atom: ob.OBAtom, cutoff: float = 2.8
) -> List[Tuple[ob.OBAtom, float]]:
    neighbors: List[Tuple[ob.OBAtom, float]] = []
    metal_idx = metal_atom.GetIdx()
    for atom in omol.atoms:
        obatom = cast(ob.OBAtom, atom.OBAtom)
        if obatom.GetIdx() == metal_idx:
            continue
        dist = metal_atom.GetDistance(obatom)
        if dist <= cutoff:
            neighbors.append((obatom, dist))
    return neighbors


def calculate_metal_penalty(omol: pybel.Molecule) -> float:
    penalty = 0.0
    for atom in omol.atoms:
        metal_atom = cast(ob.OBAtom, atom.OBAtom)
        if not metal_atom.IsMetal():
            continue
        symbol = cast(str, ob.GetSymbol(metal_atom.GetAtomicNum()))
        valence = cast(int, metal_atom.GetFormalCharge())
        if valence <= 0:
            penalty += 10
        prior_list = consts.METAL_VALENCE_AVAILABLE_PRIOR.get(symbol, [])
        minor_list = consts.METAL_VALENCE_AVAILABLE_MINOR.get(symbol, [])
        if valence not in prior_list:
            if valence in minor_list:
                penalty += 2.0
            else:
                penalty += 10.0
        neighbors = get_metal_coordination_sphere(omol, metal_atom, cutoff=2.6)
        for ligand_atom, dist in neighbors:
            ligand_charge = cast(int, ligand_atom.GetFormalCharge())
            if valence > 0:
                if ligand_charge > 0:
                    penalty += 10.0 * (ligand_charge * valence) / (dist**2)
                elif ligand_charge < 0:
                    penalty -= 2.0 * (abs(ligand_charge) * valence) / dist
    return penalty


def omol_score(omol: pybel.Molecule) -> float:
    obmol = cast(ob.OBMol, omol.OBMol)
    obmol.SetAromaticPerceived(False)
    score = calc_symmetry_penalty(omol)
    for atom_idx in range(1, obmol.NumAtoms() + 1):
        atom = cast(ob.OBAtom, obmol.GetAtom(atom_idx))
        if atom.IsMetal():
            continue
        if not atom.IsAromatic():
            score += 5
        if atom.GetSpinMultiplicity() > 0:
            score += get_deviation_score(omol, atom_idx) * 10.0
        if atom.GetFormalCharge() > 0:
            score += get_deviation_score(omol, atom_idx) * 10.0
        if atom.GetFormalCharge() < 0:
            score += (1 - get_deviation_score(omol, atom_idx)) * 10.0
    score += calculate_physchem_penalty(omol)
    score += calculate_metal_penalty(omol)
    return score


__all__ = [
    "omol_score",
]
