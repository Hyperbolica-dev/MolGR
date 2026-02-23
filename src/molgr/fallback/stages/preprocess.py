from __future__ import annotations

import itertools
from typing import List, Optional, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel


def validate_omol(
    omol: pybel.Molecule, total_charge: int = 0, total_radical_electrons: int = 0
) -> bool:
    if sum(cast(ob.OBAtom, atom.OBAtom).GetFormalCharge() for atom in omol.atoms) != total_charge:
        return False
    radical_sum = sum(cast(ob.OBAtom, atom.OBAtom).GetSpinMultiplicity() for atom in omol.atoms)
    radical_sum_singlet = sum(
        cast(ob.OBAtom, atom.OBAtom).GetSpinMultiplicity() % 2 for atom in omol.atoms
    )
    if radical_sum_singlet == total_radical_electrons:
        radical_sum = radical_sum_singlet
    return radical_sum == total_radical_electrons


def make_connections(omol: pybel.Molecule, factor: float = 1.4) -> pybel.Molecule:
    obmol = cast(ob.OBMol, omol.OBMol)
    donate_smarts = pybel.Smarts("[Nv0,Cv1,Nv3,Clv1,Clv2,Clv3,Brv1,Brv2,Brv3,Iv1,Iv2,Iv3]")
    accept_smarts = pybel.Smarts(
        "[Hv0,Bv2,Bv3,Cv0,Cv1,Cv2,Cv3,Nv1,Nv2,Ov0,Ov1,Clv0,Siv3,Pv2,Sv0,Sv1,Brv0,Iv0]"
    )
    donate_atoms: List[int] = list(itertools.chain(*donate_smarts.findall(omol)))
    accept_atoms: List[int] = list(itertools.chain(*accept_smarts.findall(omol)))
    while donate_atoms and accept_atoms:
        donate_atom_id = donate_atoms.pop(0)
        pairs = sorted(
            [
                (donate_atom_id, accept_atom_id)
                for accept_atom_id in accept_atoms
                if accept_atom_id != donate_atom_id
            ],
            key=lambda x: cast(ob.OBAtom, obmol.GetAtom(x[0])).GetDistance(obmol.GetAtom(x[1])),
        )
        if not pairs:
            continue
        for pair_1, pair_2 in pairs:
            donate_atom = cast(ob.OBAtom, obmol.GetAtom(pair_1))
            accept_atom = cast(ob.OBAtom, obmol.GetAtom(pair_2))
            distance = cast(float, donate_atom.GetDistance(accept_atom))
            if (
                distance
                < cast(
                    float,
                    ob.GetCovalentRad(donate_atom.GetAtomicNum())
                    + ob.GetCovalentRad(accept_atom.GetAtomicNum()),
                )
                * factor
                and pair_1 in donate_atoms
                and pair_2 in accept_atoms
            ):
                if obmol.GetBond(pair_1, pair_2) is None:
                    obmol.AddBond(pair_1, pair_2, 1)
                    continue
                if cast(ob.OBBond, obmol.GetBond(pair_1, pair_2)).GetBondOrder() == 0:
                    cast(ob.OBBond, obmol.GetBond(pair_1, pair_2)).SetBondOrder(1)
                    donate_atoms = list(itertools.chain(*donate_smarts.findall(omol)))
                    accept_atoms = list(itertools.chain(*accept_smarts.findall(omol)))
    return omol


def pre_clean(omol: pybel.Molecule) -> pybel.Molecule:
    obmol = cast(ob.OBMol, omol.OBMol)

    smarts = pybel.Smarts("[Cv5,Nv5,Pv5,Siv5]=,#[*]")
    while res := smarts.findall(omol):
        idxs_1 = cast(Tuple[int, int], res.pop(0))
        obbond = cast(ob.OBBond, obmol.GetBond(idxs_1[0], idxs_1[1]))
        obbond.SetBondOrder(obbond.GetBondOrder() - 1)

    smarts = pybel.Smarts("[#6]1([#6]2)([#6]3)[#7]23[#6]1")
    while res := smarts.findall(omol):
        idxs_2 = cast(Tuple[int, int, int, int, int], res.pop(0))
        bcp_n: Optional[int] = None
        bcp_c: Optional[int] = None
        for idx in idxs_2:
            indexs = set(idxs_2) - {idx}
            if all(cast(ob.OBBond, obmol.GetBond(idx, idx_2)) for idx_2 in indexs):
                if cast(ob.OBAtom, obmol.GetAtom(idx)).GetAtomicNum() == 7:
                    bcp_n = idx
                if cast(ob.OBAtom, obmol.GetAtom(idx)).GetAtomicNum() == 6:
                    bcp_c = idx
        if bcp_n is not None and bcp_c is not None:
            obmol.DeleteBond(cast(ob.OBBond, obmol.GetBond(bcp_n, bcp_c)))

    smarts = pybel.Smarts("[#6]1([#6]2)[#7]2[#6]1")
    while res := smarts.findall(omol):
        idxs_3 = cast(Tuple[int, int, int, int], res.pop(0))
        amine_n: Optional[int] = None
        butyl_c: Optional[int] = None
        for idx in idxs_3:
            indexs = set(idxs_3) - {idx}
            if all(obmol.GetBond(idx, idx_2) for idx_2 in indexs):
                if cast(ob.OBAtom, obmol.GetAtom(idx)).GetAtomicNum() == 7:
                    amine_n = idx
                if cast(ob.OBAtom, obmol.GetAtom(idx)).GetAtomicNum() == 6:
                    butyl_c = idx
        if amine_n is not None and butyl_c is not None:
            obmol.DeleteBond(obmol.GetBond(amine_n, butyl_c))

    smarts = pybel.Smarts("[Siv5]-[O,F]")
    while res := smarts.findall(omol):
        idxs = cast(Tuple[int, int], res.pop(0))
        obmol.DeleteBond(obmol.GetBond(idxs[0], idxs[1]))

    return omol


__all__ = [
    "make_connections",
    "pre_clean",
    "validate_omol",
]
