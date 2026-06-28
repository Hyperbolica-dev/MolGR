"""Early structure preparation and validation stages for fallback."""

from __future__ import annotations

import itertools
from typing import List, Optional, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.utils import smarts


def validate_omol(
    omol: pybel.Molecule, total_charge: int = 0, total_radical_electrons: int = 0
) -> bool:
    """Check whether the current molecule matches the requested charge/radical budget."""

    if sum(cast(ob.OBAtom, atom.OBAtom).GetFormalCharge() for atom in omol.atoms) != total_charge:
        return False
    radical_sum = sum(cast(ob.OBAtom, atom.OBAtom).GetSpinMultiplicity() for atom in omol.atoms)
    radical_sum_singlet = sum(
        cast(ob.OBAtom, atom.OBAtom).GetSpinMultiplicity() % 2 for atom in omol.atoms
    )
    if radical_sum_singlet == total_radical_electrons:
        radical_sum = radical_sum_singlet
    if any(
        cast(ob.OBBond, bond).GetBondOrder() > 5 or cast(ob.OBBond, bond).GetBondOrder() < 0
        for bond in ob.OBMolBondIter(cast(ob.OBMol, omol.OBMol))
    ):
        return False
    return radical_sum == total_radical_electrons


def make_connections(
    omol: pybel.Molecule,
    extra_tolerance_angstrom: float = 0.15,
) -> tuple[pybel.Molecule, bool]:
    """Reconnect obvious donor/acceptor pairs before formal valence cleanup starts."""

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False
    while True:
        donate_atoms: List[int] = list(itertools.chain(*smarts.PREPROCESS_DONATE.findall(omol)))
        accept_atoms: List[int] = list(itertools.chain(*smarts.PREPROCESS_ACCEPT.findall(omol)))
        if not donate_atoms or not accept_atoms:
            break

        changed = False
        for donate_atom_id in donate_atoms:
            pairs = sorted(
                [
                    (donate_atom_id, accept_atom_id)
                    for accept_atom_id in accept_atoms
                    if accept_atom_id != donate_atom_id
                ],
                key=lambda x: cast(ob.OBAtom, obmol.GetAtom(x[0])).GetDistance(obmol.GetAtom(x[1])),
            )
            for pair_1, pair_2 in pairs:
                donate_atom = cast(ob.OBAtom, obmol.GetAtom(pair_1))
                accept_atom = cast(ob.OBAtom, obmol.GetAtom(pair_2))
                distance = cast(float, donate_atom.GetDistance(accept_atom))
                if distance >= cast(
                    float,
                    ob.GetCovalentRad(donate_atom.GetAtomicNum())
                    + ob.GetCovalentRad(accept_atom.GetAtomicNum()),
                ) + extra_tolerance_angstrom:
                    continue

                bond = cast(Optional[ob.OBBond], obmol.GetBond(pair_1, pair_2))
                if bond is None:
                    obmol.AddBond(pair_1, pair_2, 1)
                    hit = True
                    changed = True
                    break
                if bond.GetBondOrder() == 0:
                    bond.SetBondOrder(1)
                    hit = True
                    changed = True
                    break
            if changed:
                break
        if not changed:
            break
    return omol, hit


def pre_clean(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    """Apply cheap structural fixes that simplify the later heuristic stages."""

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False

    while res := smarts.PRE_CLEAN_HYPERVALENT.findall(omol):
        idxs_1 = cast(Tuple[int, int], res.pop(0))
        obbond = cast(ob.OBBond, obmol.GetBond(idxs_1[0], idxs_1[1]))
        obbond.SetBondOrder(obbond.GetBondOrder() - 1)
        hit = True

    while res := smarts.PRE_CLEAN_BCP_RING_5.findall(omol):
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
            hit = True

    while res := smarts.PRE_CLEAN_BCP_RING_4.findall(omol):
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
            hit = True

    while res := smarts.PRE_CLEAN_SI_O_F.findall(omol):
        idxs = cast(Tuple[int, int], res.pop(0))
        obmol.DeleteBond(obmol.GetBond(idxs[0], idxs[1]))
        hit = True

    return omol, hit


__all__ = ["make_connections", "pre_clean", "validate_omol"]
