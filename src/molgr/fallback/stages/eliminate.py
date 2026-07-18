"""Charge-elimination heuristics shared by the linear and resonance cleanup paths."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.utils import consts, smarts


@dataclass(frozen=True)
class _ChargeAssignmentAction:
    atom_idx: int
    formal_charge: int
    spin_consumed: int
    charge_delta: int
    score_key: Tuple[int, ...]


@dataclass(frozen=True)
class _NegativeChargeAssignmentPattern:
    smarts: pybel.Smarts
    tier: int
    target_idx: int = 0
    requires_negative_deficit: bool = False


_NEGATIVE_CHARGE_ASSIGNMENT_PATTERNS = (
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_F, tier=10),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_O, tier=20),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_O_1, tier=21),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_CL, tier=30),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_N, tier=40),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_N_1, tier=41),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_N_2, tier=42),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_BR, tier=50),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_I, tier=60),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_S, tier=70),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_S_1, tier=71),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_SE, tier=80),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_SE_1, tier=81),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_P, tier=90),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_P_1, tier=91),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_P_2, tier=92),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_B, tier=95),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_B_1, tier=96),
    _NegativeChargeAssignmentPattern(smarts.ELIM_NEGATIVE_B_2, tier=97),
    _NegativeChargeAssignmentPattern(
        smarts.ELIM_NEGATIVE_C_V3,
        tier=100,
        requires_negative_deficit=True,
    ),
    _NegativeChargeAssignmentPattern(
        smarts.ELIM_NEGATIVE_C_LOW,
        tier=110,
        requires_negative_deficit=True,
    ),
    _NegativeChargeAssignmentPattern(
        smarts.ELIM_NEGATIVE_H,
        tier=120,
        requires_negative_deficit=True,
    ),
)


def _apply_charge_assignment_action(
    obmol: ob.OBMol,
    action: _ChargeAssignmentAction,
) -> bool:
    atom = cast(ob.OBAtom, obmol.GetAtom(action.atom_idx))
    if atom is None or action.spin_consumed <= 0:
        return False
    if atom.GetSpinMultiplicity() < action.spin_consumed:
        return False
    atom.SetSpinMultiplicity(atom.GetSpinMultiplicity() - action.spin_consumed)
    atom.SetFormalCharge(action.formal_charge)
    return True


def _atom_idx(atom: ob.OBAtom) -> int:
    return int(atom.GetIdx())


def _negative_charge_assignment_amount(atom: ob.OBAtom, given_charge: int) -> int:
    return min(atom.GetSpinMultiplicity(), max(1, abs(given_charge)))


def _positive_charge_assignment_actions(
    omol: pybel.Molecule,
    obmol: ob.OBMol,
    given_charge: int,
) -> list[_ChargeAssignmentAction]:
    if given_charge <= 0:
        return []

    actions: list[_ChargeAssignmentAction] = []
    seen: set[tuple[int, int, int]] = set()

    def append_action(atom: ob.OBAtom, *, tier: int, match_order: int, amount: int) -> None:
        if atom is None or amount <= 0:
            return
        atom_idx = _atom_idx(atom)
        key = (atom_idx, tier, amount)
        if key in seen:
            return
        seen.add(key)
        charge_after = given_charge - amount
        atomic_num = int(atom.GetAtomicNum())
        actions.append(
            _ChargeAssignmentAction(
                atom_idx=atom_idx,
                formal_charge=amount,
                spin_consumed=amount,
                charge_delta=-amount,
                score_key=(
                    tier,
                    abs(charge_after),
                    max(charge_after, 0),
                    atomic_num,
                    atom_idx,
                    match_order,
                ),
            )
        )

    for match_order, n_idxs in enumerate(
        cast(List[Tuple[int, int]], smarts.ELIM_POSITIVE_N.findall(omol))
    ):
        atom = cast(ob.OBAtom, obmol.GetAtom(n_idxs[1]))
        if atom.GetFormalCharge() == 0 and atom.GetSpinMultiplicity() >= 1:
            append_action(atom, tier=0, match_order=match_order, amount=1)

    for match_order, c_h_idxs in enumerate(
        cast(List[Tuple[int, int, int]], smarts.ELIM_POSITIVE_C_H.findall(omol))
    ):
        atom = cast(ob.OBAtom, obmol.GetAtom(c_h_idxs[0]))
        if atom.GetFormalCharge() == 0 and atom.GetSpinMultiplicity() >= 1:
            append_action(atom, tier=10, match_order=match_order, amount=1)

    for match_order, atom_iter in enumerate(ob.OBMolAtomIter(obmol)):
        atom = cast(ob.OBAtom, atom_iter)
        if atom.GetFormalCharge() == 0 and atom.GetSpinMultiplicity() >= 1:
            append_action(
                atom,
                tier=100,
                match_order=match_order,
                amount=min(atom.GetSpinMultiplicity(), given_charge),
            )

    return sorted(actions, key=lambda action: action.score_key)


def _negative_charge_assignment_actions(
    omol: pybel.Molecule,
    obmol: ob.OBMol,
    given_charge: int,
) -> list[_ChargeAssignmentAction]:
    if given_charge > 0:
        return []

    actions: list[_ChargeAssignmentAction] = []
    seen: set[tuple[int, int, int]] = set()

    def append_action(atom: ob.OBAtom, *, tier: int, match_order: int, amount: int) -> None:
        if atom is None or amount <= 0:
            return
        atom_idx = _atom_idx(atom)
        key = (atom_idx, tier, amount)
        if key in seen:
            return
        seen.add(key)
        charge_after = given_charge + amount
        atomic_num = int(atom.GetAtomicNum())
        actions.append(
            _ChargeAssignmentAction(
                atom_idx=atom_idx,
                formal_charge=-amount,
                spin_consumed=amount,
                charge_delta=amount,
                score_key=(
                    tier,
                    abs(charge_after),
                    max(charge_after, 0),
                    atomic_num,
                    atom_idx,
                    match_order,
                ),
            )
        )

    for pattern in _NEGATIVE_CHARGE_ASSIGNMENT_PATTERNS:
        if pattern.requires_negative_deficit and given_charge >= 0:
            continue
        for match_order, pattern_idxs in enumerate(
            cast(List[Tuple[int, ...]], pattern.smarts.findall(omol))
        ):
            atom = cast(ob.OBAtom, obmol.GetAtom(pattern_idxs[pattern.target_idx]))
            if atom.GetFormalCharge() == 0 and atom.GetSpinMultiplicity() >= 1:
                append_action(
                    atom,
                    tier=pattern.tier,
                    match_order=match_order,
                    amount=_negative_charge_assignment_amount(atom, given_charge),
                )

    if given_charge < 0:
        for match_order, atom_iter in enumerate(ob.OBMolAtomIter(obmol)):
            atom = cast(ob.OBAtom, atom_iter)
            if atom.GetFormalCharge() == 0 and atom.GetSpinMultiplicity() >= 1:
                append_action(
                    atom,
                    tier=1000,
                    match_order=match_order,
                    amount=_negative_charge_assignment_amount(atom, given_charge),
                )

    return sorted(actions, key=lambda action: action.score_key)


def eliminate_high_positive_charge_atoms(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    """Neutralize unstable highly positive atoms by borrowing electrons from neighbors."""

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False
    res = smarts.ELIM_HIGH_POSITIVE.findall(omol)
    while len(res):
        idxs = cast(List[Tuple[int, int]], res.pop(0))
        atom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        atom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
        if (
            -sum(cast(ob.OBAtom, atom).GetFormalCharge() for atom in ob.OBAtomAtomIter(atom1))
            > atom1.GetFormalCharge()
            or atom2.GetSpinMultiplicity() != 1
        ):
            continue
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

    def possible_carbene_atom_checker(obatom: ob.OBAtom) -> bool:
        if obatom.GetAtomicNum() not in consts.HETEROATOM and obatom.GetAtomicNum() != 6:
            return False
        if obatom.GetSpinMultiplicity() != 2:
            return False
        return obatom.GetAtomicNum() not in (8, 9, 16, 17, 35, 53)

    hit = False
    for atom in omol.atoms:
        obatom = cast(ob.OBAtom, atom.OBAtom)
        if possible_carbene_atom_checker(obatom):
            if any(
                cast(ob.OBAtom, neighbor).GetSpinMultiplicity() == 1
                for neighbor in ob.OBAtomAtomIter(obatom)
            ):
                continue
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


def assign_negative_charges_from_radicals(
    omol: pybel.Molecule, remaining_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    """Convert selected single radicals into anions and update the charge budget."""

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False

    if (
        all(cast(ob.OBAtom, atom).GetFormalCharge() == 0 for atom in ob.OBMolAtomIter(obmol))
        and sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in ob.OBMolAtomIter(obmol))
        >= 2
    ):
        radical_atoms: List[ob.OBAtom] = [
            atom
            for atom in ob.OBMolAtomIter(obmol)
            if cast(ob.OBAtom, atom).GetSpinMultiplicity() == 1
        ]
        total_radicals = sum(cast(ob.OBAtom, atom).GetSpinMultiplicity() for atom in radical_atoms)
        while total_radicals > abs(remaining_charge):
            for atom in radical_atoms:
                if atom.GetAtomicNum() in (8, 9, 17, 35, 53):
                    atom.SetFormalCharge(atom.GetFormalCharge() - atom.GetSpinMultiplicity())
                    remaining_charge += atom.GetSpinMultiplicity()
                    total_radicals -= atom.GetSpinMultiplicity()
                    atom.SetSpinMultiplicity(0)
                    radical_atoms.remove(atom)
                    hit = True
                    break
            else:
                break
        while total_radicals > abs(remaining_charge):
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 16:
                    atom.SetFormalCharge(atom.GetFormalCharge() - atom.GetSpinMultiplicity())
                    remaining_charge += atom.GetSpinMultiplicity()
                    total_radicals -= atom.GetSpinMultiplicity()
                    atom.SetSpinMultiplicity(0)
                    radical_atoms.remove(atom)
                    hit = True
                    break
            else:
                break
        while total_radicals > abs(remaining_charge):
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 7:
                    atom.SetFormalCharge(atom.GetFormalCharge() - atom.GetSpinMultiplicity())
                    remaining_charge += atom.GetSpinMultiplicity()
                    total_radicals -= atom.GetSpinMultiplicity()
                    atom.SetSpinMultiplicity(0)
                    radical_atoms.remove(atom)
                    hit = True
                    break
            else:
                break
        while total_radicals > abs(remaining_charge):
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 6 and not any(
                    _atom
                    for _atom in ob.OBAtomAtomIter(atom)
                    if cast(ob.OBAtom, _atom).GetAtomicNum() in consts.HETEROATOM
                ):
                    atom.SetFormalCharge(atom.GetFormalCharge() - atom.GetSpinMultiplicity())
                    remaining_charge += atom.GetSpinMultiplicity()
                    total_radicals -= atom.GetSpinMultiplicity()
                    atom.SetSpinMultiplicity(0)
                    radical_atoms.remove(atom)
                    hit = True
                    break
            else:
                break
        while total_radicals > abs(remaining_charge):
            for atom in radical_atoms:
                if atom.GetAtomicNum() == 6:
                    atom.SetFormalCharge(atom.GetFormalCharge() - atom.GetSpinMultiplicity())
                    remaining_charge += atom.GetSpinMultiplicity()
                    total_radicals -= atom.GetSpinMultiplicity()
                    atom.SetSpinMultiplicity(0)
                    radical_atoms.remove(atom)
                    hit = True
                    break
            else:
                break
    return omol, remaining_charge, hit


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


def eliminate_cp_like_radical_anion(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    """Convert cyclopentadienyl-like radicals into their anionic aromatic form."""

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False
    while given_charge < 0 and (res := smarts.ELIM_NEGATIVE_CP.findall(omol)):
        idxs = cast(Tuple[int, int, int, int, int], res.pop(0))
        atom = cast(ob.OBAtom, obmol.GetAtom(idxs[4]))
        if atom.GetFormalCharge() != 0:
            break
        to_add = atom.GetSpinMultiplicity()
        if to_add <= 0:
            break
        atom.SetSpinMultiplicity(atom.GetSpinMultiplicity() - to_add)
        atom.SetFormalCharge(-to_add)
        given_charge += to_add
        hit = True
    return omol, given_charge, hit


def eliminate_positive_charges(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False

    while given_charge > 0:
        actions = _positive_charge_assignment_actions(omol, obmol, given_charge)
        if not actions:
            break
        action = actions[0]
        if not _apply_charge_assignment_action(obmol, action):
            break
        given_charge += action.charge_delta
        hit = True
    return omol, given_charge, hit


def eliminate_negative_charges(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False

    while given_charge <= 0:
        actions = _negative_charge_assignment_actions(omol, obmol, given_charge)
        if not actions:
            break
        action = actions[0]
        if not _apply_charge_assignment_action(obmol, action):
            break
        given_charge += action.charge_delta
        hit = True

    return omol, given_charge, hit


__all__ = [
    "eliminate_1_3_dipole",
    "eliminate_CN_in_doubt",
    "eliminate_NNN",
    "eliminate_carboxyl",
    "eliminate_carbene_neighbor_heteroatom",
    "assign_negative_charges_from_radicals",
    "eliminate_cp_like_radical_anion",
    "eliminate_high_positive_charge_atoms",
    "eliminate_negative_charges",
    "eliminate_positive_charges",
]
