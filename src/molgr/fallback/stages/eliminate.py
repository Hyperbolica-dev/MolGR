"""Charge-elimination heuristics shared by the linear and resonance cleanup paths."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.utils import consts, smarts
from molgr.fallback.utils.electrons import (
    get_lone_pair_count,
    get_unpaired_electron_count,
    has_unresolved_two_electron_center,
    set_lone_pair_count,
    set_unpaired_electron_count,
    set_unresolved_two_electron_center,
)


@dataclass(frozen=True)
class _ChargeAssignmentAction:
    atom_idx: int
    formal_charge: int
    spin_consumed: int
    consume_unresolved_center: bool
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
    """Encode one explicitly identified electron source as formal charge.

    Radical actions consume the requested number of real unpaired electrons. An
    unresolved action instead consumes one pure deferred two-electron marker
    atomically and assigns ``+2`` or ``-2`` without first creating two radicals.
    Active lone pairs are never interchangeable inputs here.
    """

    atom = cast(ob.OBAtom, obmol.GetAtom(action.atom_idx))
    if atom is None or atom.GetFormalCharge() != 0:
        return False
    if action.consume_unresolved_center:
        if (
            action.spin_consumed != 0
            or abs(action.formal_charge) != 2
            or not has_unresolved_two_electron_center(atom)
            or get_unpaired_electron_count(atom) != 0
            or get_lone_pair_count(atom) != 0
        ):
            return False
        set_unresolved_two_electron_center(atom, False)
    else:
        if action.spin_consumed <= 0 or has_unresolved_two_electron_center(atom):
            return False
        if get_unpaired_electron_count(atom) < action.spin_consumed:
            return False
        set_unpaired_electron_count(atom, get_unpaired_electron_count(atom) - action.spin_consumed)
    atom.SetFormalCharge(action.formal_charge)
    return True


def _atom_idx(atom: ob.OBAtom) -> int:
    return int(atom.GetIdx())


def _negative_charge_assignment_amount(atom: ob.OBAtom, given_charge: int) -> int:
    """Choose how many real unpaired electrons a multivalent anion action consumes."""

    return min(get_unpaired_electron_count(atom), max(1, abs(given_charge)))


def _positive_charge_assignment_actions(
    omol: pybel.Molecule,
    obmol: ob.OBMol,
    given_charge: int,
) -> list[_ChargeAssignmentAction]:
    """Rank cation actions from radicals or a pure unresolved 2e center.

    Radical actions preserve their existing motif priority. When the remaining
    deficit is at least two, the generic tier also admits a pure unresolved center
    as an atomic ``+2`` action. Active lone pairs remain excluded.
    """

    if given_charge <= 0:
        return []

    actions: list[_ChargeAssignmentAction] = []
    seen: set[tuple[int, int, int]] = set()

    def append_action(atom: ob.OBAtom, *, tier: int, match_order: int, amount: int) -> None:
        if atom is None or amount <= 0 or has_unresolved_two_electron_center(atom):
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
                consume_unresolved_center=False,
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

    def append_unresolved_action(atom: ob.OBAtom, *, match_order: int) -> None:
        if (
            atom is None
            or given_charge < 2
            or atom.GetFormalCharge() != 0
            or not has_unresolved_two_electron_center(atom)
            or get_unpaired_electron_count(atom) != 0
            or get_lone_pair_count(atom) != 0
        ):
            return
        atom_idx = _atom_idx(atom)
        charge_after = given_charge - 2
        actions.append(
            _ChargeAssignmentAction(
                atom_idx=atom_idx,
                formal_charge=2,
                spin_consumed=0,
                consume_unresolved_center=True,
                charge_delta=-2,
                score_key=(
                    100,
                    abs(charge_after),
                    max(charge_after, 0),
                    int(atom.GetAtomicNum()),
                    atom_idx,
                    match_order,
                ),
            )
        )

    for match_order, n_idxs in enumerate(
        cast(List[Tuple[int, int]], smarts.ELIM_POSITIVE_N.findall(omol))
    ):
        atom = cast(ob.OBAtom, obmol.GetAtom(n_idxs[1]))
        if atom.GetFormalCharge() == 0 and get_unpaired_electron_count(atom) >= 1:
            append_action(atom, tier=0, match_order=match_order, amount=1)

    for match_order, c_h_idxs in enumerate(
        cast(List[Tuple[int, int, int]], smarts.ELIM_POSITIVE_C_H.findall(omol))
    ):
        atom = cast(ob.OBAtom, obmol.GetAtom(c_h_idxs[0]))
        if atom.GetFormalCharge() == 0 and get_unpaired_electron_count(atom) >= 1:
            append_action(atom, tier=10, match_order=match_order, amount=1)

    for match_order, atom_iter in enumerate(ob.OBMolAtomIter(obmol)):
        atom = cast(ob.OBAtom, atom_iter)
        if atom.GetFormalCharge() == 0 and get_unpaired_electron_count(atom) >= 1:
            append_action(
                atom,
                tier=100,
                match_order=match_order,
                amount=min(get_unpaired_electron_count(atom), given_charge),
            )
        append_unresolved_action(atom, match_order=match_order)

    return sorted(actions, key=lambda action: action.score_key)


def _negative_charge_assignment_actions(
    omol: pybel.Molecule,
    obmol: ob.OBMol,
    given_charge: int,
) -> list[_ChargeAssignmentAction]:
    """Rank anion actions from radicals or a pure unresolved 2e center.

    Radical actions preserve their existing SMARTS priority and may seed charge
    separation at zero budget. An unresolved center is admitted only for a real
    deficit of at least two and becomes ``-2``. Active lone pairs remain excluded.
    """

    if given_charge > 0:
        return []

    actions: list[_ChargeAssignmentAction] = []
    seen: set[tuple[int, int, int]] = set()

    def append_action(atom: ob.OBAtom, *, tier: int, match_order: int, amount: int) -> None:
        if atom is None or amount <= 0 or has_unresolved_two_electron_center(atom):
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
                consume_unresolved_center=False,
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

    def append_unresolved_action(atom: ob.OBAtom, *, match_order: int) -> None:
        if (
            atom is None
            or given_charge > -2
            or atom.GetFormalCharge() != 0
            or not has_unresolved_two_electron_center(atom)
            or get_unpaired_electron_count(atom) != 0
            or get_lone_pair_count(atom) != 0
        ):
            return
        atom_idx = _atom_idx(atom)
        charge_after = given_charge + 2
        actions.append(
            _ChargeAssignmentAction(
                atom_idx=atom_idx,
                formal_charge=-2,
                spin_consumed=0,
                consume_unresolved_center=True,
                charge_delta=2,
                score_key=(
                    1000,
                    abs(charge_after),
                    max(charge_after, 0),
                    int(atom.GetAtomicNum()),
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
            if atom.GetFormalCharge() == 0 and get_unpaired_electron_count(atom) >= 1:
                append_action(
                    atom,
                    tier=pattern.tier,
                    match_order=match_order,
                    amount=_negative_charge_assignment_amount(atom, given_charge),
                )

    if given_charge < 0:
        for match_order, atom_iter in enumerate(ob.OBMolAtomIter(obmol)):
            atom = cast(ob.OBAtom, atom_iter)
            if atom.GetFormalCharge() == 0 and get_unpaired_electron_count(atom) >= 1:
                append_action(
                    atom,
                    tier=1000,
                    match_order=match_order,
                    amount=_negative_charge_assignment_amount(atom, given_charge),
                )
            append_unresolved_action(atom, match_order=match_order)

    return sorted(actions, key=lambda action: action.score_key)


def eliminate_high_positive_charge_atoms(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    """Convert one neighboring monoradical into a one-electron anionic center.

    The motif only selects the site. The electronic operation consumes exactly
    one real unpaired electron on the neighbor and replaces it with
    ``formal_charge -= 1``; lone pairs and unresolved centers are not consumed.
    """

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
            or get_unpaired_electron_count(atom2) != 1
        ):
            continue
        set_unpaired_electron_count(atom2, get_unpaired_electron_count(atom2) - 1)
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
    """Convert a carboxyl-like oxygen monoradical into an ``O-`` center.

    This one-electron localization requires and consumes exactly one real
    unpaired electron, then records it through a one-unit formal charge.
    """

    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int]] = list(smarts.ELIM_CARBOXYL.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        atom_1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        if get_unpaired_electron_count(atom_1) != 1:
            continue
        set_unpaired_electron_count(atom_1, 0)
        atom_1.SetFormalCharge(atom_1.GetFormalCharge() - 1)
        given_charge += 1
        hit = True
    return omol, given_charge, hit


def eliminate_carbene_neighbor_heteroatom(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    """Close an unresolved two-electron center with a heteroatom donor.

    Raising the bond order consumes the center marker and encodes its electrons
    as ``center- / donor+``. If the donor carries an explicitly tracked active
    lone pair, one pair is consumed; a zero count means the conventional donor
    pair is represented implicitly by valence and formal charge. Real unpaired
    electrons cannot drive this closed-shell transformation.
    """

    def possible_carbene_atom_checker(obatom: ob.OBAtom) -> bool:
        """Accept only an unresolved center, never an inferred spin/lone-pair surrogate."""

        return has_unresolved_two_electron_center(obatom)

    hit = False
    for atom in omol.atoms:
        obatom = cast(ob.OBAtom, atom.OBAtom)
        if possible_carbene_atom_checker(obatom):
            if any(
                get_unpaired_electron_count(cast(ob.OBAtom, neighbor)) == 1
                for neighbor in ob.OBAtomAtomIter(obatom)
            ):
                continue
            for neighbor in ob.OBAtomAtomIter(obatom):
                if (
                    cast(ob.OBAtom, neighbor).GetAtomicNum() in consts.HETEROATOM
                    and cast(ob.OBAtom, neighbor).GetFormalCharge() == 0
                    and get_unpaired_electron_count(cast(ob.OBAtom, neighbor)) == 0
                    and not has_unresolved_two_electron_center(cast(ob.OBAtom, neighbor))
                ):
                    donor = cast(ob.OBAtom, neighbor)
                    bond = cast(ob.OBBond, obatom.GetBond(donor))
                    bond.SetBondOrder(bond.GetBondOrder() + 1)
                    set_lone_pair_count(obatom, 0)
                    set_unpaired_electron_count(obatom, 0)
                    set_unresolved_two_electron_center(obatom, False)
                    obatom.SetFormalCharge(obatom.GetFormalCharge() - 1)
                    if get_lone_pair_count(donor) > 0:
                        set_lone_pair_count(donor, get_lone_pair_count(donor) - 1)
                    donor.SetFormalCharge(donor.GetFormalCharge() + 1)
                    hit = True
                    break
    return omol, given_charge, hit


def eliminate_NNN(
    omol: pybel.Molecule, given_charge: int, positive: bool = False
) -> tuple[pybel.Molecule, int, bool]:
    """Close N-N-N motifs using explicitly classified endpoint electrons.

    In the negative branch each terminal N must supply one classified
    two-electron state: an unresolved marker, two real unpaired electrons, or
    one active lone pair. The middle N must supply exactly one real unpaired
    electron. The positive branch also requires one middle-N monoradical.
    Missing or mixed labels reject the rule before any state is changed.
    """

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False

    def has_consumable_two_electron_center(atom: ob.OBAtom) -> bool:
        """Recognize one pure unresolved, triplet, or active-singlet two-electron state.

        Mixed marker/occupancy states and a singlet accompanied by any real
        unpaired electrons are rejected before the N-N-N graph is changed.
        """

        unresolved = has_unresolved_two_electron_center(atom)
        unpaired = get_unpaired_electron_count(atom)
        lone_pairs = get_lone_pair_count(atom)
        return (
            (unresolved and unpaired == 0 and lone_pairs == 0)
            or (not unresolved and unpaired == 2 and lone_pairs == 0)
            or (not unresolved and unpaired == 0 and lone_pairs >= 1)
        )

    def consume_two_electron_center(atom: ob.OBAtom) -> None:
        """Remove exactly one prevalidated two-electron state from its own field.

        No parity conversion is performed: unresolved clears the marker, triplet
        clears two unpaired electrons, and singlet subtracts one active lone pair.
        """

        if has_unresolved_two_electron_center(atom):
            set_unresolved_two_electron_center(atom, False)
        elif get_unpaired_electron_count(atom) == 2:
            set_unpaired_electron_count(atom, 0)
        else:
            set_lone_pair_count(atom, get_lone_pair_count(atom) - 1)

    if not positive:
        while res := smarts.ELIM_NNN_NEGATIVE.findall(omol):
            idxs = cast(List[Tuple[int, int, int]], res.pop(0))
            atom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
            atom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
            atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
            bond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
            if (
                not has_consumable_two_electron_center(atom1)
                or get_unpaired_electron_count(atom2) != 1
                or not has_consumable_two_electron_center(atom3)
            ):
                break
            bond1.SetBondOrder(bond1.GetBondOrder() + 1)
            bond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
            bond2.SetBondOrder(bond2.GetBondOrder() + 1)
            consume_two_electron_center(atom1)
            atom1.SetFormalCharge(atom1.GetFormalCharge() - 1)
            set_unpaired_electron_count(atom2, get_unpaired_electron_count(atom2) - 1)
            atom2.SetFormalCharge(atom2.GetFormalCharge() + 1)
            consume_two_electron_center(atom3)
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
            if get_unpaired_electron_count(atom2) != 1:
                break
            bond1.SetBondOrder(bond1.GetBondOrder() + 1)
            atom1.SetFormalCharge(atom1.GetFormalCharge() + 1)
            set_unpaired_electron_count(atom2, get_unpaired_electron_count(atom2) - 1)
            given_charge -= 1
            hit = True
    return omol, given_charge, hit


def eliminate_1_3_dipole_postive(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    """Use positive-charge budget to close a bond in a 1,3-dipole motif.

    Forming the extra bond consumes exactly one real unpaired electron and shifts
    one unit of formal charge to the middle atom. The rule runs while
    ``given_charge >= 0``; lone pairs and unresolved centers do not satisfy it.
    """

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False
    if given_charge < 0:
        return omol, given_charge, hit

    res: List[Tuple[int, int, int]] = list(smarts.ELIM_1_3_DIPOLE_POSTIVE.findall(omol))
    while given_charge >= 0 and len(res):
        idxs = res.pop(0)
        atom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
        atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        bond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if (
            get_unpaired_electron_count(atom3) == 1
            and consts.NON_METAL_DICT[atom2.GetAtomicNum()].num_outer_electrons
            + atom2.GetTotalValence()
            == 8
        ):
            atom2.SetFormalCharge(atom2.GetFormalCharge() + 1)
            bond2.SetBondOrder(int(bond2.GetBondOrder() + 1))
            set_unpaired_electron_count(atom3, get_unpaired_electron_count(atom3) - 1)
            given_charge -= 1
            hit = True
    return omol, given_charge, hit


def eliminate_possible_cp_like_radical_anion(
    omol: pybel.Molecule,
    given_charge: int,
    total_radical_electrons: int,
) -> tuple[pybel.Molecule, int, bool]:
    """Localize an excess cyclopentadienyl-like radical as an aromatic anion.

    One real unpaired electron is consumed only while the global electron pool
    exceeds the requested radical count and pending charge budget. Diradicals,
    active lone pairs, and unresolved centers are not collapsed into this state.
    """

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False

    def available_unpaired_electrons() -> int:
        return (
            sum(
                get_unpaired_electron_count(cast(ob.OBAtom, atom))
                for atom in ob.OBMolAtomIter(obmol)
            )
            - total_radical_electrons
            - abs(given_charge)
        )

    while res := smarts.ELIM_NEGATIVE_CP.findall(omol):
        if available_unpaired_electrons() < 1:
            break
        idxs = cast(Tuple[int, int, int, int, int], res.pop(0))
        atom = cast(ob.OBAtom, obmol.GetAtom(idxs[4]))
        if atom.GetFormalCharge() != 0:
            break
        if get_unpaired_electron_count(atom) != 1:
            break
        set_unpaired_electron_count(atom, 0)
        atom.SetFormalCharge(-1)
        given_charge += 1
        hit = True
    return omol, given_charge, hit


def eliminate_positive_charges(
    omol: pybel.Molecule, given_charge: int
) -> tuple[pybel.Molecule, int, bool]:
    """Consume classified electron sources into positive formal charge.

    Radical actions can consume several explicit unpaired electrons at one atom.
    Separately, a pure unresolved center may satisfy two remaining charge units by
    becoming ``+2`` and clearing its marker; it is not treated as a diradical.
    """

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
    omol: pybel.Molecule,
    given_charge: int,
    total_radical_electrons: int | None = None,
) -> tuple[pybel.Molecule, int, bool]:
    """Consume classified electron sources into negative formal charge.

    At zero budget this may intentionally generate a charge-separated resonance
    contributor from a real radical. When the optional global radical target is
    supplied, the target count is reserved before doing so. A pure unresolved
    center is eligible only when at least two negative charge units remain and is
    consumed directly as ``-2``. Active lone pairs are not consumed.
    """

    obmol = cast(ob.OBMol, omol.OBMol)
    hit = False

    if total_radical_electrons is not None and given_charge == 0:
        real_unpaired = sum(
            get_unpaired_electron_count(cast(ob.OBAtom, atom)) for atom in ob.OBMolAtomIter(obmol)
        )
        if real_unpaired <= total_radical_electrons:
            return omol, given_charge, hit

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
    "eliminate_1_3_dipole_postive",
    "eliminate_CN_in_doubt",
    "eliminate_NNN",
    "eliminate_carboxyl",
    "eliminate_carbene_neighbor_heteroatom",
    "eliminate_possible_cp_like_radical_anion",
    "eliminate_high_positive_charge_atoms",
    "eliminate_negative_charges",
    "eliminate_positive_charges",
]
