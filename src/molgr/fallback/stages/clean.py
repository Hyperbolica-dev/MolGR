"""Radical/resonance cleanup stages for fallback."""

from __future__ import annotations

from typing import List, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.stages.fresh import assign_charge_radical_for_atom
from molgr.fallback.utils import consts, smarts
from molgr.fallback.utils.electrons import (
    get_lone_pair_count,
    get_unpaired_electron_count,
    has_unresolved_two_electron_center,
    set_lone_pair_count,
    set_unpaired_electron_count,
    set_unresolved_two_electron_center,
)


def clean_carbene_neighbor_unsaturated(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    """Move an unresolved two-electron center through an adjacent pi bond.

    ``A*-B=C`` becomes ``A=B-C``. The unresolved marker at A is consumed and
    replaced by one real unpaired electron at A plus one at C, representing the
    homolytic redistribution. Existing active lone pairs at C are preserved.
    """

    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int]] = list(smarts.CLEAN_CARBENE_NEIGHBOR_UNSAT.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        atom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        bond12 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        bond23 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if (
            has_unresolved_two_electron_center(atom1)
            and get_unpaired_electron_count(atom3) == 0
            and not has_unresolved_two_electron_center(atom3)
            and bond12.GetBondOrder() == 1
            and bond23.GetBondOrder() == 2
        ):
            bond23.SetBondOrder(int(bond23.GetBondOrder() - 1))
            bond12.SetBondOrder(int(bond12.GetBondOrder() + 1))
            set_unresolved_two_electron_center(atom1, False)
            set_lone_pair_count(atom1, 0)
            set_unpaired_electron_count(atom1, 1)
            set_unpaired_electron_count(atom3, get_unpaired_electron_count(atom3) + 1)
            hit = True
    return omol, hit


def clean_possible_1_3_dipole(
    omol: pybel.Molecule,
    given_charge: int,
    total_radical_electrons: int,
) -> tuple[pybel.Molecule, bool]:
    """Convert excess terminal radicals into a neutral 1,3-dipole.

    A conversion is available while the current real-unpaired-electron count,
    after reserving the requested system radicals and pending formal-charge
    assignments, contains at least two excess electrons. The two matched terminal
    atoms each consume one real unpaired electron. The terminal with the lower
    average explicit bond order forms the additional bond to the middle atom;
    ties are resolved by atom index. The middle atom becomes ``+1`` and the
    opposite terminal becomes ``-1``.
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

    matches: List[Tuple[int, int, int]] = list(smarts.CLEAN_POSSIBLE_1_3_DIPOLE.findall(omol))
    for idxs in matches:
        if available_unpaired_electrons() < 2:
            break

        atom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        atom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
        atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        bond12 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        bond23 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if not all((atom1, atom2, atom3, bond12, bond23)):
            continue
        if any(atom.GetFormalCharge() != 0 for atom in (atom1, atom2, atom3)):
            continue
        if (
            get_unpaired_electron_count(atom1) < 1
            or get_unpaired_electron_count(atom3) < 1
            or has_unresolved_two_electron_center(atom1)
            or has_unresolved_two_electron_center(atom3)
        ):
            continue

        degree1 = int(atom1.GetExplicitDegree())
        degree3 = int(atom3.GetExplicitDegree())
        if degree1 <= 0 or degree3 <= 0:
            continue
        valence1 = int(atom1.GetExplicitValence())
        valence3 = int(atom3.GetExplicitValence())
        left_average_key = valence1 * degree3
        right_average_key = valence3 * degree1
        add_left_bond = left_average_key < right_average_key or (
            left_average_key == right_average_key and atom1.GetIdx() < atom3.GetIdx()
        )
        bond_to_increase = bond12 if add_left_bond else bond23
        negative_atom = atom3 if add_left_bond else atom1
        if bond_to_increase.GetBondOrder() >= 3:
            continue

        bond_to_increase.SetBondOrder(bond_to_increase.GetBondOrder() + 1)
        atom2.SetFormalCharge(1)
        negative_atom.SetFormalCharge(-1)
        set_unpaired_electron_count(atom1, get_unpaired_electron_count(atom1) - 1)
        set_unpaired_electron_count(atom3, get_unpaired_electron_count(atom3) - 1)
        hit = True

    return omol, hit


def clean_neighbor_radicals(
    omol: pybel.Molecule,
    given_charge: int,
    total_radical_electrons: int,
) -> tuple[pybel.Molecule, bool]:
    """Resolve adjacent radical-compatible electron states.

    Each endpoint contributes either its real unpaired electrons or the two
    electrons represented by an unresolved center. A matched pair consumes one
    electron from each endpoint and increases the bond by one order. The
    operation is limited by the excess-electron budget, so two unresolved
    centers can consume one pair and become a double bond while their remaining
    electrons are classified on the next pass. A real radical adjacent to an
    unresolved center therefore preserves the real-radical count after one
    pair is consumed, while two real radicals require two excess electrons.
    """

    obmol = cast(ob.OBMol, omol.OBMol)

    def available_electron_budget() -> int:
        return (
            sum(
                get_unpaired_electron_count(cast(ob.OBAtom, atom))
                + (2 if has_unresolved_two_electron_center(cast(ob.OBAtom, atom)) else 0)
                for atom in ob.OBMolAtomIter(obmol)
            )
            - total_radical_electrons
            - abs(given_charge)
        )

    hit = False
    for bond in list(ob.OBMolBondIter(obmol)):
        begin_atom = cast(ob.OBAtom, bond.GetBeginAtom())
        end_atom = cast(ob.OBAtom, bond.GetEndAtom())
        begin_unpaired = get_unpaired_electron_count(begin_atom)
        end_unpaired = get_unpaired_electron_count(end_atom)
        begin_unresolved = has_unresolved_two_electron_center(begin_atom)
        end_unresolved = has_unresolved_two_electron_center(end_atom)

        begin_capacity = 2 if begin_unresolved else begin_unpaired
        end_capacity = 2 if end_unresolved else end_unpaired
        if begin_capacity <= 0 or end_capacity <= 0:
            continue

        available = available_electron_budget()
        if available < 2:
            continue
        bond_to_add = min(begin_capacity, end_capacity, available // 2)
        if bond_to_add <= 0:
            continue

        bond.SetBondOrder(bond.GetBondOrder() + bond_to_add)
        for atom, unresolved, capacity in (
            (begin_atom, begin_unresolved, begin_capacity),
            (end_atom, end_unresolved, end_capacity),
        ):
            remaining = capacity - bond_to_add
            if unresolved:
                set_unresolved_two_electron_center(atom, False)
                set_lone_pair_count(atom, 0)
                set_unpaired_electron_count(atom, remaining)
                # A partial consumption leaves one real radical at this
                # endpoint. Re-running the local classifier would see the
                # post-bond two-electron deficit and mark it unresolved again.
                if remaining == 0:
                    assign_charge_radical_for_atom(atom)
            else:
                set_unpaired_electron_count(atom, remaining)
                assign_charge_radical_for_atom(atom)
        hit = True
    return omol, hit


def clean_1_4_radicals(
    omol: pybel.Molecule,
    given_charge: int,
    total_radical_electrons: int,
) -> tuple[pybel.Molecule, bool]:
    """Resolve separated terminal radicals or unresolved electron centers.

    For ``A-B=C-D``, consume one eligible electron state at each terminal and
    shift the middle pi bond to the two outer bonds. Real radicals obey the
    same global excess-electron budget as :func:`clean_neighbor_radicals`;
    unresolved two-electron centers are consumed directly.
    """

    obmol = cast(ob.OBMol, omol.OBMol)

    def available_unpaired_electrons() -> int:
        return (
            sum(
                get_unpaired_electron_count(cast(ob.OBAtom, atom))
                for atom in ob.OBMolAtomIter(obmol)
            )
            - total_radical_electrons
            - abs(given_charge)
        )

    hit = False
    for idxs in list(smarts.CLEAN_1_4_RADICALS.findall(omol)):
        atom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        atom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
        atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        atom4 = cast(ob.OBAtom, obmol.GetAtom(idxs[3]))
        bond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        bond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        bond3 = cast(ob.OBBond, obmol.GetBond(idxs[2], idxs[3]))
        if not all((atom1, atom2, atom3, atom4, bond1, bond2, bond3)):
            continue
        if bond1.GetBondOrder() != 1 or bond2.GetBondOrder() != 2 or bond3.GetBondOrder() != 1:
            continue

        endpoints = (atom1, atom4)
        endpoint_radicals = [get_unpaired_electron_count(atom) for atom in endpoints]
        endpoint_unresolved = [has_unresolved_two_electron_center(atom) for atom in endpoints]
        if not all(
            radical > 0 or unresolved
            for radical, unresolved in zip(endpoint_radicals, endpoint_unresolved)
        ):
            continue
        if (
            sum(radical > 0 for radical in endpoint_radicals) == 2
            and available_unpaired_electrons() < 2
        ):
            continue

        bond1.SetBondOrder(2)
        bond2.SetBondOrder(1)
        bond3.SetBondOrder(2)
        for atom, radical, unresolved in zip(endpoints, endpoint_radicals, endpoint_unresolved):
            if radical > 0:
                set_unpaired_electron_count(atom, radical - 1)
            if unresolved:
                set_unresolved_two_electron_center(atom, False)
                set_lone_pair_count(atom, 0)
            assign_charge_radical_for_atom(atom)
        hit = True
    return omol, hit


def clean_1_6_radicals(
    omol: pybel.Molecule,
    given_charge: int,
    total_radical_electrons: int,
) -> tuple[pybel.Molecule, bool]:
    """Resolve terminal electron states across an alternating six-atom path.

    For ``A-B=C-D=E-F``, consume one eligible state at A and F and shift the
    two existing pi bonds outward, producing ``A=B-C=D-E=F``. Endpoint and
    global-electron conditions are identical to :func:`clean_1_4_radicals`.
    """

    obmol = cast(ob.OBMol, omol.OBMol)

    def available_unpaired_electrons() -> int:
        return (
            sum(
                get_unpaired_electron_count(cast(ob.OBAtom, atom))
                for atom in ob.OBMolAtomIter(obmol)
            )
            - total_radical_electrons
            - abs(given_charge)
        )

    hit = False
    for idxs in list(smarts.CLEAN_1_6_RADICALS.findall(omol)):
        atoms = tuple(cast(ob.OBAtom, obmol.GetAtom(idx)) for idx in idxs)
        bonds = tuple(
            cast(ob.OBBond, obmol.GetBond(begin_idx, end_idx))
            for begin_idx, end_idx in zip(idxs, idxs[1:])
        )
        if not all((*atoms, *bonds)):
            continue
        if tuple(bond.GetBondOrder() for bond in bonds) != (1, 2, 1, 2, 1):
            continue

        endpoints = (atoms[0], atoms[-1])
        endpoint_radicals = [get_unpaired_electron_count(atom) for atom in endpoints]
        endpoint_unresolved = [has_unresolved_two_electron_center(atom) for atom in endpoints]
        if not all(
            radical > 0 or unresolved
            for radical, unresolved in zip(endpoint_radicals, endpoint_unresolved)
        ):
            continue
        if (
            sum(radical > 0 for radical in endpoint_radicals) == 2
            and available_unpaired_electrons() < 2
        ):
            continue

        for bond, order in zip(bonds, (2, 1, 2, 1, 2)):
            bond.SetBondOrder(order)
        for atom, radical, unresolved in zip(endpoints, endpoint_radicals, endpoint_unresolved):
            if radical > 0:
                set_unpaired_electron_count(atom, radical - 1)
            if unresolved:
                set_unresolved_two_electron_center(atom, False)
                set_lone_pair_count(atom, 0)
            assign_charge_radical_for_atom(atom)
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
            obatom1.GetFormalCharge() == -1
            and obatom4.GetFormalCharge() == 1
            and obbond1.GetBondOrder() == 1
            and obbond2.GetBondOrder() == 2
            and consts.NON_METAL_DICT[obatom1.GetAtomicNum()].default_valence
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
        obatom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
        obatom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[-1]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        obbond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        if (
            obatom1.GetFormalCharge() == -1
            and obatom2.GetFormalCharge() == 1
            and obatom3.GetFormalCharge() == 0
            and obbond1.GetBondOrder() == 2
            and obbond2.GetBondOrder() == 2
        ):
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
        obatom5_room = consts.NON_METAL_DICT[obatom5.GetAtomicNum()].default_valence - (
            obatom5.GetTotalValence() + 1
        )
        if (
            obbond4.GetBondOrder() == 1
            and obbond3.GetBondOrder() == 2
            and obbond2.GetBondOrder() == 1
            and obbond1.GetBondOrder() == 2
            and obatom1.GetFormalCharge() == 0
            and obatom5.GetFormalCharge() == -1
            and obatom5_room >= 0
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
    """Shift the rule-3 charge/pi pattern, then rebuild its two endpoints.

    Internal atoms exchange one bond-order unit on each side and keep the same
    total valence. Only the two charged endpoints change total valence, so only
    they are allowed to replace unpaired, active-lone-pair, or unresolved fields.
    """

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
        if (
            obatom1.GetFormalCharge() == 1
            and obatom5.GetFormalCharge() == -1
            and obbond1.GetBondOrder() == 2
            and obbond2.GetBondOrder() == 1
            and obbond3.GetBondOrder() == 2
            and obbond4.GetBondOrder() == 1
        ):
            obbond1.SetBondOrder(obbond1.GetBondOrder() - 1)
            obbond2.SetBondOrder(obbond2.GetBondOrder() + 1)
            obbond3.SetBondOrder(obbond3.GetBondOrder() - 1)
            obbond4.SetBondOrder(obbond4.GetBondOrder() + 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() - 1)
            obatom5.SetFormalCharge(obatom5.GetFormalCharge() + 1)
            assign_charge_radical_for_atom(obatom1)
            assign_charge_radical_for_atom(obatom5)
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
        if (
            obatom1.GetFormalCharge() == 1
            and obatom3.GetFormalCharge() == -1
            and obbond1.GetBondOrder() == 2
            and obbond2.GetBondOrder() == 1
        ):
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
        if (
            obatom1.GetFormalCharge() == 0
            and obatom3.GetFormalCharge() == -1
            and obbond1.GetBondOrder() == 2
            and obbond2.GetBondOrder() == 2
        ):
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
        obatom1_room = consts.NON_METAL_DICT[obatom1.GetAtomicNum()].default_valence - (
            obatom1.GetTotalValence() + 1
        )
        if (
            obbond1.GetBondOrder() == 1
            and obbond2.GetBondOrder() == 2
            and obatom1.GetFormalCharge() == -1
            and obatom3.GetFormalCharge() == 0
            and obatom1_room >= 0
        ):
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() + 1)
            obatom3.SetFormalCharge(obatom3.GetFormalCharge() - 1)
            hit = True
    return omol, hit


def clean_resonances_8(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    """Shift the rule-8 charge/pi pattern and rebuild both endpoint states.

    Endpoint refresh may produce unpaired electrons, active lone pairs, or a new
    unresolved center from the changed local valence; no old field is blindly
    carried across the bond/charge change.
    """

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
            obatom1.GetFormalCharge() == -1
            and obatom5.GetFormalCharge() == 0
            and obbond1.GetBondOrder() == 1
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
            assign_charge_radical_for_atom(obatom1)
            assign_charge_radical_for_atom(obatom5)
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
            obatom1.GetFormalCharge() > 0
            and obatom2.GetFormalCharge() < 0
            and consts.NON_METAL_DICT[obatom1.GetAtomicNum()].default_valence
            - obatom1.GetTotalValence()
            >= 1
            and consts.NON_METAL_DICT[obatom2.GetAtomicNum()].default_valence
            - obatom2.GetTotalValence()
            >= 1
            and obbond1.GetBondOrder() in (1, 2)
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
    """Consume two terminal monoradicals while redistributing one middle pi bond.

    Each terminal must carry exactly one real unpaired electron. Both are consumed
    to form the two new bond-order units; lone pairs and unresolved centers cannot
    substitute for either radical.
    """

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
        if get_unpaired_electron_count(obatom1) == 1 and get_unpaired_electron_count(obatom3) == 1:
            if not (
                obbond1.GetBondOrder() == 1
                and obbond2.GetBondOrder() in (2, 3)
                and obbond3.GetBondOrder() == 1
            ):
                continue
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond3.SetBondOrder(obbond3.GetBondOrder() + 1)
            set_unpaired_electron_count(obatom1, get_unpaired_electron_count(obatom1) - 1)
            set_unpaired_electron_count(obatom3, get_unpaired_electron_count(obatom3) - 1)
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
        if (
            obatom1.GetFormalCharge() == 0
            and obatom2.GetFormalCharge() == 1
            and consts.NON_METAL_DICT[obatom2.GetAtomicNum()].default_valence
            - obatom2.GetTotalValence()
            >= 1
            and (obbond1.GetBondOrder() == 1 or obbond1.GetBondOrder() == 2)
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
        if (
            obatom1.GetFormalCharge() == 0
            and obatom4.GetFormalCharge() == 1
            and consts.NON_METAL_DICT[obatom4.GetAtomicNum()].default_valence
            - obatom4.GetTotalValence()
            >= 1
            and (
                obbond1.GetBondOrder() == 1
                and obbond2.GetBondOrder() == 2
                and obbond3.GetBondOrder() == 1
            )
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
            obatom1.GetFormalCharge() == -1
            and obatom3.GetFormalCharge() == 0
            and consts.NON_METAL_DICT[obatom1.GetAtomicNum()].default_valence
            - obatom1.GetTotalValence()
            >= 1
            and obbond1.GetBondOrder() == 1
            and obbond2.GetBondOrder() == 2
        ):
            obbond1.SetBondOrder(obbond1.GetBondOrder() + 1)
            obbond2.SetBondOrder(obbond2.GetBondOrder() - 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() + 1)
            obatom3.SetFormalCharge(obatom3.GetFormalCharge() - 1)
            hit = True
    return omol, hit


def clean_resonances_14(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    """Neutralize the rule-14 triple-bond ion pair and rebuild both endpoints.

    Lowering the bond and neutralizing charges can expose new unpaired/lone-pair
    occupancy. Endpoint-local assignment replaces stale explicit fields and may
    create an unresolved C/N/P two-electron center.
    """

    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int]] = list(smarts.CLEAN_RESONANCE_14.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        obatom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        obatom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
        obbond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        if (
            obatom1.GetFormalCharge() == -1
            and obatom2.GetFormalCharge() == 1
            and obbond1.GetBondOrder() == 3
        ):
            obbond1.SetBondOrder(obbond1.GetBondOrder() - 1)
            obatom1.SetFormalCharge(obatom1.GetFormalCharge() + 1)
            obatom2.SetFormalCharge(obatom2.GetFormalCharge() - 1)
            assign_charge_radical_for_atom(obatom1)
            assign_charge_radical_for_atom(obatom2)
            hit = True
    return omol, hit


def clean_resonances_16(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    """Neutralize the rule-16 conjugated ion pair and rebuild both endpoints.

    The input requires zero explicit unpaired and active-lone-pair counts at both
    charged endpoints. After charge/pi migration, endpoint-local assignment
    derives any newly exposed electron state, including an unresolved center.
    """

    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int, int, int]] = list(smarts.CLEAN_RESONANCE_16.findall(omol))
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
            obatom1.GetFormalCharge() == -1
            and get_unpaired_electron_count(obatom1) == 0
            and get_lone_pair_count(obatom1) == 0
            and obatom5.GetFormalCharge() == 1
            and get_unpaired_electron_count(obatom5) == 0
            and get_lone_pair_count(obatom5) == 0
            and obbond1.GetBondOrder() == 1
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
            assign_charge_radical_for_atom(obatom1)
            assign_charge_radical_for_atom(obatom5)
            hit = True
    return omol, hit


def clean_resonances_17(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    """Replace a ring allene with its charge-shifted conjugated resonance form.

    ``A(-)-B=C=D`` becomes ``A=B-C(-)=D`` when B, C, and D are ring atoms.
    The transformation conserves total charge and each atom's electron count;
    only the location of one negative charge and one pi bond changes.
    """

    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int, int]] = list(smarts.CLEAN_RESONANCE_17.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        atom1 = cast(ob.OBAtom, obmol.GetAtom(idxs[0]))
        atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        bond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        bond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        bond3 = cast(ob.OBBond, obmol.GetBond(idxs[2], idxs[3]))
        atom1_room = consts.NON_METAL_DICT[atom1.GetAtomicNum()].default_valence - (
            atom1.GetTotalValence() + 1
        )
        if (
            atom1.GetFormalCharge() == -1
            and atom3.GetFormalCharge() == 0
            and bond1.GetBondOrder() == 1
            and bond2.GetBondOrder() == 2
            and bond3.GetBondOrder() == 2
            and atom1_room >= 0
        ):
            bond1.SetBondOrder(2)
            bond2.SetBondOrder(1)
            atom1.SetFormalCharge(0)
            atom3.SetFormalCharge(-1)
            hit = True
    return omol, hit


def clean_resonances_18(omol: pybel.Molecule) -> tuple[pybel.Molecule, bool]:
    """Convert an unresolved terminal diazene fragment to a charge-separated azide.

    ``[*]-N=N-[N]`` becomes ``[*]-N=[N+]=[N-]`` when the terminal nitrogen is
    an unresolved two-electron center. The transformation consumes that deferred
    pair as the additional N-N pi bond and does not alter the system radical
    budget.
    """

    obmol = cast(ob.OBMol, omol.OBMol)
    res: List[Tuple[int, int, int, int]] = list(smarts.CLEAN_RESONANCE_18.findall(omol))
    hit = False
    while len(res):
        idxs = res.pop(0)
        atom2 = cast(ob.OBAtom, obmol.GetAtom(idxs[1]))
        atom3 = cast(ob.OBAtom, obmol.GetAtom(idxs[2]))
        atom4 = cast(ob.OBAtom, obmol.GetAtom(idxs[3]))
        bond1 = cast(ob.OBBond, obmol.GetBond(idxs[0], idxs[1]))
        bond2 = cast(ob.OBBond, obmol.GetBond(idxs[1], idxs[2]))
        bond3 = cast(ob.OBBond, obmol.GetBond(idxs[2], idxs[3]))
        if not all((atom2, atom3, atom4, bond1, bond2, bond3)):
            continue
        if (
            atom2.GetFormalCharge() != 0
            or atom3.GetFormalCharge() != 0
            or atom4.GetFormalCharge() != 0
            or get_unpaired_electron_count(atom2) != 0
            or get_unpaired_electron_count(atom3) != 0
            or get_unpaired_electron_count(atom4) != 0
            or has_unresolved_two_electron_center(atom2)
            or has_unresolved_two_electron_center(atom3)
            or not has_unresolved_two_electron_center(atom4)
            or bond1.GetBondOrder() != 1
            or bond2.GetBondOrder() != 2
            or bond3.GetBondOrder() != 1
        ):
            continue

        bond3.SetBondOrder(bond3.GetBondOrder() + 1)
        atom3.SetFormalCharge(1)
        atom4.SetFormalCharge(-1)
        set_unresolved_two_electron_center(atom4, False)
        set_lone_pair_count(atom4, 0)
        set_unpaired_electron_count(atom4, 0)
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
    omol, stage_hit = clean_resonances_14(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_16(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_17(omol)
    hit = stage_hit or hit
    omol, stage_hit = clean_resonances_18(omol)
    hit = stage_hit or hit
    return omol, hit


__all__ = [
    "clean_carbene_neighbor_unsaturated",
    "clean_1_4_radicals",
    "clean_1_6_radicals",
    "clean_neighbor_radicals",
    "clean_possible_1_3_dipole",
    "clean_resonances",
    "clean_resonances_18",
]
