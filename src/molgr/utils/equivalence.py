"""
Author: TMJ
Date: 2026-02-22 20:18:39
LastEditors: TMJ
LastEditTime: 2026-04-19 20:48:18
Description: 用于判定对于同一个分子坐标输入下不同的分子图重建算法结果的一致性，只要和标准答案在共振异构级别上有一个一致，就认为是一致的。
"""

from __future__ import annotations

from collections import Counter
from contextlib import suppress
from dataclasses import dataclass
from enum import Enum
from functools import partial
from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import ResonanceMolSupplier, inchi

from .radical_resonance import enumerate_resonance_radical


# ================================
# Enum
# ================================


class EquivalenceMethod(str, Enum):
    IDEAL = "ideal"
    ISOMORPHIC = "isomorphic"
    CARBENE_ZWITTERION = "carbene_zwitterion"
    COORDINATION_STRIPPED = "coordination_stripped"
    RESONANCE = "resonance"
    INCHI_CONNECTIVITY = "inchi_connectivity"


# ================================
# Dataclasses
# ================================


@dataclass
class PropertyCheck:
    mol1: int | str
    mol2: int | str
    passed: bool


@dataclass
class CanonicalSmilesDetail:
    mol1: str
    mol2: str
    use_chirality: bool


@dataclass
class IsomorphicDetail:
    mol1_in_mol2: bool
    mol2_in_mol1: bool


@dataclass
class ResonanceDetail:
    max_resonance: int
    resonance_flags: int
    mol2_resonance_count: int
    hit_smiles: Optional[str] = None


@dataclass
class CarbeneZwitterionDetail:
    mol1_normalized: str
    mol2_normalized: str


@dataclass
class CoordinationStrippedDetail:
    mol1_stripped: str
    mol2_stripped: str


@dataclass
class EquivalenceChecks:
    formal_charge: PropertyCheck
    radical_electrons: PropertyCheck
    num_atoms: PropertyCheck
    formula: PropertyCheck


@dataclass
class EquivalenceInfo:
    equivalent: bool = False
    method: Optional[EquivalenceMethod] = None
    reason: str = ""

    checks: Optional[EquivalenceChecks] = None

    canonical_smiles: Optional[CanonicalSmilesDetail] = None
    isomorphic: Optional[IsomorphicDetail] = None
    carbene_zwitterion: Optional[CarbeneZwitterionDetail] = None
    coordination_stripped: Optional[CoordinationStrippedDetail] = None
    resonance: Optional[ResonanceDetail] = None


_NON_METAL_ATOMIC_NUMBERS = frozenset(
    {
        1,
        5,
        6,
        7,
        8,
        9,
        14,
        15,
        16,
        17,
        33,
        34,
        35,
        51,
        52,
        53,
    }
)

_SULFIMIDE_RESONANCE_PATTERN = Chem.MolFromSmarts("[N:1]=[S:2](=[O:3])([O-:4])[!#1:5]")
_THIOSEMICARBAZONE_RESONANCE_PATTERNS = tuple(
    pattern
    for pattern in (
        Chem.MolFromSmarts("[N:1]-[N:2]=[C:3]([S:4])[N:5]"),
        Chem.MolFromSmarts("[N:1]=[N:2]-[C:3](=[S:4])[N:5]"),
    )
    if pattern is not None
)


def _canon_smiles(m: Chem.Mol, use_chirality: bool) -> str:
    if Chem.SanitizeMol(m) != Chem.SanitizeFlags.SANITIZE_NONE:
        raise ValueError("Molecule is not sanitized.")
    return Chem.MolToSmiles(m, canonical=True, isomericSmiles=use_chirality)


def _inchi_connectivity_key(m: Chem.Mol) -> str | None:
    try:
        key = inchi.MolToInchiKey(m)
    except Exception:  # noqa: BLE001
        return None
    if not key:
        return None
    return key.split("-")[0]


def _is_metal_atom(atom: Chem.Atom) -> bool:
    return int(atom.GetAtomicNum()) not in _NON_METAL_ATOMIC_NUMBERS


def _safe_smiles(mol: Chem.Mol, *, use_chirality: bool) -> str:
    clone = Chem.Mol(mol)
    clone.UpdatePropertyCache(strict=False)
    with suppress(Exception):
        Chem.SanitizeMol(clone)
    return Chem.MolToSmiles(clone, canonical=True, isomericSmiles=use_chirality)


def _best_effort_sanitize(mol: Chem.Mol) -> str | None:
    mol.UpdatePropertyCache(strict=False)
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:  # noqa: BLE001
        mol.UpdatePropertyCache(strict=False)
        return f"{type(exc).__name__}: {exc}"
    return None


def _remove_hs_without_sanitize(mol: Chem.Mol) -> Chem.Mol:
    try:
        stripped = Chem.RemoveHs(mol, sanitize=False)
    except TypeError:
        stripped = Chem.RemoveHs(mol)
    stripped.UpdatePropertyCache(strict=False)
    return stripped


def _prepare_equivalence_mol(mol: Chem.Mol) -> tuple[Chem.Mol, list[str]]:
    prepared = Chem.Mol(mol)
    prep_errors: list[str] = []

    sanitize_error = _best_effort_sanitize(prepared)
    if sanitize_error is not None:
        prep_errors.append(f"sanitize failed: {sanitize_error}")

    try:
        prepared = _remove_hs_without_sanitize(prepared)
    except Exception as exc:  # noqa: BLE001
        prep_errors.append(f"RemoveHs failed: {type(exc).__name__}: {exc}")
        return prepared, prep_errors

    post_removehs_error = _best_effort_sanitize(prepared)
    if post_removehs_error is not None:
        prep_errors.append(f"post-RemoveHs sanitize failed: {post_removehs_error}")

    return prepared, prep_errors


def _iter_resonance_structures(
    mol: Chem.Mol,
    *,
    max_resonance: int,
    resonance_flags: Chem.ResonanceFlags,
):
    for resonance_mol in ResonanceMolSupplier(
        mol,
        maxStructs=max_resonance,
        flags=resonance_flags,
    ):
        if resonance_mol is not None:
            yield resonance_mol


def _formula_key_without_hydrogen(mol: Chem.Mol) -> str:
    counts = Counter(atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1)
    if not counts:
        return ""

    ordered_symbols: list[str] = []
    if "C" in counts:
        ordered_symbols.append("C")
    ordered_symbols.extend(symbol for symbol in sorted(counts) if symbol != "C")
    return "".join(
        f"{symbol}{counts[symbol] if counts[symbol] != 1 else ''}" for symbol in ordered_symbols
    )


def _strip_metal_coordination_bonds(mol: Chem.Mol) -> tuple[Chem.Mol, bool]:
    rw_mol = Chem.RWMol(Chem.Mol(mol))
    bonds_to_remove: list[tuple[int, int]] = []
    for bond in rw_mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if _is_metal_atom(begin_atom) == _is_metal_atom(end_atom):
            continue
        bonds_to_remove.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    for begin_atom_idx, end_atom_idx in bonds_to_remove:
        rw_mol.RemoveBond(begin_atom_idx, end_atom_idx)
    stripped = rw_mol.GetMol()
    stripped.UpdatePropertyCache(strict=False)
    with suppress(Exception):
        Chem.SanitizeMol(stripped)
    return stripped, bool(bonds_to_remove)


def _is_carbene_zwitterion_candidate_ring(mol: Chem.Mol, ring_atom_ids: tuple[int, ...]) -> bool:
    ring_atom_set = set(ring_atom_ids)
    ring_atoms = [mol.GetAtomWithIdx(atom_idx) for atom_idx in ring_atom_ids]
    nitrogen_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
    carbon_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
    if nitrogen_count != 2 or carbon_count != 3:
        return False

    has_carbene_like_carbon = any(
        atom.GetAtomicNum() == 6 and atom.GetNumRadicalElectrons() > 0 for atom in ring_atoms
    )
    has_charge_separated_ring = any(
        atom.GetFormalCharge() != 0 or atom.GetNumRadicalElectrons() != 0 for atom in ring_atoms
    )
    has_metal_bound_carbene_carbon = any(
        atom.GetAtomicNum() == 6
        and any(
            bond.GetBondType() == Chem.BondType.DATIVE
            and bond.GetOtherAtom(atom).GetIdx() not in ring_atom_set
            and _is_metal_atom(bond.GetOtherAtom(atom))
            for bond in atom.GetBonds()
        )
        for atom in ring_atoms
    )

    return has_carbene_like_carbon or has_charge_separated_ring or has_metal_bound_carbene_carbon


def _normalize_carbene_zwitterion_forms(mol: Chem.Mol) -> Chem.Mol:
    normalized = Chem.Mol(mol)
    rw_mol = Chem.RWMol(normalized)
    ring_info = normalized.GetRingInfo()

    for ring_atom_ids in ring_info.AtomRings():
        if len(ring_atom_ids) != 5:
            continue
        ring_atom_ids_tuple = tuple(int(atom_idx) for atom_idx in ring_atom_ids)
        if not _is_carbene_zwitterion_candidate_ring(normalized, ring_atom_ids_tuple):
            continue

        ring_atom_set = set(ring_atom_ids_tuple)
        for atom_idx in ring_atom_ids_tuple:
            atom = rw_mol.GetAtomWithIdx(atom_idx)
            atom.SetFormalCharge(0)
            atom.SetNumRadicalElectrons(0)
            atom.SetIsAromatic(False)
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(False)
        for bond in normalized.GetBonds():
            begin_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            if begin_atom_idx not in ring_atom_set or end_atom_idx not in ring_atom_set:
                continue
            ring_bond = rw_mol.GetBondBetweenAtoms(begin_atom_idx, end_atom_idx)
            if ring_bond is None:
                continue
            ring_bond.SetBondType(Chem.BondType.SINGLE)
            ring_bond.SetIsAromatic(False)

    normalized_mol = rw_mol.GetMol()
    normalized_mol.UpdatePropertyCache(strict=False)
    return normalized_mol


def _carbene_zwitterion_normalized_smiles(
    mol: Chem.Mol,
    *,
    use_chirality: bool,
) -> str:
    normalized = _normalize_carbene_zwitterion_forms(mol)
    return Chem.MolToSmiles(
        normalized,
        canonical=True,
        isomericSmiles=use_chirality,
    )


def _normalize_special_resonance_forms(mol: Chem.Mol) -> Chem.Mol:
    normalized = Chem.Mol(mol)
    rw_mol = Chem.RWMol(normalized)

    if _SULFIMIDE_RESONANCE_PATTERN is not None:
        base = rw_mol.GetMol()
        for match in base.GetSubstructMatches(_SULFIMIDE_RESONANCE_PATTERN):
            n_idx, s_idx, _, o_minus_idx, _ = (int(atom_idx) for atom_idx in match)
            n_atom = rw_mol.GetAtomWithIdx(n_idx)
            o_minus_atom = rw_mol.GetAtomWithIdx(o_minus_idx)
            ns_bond = rw_mol.GetBondBetweenAtoms(n_idx, s_idx)
            so_bond = rw_mol.GetBondBetweenAtoms(s_idx, o_minus_idx)
            if ns_bond is None or so_bond is None:
                continue
            ns_bond.SetBondType(Chem.BondType.SINGLE)
            so_bond.SetBondType(Chem.BondType.DOUBLE)
            n_atom.SetFormalCharge(n_atom.GetFormalCharge() - 1)
            o_minus_atom.SetFormalCharge(0)

    def collapse_thiosemicarbazone_match(match: tuple[int, ...]) -> None:
        n1_idx, n2_idx, c_idx, s_idx, n3_idx = (int(atom_idx) for atom_idx in match)
        bond_pairs = (
            (n1_idx, n2_idx),
            (n2_idx, c_idx),
            (c_idx, s_idx),
            (c_idx, n3_idx),
        )
        for begin_idx, end_idx in bond_pairs:
            bond = rw_mol.GetBondBetweenAtoms(begin_idx, end_idx)
            if bond is not None:
                bond.SetBondType(Chem.BondType.SINGLE)
                bond.SetIsAromatic(False)
        for atom_idx in (n1_idx, n2_idx, c_idx, s_idx, n3_idx):
            atom = rw_mol.GetAtomWithIdx(atom_idx)
            atom.SetFormalCharge(0)
            atom.SetNumRadicalElectrons(0)
            atom.SetIsAromatic(False)

    base = rw_mol.GetMol()
    for pattern in _THIOSEMICARBAZONE_RESONANCE_PATTERNS:
        for match in base.GetSubstructMatches(pattern):
            collapse_thiosemicarbazone_match(match)

    normalized_mol = rw_mol.GetMol()
    normalized_mol.UpdatePropertyCache(strict=False)
    with suppress(Exception):
        Chem.SanitizeMol(normalized_mol)
    return normalized_mol


def _special_resonance_normalized_smiles(
    mol: Chem.Mol,
    *,
    use_chirality: bool,
) -> str:
    normalized = _normalize_special_resonance_forms(mol)
    return Chem.MolToSmiles(
        normalized,
        canonical=True,
        isomericSmiles=use_chirality,
    )


def _apply_exception_fallback(
    info: EquivalenceInfo,
    m1: Chem.Mol,
    m2: Chem.Mol,
    phase_errors: list[str],
) -> Tuple[bool, EquivalenceInfo]:
    inchi_key1 = _inchi_connectivity_key(m1)
    inchi_key2 = _inchi_connectivity_key(m2)
    joined_errors = "; ".join(phase_errors)

    if inchi_key1 is not None and inchi_key1 == inchi_key2:
        info.equivalent = True
        info.method = EquivalenceMethod.INCHI_CONNECTIVITY
        info.reason = (
            "Equivalent: InChI connectivity key fallback matched after equivalence-phase exception(s): "
            f"{joined_errors}"
        )
        return True, info

    info.reason = (
        "Not equivalent: equivalence-phase exception(s) prevented full comparison and "
        f"InChI connectivity fallback did not match: {joined_errors}"
    )
    return False, info


def _check_coordination_stripped_equivalence(
    info: EquivalenceInfo,
    m1: Chem.Mol,
    m2: Chem.Mol,
    *,
    use_chirality: bool,
    max_resonance: int,
    resonance_flags: Chem.ResonanceFlags,
    phase_errors: list[str],
    coordination_bonds_already_stripped: bool,
) -> Tuple[bool, EquivalenceInfo]:
    if coordination_bonds_already_stripped:
        return False, info

    try:
        stripped_m1, stripped_m1_changed = _strip_metal_coordination_bonds(m1)
        stripped_m2, stripped_m2_changed = _strip_metal_coordination_bonds(m2)
        if not stripped_m1_changed and not stripped_m2_changed:
            return False, info

        info.coordination_stripped = CoordinationStrippedDetail(
            mol1_stripped=_safe_smiles(stripped_m1, use_chirality=use_chirality),
            mol2_stripped=_safe_smiles(stripped_m2, use_chirality=use_chirality),
        )

        stripped_equivalent, stripped_info = check_equivalence(
            stripped_m1,
            stripped_m2,
            use_chirality=use_chirality,
            max_resonance=max_resonance,
            resonance_flags=resonance_flags,
            _coordination_bonds_already_stripped=True,
        )
        if stripped_equivalent:
            info.equivalent = True
            info.method = EquivalenceMethod.COORDINATION_STRIPPED
            info.reason = (
                "Equivalent after stripping metal-ligand coordination bonds: "
                f"{stripped_info.reason}"
            )
            return True, info
    except Exception as exc:  # noqa: BLE001
        phase_errors.append(f"coordination-stripped check failed: {type(exc).__name__}: {exc}")

    return False, info


def check_equivalence(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    use_chirality: bool = True,
    max_resonance: int = 50,
    resonance_flags: Chem.ResonanceFlags = Chem.ResonanceFlags.UNCONSTRAINED_CATIONS
    | Chem.ResonanceFlags.UNCONSTRAINED_ANIONS,
    _coordination_bonds_already_stripped: bool = False,
) -> Tuple[bool, EquivalenceInfo]:
    """
    Judge whether two molecules are equivalent.

    1) ideally equivalence: when canonical SMILES are equal.

    2) isomorphic equivalence: when two molecules are isomorphic, which means they are the subgraph of each other.

    3) resonance equivalence: when two molecules can be transformed into each other by resonance.

    Any of the above three methods is True, then the two molecules are equivalent.

    Parameters:
        mol1 (Chem.Mol): The first molecule.
        mol2 (Chem.Mol): The second molecule.
        use_chirality (bool, optional): Whether to use chirality. Defaults to True.
        max_resonance (int, optional): The maximum number of resonance structures to generate. Defaults to 100.
        resonance_flags (Chem.ResonanceFlags, optional): The resonance flags. Defaults to Chem.ResonanceFlags.UNCONSTRAINED_CATIONS | Chem.ResonanceFlags.UNCONSTRAINED_ANIONS.

    Returns:
        (Tuple[bool, EquivalenceInfo]): Whether the two molecules are equivalent and the reason.
    """
    original_m1 = Chem.Mol(mol1)
    original_m2 = Chem.Mol(mol2)
    m1, prep_errors_1 = _prepare_equivalence_mol(mol1)
    m2, prep_errors_2 = _prepare_equivalence_mol(mol2)

    fc1 = sum(m1.GetAtomWithIdx(atom_id).GetFormalCharge() for atom_id in range(m1.GetNumAtoms()))
    fc2 = sum(m2.GetAtomWithIdx(atom_id).GetFormalCharge() for atom_id in range(m2.GetNumAtoms()))

    rad1 = sum(
        m1.GetAtomWithIdx(atom_id).GetNumRadicalElectrons() for atom_id in range(m1.GetNumAtoms())
    )
    rad2 = sum(
        m2.GetAtomWithIdx(atom_id).GetNumRadicalElectrons() for atom_id in range(m2.GetNumAtoms())
    )

    n1 = m1.GetNumAtoms()
    n2 = m2.GetNumAtoms()

    f1 = _formula_key_without_hydrogen(m1)
    f2 = _formula_key_without_hydrogen(m2)

    checks = EquivalenceChecks(
        formal_charge=PropertyCheck(fc1, fc2, fc1 == fc2),
        radical_electrons=PropertyCheck(rad1, rad2, rad1 == rad2),
        num_atoms=PropertyCheck(n1, n2, n1 == n2),
        formula=PropertyCheck(f1, f2, f1 == f2),
    )

    info = EquivalenceInfo(checks=checks)
    phase_errors: list[str] = [*prep_errors_1, *prep_errors_2]

    # Early exits
    if not checks.formal_charge.passed:
        info.reason = "Not equivalent: total formal charges differ."
        return False, info

    if not checks.num_atoms.passed:
        inchi_key1 = _inchi_connectivity_key(m1)
        inchi_key2 = _inchi_connectivity_key(m2)
        if inchi_key1 is not None and inchi_key1 == inchi_key2:
            info.equivalent = True
            info.method = EquivalenceMethod.INCHI_CONNECTIVITY
            info.reason = "Equivalent: InChI connectivity key matches despite explicit-hydrogen atom-count mismatch."
            return True, info
        info.reason = "Not equivalent: number of atoms differs."
        return False, info

    if not checks.formula.passed:
        info.reason = "Not equivalent: molecular formulas differ."
        return False, info

    if not checks.radical_electrons.passed:
        try:
            normalized_smiles_1 = _carbene_zwitterion_normalized_smiles(
                m1,
                use_chirality=use_chirality,
            )
            normalized_smiles_2 = _carbene_zwitterion_normalized_smiles(
                m2,
                use_chirality=use_chirality,
            )
            info.carbene_zwitterion = CarbeneZwitterionDetail(
                mol1_normalized=normalized_smiles_1,
                mol2_normalized=normalized_smiles_2,
            )
            if normalized_smiles_1 == normalized_smiles_2:
                info.equivalent = True
                info.method = EquivalenceMethod.CARBENE_ZWITTERION
                info.reason = (
                    "Equivalent: carbene/zwitterion normalization matched despite "
                    "different local charge or radical representation."
                )
                return True, info
        except Exception as exc:  # noqa: BLE001
            phase_errors.append(f"carbene/zwitterion check failed: {type(exc).__name__}: {exc}")
        stripped_equivalent, info = _check_coordination_stripped_equivalence(
            info,
            original_m1,
            original_m2,
            use_chirality=use_chirality,
            max_resonance=max_resonance,
            resonance_flags=resonance_flags,
            phase_errors=phase_errors,
            coordination_bonds_already_stripped=_coordination_bonds_already_stripped,
        )
        if stripped_equivalent:
            return True, info
        info.reason = "Not equivalent: total radical electron counts differ."
        return False, info

    # 1) Ideal equivalence
    try:
        s1 = _canon_smiles(m1, use_chirality)
        s2 = _canon_smiles(m2, use_chirality)

        info.canonical_smiles = CanonicalSmilesDetail(s1, s2, use_chirality)

        if s1 == s2:
            info.equivalent = True
            info.method = EquivalenceMethod.IDEAL
            info.reason = "Equivalent: canonical SMILES are identical."
            return True, info
    except Exception as exc:  # noqa: BLE001
        phase_errors.append(f"ideal check failed: {type(exc).__name__}: {exc}")

    # 2) Isomorphic equivalence
    m1_in_m2 = False
    m2_in_m1 = False
    try:
        m1_in_m2 = m2.HasSubstructMatch(m1, useChirality=use_chirality)
        m2_in_m1 = m1.HasSubstructMatch(m2, useChirality=use_chirality)
    except Exception as exc:  # noqa: BLE001
        phase_errors.append(f"isomorphic check failed: {type(exc).__name__}: {exc}")

    info.isomorphic = IsomorphicDetail(m1_in_m2, m2_in_m1)

    if m1_in_m2 and m2_in_m1:
        info.equivalent = True
        info.method = EquivalenceMethod.ISOMORPHIC
        info.reason = "Equivalent: molecules are mutually substructure-matching (graph isomorphic)."
        return True, info

    # 3) Resonance equivalence
    hit_smiles = None
    resonance_count = 0
    try:
        canon = partial(_canon_smiles, use_chirality=use_chirality)

        mh1 = Chem.AddHs(m1)
        mh2 = Chem.AddHs(m2)
        try:
            Chem.Kekulize(mh1, clearAromaticFlags=True)
            Chem.Kekulize(mh2, clearAromaticFlags=True)
        except Chem.rdchem.KekulizeException:
            pass

        res2_set = {
            canon(rm_radical)
            for rm_charge in _iter_resonance_structures(
                mh2,
                max_resonance=max_resonance,
                resonance_flags=resonance_flags,
            )
            for rm_radical in enumerate_resonance_radical(rm_charge, depth=3)
            if rm_radical is not None
        }
        resonance_count = len(res2_set)

        for rm_charge in _iter_resonance_structures(
            mh1,
            max_resonance=max_resonance,
            resonance_flags=resonance_flags,
        ):
            for rm_radical in enumerate_resonance_radical(rm_charge, depth=3):
                if rm_radical is None:
                    continue
                s = canon(rm_radical)
                if s in res2_set:
                    hit_smiles = s
                    break
    except Exception as exc:  # noqa: BLE001
        phase_errors.append(f"resonance check failed: {type(exc).__name__}: {exc}")

    info.resonance = ResonanceDetail(
        max_resonance=max_resonance,
        resonance_flags=int(resonance_flags),
        mol2_resonance_count=resonance_count,
        hit_smiles=hit_smiles,
    )

    if hit_smiles is not None:
        info.equivalent = True
        info.method = EquivalenceMethod.RESONANCE
        info.reason = "Equivalent: at least one resonance structure matches in canonical SMILES."
        return True, info

    try:
        normalized_smiles_1 = _special_resonance_normalized_smiles(
            m1,
            use_chirality=use_chirality,
        )
        normalized_smiles_2 = _special_resonance_normalized_smiles(
            m2,
            use_chirality=use_chirality,
        )
        if normalized_smiles_1 == normalized_smiles_2:
            info.resonance.hit_smiles = normalized_smiles_1
            info.equivalent = True
            info.method = EquivalenceMethod.RESONANCE
            info.reason = "Equivalent: special resonance normalization matched canonical SMILES."
            return True, info
    except Exception as exc:  # noqa: BLE001
        phase_errors.append(f"special resonance normalization failed: {type(exc).__name__}: {exc}")

    stripped_equivalent, info = _check_coordination_stripped_equivalence(
        info,
        original_m1,
        original_m2,
        use_chirality=use_chirality,
        max_resonance=max_resonance,
        resonance_flags=resonance_flags,
        phase_errors=phase_errors,
        coordination_bonds_already_stripped=_coordination_bonds_already_stripped,
    )
    if stripped_equivalent:
        return True, info

    if phase_errors:
        return _apply_exception_fallback(info, m1, m2, phase_errors)

    info.reason = "Not equivalent: none of ideal, isomorphic, or resonance checks matched."
    return False, info
