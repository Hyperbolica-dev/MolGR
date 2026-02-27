"""
Author: TMJ
Date: 2026-02-22 20:18:39
LastEditors: TMJ
LastEditTime: 2026-02-25 19:18:00
Description: 用于判定对于同一个分子坐标输入下不同的分子图重建算法结果的一致性，只要和标准答案在共振异构级别上有一个一致，就认为是一致的。
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from functools import partial
from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import ResonanceMolSupplier, inchi
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from .radical_resonance import enumerate_resonance_radical


# ================================
# Enum
# ================================


class EquivalenceMethod(str, Enum):
    IDEAL = "ideal"
    ISOMORPHIC = "isomorphic"
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
    resonance: Optional[ResonanceDetail] = None


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


def check_equivalence(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    use_chirality: bool = True,
    max_resonance: int = 50,
    resonance_flags: Chem.ResonanceFlags = Chem.ResonanceFlags.UNCONSTRAINED_CATIONS
    | Chem.ResonanceFlags.UNCONSTRAINED_ANIONS,
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

    m1 = Chem.RemoveHs(Chem.Mol(mol1))
    m2 = Chem.RemoveHs(Chem.Mol(mol2))

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

    f1 = CalcMolFormula(m1)
    f2 = CalcMolFormula(m2)

    checks = EquivalenceChecks(
        formal_charge=PropertyCheck(fc1, fc2, fc1 == fc2),
        radical_electrons=PropertyCheck(rad1, rad2, rad1 == rad2),
        num_atoms=PropertyCheck(n1, n2, n1 == n2),
        formula=PropertyCheck(f1, f2, f1 == f2),
    )

    info = EquivalenceInfo(checks=checks)
    phase_errors: list[str] = []

    # Early exits
    if not checks.formal_charge.passed:
        info.reason = "Not equivalent: total formal charges differ."
        return False, info

    if not checks.radical_electrons.passed:
        info.reason = "Not equivalent: total radical electron counts differ."
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
    m1_in_m2 = m2.HasSubstructMatch(m1, useChirality=use_chirality)
    m2_in_m1 = m1.HasSubstructMatch(m2, useChirality=use_chirality)

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
            for rm_charge in ResonanceMolSupplier(
                mh2, maxStructs=max_resonance, flags=resonance_flags
            )
            for rm_radical in enumerate_resonance_radical(rm_charge, depth=3)
            if rm_radical is not None
        }
        resonance_count = len(res2_set)

        for rm_charge in ResonanceMolSupplier(mh1, maxStructs=max_resonance, flags=resonance_flags):
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

    if phase_errors:
        return _apply_exception_fallback(info, m1, m2, phase_errors)

    info.reason = "Not equivalent: none of ideal, isomorphic, or resonance checks matched."
    return False, info
