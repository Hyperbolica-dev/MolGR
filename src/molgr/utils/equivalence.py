# pyright: reportCallIssue=false
from __future__ import annotations

from collections import Counter
from contextlib import suppress
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import ResonanceMolSupplier, inchi


pt = Chem.GetPeriodicTable()


class EquivalenceMethod(str, Enum):
    IDEAL = "ideal"
    INCHI_KEY = "inchi_key"
    CARBENE_ZWITTERION = "carbene_zwitterion"
    RESONANCE = "resonance"


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
class CarbeneZwitterionDetail:
    mol1_normalized: str
    mol2_normalized: str


@dataclass
class ResonanceDetail:
    max_resonance: int
    resonance_flags: int
    mol1_resonance_count: int
    mol2_resonance_count: int
    hit_smiles: Optional[str] = None


@dataclass
class EquivalenceChecks:
    formal_charge: PropertyCheck
    radical_electrons: PropertyCheck
    num_atoms: PropertyCheck
    heavy_atom_formula: PropertyCheck
    explicit_h_formula: PropertyCheck


@dataclass
class EquivalenceInfo:
    equivalent: bool = False
    method: Optional[EquivalenceMethod] = None
    reason: str = ""
    checks: Optional[EquivalenceChecks] = None
    canonical_smiles: Optional[CanonicalSmilesDetail] = None
    carbene_zwitterion: Optional[CarbeneZwitterionDetail] = None
    resonance: Optional[ResonanceDetail] = None


_NON_METAL_ATOMIC_NUMBERS = frozenset(
    {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53}
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


def _is_metal_atom(atom: Chem.Atom) -> bool:
    return int(atom.GetAtomicNum()) not in _NON_METAL_ATOMIC_NUMBERS


def _total_formal_charge(mol: Chem.Mol) -> int:
    return sum(int(atom.GetFormalCharge()) for atom in mol.GetAtoms())


def _total_radical_electrons(mol: Chem.Mol) -> int:
    return sum(int(atom.GetNumRadicalElectrons()) for atom in mol.GetAtoms())


def _formula_key(mol: Chem.Mol, *, include_hydrogen: bool) -> str:
    counts = Counter(
        atom.GetSymbol()
        for atom in mol.GetAtoms()
        if include_hydrogen or atom.GetAtomicNum() != 1
    )
    if not counts:
        return ""
    ordered_symbols: list[str] = []
    if "C" in counts:
        ordered_symbols.append("C")
    if include_hydrogen and "H" in counts:
        ordered_symbols.append("H")
    ordered_symbols.extend(symbol for symbol in sorted(counts) if symbol not in ordered_symbols)
    return "".join(f"{symbol}{counts[symbol] if counts[symbol] != 1 else ''}" for symbol in ordered_symbols)


def _safe_copy(mol: Chem.Mol) -> Chem.Mol:
    clone = Chem.Mol(mol)
    clone.UpdatePropertyCache(strict=False)
    with suppress(Exception):
        Chem.SanitizeMol(clone)
    return clone


def _clear_metal_radicals(mol: Chem.Mol) -> Chem.Mol:
    rw_mol = Chem.RWMol(mol)
    for atom in rw_mol.GetAtoms():
        if _is_metal_atom(atom):
            atom.SetNumRadicalElectrons(0)
    normalized = rw_mol.GetMol()
    normalized.UpdatePropertyCache(strict=False)
    return normalized


def _safe_smiles(mol: Chem.Mol, *, use_chirality: bool) -> str:
    clone = _safe_copy(mol)
    return Chem.MolToSmiles(clone, canonical=True, isomericSmiles=use_chirality)


def _canon_smiles(mol: Chem.Mol, use_chirality: bool) -> str:
    return Chem.CanonSmiles(_safe_smiles(mol, use_chirality=use_chirality))


def _inchi_key(mol: Chem.Mol) -> str | None:
    try:
        key = inchi.MolToInchiKey(_safe_copy(mol))
    except Exception:  # noqa: BLE001
        return None
    if not key:
        return None
    return key


def _remove_hs_without_sanitize(mol: Chem.Mol) -> Chem.Mol:
    try:
        stripped = Chem.RemoveHs(mol, sanitize=False)
    except TypeError:
        stripped = Chem.RemoveHs(mol)
    stripped.UpdatePropertyCache(strict=False)
    return stripped


def _add_hs_without_sanitize(mol: Chem.Mol) -> Chem.Mol:
    added = Chem.AddHs(mol)
    added.UpdatePropertyCache(strict=False)
    return added


def _standardize_metal_bonds(mol: Chem.Mol) -> Chem.Mol:
    rw_mol = Chem.RWMol(_safe_copy(_clear_metal_radicals(mol)))
    for bond in list(rw_mol.GetBonds()):
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        begin_is_metal = _is_metal_atom(begin_atom)
        end_is_metal = _is_metal_atom(end_atom)
        if begin_is_metal == end_is_metal:
            continue
        if bond.GetBondType() == Chem.BondType.DATIVE:
            continue

        metal_atom = begin_atom if begin_is_metal else end_atom
        ligand_atom = end_atom if begin_is_metal else begin_atom
        bond_order = max(1, int(round(bond.GetBondTypeAsDouble())))

        metal_atom.SetFormalCharge(metal_atom.GetFormalCharge() + bond_order)
        ligand_atom.SetFormalCharge(ligand_atom.GetFormalCharge() - bond_order)

        bond.SetBondType(Chem.BondType.DATIVE)
        bond.SetBondDir(Chem.BondDir.NONE)
        bond.SetIsAromatic(False)

    standardized = rw_mol.GetMol()
    standardized.UpdatePropertyCache(strict=False)
    with suppress(Exception):
        Chem.SanitizeMol(standardized)
    return standardized


def _remove_coordination_bonds(mol: Chem.Mol) -> Chem.Mol:
    rw_mol = Chem.RWMol(_safe_copy(mol))
    bonds_to_remove: list[tuple[int, int]] = []
    for bond in rw_mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if _is_metal_atom(begin_atom) == _is_metal_atom(end_atom):
            continue
        if bond.GetBondType() == Chem.BondType.DATIVE:
            bonds_to_remove.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    for begin_idx, end_idx in bonds_to_remove:
        rw_mol.RemoveBond(begin_idx, end_idx)
    stripped = rw_mol.GetMol()
    stripped.UpdatePropertyCache(strict=False)
    with suppress(Exception):
        Chem.SanitizeMol(stripped)
    return stripped


def _remove_metal_atoms(mol: Chem.Mol) -> Chem.Mol:
    rw_mol = Chem.RWMol(_safe_copy(mol))
    metal_indices = [atom.GetIdx() for atom in rw_mol.GetAtoms() if _is_metal_atom(atom)]
    for atom_idx in sorted(metal_indices, reverse=True):
        rw_mol.RemoveAtom(int(atom_idx))
    stripped = rw_mol.GetMol()
    stripped.UpdatePropertyCache(strict=False)
    with suppress(Exception):
        Chem.SanitizeMol(stripped)
    return stripped


def _prepare_organic_mol(mol: Chem.Mol) -> Chem.Mol:
    standardized = _standardize_metal_bonds(mol)
    standardized = _remove_coordination_bonds(standardized)
    standardized = _remove_metal_atoms(standardized)
    standardized = _add_hs_without_sanitize(standardized)
    standardized.UpdatePropertyCache(strict=False)
    with suppress(Exception):
        Chem.SanitizeMol(standardized)
    return standardized


def _metal_signature_key(mol: Chem.Mol) -> tuple[tuple[int, int], ...]:
    signature = []
    for atom in mol.GetAtoms():
        if not _is_metal_atom(atom):
            continue
        signature.append((int(atom.GetAtomicNum()), int(atom.GetFormalCharge())))
    return tuple(sorted(signature))


def _carbene_zwitterion_normalized_smiles(mol: Chem.Mol, *, use_chirality: bool) -> str:
    normalized = _safe_copy(mol)
    rw_mol = Chem.RWMol(normalized)
    with suppress(Exception):
        Chem.Kekulize(rw_mol)

    smarts = Chem.MolFromSmarts("[*-]=[*+]")
    if smarts is not None:
        for match in rw_mol.GetSubstructMatches(smarts):
            begin_idx, end_idx = (int(atom_idx) for atom_idx in match)
            bond = rw_mol.GetBondBetweenAtoms(begin_idx, end_idx)
            atom_begin = rw_mol.GetAtomWithIdx(begin_idx)
            atom_end = rw_mol.GetAtomWithIdx(end_idx)
            if bond is None:
                continue
            if (
                pt.GetDefaultValence(atom_begin.GetAtomicNum()) - atom_begin.GetTotalValence() == 1
                and pt.GetDefaultValence(atom_end.GetAtomicNum()) - atom_end.GetTotalValence() == -1
            ):
                bond.SetBondType(Chem.BondType.SINGLE)
                atom_begin.SetFormalCharge(atom_begin.GetFormalCharge() + 1)
                atom_end.SetFormalCharge(atom_end.GetFormalCharge() - 1)
                atom_begin.SetNumRadicalElectrons(2)
                atom_end.SetNumRadicalElectrons(0)

    normalized_mol = rw_mol.GetMol()
    normalized_mol.UpdatePropertyCache(strict=False)
    with suppress(Exception):
        Chem.SanitizeMol(normalized_mol)
    return _canon_smiles(normalized_mol, use_chirality)


def _normalize_special_resonance_forms(mol: Chem.Mol) -> Chem.Mol:
    normalized = _safe_copy(mol)
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
        for begin_idx, end_idx in (
            (n1_idx, n2_idx),
            (n2_idx, c_idx),
            (c_idx, s_idx),
            (c_idx, n3_idx),
        ):
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


def _resonance_match(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    *,
    use_chirality: bool,
    max_resonance: int,
    resonance_flags: Chem.ResonanceFlags,
) -> tuple[bool, int, int, str | None]:
    reference = _normalize_special_resonance_forms(mol2)
    try:
        reference_smiles = _canon_smiles(reference, use_chirality)
    except Exception:  # noqa: BLE001
        return False, 0, 0, None

    mol1_count = 0
    try:
        source = _normalize_special_resonance_forms(mol1)
        source = _add_hs_without_sanitize(source)
        with suppress(Exception):
            Chem.Kekulize(source, clearAromaticFlags=True)
        for resonance_mol in ResonanceMolSupplier(
            source,
            maxStructs=max_resonance,
            flags=resonance_flags,
        ):
            if resonance_mol is None:
                continue
            mol1_count += 1
            with suppress(Exception):
                if _canon_smiles(resonance_mol, use_chirality) == reference_smiles:
                    return True, mol1_count, 0, reference_smiles
    except Exception:  # noqa: BLE001
        pass

    mol2_count = 0
    try:
        source = _normalize_special_resonance_forms(mol2)
        source = _add_hs_without_sanitize(source)
        with suppress(Exception):
            Chem.Kekulize(source, clearAromaticFlags=True)
        for resonance_mol in ResonanceMolSupplier(
            source,
            maxStructs=max_resonance,
            flags=resonance_flags,
        ):
            if resonance_mol is None:
                continue
            mol2_count += 1
            with suppress(Exception):
                if _canon_smiles(resonance_mol, use_chirality) == reference_smiles:
                    return True, mol1_count, mol2_count, reference_smiles
    except Exception:  # noqa: BLE001
        pass

    return False, mol1_count, mol2_count, None


def check_equivalence(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    use_chirality: bool = True,
    max_resonance: int = 50,
    resonance_flags: Chem.ResonanceFlags = Chem.ResonanceFlags.UNCONSTRAINED_CATIONS
    | Chem.ResonanceFlags.UNCONSTRAINED_ANIONS,
) -> Tuple[bool, EquivalenceInfo]:
    standardized_1 = _standardize_metal_bonds(mol1)
    standardized_2 = _standardize_metal_bonds(mol2)

    organic_1 = _prepare_organic_mol(standardized_1)
    organic_2 = _prepare_organic_mol(standardized_2)

    checks = EquivalenceChecks(
        formal_charge=PropertyCheck(
            _total_formal_charge(organic_1),
            _total_formal_charge(organic_2),
            _total_formal_charge(organic_1) == _total_formal_charge(organic_2),
        ),
        radical_electrons=PropertyCheck(
            _total_radical_electrons(organic_1),
            _total_radical_electrons(organic_2),
            _total_radical_electrons(organic_1) == _total_radical_electrons(organic_2),
        ),
        num_atoms=PropertyCheck(
            organic_1.GetNumAtoms(),
            organic_2.GetNumAtoms(),
            organic_1.GetNumAtoms() == organic_2.GetNumAtoms(),
        ),
        heavy_atom_formula=PropertyCheck(
            _formula_key(organic_1, include_hydrogen=False),
            _formula_key(organic_2, include_hydrogen=False),
            _formula_key(organic_1, include_hydrogen=False) == _formula_key(organic_2, include_hydrogen=False),
        ),
        explicit_h_formula=PropertyCheck(
            _formula_key(organic_1, include_hydrogen=True),
            _formula_key(organic_2, include_hydrogen=True),
            _formula_key(organic_1, include_hydrogen=True) == _formula_key(organic_2, include_hydrogen=True),
        ),
    )
    info = EquivalenceInfo(checks=checks)

    if not checks.heavy_atom_formula.passed:
        info.reason = "Not equivalent: heavy-atom element counts differ."
        return False, info

    if _metal_signature_key(standardized_1) != _metal_signature_key(standardized_2):
        info.reason = "Not equivalent: metal valence assignment differs."
        return False, info

    if not checks.explicit_h_formula.passed:
        info.reason = "Not equivalent: explicit-hydrogen element counts differ."
        return False, info

    if not checks.num_atoms.passed:
        info.reason = "Not equivalent: explicit-hydrogen atom counts differ."
        return False, info

    try:
        smiles_1 = _canon_smiles(organic_1, use_chirality)
        smiles_2 = _canon_smiles(organic_2, use_chirality)
        info.canonical_smiles = CanonicalSmilesDetail(smiles_1, smiles_2, use_chirality)
        if smiles_1 == smiles_2:
            info.equivalent = True
            info.method = EquivalenceMethod.IDEAL
            info.reason = "Equivalent: canonical SMILES are identical after standardization."
            return True, info
        if use_chirality and _canon_smiles(organic_1, False) == _canon_smiles(organic_2, False):
            info.reason = "Not equivalent: stereochemistry differs."
            return False, info
    except Exception as exc:  # noqa: BLE001
        info.reason = f"Not equivalent: canonical SMILES comparison failed: {type(exc).__name__}: {exc}"

    inchi_key_1 = _inchi_key(organic_1)
    inchi_key_2 = _inchi_key(organic_2)
    if inchi_key_1 is not None and inchi_key_1 == inchi_key_2:
        info.equivalent = True
        info.method = EquivalenceMethod.INCHI_KEY
        info.reason = "Equivalent: full InChIKey matches after standardization."
        return True, info

    try:
        normalized_1 = _carbene_zwitterion_normalized_smiles(organic_1, use_chirality=use_chirality)
        normalized_2 = _carbene_zwitterion_normalized_smiles(organic_2, use_chirality=use_chirality)
        info.carbene_zwitterion = CarbeneZwitterionDetail(
            mol1_normalized=normalized_1,
            mol2_normalized=normalized_2,
        )
        if normalized_1 == normalized_2:
            info.equivalent = True
            info.method = EquivalenceMethod.CARBENE_ZWITTERION
            info.reason = "Equivalent: carbene/zwitterion normalization matched."
            return True, info
    except Exception:
        pass

    resonance_matched, mol1_count, mol2_count, hit_smiles = _resonance_match(
        organic_1,
        organic_2,
        use_chirality=use_chirality,
        max_resonance=max_resonance,
        resonance_flags=resonance_flags,
    )
    if resonance_matched:
        info.equivalent = True
        info.method = EquivalenceMethod.RESONANCE
        info.resonance = ResonanceDetail(
            max_resonance=max_resonance,
            resonance_flags=int(resonance_flags),
            mol1_resonance_count=mol1_count,
            mol2_resonance_count=mol2_count,
            hit_smiles=hit_smiles,
        )
        info.reason = "Equivalent: resonance normalization matched."
        return True, info

    info.reason = "Not equivalent: no standardized comparison path matched."
    return False, info
