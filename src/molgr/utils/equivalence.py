# pyright: reportCallIssue=false
"""
Author: TMJ
Date: 2026-02-22 20:18:39
LastEditors: TMJ
LastEditTime: 2026-04-28 22:58:17
Description: 用于判定对于同一个分子坐标输入下不同的分子图重建算法结果的一致性，只要和标准答案在共振异构级别上有一个一致，就认为是一致的。
"""

from __future__ import annotations

from collections import Counter
from contextlib import suppress
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import ResonanceMolSupplier, inchi, rdMolHash


pt = Chem.GetPeriodicTable()


# ================================
# Enum
# ================================


class EquivalenceMethod(str, Enum):
    IDEAL = "ideal"
    ISOMORPHIC = "isomorphic"
    TOPOLOGICAL_ISOMORPHIC = "topological_isomorphic"
    CARBENE_ZWITTERION = "carbene_zwitterion"
    METAL_CARBENE_VALENCE = "metal_carbene_valence"
    COORDINATION_STRIPPED = "coordination_stripped"
    RESONANCE = "resonance"
    INCHI_CONNECTIVITY = "inchi_connectivity"
    MOLHASH = "molhash"


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
class MetalCarbeneValenceDetail:
    mol1_normalized: str
    mol2_normalized: str
    mol1_transformed_count: int
    mol2_transformed_count: int


@dataclass(frozen=True)
class _MappedResonanceChange:
    mol1_search_atom_indices: frozenset[int]
    mol2_search_atom_indices: frozenset[int]
    changed_atom_count: int
    changed_bond_count: int


@dataclass(frozen=True)
class _MappedResonanceMatch:
    mol1_resonance_count: int
    mol2_resonance_count: int
    changed_atom_count: int
    changed_bond_count: int
    hit_smiles: Optional[str]


@dataclass
class CoordinationStrippedDetail:
    mol1_stripped: str
    mol2_stripped: str


@dataclass
class MolHashDetail:
    mol1_arthor_substructure_order: str
    mol2_arthor_substructure_order: str
    mol1_mesomer: str
    mol2_mesomer: str


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
    metal_carbene_valence: Optional[MetalCarbeneValenceDetail] = None
    coordination_stripped: Optional[CoordinationStrippedDetail] = None
    resonance: Optional[ResonanceDetail] = None
    molhash: Optional[MolHashDetail] = None


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
_METAL_CARBENE_VALENCE_TRANSFORMED_ATOMIC_NUMS = frozenset({6})
_METAL_CARBENE_VALENCE_BOND_TYPES = frozenset(
    {
        Chem.BondType.SINGLE,
        Chem.BondType.DOUBLE,
        Chem.BondType.DATIVE,
    }
)


def _canon_smiles(m: Chem.Mol, use_chirality: bool) -> str:
    if Chem.SanitizeMol(m) != Chem.SanitizeFlags.SANITIZE_NONE:
        raise ValueError("Molecule is not sanitized.")
    return Chem.CanonSmiles(Chem.MolToSmiles(m, canonical=True, isomericSmiles=use_chirality))


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


def _has_metal_atom(mol: Chem.Mol) -> bool:
    return any(_is_metal_atom(atom) for atom in mol.GetAtoms())


def _safe_smiles(mol: Chem.Mol, *, use_chirality: bool) -> str:
    clone = Chem.Mol(mol)
    clone.UpdatePropertyCache(strict=False)
    with suppress(Exception):
        Chem.SanitizeMol(clone)
    return Chem.CanonSmiles(Chem.MolToSmiles(clone, canonical=True, isomericSmiles=use_chirality))


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
    counts = Counter(
        mol.GetAtomWithIdx(atom_idx).GetSymbol()
        for atom_idx in range(mol.GetNumAtoms())
        if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() != 1
    )
    if not counts:
        return ""

    ordered_symbols: list[str] = []
    if "C" in counts:
        ordered_symbols.append("C")
    ordered_symbols.extend(symbol for symbol in sorted(counts) if symbol != "C")
    return "".join(
        f"{symbol}{counts[symbol] if counts[symbol] != 1 else ''}" for symbol in ordered_symbols
    )


def _total_formal_charge(mol: Chem.Mol) -> int:
    return sum(
        mol.GetAtomWithIdx(atom_idx).GetFormalCharge() for atom_idx in range(mol.GetNumAtoms())
    )


def _total_radical_electrons(mol: Chem.Mol) -> int:
    return sum(
        mol.GetAtomWithIdx(atom_idx).GetNumRadicalElectrons()
        for atom_idx in range(mol.GetNumAtoms())
    )


def _strip_metal_coordination_bonds(mol: Chem.Mol) -> tuple[Chem.Mol, bool]:
    rw_mol = Chem.RWMol(Chem.Mol(mol))
    bonds_to_remove: list[tuple[int, int]] = []
    for bond_idx in range(rw_mol.GetNumBonds()):
        bond = rw_mol.GetBondWithIdx(bond_idx)
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


def _simplified_topology_mol(mol: Chem.Mol) -> Chem.Mol:
    simplified_rw = Chem.RWMol(Chem.Mol(mol))
    for atom in simplified_rw.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNumRadicalElectrons(0)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
        atom.SetIsAromatic(False)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in simplified_rw.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetBondDir(Chem.BondDir.NONE)
        bond.SetStereo(Chem.BondStereo.STEREONONE)
        bond.SetIsAromatic(False)
    simplified = simplified_rw.GetMol()
    simplified.UpdatePropertyCache(strict=False)
    return simplified


def _full_simplified_topology_matches(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
) -> tuple[tuple[int, ...], ...]:
    if mol1.GetNumAtoms() != mol2.GetNumAtoms():
        return ()
    simplified_mol1 = _simplified_topology_mol(mol1)
    simplified_mol2 = _simplified_topology_mol(mol2)

    try:
        topology_smiles_1 = Chem.MolToSmiles(
            simplified_mol1,
            canonical=True,
            isomericSmiles=False,
        )
        topology_smiles_2 = Chem.MolToSmiles(
            simplified_mol2,
            canonical=True,
            isomericSmiles=False,
        )
    except Exception:  # noqa: BLE001
        return ()
    if topology_smiles_1 != topology_smiles_2:
        return ()

    ranks1 = tuple(int(rank) for rank in Chem.CanonicalRankAtoms(simplified_mol1, breakTies=True))
    ranks2 = tuple(int(rank) for rank in Chem.CanonicalRankAtoms(simplified_mol2, breakTies=True))
    if len(set(ranks1)) != mol1.GetNumAtoms() or len(set(ranks2)) != mol2.GetNumAtoms():
        return ()
    if set(ranks1) != set(ranks2):
        return ()

    mol2_idx_by_rank = {rank: atom_idx for atom_idx, rank in enumerate(ranks2)}
    match = tuple(mol2_idx_by_rank[rank] for rank in ranks1)
    if not _simplified_topology_mapping_is_valid(simplified_mol1, simplified_mol2, match):
        return ()
    return (match,)


def _simplified_topology_mapping_is_valid(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    mol1_to_mol2: tuple[int, ...],
) -> bool:
    if len(mol1_to_mol2) != mol1.GetNumAtoms():
        return False
    if len(set(mol1_to_mol2)) != mol2.GetNumAtoms():
        return False
    if mol1.GetNumBonds() != mol2.GetNumBonds():
        return False

    for atom1 in mol1.GetAtoms():
        atom1_idx = int(atom1.GetIdx())
        atom2 = mol2.GetAtomWithIdx(mol1_to_mol2[atom1_idx])
        if atom1.GetAtomicNum() != atom2.GetAtomicNum():
            return False

    for bond1 in mol1.GetBonds():
        begin2_idx = mol1_to_mol2[int(bond1.GetBeginAtomIdx())]
        end2_idx = mol1_to_mol2[int(bond1.GetEndAtomIdx())]
        if mol2.GetBondBetweenAtoms(begin2_idx, end2_idx) is None:
            return False

    return True


def _simplified_topology_smiles(mol: Chem.Mol) -> str | None:
    simplified_mol = _simplified_topology_mol(mol)
    try:
        return Chem.MolToSmiles(
            simplified_mol,
            canonical=True,
            isomericSmiles=False,
        )
    except Exception:  # noqa: BLE001
        return None


def _has_same_simplified_topology(mol1: Chem.Mol, mol2: Chem.Mol) -> bool:
    if mol1.GetNumAtoms() != mol2.GetNumAtoms():
        return False
    if mol1.GetNumBonds() != mol2.GetNumBonds():
        return False
    topology_smiles_1 = _simplified_topology_smiles(mol1)
    if topology_smiles_1 is None:
        return False
    return topology_smiles_1 == _simplified_topology_smiles(mol2)


def _fragment_formula_topology_counts(
    mol: Chem.Mol,
) -> tuple[Counter[str], dict[str, Counter[str]]] | None:
    formula_counts: Counter[str] = Counter()
    topology_counts_by_formula: dict[str, Counter[str]] = {}

    try:
        fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    except Exception:  # noqa: BLE001
        return None

    for fragment in fragments:
        fragment.UpdatePropertyCache(strict=False)
        with suppress(Exception):
            Chem.SanitizeMol(fragment)
        formula = _formula_key_without_hydrogen(fragment)
        topology_smiles = _simplified_topology_smiles(fragment)
        if topology_smiles is None:
            return None
        formula_counts[formula] += 1
        topology_counts_by_formula.setdefault(formula, Counter())[topology_smiles] += 1

    return formula_counts, topology_counts_by_formula


def _counter_summary(counter: Counter[str]) -> str:
    return ".".join(f"{key}:{counter[key]}" for key in sorted(counter))


def _atom_resonance_signature(atom: Chem.Atom) -> tuple[int, int, bool]:
    return (
        atom.GetFormalCharge(),
        atom.GetNumRadicalElectrons(),
        atom.GetIsAromatic(),
    )


def _bond_resonance_signature(bond: Chem.Bond) -> tuple[Chem.BondType, bool]:
    return bond.GetBondType(), bond.GetIsAromatic()


def _is_resonance_domain_bond(mol: Chem.Mol, begin_atom_idx: int, end_atom_idx: int) -> bool:
    bond = mol.GetBondBetweenAtoms(begin_atom_idx, end_atom_idx)
    if bond is None:
        return False
    return (
        bond.GetIsConjugated()
        or bond.GetIsAromatic()
        or bond.GetBondType()
        in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE, Chem.BondType.AROMATIC)
    )


def _mapped_resonance_change(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    mol1_to_mol2: tuple[int, ...],
) -> _MappedResonanceChange | None:
    changed_atoms: set[int] = set()
    changed_bonds: set[tuple[int, int]] = set()
    domain_edges: set[tuple[int, int]] = set()
    domain_adjacency: dict[int, set[int]] = {idx: set() for idx in range(mol1.GetNumAtoms())}

    for atom1 in mol1.GetAtoms():
        atom1_idx = int(atom1.GetIdx())
        atom2 = mol2.GetAtomWithIdx(mol1_to_mol2[atom1_idx])
        if _atom_resonance_signature(atom1) != _atom_resonance_signature(atom2):
            changed_atoms.add(atom1_idx)

    for bond1 in mol1.GetBonds():
        begin1_idx = int(bond1.GetBeginAtomIdx())
        end1_idx = int(bond1.GetEndAtomIdx())
        begin2_idx = mol1_to_mol2[begin1_idx]
        end2_idx = mol1_to_mol2[end1_idx]
        bond2 = mol2.GetBondBetweenAtoms(begin2_idx, end2_idx)
        if bond2 is None:
            return None

        edge = (begin1_idx, end1_idx) if begin1_idx <= end1_idx else (end1_idx, begin1_idx)
        is_domain_edge = _is_resonance_domain_bond(
            mol1,
            begin1_idx,
            end1_idx,
        ) or _is_resonance_domain_bond(mol2, begin2_idx, end2_idx)
        if is_domain_edge:
            domain_edges.add(edge)
            domain_adjacency[begin1_idx].add(end1_idx)
            domain_adjacency[end1_idx].add(begin1_idx)

        if _bond_resonance_signature(bond1) != _bond_resonance_signature(bond2):
            changed_bonds.add(edge)
            changed_atoms.update(edge)

    if not changed_atoms and not changed_bonds:
        return None
    if any(edge not in domain_edges for edge in changed_bonds):
        return None

    start_atom_idx = next(iter(changed_atoms))
    seen = {start_atom_idx}
    stack = [start_atom_idx]
    while stack:
        atom_idx = stack.pop()
        for neighbor_idx in domain_adjacency[atom_idx]:
            if neighbor_idx in seen:
                continue
            seen.add(neighbor_idx)
            stack.append(neighbor_idx)

    if not changed_atoms <= seen:
        return None

    search_atom_indices = set(changed_atoms)
    for atom_idx in tuple(changed_atoms):
        search_atom_indices.update(
            int(neighbor.GetIdx()) for neighbor in mol1.GetAtomWithIdx(atom_idx).GetNeighbors()
        )
    return _MappedResonanceChange(
        mol1_search_atom_indices=frozenset(search_atom_indices),
        mol2_search_atom_indices=frozenset(mol1_to_mol2[idx] for idx in search_atom_indices),
        changed_atom_count=len(changed_atoms),
        changed_bond_count=len(changed_bonds),
    )


def _fragment_canonical_smiles(
    mol: Chem.Mol,
    atom_indices: frozenset[int],
    *,
    use_chirality: bool,
) -> str | None:
    if not atom_indices:
        return None
    smiles = Chem.MolFragmentToSmiles(
        mol,
        atomsToUse=sorted(atom_indices),
        canonical=True,
        isomericSmiles=use_chirality,
    )
    fragment_mol = Chem.MolFromSmiles(smiles)
    if fragment_mol is None:
        return None
    return _canon_smiles(Chem.AddHs(fragment_mol), use_chirality)


def _directed_resonance_fragment_match(
    source_mol: Chem.Mol,
    target_mol: Chem.Mol,
    source_atom_indices: frozenset[int],
    target_atom_indices: frozenset[int],
    *,
    use_chirality: bool,
    max_resonance: int,
    resonance_flags: Chem.ResonanceFlags,
) -> tuple[bool, int, str | None]:
    source_fragment_smiles = Chem.MolFragmentToSmiles(
        source_mol,
        atomsToUse=sorted(source_atom_indices),
        canonical=True,
        isomericSmiles=use_chirality,
    )
    source_fragment_mol = Chem.MolFromSmiles(source_fragment_smiles)
    if source_fragment_mol is None:
        return False, 0, None

    target_smiles = _fragment_canonical_smiles(
        target_mol,
        target_atom_indices,
        use_chirality=use_chirality,
    )
    if target_smiles is None:
        return False, 0, None

    target_charge = _total_formal_charge(source_fragment_mol)
    source_fragment_mol = Chem.AddHs(source_fragment_mol)
    with suppress(Exception):
        Chem.Kekulize(source_fragment_mol, clearAromaticFlags=True)

    resonance_count = 0
    try:
        resonance_iter = _iter_resonance_structures(
            source_fragment_mol,
            max_resonance=max_resonance,
            resonance_flags=resonance_flags,
        )
        for resonance_mol in resonance_iter:
            if resonance_mol is None:
                continue
            resonance_count += 1
            if _total_formal_charge(resonance_mol) != target_charge:
                continue
            with suppress(Exception):
                resonance_smiles = _canon_smiles(resonance_mol, use_chirality)
                if resonance_smiles == target_smiles:
                    return True, resonance_count, resonance_smiles
    except Exception:  # noqa: BLE001
        pass

    return False, resonance_count, None


def _mapped_resonance_fragment_match(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    *,
    use_chirality: bool,
    max_resonance: int,
    resonance_flags: Chem.ResonanceFlags,
    coordination_bonds_already_stripped: bool,
) -> _MappedResonanceMatch | None:
    if coordination_bonds_already_stripped:
        stripped_m1 = Chem.Mol(mol1)
        stripped_m2 = Chem.Mol(mol2)
    else:
        stripped_m1, _ = _strip_metal_coordination_bonds(mol1)
        stripped_m2, _ = _strip_metal_coordination_bonds(mol2)

    if _total_formal_charge(stripped_m1) != _total_formal_charge(stripped_m2):
        return None
    if _total_radical_electrons(stripped_m1) != _total_radical_electrons(stripped_m2):
        return None
    if _formula_key_without_hydrogen(stripped_m1) != _formula_key_without_hydrogen(stripped_m2):
        return None
    if stripped_m1.GetNumAtoms() != stripped_m2.GetNumAtoms():
        return None

    with suppress(Exception):
        if _safe_smiles(stripped_m1, use_chirality=use_chirality) == _safe_smiles(
            stripped_m2,
            use_chirality=use_chirality,
        ):
            return None

    resonance_change = None
    for match in _full_simplified_topology_matches(stripped_m1, stripped_m2):
        resonance_change = _mapped_resonance_change(
            stripped_m1,
            stripped_m2,
            match,
        )
        if resonance_change is not None:
            break
    if resonance_change is None:
        return None

    mol1_to_mol2_matched, mol1_count, mol1_hit_smiles = _directed_resonance_fragment_match(
        stripped_m1,
        stripped_m2,
        resonance_change.mol1_search_atom_indices,
        resonance_change.mol2_search_atom_indices,
        use_chirality=use_chirality,
        max_resonance=max_resonance,
        resonance_flags=resonance_flags,
    )
    mol2_to_mol1_matched, mol2_count, mol2_hit_smiles = _directed_resonance_fragment_match(
        stripped_m2,
        stripped_m1,
        resonance_change.mol2_search_atom_indices,
        resonance_change.mol1_search_atom_indices,
        use_chirality=use_chirality,
        max_resonance=max_resonance,
        resonance_flags=resonance_flags,
    )
    if not mol1_to_mol2_matched and not mol2_to_mol1_matched:
        return None

    return _MappedResonanceMatch(
        mol1_resonance_count=mol1_count,
        mol2_resonance_count=mol2_count,
        changed_atom_count=resonance_change.changed_atom_count,
        changed_bond_count=resonance_change.changed_bond_count,
        hit_smiles=mol1_hit_smiles if mol1_to_mol2_matched else mol2_hit_smiles,
    )


def _normalize_carbene_zwitterion_forms(mol: Chem.Mol) -> Chem.Mol:
    normalized = Chem.Mol(mol)
    rw_mol = Chem.RWMol(normalized)
    Chem.Kekulize(rw_mol)

    smarts = Chem.MolFromSmarts("[*-]=[*+]")
    for match in rw_mol.GetSubstructMatches(smarts):
        begin_idx, end_idx = (int(atom_idx) for atom_idx in match)
        bond = rw_mol.GetBondBetweenAtoms(begin_idx, end_idx)
        atom_begin = rw_mol.GetAtomWithIdx(begin_idx)
        atom_end = rw_mol.GetAtomWithIdx(end_idx)
        if (
            pt.GetDefaultValence(atom_begin.GetAtomicNum()) - atom_begin.GetTotalValence() == 1
        ) and (pt.GetDefaultValence(atom_end.GetAtomicNum()) - atom_end.GetTotalValence() == -1):
            bond.SetBondType(Chem.BondType.SINGLE)
            atom_begin.SetFormalCharge(atom_begin.GetFormalCharge() + 1)
            atom_end.SetFormalCharge(atom_end.GetFormalCharge() - 1)
            atom_begin.SetNumRadicalElectrons(2)
            atom_end.SetNumRadicalElectrons(0)

    normalized_mol = rw_mol.GetMol()
    normalized_mol.UpdatePropertyCache(strict=False)
    Chem.Cleanup(normalized_mol)
    for atom_idx in range(normalized_mol.GetNumAtoms()):
        normalized_mol.GetAtomWithIdx(atom_idx).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        normalized_mol.GetAtomWithIdx(atom_idx).SetIsAromatic(False)
    return normalized_mol


def _carbene_zwitterion_normalized_smiles(
    mol: Chem.Mol,
    *,
    use_chirality: bool,
) -> str:
    normalized = _normalize_carbene_zwitterion_forms(mol)
    return Chem.CanonSmiles(
        Chem.MolToSmiles(
            normalized,
            canonical=True,
            isomericSmiles=use_chirality,
        )
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


def _bond_order_for_valence(bond: Chem.Bond) -> float:
    if bond.GetBondType() == Chem.BondType.DATIVE:
        return 1.0
    return float(bond.GetBondTypeAsDouble())


def _neutral_carbenic_carbon_radical_count(atom: Chem.Atom) -> int:
    valence_order = float(atom.GetTotalNumHs())
    valence_order += sum(_bond_order_for_valence(bond) for bond in atom.GetBonds())
    return max(0, int(round(4.0 - valence_order)))


def _is_metal_carbene_valence_form(carbon_atom: Chem.Atom, bond: Chem.Bond) -> bool:
    if carbon_atom.GetAtomicNum() not in _METAL_CARBENE_VALENCE_TRANSFORMED_ATOMIC_NUMS:
        return False
    bond_type = bond.GetBondType()
    if bond_type not in _METAL_CARBENE_VALENCE_BOND_TYPES:
        return False
    if carbon_atom.GetFormalCharge() == -2:
        return True
    if carbon_atom.GetFormalCharge() != 0:
        return False
    return bond_type == Chem.BondType.DOUBLE or (
        bond_type == Chem.BondType.DATIVE and carbon_atom.GetNumRadicalElectrons() > 0
    )


def _has_metal_carbene_valence_candidate(mol: Chem.Mol) -> bool:
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if _is_metal_atom(begin_atom) == _is_metal_atom(end_atom):
            continue
        organic_atom = end_atom if _is_metal_atom(begin_atom) else begin_atom
        if (
            organic_atom.GetAtomicNum() in _METAL_CARBENE_VALENCE_TRANSFORMED_ATOMIC_NUMS
            and organic_atom.GetFormalCharge() == -2
            and bond.GetBondType() in _METAL_CARBENE_VALENCE_BOND_TYPES
        ):
            return True
    return False


def _fold_explicit_metal_hydrogens(mol: Chem.Mol) -> Chem.Mol:
    normalized = Chem.Mol(mol)
    rw_mol = Chem.RWMol(normalized)
    metal_hydrogen_counts: Counter[int] = Counter()
    hydrogens_to_remove: list[int] = []

    for atom in rw_mol.GetAtoms():
        if (
            atom.GetAtomicNum() != 1
            or atom.GetDegree() != 1
            or atom.GetFormalCharge() != 0
            or atom.GetNumRadicalElectrons() != 0
            or atom.GetIsotope() != 0
        ):
            continue
        bond = atom.GetBonds()[0]
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        neighbor = atom.GetNeighbors()[0]
        if not _is_metal_atom(neighbor):
            continue
        metal_hydrogen_counts[int(neighbor.GetIdx())] += 1
        hydrogens_to_remove.append(int(atom.GetIdx()))

    if not hydrogens_to_remove:
        return normalized

    for metal_idx, hydrogen_count in metal_hydrogen_counts.items():
        metal_atom = rw_mol.GetAtomWithIdx(metal_idx)
        metal_atom.SetNumExplicitHs(metal_atom.GetNumExplicitHs() + hydrogen_count)
        metal_atom.SetNoImplicit(True)

    for atom_idx in sorted(hydrogens_to_remove, reverse=True):
        rw_mol.RemoveAtom(atom_idx)

    folded = rw_mol.GetMol()
    folded.UpdatePropertyCache(strict=False)
    with suppress(Exception):
        Chem.SanitizeMol(folded)
    return folded


def _normalize_metal_carbene_valence_forms(
    mol: Chem.Mol,
) -> tuple[Chem.Mol, int]:
    normalized = Chem.Mol(mol)
    rw_mol = Chem.RWMol(normalized)
    transformed_count = 0
    normalized_carbon_indices: set[int] = set()

    for bond in list(rw_mol.GetBonds()):
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if _is_metal_atom(begin_atom) == _is_metal_atom(end_atom):
            continue

        metal_atom = begin_atom if _is_metal_atom(begin_atom) else end_atom
        organic_atom = end_atom if metal_atom.GetIdx() == begin_atom.GetIdx() else begin_atom
        if not _is_metal_carbene_valence_form(organic_atom, bond):
            continue

        if organic_atom.GetFormalCharge() == -2:
            metal_atom.SetFormalCharge(metal_atom.GetFormalCharge() - 2)
            transformed_count += 1

        bond.SetBondType(Chem.BondType.DOUBLE)
        bond.SetIsAromatic(False)
        organic_atom.SetFormalCharge(0)
        organic_atom.SetIsAromatic(False)
        organic_atom.SetNoImplicit(True)
        normalized_carbon_indices.add(int(organic_atom.GetIdx()))

    rw_mol.UpdatePropertyCache(strict=False)
    for atom_idx in normalized_carbon_indices:
        carbon_atom = rw_mol.GetAtomWithIdx(atom_idx)
        carbon_atom.SetNumRadicalElectrons(_neutral_carbenic_carbon_radical_count(carbon_atom))

    normalized_mol = rw_mol.GetMol()
    normalized_mol.UpdatePropertyCache(strict=False)
    with suppress(Exception):
        Chem.SanitizeMol(normalized_mol)
    return normalized_mol, transformed_count


def _metal_carbene_valence_normalized_smiles(
    mol: Chem.Mol,
    *,
    use_chirality: bool,
) -> tuple[str, int]:
    mol = _fold_explicit_metal_hydrogens(mol)
    normalized, transformed_count = _normalize_metal_carbene_valence_forms(mol)
    normalized = _normalize_carbene_zwitterion_forms(normalized)
    normalized = _fold_explicit_metal_hydrogens(normalized)
    normalized, second_transformed_count = _normalize_metal_carbene_valence_forms(normalized)
    transformed_count += second_transformed_count
    return (
        Chem.CanonSmiles(
            Chem.MolToSmiles(
                normalized,
                canonical=True,
                isomericSmiles=use_chirality,
            )
        ),
        transformed_count,
    )


def _check_metal_carbene_valence_equivalence(
    info: EquivalenceInfo,
    m1: Chem.Mol,
    m2: Chem.Mol,
    *,
    use_chirality: bool,
    phase_errors: list[str],
) -> Tuple[bool, EquivalenceInfo]:
    try:
        normalized_smiles_1, transformed_count_1 = _metal_carbene_valence_normalized_smiles(
            m1,
            use_chirality=use_chirality,
        )
        normalized_smiles_2, transformed_count_2 = _metal_carbene_valence_normalized_smiles(
            m2,
            use_chirality=use_chirality,
        )
        info.metal_carbene_valence = MetalCarbeneValenceDetail(
            mol1_normalized=normalized_smiles_1,
            mol2_normalized=normalized_smiles_2,
            mol1_transformed_count=transformed_count_1,
            mol2_transformed_count=transformed_count_2,
        )
        if (
            normalized_smiles_1 == normalized_smiles_2
            and transformed_count_1 != transformed_count_2
            and (transformed_count_1 > 0 or transformed_count_2 > 0)
        ):
            info.equivalent = True
            info.method = EquivalenceMethod.METAL_CARBENE_VALENCE
            info.reason = (
                "Equivalent: metal-carbene valence normalization matched a subjective "
                "M=C versus C2-/metal-valence assignment."
            )
            return True, info
    except Exception as exc:  # noqa: BLE001
        phase_errors.append(f"metal-carbene valence check failed: {type(exc).__name__}: {exc}")

    return False, info


def _special_resonance_normalized_smiles(
    mol: Chem.Mol,
    *,
    use_chirality: bool,
) -> str:
    normalized = _normalize_special_resonance_forms(mol)
    return Chem.CanonSmiles(
        Chem.MolToSmiles(
            normalized,
            canonical=True,
            isomericSmiles=use_chirality,
        )
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


def _check_inchi_connectivity_equivalence(
    info: EquivalenceInfo,
    m1: Chem.Mol,
    m2: Chem.Mol,
    *,
    reason: str,
) -> Tuple[bool, EquivalenceInfo]:
    inchi_key1 = _inchi_connectivity_key(m1)
    inchi_key2 = _inchi_connectivity_key(m2)

    if inchi_key1 is not None and inchi_key1 == inchi_key2:
        info.equivalent = True
        info.method = EquivalenceMethod.INCHI_CONNECTIVITY
        info.reason = reason
        return True, info

    return False, info


def _molhash(mol: Chem.Mol, hash_function: rdMolHash.HashFunction) -> str | None:
    try:
        return rdMolHash.MolHash(mol, hash_function)
    except Exception:  # noqa: BLE001
        return None


def _check_molhash_equivalence(
    info: EquivalenceInfo,
    m1: Chem.Mol,
    m2: Chem.Mol,
) -> Tuple[bool, EquivalenceInfo]:
    if _has_metal_atom(m1) or _has_metal_atom(m2):
        return False, info
    if not _has_same_simplified_topology(m1, m2):
        return False, info

    arthor_1 = _molhash(m1, rdMolHash.HashFunction.ArthorSubstructureOrder)
    arthor_2 = _molhash(m2, rdMolHash.HashFunction.ArthorSubstructureOrder)
    mesomer_1 = _molhash(m1, rdMolHash.HashFunction.Mesomer)
    mesomer_2 = _molhash(m2, rdMolHash.HashFunction.Mesomer)

    if arthor_1 is None or arthor_2 is None or mesomer_1 is None or mesomer_2 is None:
        return False, info

    info.molhash = MolHashDetail(
        mol1_arthor_substructure_order=arthor_1,
        mol2_arthor_substructure_order=arthor_2,
        mol1_mesomer=mesomer_1,
        mol2_mesomer=mesomer_2,
    )

    if arthor_1 == arthor_2 and mesomer_1 == mesomer_2:
        info.equivalent = True
        info.method = EquivalenceMethod.MOLHASH
        info.reason = (
            "Equivalent: RDKit MolHash ArthorSubstructureOrder and Mesomer hashes match "
            "for molecules with the same simplified topology."
        )
        return True, info

    return False, info


def _check_mapped_resonance_equivalence(
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
    try:
        resonance_match = _mapped_resonance_fragment_match(
            m1,
            m2,
            use_chirality=use_chirality,
            max_resonance=max_resonance,
            resonance_flags=resonance_flags,
            coordination_bonds_already_stripped=coordination_bonds_already_stripped,
        )
        if resonance_match is None:
            return False, info

        info.resonance = ResonanceDetail(
            max_resonance=max_resonance,
            resonance_flags=int(resonance_flags),
            mol2_resonance_count=resonance_match.mol2_resonance_count,
            hit_smiles=resonance_match.hit_smiles,
        )

        info.equivalent = True
        info.method = EquivalenceMethod.RESONANCE
        info.reason = "Equivalent: a resonance move on the mapped difference fragment matched."
        return True, info
    except Exception as exc:  # noqa: BLE001
        phase_errors.append(f"mapped resonance check failed: {type(exc).__name__}: {exc}")

    return False, info


def _check_stripped_fragment_topology_equivalence(
    info: EquivalenceInfo,
    m1: Chem.Mol,
    m2: Chem.Mol,
    *,
    use_chirality: bool,
    phase_errors: list[str],
    coordination_bonds_already_stripped: bool,
) -> Tuple[bool | None, EquivalenceInfo]:
    if coordination_bonds_already_stripped:
        return None, info

    if not (_has_metal_atom(m1) or _has_metal_atom(m2)):
        return None, info

    try:
        stripped_m1, _ = _strip_metal_coordination_bonds(m1)
        stripped_m2, _ = _strip_metal_coordination_bonds(m2)

        info.coordination_stripped = CoordinationStrippedDetail(
            mol1_stripped=_safe_smiles(stripped_m1, use_chirality=use_chirality),
            mol2_stripped=_safe_smiles(stripped_m2, use_chirality=use_chirality),
        )

        fragment_counts_1 = _fragment_formula_topology_counts(stripped_m1)
        fragment_counts_2 = _fragment_formula_topology_counts(stripped_m2)
        if fragment_counts_1 is None or fragment_counts_2 is None:
            return None, info

        formula_counts_1, topology_counts_by_formula_1 = fragment_counts_1
        formula_counts_2, topology_counts_by_formula_2 = fragment_counts_2
        if formula_counts_1 != formula_counts_2:
            info.reason = (
                "Not equivalent after stripping coordination bonds: fragment element "
                f"composition differs ({_counter_summary(formula_counts_1)} vs "
                f"{_counter_summary(formula_counts_2)})."
            )
            return False, info

        if topology_counts_by_formula_1 != topology_counts_by_formula_2:
            info.reason = (
                "Not equivalent after stripping coordination bonds: matched fragment "
                "simplified topologies differ."
            )
            return False, info

        info.equivalent = True
        info.method = EquivalenceMethod.COORDINATION_STRIPPED
        info.reason = (
            "Equivalent: after stripping metal-ligand coordination bonds, fragments "
            "matched by element composition and simplified topology."
        )
        return True, info
    except Exception as exc:  # noqa: BLE001
        phase_errors.append(
            f"coordination-stripped fragment topology check failed: {type(exc).__name__}: {exc}"
        )

    return None, info


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
    if not use_chirality:
        m1 = _fold_explicit_metal_hydrogens(m1)
        m2 = _fold_explicit_metal_hydrogens(m2)

    fc1 = _total_formal_charge(m1)
    fc2 = _total_formal_charge(m2)

    rad1 = _total_radical_electrons(m1)
    rad2 = _total_radical_electrons(m2)

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

    if checks.formula.passed and (
        _has_metal_carbene_valence_candidate(m1) or _has_metal_carbene_valence_candidate(m2)
    ):
        metal_carbene_valence_equivalent, info = _check_metal_carbene_valence_equivalence(
            info,
            m1,
            m2,
            use_chirality=use_chirality,
            phase_errors=phase_errors,
        )
        if metal_carbene_valence_equivalent:
            return True, info

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

        fragment_topology_decision, info = _check_stripped_fragment_topology_equivalence(
            info,
            m1,
            m2,
            use_chirality=use_chirality,
            phase_errors=phase_errors,
            coordination_bonds_already_stripped=_coordination_bonds_already_stripped,
        )
        if fragment_topology_decision is not None:
            return fragment_topology_decision, info

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

    fragment_topology_decision, info = _check_stripped_fragment_topology_equivalence(
        info,
        m1,
        m2,
        use_chirality=use_chirality,
        phase_errors=phase_errors,
        coordination_bonds_already_stripped=_coordination_bonds_already_stripped,
    )
    if fragment_topology_decision is not None:
        return fragment_topology_decision, info

    try:
        if (
            not _coordination_bonds_already_stripped
            and (_has_metal_atom(m1) or _has_metal_atom(m2))
            and _has_same_simplified_topology(m1, m2)
        ):
            info.equivalent = True
            info.method = EquivalenceMethod.TOPOLOGICAL_ISOMORPHIC
            info.reason = (
                "Equivalent: metal-complex simplified topology matched after ignoring bond order, "
                "formal charge placement, radical labels, aromaticity, and chirality."
            )
            return True, info
    except Exception as exc:  # noqa: BLE001
        phase_errors.append(f"topological isomorphic check failed: {type(exc).__name__}: {exc}")

    # 3) Resonance equivalence
    mapped_resonance_equivalent, info = _check_mapped_resonance_equivalence(
        info,
        original_m1,
        original_m2,
        use_chirality=use_chirality,
        max_resonance=max_resonance,
        resonance_flags=resonance_flags,
        phase_errors=phase_errors,
        coordination_bonds_already_stripped=_coordination_bonds_already_stripped,
    )
    if mapped_resonance_equivalent:
        return True, info

    info.resonance = ResonanceDetail(
        max_resonance=max_resonance,
        resonance_flags=int(resonance_flags),
        mol2_resonance_count=0,
    )

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

    inchi_equivalent, info = _check_inchi_connectivity_equivalence(
        info,
        m1,
        m2,
        reason=(
            "Equivalent: InChI connectivity key matches despite different canonical SMILES, "
            "isomorphism, or resonance representation."
        ),
    )
    if inchi_equivalent:
        return True, info

    molhash_equivalent, info = _check_molhash_equivalence(info, m1, m2)
    if molhash_equivalent:
        return True, info

    info.reason = "Not equivalent: none of ideal, isomorphic, or resonance checks matched."
    return False, info
