# pyright: reportCallIssue=false
from __future__ import annotations

from collections import Counter
from contextlib import suppress
from dataclasses import dataclass, field
from enum import Enum
from functools import lru_cache
from typing import Optional, Tuple

from rdkit import Chem, rdBase
from rdkit.Chem import ResonanceMolSupplier, inchi


pt = Chem.GetPeriodicTable()


class EquivalenceMethod(str, Enum):
    IDEAL = "ideal"
    INCHI_KEY = "inchi_key"
    CARBENE_ZWITTERION = "carbene_zwitterion"
    RESONANCE = "resonance"


class EquivalenceDecision(str, Enum):
    EQUIVALENT = "equivalent"
    NOT_EQUIVALENT = "not_equivalent"
    INCONCLUSIVE = "inconclusive"


class EquivalenceRelation(str, Enum):
    NORMALIZED_GRAPH_IDENTITY = "normalized_graph_identity"
    IDENTIFIER_EQUIVALENCE = "identifier_equivalence"
    CARBENE_ZWITTERION_EQUIVALENCE = "carbene_zwitterion_equivalence"
    RESONANCE_EQUIVALENCE = "resonance_equivalence"
    NONE = "none"


class InvariantStatus(str, Enum):
    PASSED = "passed"
    FAILED = "failed"
    NOT_EVALUATED = "not_evaluated"


@dataclass
class InvariantResult:
    status: InvariantStatus
    mol1: object | None = None
    mol2: object | None = None
    reason: str = ""


@dataclass
class BoundedSearchMetadata:
    attempted: bool = False
    limit: int | None = None
    mol1_count: int = 0
    mol2_count: int = 0
    limit_reached: bool = False
    exhaustive: bool | None = None


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


UNKNOWN_ELECTRON_METADATA = "unknown"


@dataclass
class RawMoleculeDiagnostics:
    """Read-only diagnostics captured before equivalence normalization.

    The three MolGR-specific fields deliberately use ``"unknown"`` when the
    source molecule does not carry the corresponding property.  A missing
    property is not evidence for a zero-valued electronic state, especially
    for ordinary reference RDKit molecules.
    """

    formula: str
    hydrogen_count: int
    explicit_hydrogen_count: int
    formal_charge: int
    rdkit_radical_electrons: int
    metal_formal_state: tuple[tuple[int, str, int], ...]
    molgr_metal_unpaired_electrons: object = UNKNOWN_ELECTRON_METADATA
    active_lone_pair_properties: object = UNKNOWN_ELECTRON_METADATA
    unresolved_two_electron_center_properties: object = UNKNOWN_ELECTRON_METADATA


@dataclass
class EquivalenceChecks:
    formal_charge: PropertyCheck
    radical_electrons: PropertyCheck
    num_atoms: PropertyCheck
    heavy_atom_formula: PropertyCheck
    explicit_h_formula: PropertyCheck


@dataclass
class EquivalenceResult:
    decision: EquivalenceDecision = EquivalenceDecision.INCONCLUSIVE
    relation: EquivalenceRelation = EquivalenceRelation.NONE
    equivalent: bool = False
    method: Optional[EquivalenceMethod] = None
    reason: str = ""
    invariants: dict[str, InvariantResult] = field(default_factory=dict)
    contradictions: list[str] = field(default_factory=list)
    bounded_search: BoundedSearchMetadata | None = None
    checks: Optional[EquivalenceChecks] = None
    canonical_smiles: Optional[CanonicalSmilesDetail] = None
    carbene_zwitterion: Optional[CarbeneZwitterionDetail] = None
    resonance: Optional[ResonanceDetail] = None
    # ``mol1``/``mol2`` are the evaluator's positional inputs.  The candidate
    # and reference aliases make the same evidence convenient for reviewer
    # callers without changing the historical positional API.
    raw_mol1: Optional[RawMoleculeDiagnostics] = None
    raw_mol2: Optional[RawMoleculeDiagnostics] = None
    raw_diagnostics: dict[str, RawMoleculeDiagnostics] = field(default_factory=dict)
    normalization_electronic_effects: dict[str, dict[str, int]] = field(default_factory=dict)


# Compatibility name retained for callers that imported the former detail type.
EquivalenceInfo = EquivalenceResult


_NON_METAL_ATOMIC_NUMBERS = frozenset({1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53})

_SULFIMIDE_RESONANCE_PATTERN = Chem.MolFromSmarts("[N:1]=[S:2](=[O:3])([O-:4])[!#1:5]")
_CARBENE_ZWITTERION_PATTERN = Chem.MolFromSmarts("[*-]=[*+]")
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
        atom.GetSymbol() for atom in mol.GetAtoms() if include_hydrogen or atom.GetAtomicNum() != 1
    )
    if not counts:
        return ""
    ordered_symbols: list[str] = []
    if "C" in counts:
        ordered_symbols.append("C")
    if include_hydrogen and "H" in counts:
        ordered_symbols.append("H")
    ordered_symbols.extend(symbol for symbol in sorted(counts) if symbol not in ordered_symbols)
    return "".join(
        f"{symbol}{counts[symbol] if counts[symbol] != 1 else ''}" for symbol in ordered_symbols
    )


def _raw_hydrogen_count(mol: Chem.Mol) -> tuple[int, int]:
    """Return explicit and total hydrogen counts without normalizing ``mol``."""

    explicit = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    total = explicit
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        with suppress(Exception):
            # Do not count explicit hydrogen neighbors here: they are already
            # included in ``explicit``.  ``includeNeighbors=True`` would count
            # those neighbors a second time for molecules passed through
            # ``Chem.AddHs``.
            total += int(atom.GetTotalNumHs(includeNeighbors=False))
    return explicit, total


def _raw_formula_key(mol: Chem.Mol, hydrogen_count: int) -> str:
    counts = Counter(
        atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1
    )
    if hydrogen_count:
        counts["H"] = hydrogen_count
    ordered_symbols: list[str] = []
    if "C" in counts:
        ordered_symbols.append("C")
    if "H" in counts:
        ordered_symbols.append("H")
    ordered_symbols.extend(symbol for symbol in sorted(counts) if symbol not in ordered_symbols)
    return "".join(
        f"{symbol}{counts[symbol] if counts[symbol] != 1 else ''}" for symbol in ordered_symbols
    )


def _raw_atom_property_entries(
    mol: Chem.Mol,
    property_name: str,
    *,
    metals_only: bool = False,
    boolean: bool = False,
) -> object:
    """Collect only present atom properties, preserving missing metadata."""

    entries: list[tuple[int, object]] = []
    for atom in mol.GetAtoms():
        if metals_only and not _is_metal_atom(atom):
            continue
        if not atom.HasProp(property_name):
            continue
        try:
            value: object = atom.GetBoolProp(property_name) if boolean else atom.GetIntProp(property_name)
        except (RuntimeError, TypeError, ValueError):
            try:
                raw_value = atom.GetProp(property_name)
                value = raw_value.lower() in {"1", "true"} if boolean else int(raw_value)
            except (RuntimeError, TypeError, ValueError):
                value = UNKNOWN_ELECTRON_METADATA
        entries.append((int(atom.GetIdx()), value))
    return tuple(entries) if entries else UNKNOWN_ELECTRON_METADATA


def _raw_molecule_diagnostics(mol: Chem.Mol) -> RawMoleculeDiagnostics:
    """Capture source-state evidence before metal or octet normalization."""

    explicit_h, total_h = _raw_hydrogen_count(mol)
    metal_formal_state = tuple(
        sorted(
            (
                int(atom.GetIdx()),
                atom.GetSymbol(),
                int(atom.GetFormalCharge()),
            )
            for atom in mol.GetAtoms()
            if _is_metal_atom(atom)
        )
    )
    return RawMoleculeDiagnostics(
        formula=_raw_formula_key(mol, total_h),
        hydrogen_count=total_h,
        explicit_hydrogen_count=explicit_h,
        formal_charge=_total_formal_charge(mol),
        rdkit_radical_electrons=_total_radical_electrons(mol),
        metal_formal_state=metal_formal_state,
        molgr_metal_unpaired_electrons=_raw_atom_property_entries(
            mol,
            "MOLGR_METAL_UNPAIRED_ELECTRONS",
            metals_only=True,
        ),
        active_lone_pair_properties=_raw_atom_property_entries(
            mol,
            "MOLGR_LONE_PAIR_COUNT",
        ),
        unresolved_two_electron_center_properties=_raw_atom_property_entries(
            mol,
            "MOLGR_UNRESOLVED_TWO_ELECTRON_CENTER",
            boolean=True,
        ),
    )


def _safe_copy(mol: Chem.Mol) -> Chem.Mol:
    clone = Chem.Mol(mol)
    clone.UpdatePropertyCache(strict=False)
    _sanitize_if_possible(clone)
    return clone


def _sanitize_if_possible(mol: Chem.Mol) -> None:
    try:
        Chem.SanitizeMol(mol)
    except TimeoutError:
        raise
    except Exception:  # noqa: BLE001
        pass


def _kekulize_if_possible(mol: Chem.Mol, *, clear_aromatic_flags: bool = False) -> None:
    try:
        Chem.Kekulize(mol, clearAromaticFlags=clear_aromatic_flags)
    except TimeoutError:
        raise
    except Exception:  # noqa: BLE001
        pass


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
    if not use_chirality:
        Chem.RemoveStereochemistry(clone)
        clone.UpdatePropertyCache(strict=False)
    return Chem.MolToSmiles(clone, canonical=True, isomericSmiles=use_chirality)


def _canon_smiles(mol: Chem.Mol, use_chirality: bool) -> str:
    # MolToSmiles(canonical=True) already canonicalizes the molecular graph.
    return _safe_smiles(mol, use_chirality=use_chirality)


def _resonance_form_smiles(mol: Chem.Mol, use_chirality: bool) -> str:
    try:
        # ResonanceMolSupplier returns sanitized molecules. Avoid copying and
        # sanitizing every form in the hot loop, but retain the safe fallback.
        if not use_chirality:
            # ``isomericSmiles=False`` suppresses most stereo output, but an
            # explicit copy also removes directional bond/stereo state before
            # resonance forms are cached and compared.
            achiral = Chem.Mol(mol)
            Chem.RemoveStereochemistry(achiral)
            for bond in achiral.GetBonds():
                bond.SetBondDir(Chem.BondDir.NONE)
                bond.SetStereo(Chem.BondStereo.STEREONONE)
            achiral.UpdatePropertyCache(strict=False)
            mol = achiral
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=use_chirality)
    except TimeoutError:
        raise
    except Exception:  # noqa: BLE001
        return _canon_smiles(mol, use_chirality)


def _inchi_key(mol: Chem.Mol, *, use_chirality: bool = True) -> str | None:
    try:
        clone = _safe_copy(mol)
        if not use_chirality:
            Chem.RemoveStereochemistry(clone)
            clone.UpdatePropertyCache(strict=False)
        key = inchi.MolToInchiKey(clone)
    except TimeoutError:
        raise
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
    rw_mol = Chem.RWMol(_safe_copy(_clear_metal_radicals(_add_hs_without_sanitize(mol))))
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

        if bond_order == 2:
            # Return the metal-ligand double-bond electron pair to the ligand.
            # Once the metal is removed, this is the carbene-like two-electron
            # center represented by an explicit pair of radical electrons. A
            # carbene coordination bond does not change the metal valence.
            ligand_atom.SetNumRadicalElectrons(ligand_atom.GetNumRadicalElectrons() + 2)
        else:
            metal_atom.SetFormalCharge(metal_atom.GetFormalCharge() + bond_order)
            ligand_atom.SetFormalCharge(ligand_atom.GetFormalCharge() - bond_order)

        bond.SetBondType(Chem.BondType.DATIVE)
        bond.SetBondDir(Chem.BondDir.NONE)
        bond.SetIsAromatic(False)

    standardized = rw_mol.GetMol()
    standardized.UpdatePropertyCache(strict=False)
    _sanitize_if_possible(standardized)
    return standardized


def _remove_coordination_bonds(mol: Chem.Mol) -> Chem.Mol:
    rw_mol = Chem.RWMol(Chem.Mol(mol))
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
    _sanitize_if_possible(stripped)
    return stripped


def _remove_metal_atoms(mol: Chem.Mol) -> Chem.Mol:
    rw_mol = Chem.RWMol(_safe_copy(mol))
    metal_indices = [atom.GetIdx() for atom in rw_mol.GetAtoms() if _is_metal_atom(atom)]
    for atom_idx in sorted(metal_indices, reverse=True):
        rw_mol.RemoveAtom(int(atom_idx))
    stripped = rw_mol.GetMol()
    stripped.UpdatePropertyCache(strict=False)
    _sanitize_if_possible(stripped)
    return stripped


def _prepare_organic_mol(mol: Chem.Mol, *, already_standardized: bool = False) -> Chem.Mol:
    standardized = mol if already_standardized else _standardize_metal_bonds(mol)
    standardized = _remove_coordination_bonds(standardized)
    standardized = _remove_metal_atoms(standardized)
    standardized.UpdatePropertyCache(strict=False)
    _sanitize_if_possible(standardized)
    return standardized


# These elements are allowed to appear in tmQMg in hypervalent, multiple-bond
# forms.  Equivalence comparison uses their charge-separated octet form only;
# reconstruction output and stored reference graphs are left untouched.
_OCTET_NORMALIZED_ATOMIC_NUMS = frozenset({7, 8, 9, 15, 16, 17, 33, 34, 35, 53})
_OCTET_VALENCE_LIMITS = {
    9: 1,
    17: 1,
    35: 1,
    53: 1,
}


def _bond_type_for_order(order: int) -> Chem.BondType:
    return {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
        4: Chem.BondType.QUADRUPLE,
    }[order]


def _normalize_nonmetal_octet(mol: Chem.Mol) -> Chem.Mol:
    """Canonicalize supported hypervalent nonmetals to charge-separated octets."""

    rw_mol = Chem.RWMol(_safe_copy(mol))
    # A single deterministic pass is insufficient when reducing one bond
    # exposes another overvalent atom, so process until no excess remains.
    changed = True
    while changed:
        changed = False
        for atom in rw_mol.GetAtoms():
            atomic_num = int(atom.GetAtomicNum())
            if atomic_num not in _OCTET_NORMALIZED_ATOMIC_NUMS:
                continue
            limit = _OCTET_VALENCE_LIMITS.get(atomic_num, 4)
            bond_order_sum = sum(
                int(round(bond.GetBondTypeAsDouble()))
                for bond in atom.GetBonds()
                if bond.GetBondType() != Chem.BondType.DATIVE
            )
            if bond_order_sum <= limit:
                continue

            reducible = sorted(
                (
                    bond
                    for bond in atom.GetBonds()
                    if int(round(bond.GetBondTypeAsDouble())) > 1
                    and not bond.GetIsAromatic()
                    and bond.GetBondType() != Chem.BondType.DATIVE
                ),
                key=lambda bond: (
                    -int(round(bond.GetBondTypeAsDouble())),
                    int(bond.GetOtherAtomIdx(atom.GetIdx())),
                ),
            )
            if not reducible:
                continue
            bond = reducible[0]
            order = int(round(bond.GetBondTypeAsDouble()))
            other = rw_mol.GetAtomWithIdx(bond.GetOtherAtomIdx(atom.GetIdx()))
            bond.SetBondType(_bond_type_for_order(order - 1))
            atom.SetFormalCharge(atom.GetFormalCharge() + 1)
            other.SetFormalCharge(other.GetFormalCharge() - 1)
            changed = True
            break

    normalized = rw_mol.GetMol()
    normalized.UpdatePropertyCache(strict=False)
    return normalized


def _prepare_resonance_source(mol: Chem.Mol) -> Chem.Mol:
    source = _add_hs_without_sanitize(mol)
    _kekulize_if_possible(source, clear_aromatic_flags=True)
    return source


def _charge_shift_resonance_forms(mol: Chem.Mol) -> tuple[Chem.Mol, ...]:
    """Generate the generic ``[*-]-[*]=[*]`` charge-shift form.

    RDKit's resonance supplier can omit this elementary migration when a
    molecule has many other resonance states.  Keep the transformation local
    and valence-safe: the negative endpoint gains a pi bond and is neutralized,
    while the terminal atom receives the negative charge as its pi bond is
    reduced to a single bond.
    """

    forms: list[Chem.Mol] = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE or bond.GetIsAromatic():
            continue
        for negative_idx, center_idx in (
            (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            (bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()),
        ):
            negative = mol.GetAtomWithIdx(negative_idx)
            center = mol.GetAtomWithIdx(center_idx)
            if (
                negative.GetFormalCharge() != -1
                or center.GetFormalCharge() != 0
                or negative.GetNumRadicalElectrons() != 0
                or center.GetNumRadicalElectrons() != 0
            ):
                continue
            default_valence = pt.GetDefaultValence(negative.GetAtomicNum())
            if default_valence <= 0 or negative.GetTotalValence() + 1 > default_valence:
                continue
            for pi_bond in center.GetBonds():
                terminal = pi_bond.GetOtherAtom(center)
                if terminal.GetIdx() == negative_idx:
                    continue
                if (
                    pi_bond.GetBondType() != Chem.BondType.DOUBLE
                    or pi_bond.GetIsAromatic()
                    or terminal.GetFormalCharge() != 0
                    or terminal.GetNumRadicalElectrons() != 0
                ):
                    continue
                shifted = Chem.RWMol(mol)
                shifted.GetBondBetweenAtoms(negative_idx, center_idx).SetBondType(
                    Chem.BondType.DOUBLE
                )
                shifted.GetBondBetweenAtoms(center_idx, terminal.GetIdx()).SetBondType(
                    Chem.BondType.SINGLE
                )
                shifted.GetAtomWithIdx(negative_idx).SetFormalCharge(0)
                shifted.GetAtomWithIdx(terminal.GetIdx()).SetFormalCharge(-1)
                form = shifted.GetMol()
                form.UpdatePropertyCache(strict=False)
                _sanitize_if_possible(form)
                forms.append(form)
    return tuple(forms)


@lru_cache(maxsize=128)
def _cached_resonance_smiles(
    source_smiles: str,
    *,
    use_chirality: bool,
    max_resonance: int,
    resonance_flags: int,
) -> tuple[str, ...]:
    source_mol = Chem.MolFromSmiles(source_smiles)
    if source_mol is None:
        return ()
    source = _prepare_resonance_source(source_mol)
    seen: list[str] = []

    def add_form(form: Chem.Mol) -> None:
        if len(seen) >= max_resonance:
            return
        smiles = _resonance_form_smiles(form, use_chirality)
        if smiles not in seen:
            seen.append(smiles)

    supplier_forms: list[Chem.Mol] = []
    try:
        supplier = ResonanceMolSupplier(
            source,
            maxStructs=max_resonance,
            flags=Chem.ResonanceFlags(resonance_flags),
        )
        for resonance_mol in supplier:
            if resonance_mol is None:
                continue
            supplier_forms.append(resonance_mol)
    except TimeoutError:
        raise
    except Exception:  # noqa: BLE001
        supplier_forms = []

    # Reserve the bounded result set for explicit charge-shift forms before
    # adding RDKit's broader enumeration.  Walk the local transformation to a
    # bounded closure because a charge may need to cross more than one
    # conjugated bond before reaching the reference representation.
    charge_shift_queue = [source, *supplier_forms]
    queued_smiles: set[str] = set()
    queue_index = 0
    while queue_index < len(charge_shift_queue) and len(seen) < max_resonance:
        resonance_form = charge_shift_queue[queue_index]
        queue_index += 1
        form_smiles = _resonance_form_smiles(resonance_form, use_chirality)
        if form_smiles in queued_smiles:
            continue
        queued_smiles.add(form_smiles)
        for shifted_form in _charge_shift_resonance_forms(resonance_form):
            shifted_smiles = _resonance_form_smiles(shifted_form, use_chirality)
            if shifted_smiles in queued_smiles:
                continue
            queued_smiles.add(shifted_smiles)
            charge_shift_queue.append(shifted_form)
            add_form(shifted_form)
            if len(seen) >= max_resonance:
                return tuple(dict.fromkeys(seen))
    for resonance_form in supplier_forms:
        add_form(resonance_form)
    return tuple(dict.fromkeys(seen))


def _resonance_topology_key(mol: Chem.Mol) -> str:
    rw_mol = Chem.RWMol(Chem.Mol(mol))
    for atom in rw_mol.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNumRadicalElectrons(0)
        atom.SetIsAromatic(False)
        atom.SetNoImplicit(True)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in rw_mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
        bond.SetBondDir(Chem.BondDir.NONE)
        bond.SetStereo(Chem.BondStereo.STEREONONE)
    topology = rw_mol.GetMol()
    topology.UpdatePropertyCache(strict=False)
    return Chem.MolToSmiles(topology, canonical=True, isomericSmiles=False)


def _has_defined_stereochemistry(mol: Chem.Mol) -> bool:
    if any(atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for atom in mol.GetAtoms()):
        return True
    return any(bond.GetStereo() != Chem.BondStereo.STEREONONE for bond in mol.GetBonds())


def _component_electron_signature(mol: Chem.Mol) -> tuple[tuple[str, int, int], ...]:
    signature = []
    for fragment in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False):
        signature.append(
            (
                _resonance_topology_key(fragment),
                _total_formal_charge(fragment),
                _total_radical_electrons(fragment),
            )
        )
    return tuple(sorted(signature))


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
    _kekulize_if_possible(rw_mol)

    if _CARBENE_ZWITTERION_PATTERN is not None:
        for match in rw_mol.GetSubstructMatches(_CARBENE_ZWITTERION_PATTERN):
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
    _sanitize_if_possible(normalized_mol)
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
    _sanitize_if_possible(normalized_mol)
    return normalized_mol


def _resonance_match(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    *,
    use_chirality: bool,
    max_resonance: int,
    resonance_flags: Chem.ResonanceFlags,
) -> tuple[bool, int, int, str | None]:
    def match_sets(mol1_smiles: set[str], mol2_smiles: set[str], target_2_smiles: str):
        if target_2_smiles in mol1_smiles:
            return True, len(mol1_smiles), 0, target_2_smiles
        intersection = mol1_smiles & mol2_smiles
        if intersection:
            return True, len(mol1_smiles), len(mol2_smiles), next(iter(intersection))
        return False, len(mol1_smiles), len(mol2_smiles), None

    def local_charge_shift_smiles(mol: Chem.Mol) -> set[str]:
        source = _prepare_resonance_source(mol)
        queue = [source]
        smiles: set[str] = set()
        queue_index = 0
        while queue_index < len(queue) and len(smiles) < max_resonance:
            form = queue[queue_index]
            queue_index += 1
            form_smiles = _resonance_form_smiles(form, use_chirality)
            if form_smiles in smiles:
                continue
            smiles.add(form_smiles)
            for shifted_form in _charge_shift_resonance_forms(form):
                shifted_smiles = _resonance_form_smiles(shifted_form, use_chirality)
                if shifted_smiles not in smiles:
                    queue.append(shifted_form)
        return smiles

    def cached_resonance_sets(target_1_smiles: str, target_2_smiles: str):
        try:
            mol1_smiles = set(
                _cached_resonance_smiles(
                    target_1_smiles,
                    use_chirality=use_chirality,
                    max_resonance=max_resonance,
                    resonance_flags=int(resonance_flags),
                )
            )
            mol2_smiles = set(
                _cached_resonance_smiles(
                    target_2_smiles,
                    use_chirality=use_chirality,
                    max_resonance=max_resonance,
                    resonance_flags=int(resonance_flags),
                )
            )
        except TimeoutError:
            raise
        return match_sets(mol1_smiles, mol2_smiles, target_2_smiles)

    try:
        # Preserve formal charges and bond orders for the generic local
        # [*-]-[*]=[*] migration.  The special normalizer below intentionally
        # collapses thiosemicarbazone-like motifs and would erase this input.
        raw_target_1 = _canon_smiles(mol1, use_chirality)
        raw_target_2 = _canon_smiles(mol2, use_chirality)
        raw_result = match_sets(
            local_charge_shift_smiles(mol1),
            local_charge_shift_smiles(mol2),
            raw_target_2,
        )
        if raw_result[0]:
            return raw_result

        normalized_1 = _normalize_special_resonance_forms(mol1)
        normalized_2 = _normalize_special_resonance_forms(mol2)
        normalized_target_1 = _canon_smiles(normalized_1, use_chirality)
        normalized_target_2 = _canon_smiles(normalized_2, use_chirality)
    except TimeoutError:
        raise
    except Exception:  # noqa: BLE001
        return False, 0, 0, None

    if (normalized_target_1, normalized_target_2) == (raw_target_1, raw_target_2):
        return cached_resonance_sets(raw_target_1, raw_target_2)
    return cached_resonance_sets(normalized_target_1, normalized_target_2)


def evaluate_equivalence(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    use_chirality: bool = True,
    max_resonance: int = 50,
    resonance_flags: Chem.ResonanceFlags = Chem.ResonanceFlags.UNCONSTRAINED_CATIONS
    | Chem.ResonanceFlags.UNCONSTRAINED_ANIONS,
) -> EquivalenceResult:
    """Evaluate molecular equivalence and retain the proof and invariant state.

    The result separates the comparison decision from the relation that supplied
    positive evidence.  Identifier equality is evidence only after the shared
    graph/electronic invariants have passed.
    """

    # Capture the source molecules before any graph normalization.  In
    # particular, _standardize_metal_bonds intentionally creates bookkeeping
    # radicals for metal double bonds; those generated labels must never be
    # presented as the source molecule's physical open-shell state.
    raw_mol1 = _raw_molecule_diagnostics(mol1)
    raw_mol2 = _raw_molecule_diagnostics(mol2)

    standardized_1 = _standardize_metal_bonds(mol1)
    standardized_2 = _standardize_metal_bonds(mol2)
    normalization_electronic_effects = {
        "mol1": {
            "raw_formal_charge": raw_mol1.formal_charge,
            "standardized_formal_charge": _total_formal_charge(standardized_1),
            "raw_rdkit_radical_electrons": raw_mol1.rdkit_radical_electrons,
            "standardized_rdkit_radical_electrons": _total_radical_electrons(standardized_1),
        },
        "mol2": {
            "raw_formal_charge": raw_mol2.formal_charge,
            "standardized_formal_charge": _total_formal_charge(standardized_2),
            "raw_rdkit_radical_electrons": raw_mol2.rdkit_radical_electrons,
            "standardized_rdkit_radical_electrons": _total_radical_electrons(standardized_2),
        },
    }
    normalization_electronic_effects["candidate"] = normalization_electronic_effects["mol1"]
    normalization_electronic_effects["reference"] = normalization_electronic_effects["mol2"]

    prepared_organic_1 = _prepare_organic_mol(standardized_1, already_standardized=True)
    prepared_organic_2 = _prepare_organic_mol(standardized_2, already_standardized=True)
    organic_1 = _normalize_nonmetal_octet(prepared_organic_1)
    organic_2 = _normalize_nonmetal_octet(prepared_organic_2)

    formal_charge_1 = _total_formal_charge(organic_1)
    formal_charge_2 = _total_formal_charge(organic_2)
    radical_electrons_1 = _total_radical_electrons(organic_1)
    radical_electrons_2 = _total_radical_electrons(organic_2)
    num_atoms_1 = organic_1.GetNumAtoms()
    num_atoms_2 = organic_2.GetNumAtoms()
    heavy_atom_formula_1 = _formula_key(organic_1, include_hydrogen=False)
    heavy_atom_formula_2 = _formula_key(organic_2, include_hydrogen=False)
    explicit_h_formula_1 = _formula_key(organic_1, include_hydrogen=True)
    explicit_h_formula_2 = _formula_key(organic_2, include_hydrogen=True)

    checks = EquivalenceChecks(
        formal_charge=PropertyCheck(
            formal_charge_1,
            formal_charge_2,
            formal_charge_1 == formal_charge_2,
        ),
        radical_electrons=PropertyCheck(
            radical_electrons_1,
            radical_electrons_2,
            radical_electrons_1 == radical_electrons_2,
        ),
        num_atoms=PropertyCheck(
            num_atoms_1,
            num_atoms_2,
            num_atoms_1 == num_atoms_2,
        ),
        heavy_atom_formula=PropertyCheck(
            heavy_atom_formula_1,
            heavy_atom_formula_2,
            heavy_atom_formula_1 == heavy_atom_formula_2,
        ),
        explicit_h_formula=PropertyCheck(
            explicit_h_formula_1,
            explicit_h_formula_2,
            explicit_h_formula_1 == explicit_h_formula_2,
        ),
    )
    invariants = {
        "heavy_atom_formula": InvariantResult(
            InvariantStatus.PASSED if checks.heavy_atom_formula.passed else InvariantStatus.FAILED,
            checks.heavy_atom_formula.mol1,
            checks.heavy_atom_formula.mol2,
        ),
        "explicit_h_formula": InvariantResult(
            InvariantStatus.PASSED if checks.explicit_h_formula.passed else InvariantStatus.FAILED,
            checks.explicit_h_formula.mol1,
            checks.explicit_h_formula.mol2,
        ),
        "atom_count": InvariantResult(
            InvariantStatus.PASSED if checks.num_atoms.passed else InvariantStatus.FAILED,
            checks.num_atoms.mol1,
            checks.num_atoms.mol2,
        ),
        "formal_charge": InvariantResult(
            InvariantStatus.PASSED if checks.formal_charge.passed else InvariantStatus.FAILED,
            checks.formal_charge.mol1,
            checks.formal_charge.mol2,
        ),
        "radical_electrons": InvariantResult(
            InvariantStatus.PASSED if checks.radical_electrons.passed else InvariantStatus.FAILED,
            checks.radical_electrons.mol1,
            checks.radical_electrons.mol2,
        ),
        "metal_state": InvariantResult(InvariantStatus.NOT_EVALUATED),
        "nonmetal_connectivity": InvariantResult(InvariantStatus.NOT_EVALUATED),
        "component_electrons": InvariantResult(InvariantStatus.NOT_EVALUATED),
        "stereochemistry": InvariantResult(InvariantStatus.NOT_EVALUATED),
    }
    info = EquivalenceResult(
        checks=checks,
        invariants=invariants,
        contradictions=[],
        bounded_search=BoundedSearchMetadata(limit=max_resonance),
        raw_mol1=raw_mol1,
        raw_mol2=raw_mol2,
        raw_diagnostics={
            "mol1": raw_mol1,
            "mol2": raw_mol2,
            # Positional aliases are intentional: callers that label the
            # inputs Candidate/Reference can consume the same snapshot without
            # a second normalization pass.
            "candidate": raw_mol1,
            "reference": raw_mol2,
        },
        normalization_electronic_effects=normalization_electronic_effects,
    )

    def finish(
        decision: EquivalenceDecision,
        reason: str,
        *,
        relation: EquivalenceRelation = EquivalenceRelation.NONE,
        method: EquivalenceMethod | None = None,
    ) -> EquivalenceResult:
        info.decision = decision
        info.relation = relation
        info.equivalent = decision == EquivalenceDecision.EQUIVALENT
        info.method = method
        info.reason = reason
        return info

    if not checks.heavy_atom_formula.passed:
        return finish(
            EquivalenceDecision.NOT_EQUIVALENT,
            "Not equivalent: heavy-atom element counts differ.",
        )

    metal_signature_1 = _metal_signature_key(standardized_1)
    metal_signature_2 = _metal_signature_key(standardized_2)
    invariants["metal_state"] = InvariantResult(
        InvariantStatus.PASSED
        if metal_signature_1 == metal_signature_2
        else InvariantStatus.FAILED,
        metal_signature_1,
        metal_signature_2,
    )
    if metal_signature_1 != metal_signature_2:
        return finish(
            EquivalenceDecision.NOT_EQUIVALENT,
            "Not equivalent: metal valence assignment differs.",
        )

    if not checks.explicit_h_formula.passed:
        return finish(
            EquivalenceDecision.NOT_EQUIVALENT,
            "Not equivalent: explicit-hydrogen element counts differ.",
        )

    if not checks.num_atoms.passed:
        return finish(
            EquivalenceDecision.NOT_EQUIVALENT,
            "Not equivalent: explicit-hydrogen atom counts differ.",
        )

    try:
        smiles_1 = _canon_smiles(organic_1, use_chirality)
        smiles_2 = _canon_smiles(organic_2, use_chirality)
        info.canonical_smiles = CanonicalSmilesDetail(smiles_1, smiles_2, use_chirality)
        topology_1 = _resonance_topology_key(organic_1)
        topology_2 = _resonance_topology_key(organic_2)
        topology_matches = topology_1 == topology_2
        invariants["nonmetal_connectivity"] = InvariantResult(
            InvariantStatus.PASSED if topology_matches else InvariantStatus.FAILED,
            topology_1,
            topology_2,
        )
        if not topology_matches:
            return finish(
                EquivalenceDecision.NOT_EQUIVALENT,
                "Not equivalent: non-metal connectivity differs.",
            )

        component_signature_1 = _component_electron_signature(organic_1)
        component_signature_2 = _component_electron_signature(organic_2)
        component_electrons_match = component_signature_1 == component_signature_2
        invariants["component_electrons"] = InvariantResult(
            InvariantStatus.PASSED if component_electrons_match else InvariantStatus.FAILED,
            component_signature_1,
            component_signature_2,
        )

        achiral_smiles_match = _canon_smiles(organic_1, False) == _canon_smiles(organic_2, False)
        stereo_matches = not use_chirality or not achiral_smiles_match or smiles_1 == smiles_2
        invariants["stereochemistry"] = InvariantResult(
            InvariantStatus.PASSED if stereo_matches else InvariantStatus.FAILED,
            smiles_1,
            smiles_2,
        )
    except TimeoutError:
        raise
    except Exception as exc:  # noqa: BLE001
        return finish(
            EquivalenceDecision.INCONCLUSIVE,
            f"Inconclusive: canonical comparison failed: {type(exc).__name__}: {exc}",
        )

    prepared_inchi_source_1 = prepared_organic_1
    prepared_inchi_source_2 = prepared_organic_2
    if not any(_is_metal_atom(atom) for atom in mol1.GetAtoms()) and not any(
        _is_metal_atom(atom) for atom in mol2.GetAtoms()
    ):
        prepared_inchi_source_1 = mol1
        prepared_inchi_source_2 = mol2
    prepared_inchi_key_1 = _inchi_key(prepared_inchi_source_1, use_chirality=use_chirality)
    prepared_inchi_key_2 = _inchi_key(prepared_inchi_source_2, use_chirality=use_chirality)
    full_inchi_key_1 = _inchi_key(organic_1, use_chirality=use_chirality)
    full_inchi_key_2 = _inchi_key(organic_2, use_chirality=use_chirality)
    prepared_identifier_match = (
        prepared_inchi_key_1 is not None and prepared_inchi_key_1 == prepared_inchi_key_2
    )
    full_identifier_match = full_inchi_key_1 is not None and full_inchi_key_1 == full_inchi_key_2
    identifier_match = prepared_identifier_match or full_identifier_match

    try:
        normalized_1 = _carbene_zwitterion_normalized_smiles(organic_1, use_chirality=use_chirality)
        normalized_2 = _carbene_zwitterion_normalized_smiles(organic_2, use_chirality=use_chirality)
        info.carbene_zwitterion = CarbeneZwitterionDetail(
            mol1_normalized=normalized_1,
            mol2_normalized=normalized_2,
        )
        carbene_match = normalized_1 == normalized_2
    except TimeoutError:
        raise
    except Exception:
        carbene_match = False

    if not component_electrons_match:
        if identifier_match:
            info.contradictions.append("identifier_match_despite_component_electron_mismatch")
        if carbene_match:
            info.contradictions.append(
                "carbene_normalization_match_despite_component_electron_mismatch"
            )
        return finish(
            EquivalenceDecision.NOT_EQUIVALENT,
            "Not equivalent: charge or radical count differs within a component.",
        )

    if use_chirality and not stereo_matches:
        if identifier_match:
            info.contradictions.append("identifier_match_despite_stereochemistry_mismatch")
        return finish(
            EquivalenceDecision.NOT_EQUIVALENT,
            "Not equivalent: stereochemistry differs.",
        )

    if smiles_1 == smiles_2:
        return finish(
            EquivalenceDecision.EQUIVALENT,
            "Equivalent: canonical SMILES are identical after standardization.",
            relation=EquivalenceRelation.NORMALIZED_GRAPH_IDENTITY,
            method=EquivalenceMethod.IDEAL,
        )

    if carbene_match:
        return finish(
            EquivalenceDecision.EQUIVALENT,
            "Equivalent: carbene/zwitterion normalization matched.",
            relation=EquivalenceRelation.CARBENE_ZWITTERION_EQUIVALENCE,
            method=EquivalenceMethod.CARBENE_ZWITTERION,
        )

    # A matching topology plus matching component electron totals is the
    # evaluator's existing stronger radical-resonance evidence.  It is not an
    # InChI-only result, so it remains valid even when identifier generation is
    # unavailable.  Closed-shell cases still require the bounded enumerator.
    if (
        component_electrons_match
        and checks.formal_charge.passed
        and checks.radical_electrons.passed
        and radical_electrons_1 > 0
        and (
            not use_chirality
            or (
                not _has_defined_stereochemistry(organic_1)
                and not _has_defined_stereochemistry(organic_2)
            )
        )
    ):
        info.resonance = ResonanceDetail(
            max_resonance=max_resonance,
            resonance_flags=int(resonance_flags),
            mol1_resonance_count=0,
            mol2_resonance_count=0,
        )
        return finish(
            EquivalenceDecision.EQUIVALENT,
            "Equivalent: resonance topology and component electron counts match.",
            relation=EquivalenceRelation.RESONANCE_EQUIVALENCE,
            method=EquivalenceMethod.RESONANCE,
        )

    assert info.bounded_search is not None
    info.bounded_search.attempted = True
    resonance_matched, mol1_count, mol2_count, hit_smiles = _resonance_match(
        organic_1,
        organic_2,
        use_chirality=use_chirality,
        max_resonance=max_resonance,
        resonance_flags=resonance_flags,
    )
    info.bounded_search.mol1_count = mol1_count
    info.bounded_search.mol2_count = mol2_count
    info.bounded_search.limit_reached = not resonance_matched and (
        mol1_count >= max_resonance or mol2_count >= max_resonance
    )
    info.bounded_search.exhaustive = not info.bounded_search.limit_reached
    if resonance_matched:
        info.resonance = ResonanceDetail(
            max_resonance=max_resonance,
            resonance_flags=int(resonance_flags),
            mol1_resonance_count=mol1_count,
            mol2_resonance_count=mol2_count,
            hit_smiles=hit_smiles,
        )
        return finish(
            EquivalenceDecision.EQUIVALENT,
            "Equivalent: resonance normalization matched.",
            relation=EquivalenceRelation.RESONANCE_EQUIVALENCE,
            method=EquivalenceMethod.RESONANCE,
        )

    if identifier_match:
        reason = (
            "Inconclusive: identifier_equivalence (prepared organic InChIKey agreement) "
            "has no independent stronger structural evidence."
            if prepared_identifier_match
            else "Inconclusive: identifier_equivalence (full InChIKey agreement) "
            "has no independent stronger structural evidence."
        )
        return finish(
            EquivalenceDecision.INCONCLUSIVE,
            reason,
            relation=EquivalenceRelation.IDENTIFIER_EQUIVALENCE,
            method=EquivalenceMethod.INCHI_KEY,
        )

    if info.bounded_search.limit_reached:
        return finish(
            EquivalenceDecision.INCONCLUSIVE,
            "Inconclusive: bounded resonance search reached its configured limit.",
        )

    return finish(
        EquivalenceDecision.NOT_EQUIVALENT,
        "Not equivalent: no standardized comparison path matched.",
    )


def check_equivalence(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    use_chirality: bool = True,
    max_resonance: int = 50,
    resonance_flags: Chem.ResonanceFlags = Chem.ResonanceFlags.UNCONSTRAINED_CATIONS
    | Chem.ResonanceFlags.UNCONSTRAINED_ANIONS,
) -> Tuple[bool, EquivalenceInfo]:
    with rdBase.BlockLogs():
        result = evaluate_equivalence(
            mol1,
            mol2,
            use_chirality=use_chirality,
            max_resonance=max_resonance,
            resonance_flags=resonance_flags,
        )
    return result.decision == EquivalenceDecision.EQUIVALENT, result
