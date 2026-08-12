from __future__ import annotations

import gzip
import hashlib
import heapq
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


EXPECTED_FILENAME = "20200415_radical_database.sdf.gz"
EXPECTED_RECORD_COUNT = 289_639
STRATUM_POPULATIONS = {
    "closed_shell": 43_276,
    "c_radical": 205_572,
    "n_radical": 25_434,
    "o_radical": 15_355,
}
REQUIRED_EDGE_SMILES = frozenset({"[H]", "[O-]"})


@dataclass(frozen=True)
class BDECase:
    case_idx: int
    case_id: str
    source_record_index: int
    xyz: str
    total_charge: int
    spin_multiplicity: int
    reference_mol: Any
    reference_smiles: str
    parent_id: str | None
    radical_site: int | None
    source_metadata: dict[str, Any]
    stratum: str

    def to_method_case(self) -> dict[str, Any]:
        return {
            "case_idx": self.case_idx,
            "case_id": self.case_id,
            "input_smiles": self.reference_smiles,
            "xyz_block": self.xyz,
            "total_charge": self.total_charge,
            "total_radical_electrons": self.spin_multiplicity - 1,
            "ground_truth_rdmol": self.reference_mol,
            "ground_truth_smiles": self.reference_smiles,
        }


@dataclass
class LoadDiagnostics:
    scanned_records: int = 0
    eligible_records: int = 0
    selected_records: int = 0
    strata_seen: dict[str, int] = field(default_factory=dict)
    failures: list[dict[str, str | int]] = field(default_factory=list)


def _score(seed: int, record_index: int, case_id: str) -> int:
    payload = f"{seed}:{record_index}:{case_id}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def _formal_charge(mol: Chem.Mol) -> int:
    return sum(int(atom.GetFormalCharge()) for atom in mol.GetAtoms())


def _radical_atoms(mol: Chem.Mol) -> list[tuple[int, int, str]]:
    return [
        (atom.GetIdx(), int(atom.GetNumRadicalElectrons()), atom.GetSymbol())
        for atom in mol.GetAtoms()
        if atom.GetNumRadicalElectrons()
    ]


def _derive_multiplicity(radical_atoms: list[tuple[int, int, str]]) -> int:
    radical_electrons = sum(count for _, count, _ in radical_atoms)
    if radical_electrons == 0:
        return 1
    if radical_electrons == 1:
        return 2
    raise ValueError(
        "spin multiplicity is not uniquely supported for structures with "
        f"{radical_electrons} encoded radical electrons"
    )


def _stratum(total_charge: int, radical_atoms: list[tuple[int, int, str]]) -> str:
    if total_charge != 0:
        return "charged"
    radical_electrons = sum(count for _, count, _ in radical_atoms)
    if radical_electrons == 0:
        return "closed_shell"
    if radical_electrons != 1 or len(radical_atoms) != 1:
        return "unsupported_multiradical"
    element = radical_atoms[0][2].lower()
    if element in {"c", "n", "o"}:
        return f"{element}_radical"
    return "other_radical"


def _metadata(
    mol: Chem.Mol,
    *,
    source_isomeric_smiles: str,
    sdf_isomeric_smiles: str,
) -> dict[str, Any]:
    properties = {
        name: mol.GetProp(name)
        for name in mol.GetPropNames(includePrivate=False, includeComputed=False)
    }
    return {
        "name": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "properties": properties,
        "source_isomeric_smiles": source_isomeric_smiles,
        "sdf_isomeric_smiles": sdf_isomeric_smiles,
        "stereo_agrees_with_smiles_property": source_isomeric_smiles == sdf_isomeric_smiles,
    }


def _without_hs_for_graph(mol: Chem.Mol) -> Chem.Mol:
    if all(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()):
        return Chem.Mol(mol)
    return Chem.RemoveAllHs(mol)


def _case_from_mol(mol: Chem.Mol, record_index: int) -> BDECase:
    if mol.GetNumConformers() != 1:
        raise ValueError(f"expected exactly one conformer, found {mol.GetNumConformers()}")
    if not mol.HasProp("SMILES"):
        raise ValueError("missing required SDF property: SMILES")
    reference_smiles = mol.GetProp("SMILES")
    source_smiles_mol = Chem.MolFromSmiles(reference_smiles)
    if source_smiles_mol is None:
        raise ValueError(f"invalid SMILES property: {reference_smiles!r}")
    source_isomeric_smiles = Chem.MolToSmiles(
        source_smiles_mol,
        canonical=True,
        isomericSmiles=True,
    )
    sdf_graph_mol = _without_hs_for_graph(mol)
    sdf_isomeric_smiles = Chem.MolToSmiles(
        sdf_graph_mol,
        canonical=True,
        isomericSmiles=True,
    )
    source_graph_smiles = Chem.MolToSmiles(
        source_smiles_mol,
        canonical=True,
        isomericSmiles=False,
    )
    sdf_graph_smiles = Chem.MolToSmiles(
        sdf_graph_mol,
        canonical=True,
        isomericSmiles=False,
    )
    if source_graph_smiles != sdf_graph_smiles:
        raise ValueError(
            "SDF reference graph disagrees with SMILES property: "
            f"{sdf_graph_smiles!r} != {source_graph_smiles!r}"
        )
    total_charge = _formal_charge(mol)
    radical_atoms = _radical_atoms(mol)
    spin_multiplicity = _derive_multiplicity(radical_atoms)
    stratum = _stratum(total_charge, radical_atoms)
    name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
    case_id = name or f"sdf-record-{record_index}"
    radical_site = radical_atoms[0][0] if len(radical_atoms) == 1 else None
    return BDECase(
        case_idx=0,
        case_id=case_id,
        source_record_index=record_index,
        xyz=Chem.MolToXYZBlock(mol),
        total_charge=total_charge,
        spin_multiplicity=spin_multiplicity,
        reference_mol=mol,
        reference_smiles=source_isomeric_smiles,
        parent_id=None,
        radical_site=radical_site,
        source_metadata=_metadata(
            mol,
            source_isomeric_smiles=source_isomeric_smiles,
            sdf_isomeric_smiles=sdf_isomeric_smiles,
        ),
        stratum=stratum,
    )


def _iter_sdf(path: Path) -> Iterator[tuple[int, Chem.Mol | None]]:
    with gzip.open(path, "rb") as stream:
        supplier = Chem.ForwardSDMolSupplier(
            stream,
            sanitize=True,
            removeHs=False,
            strictParsing=True,
        )
        yield from enumerate(supplier, start=1)


def _quota(limit: int) -> dict[str, int]:
    if limit <= 0:
        return dict.fromkeys(STRATUM_POPULATIONS, 0)
    total_weight = sum(STRATUM_POPULATIONS.values())
    quotas = {name: limit * weight // total_weight for name, weight in STRATUM_POPULATIONS.items()}
    remaining = limit - sum(quotas.values())
    remainder_order = sorted(
        STRATUM_POPULATIONS,
        key=lambda name: (-(limit * STRATUM_POPULATIONS[name] % total_weight), name),
    )
    for name in remainder_order:
        if remaining == 0:
            break
        quotas[name] += 1
        remaining -= 1
    return quotas


def load_bde_cases(
    input_path: Path,
    *,
    limit: int = 100,
    start: int | None = None,
    end: int | None = None,
    seed: int = 0,
) -> tuple[list[BDECase], LoadDiagnostics]:
    if not input_path.is_file():
        raise FileNotFoundError(f"BDE-db SDF does not exist: {input_path}")
    if limit < 0:
        raise ValueError("limit must be non-negative")
    if start is not None and start < 1:
        raise ValueError("start must be at least 1")
    if end is not None and end < 1:
        raise ValueError("end must be at least 1")
    if start is not None and end is not None and start > end:
        raise ValueError("start cannot be greater than end")

    diagnostics = LoadDiagnostics()
    candidates: dict[str, list[tuple[int, int, BDECase]]] = {}
    core_limit = max(0, limit - min(limit, len(REQUIRED_EDGE_SMILES)))
    quotas = _quota(core_limit)
    edge_cases: dict[str, BDECase] = {}
    for record_index, mol in _iter_sdf(input_path):
        diagnostics.scanned_records += 1
        if start is not None and record_index < start:
            continue
        if end is not None and record_index > end:
            break
        if mol is None:
            diagnostics.failures.append(
                {"record_index": record_index, "error": "RDKit failed to parse SDF record"}
            )
            continue
        try:
            case = _case_from_mol(mol, record_index)
        except Exception as exc:
            diagnostics.failures.append(
                {"record_index": record_index, "error": f"{type(exc).__name__}: {exc}"}
            )
            continue
        diagnostics.eligible_records += 1
        diagnostics.strata_seen[case.stratum] = diagnostics.strata_seen.get(case.stratum, 0) + 1
        if case.reference_smiles in REQUIRED_EDGE_SMILES:
            edge_cases[case.reference_smiles] = case
            continue
        capacity = quotas.get(case.stratum, 0) + len(REQUIRED_EDGE_SMILES)
        if capacity == 0:
            continue
        score = _score(seed, record_index, case.case_id)
        heap = candidates.setdefault(case.stratum, [])
        item = (-score, record_index, case)
        if len(heap) < capacity:
            heapq.heappush(heap, item)
        elif score < -heap[0][0]:
            heapq.heapreplace(heap, item)

    ordered_candidates = {
        name: [item[2] for item in sorted(heap, key=lambda item: (-item[0], item[1]))]
        for name, heap in candidates.items()
    }
    selected: list[BDECase] = []
    selected_keys: set[tuple[str, int]] = set()
    for smiles in sorted(REQUIRED_EDGE_SMILES):
        case = edge_cases.get(smiles)
        if case is None or len(selected) == limit:
            continue
        selected.append(case)
        selected_keys.add((case.case_id, case.source_record_index))
    remaining_limit = limit - len(selected)
    for name, count in _quota(remaining_limit).items():
        for case in ordered_candidates.get(name, [])[:count]:
            selected.append(case)
            selected_keys.add((case.case_id, case.source_record_index))
    if len(selected) < limit:
        represented = {case.stratum for case in selected}
        for name in STRATUM_POPULATIONS:
            if len(selected) == limit:
                break
            if name in represented or not ordered_candidates.get(name):
                continue
            case = ordered_candidates[name][0]
            selected.append(case)
            selected_keys.add((case.case_id, case.source_record_index))
            represented.add(name)
    if len(selected) < limit:
        remainder = sorted(
            (
                _score(seed, case.source_record_index, case.case_id),
                case.source_record_index,
                case,
            )
            for cases in ordered_candidates.values()
            for case in cases
            if (case.case_id, case.source_record_index) not in selected_keys
        )
        selected.extend(item[2] for item in remainder[: limit - len(selected)])
    selected.sort(key=lambda case: (case.source_record_index, case.case_id))
    indexed = [
        BDECase(
            case_idx=index,
            case_id=case.case_id,
            source_record_index=case.source_record_index,
            xyz=case.xyz,
            total_charge=case.total_charge,
            spin_multiplicity=case.spin_multiplicity,
            reference_mol=case.reference_mol,
            reference_smiles=case.reference_smiles,
            parent_id=case.parent_id,
            radical_site=case.radical_site,
            source_metadata=case.source_metadata,
            stratum=case.stratum,
        )
        for index, case in enumerate(selected[:limit], start=1)
    ]
    diagnostics.selected_records = len(indexed)
    return indexed, diagnostics
