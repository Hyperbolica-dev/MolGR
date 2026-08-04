#!/usr/bin/env python3
"""Build review-approved molecule graph fixtures."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import threading
from pathlib import Path
from typing import Any, Mapping

from rdkit import Chem

from molgr.interface import xyz_to_rdmol
from molgr.utils.converter import (
    METAL_UNPAIRED_ELECTRONS_PROP,
    get_atom_unpaired_electrons,
)


APP_DIR = Path(__file__).resolve().parent
REPO_ROOT = APP_DIR.parents[1]
DEFAULT_FIXTURES_DIR = Path(
    os.environ.get(
        "MOLGR_REVIEW_FIXTURES_DIR",
        REPO_ROOT / "tests/data/reviewed",
    )
)
MANIFEST_NAME = "manifest.json"
ACCEPT_BOTH_STATUS = "accept_both"
APPROVED_GRAPH_STATUSES = {"accept_candidate", ACCEPT_BOTH_STATUS}
REFERENCE_GRAPH_STATUS = "accept_reference"
MANUAL_REFERENCE_STATUS = "manual_reference"
ANNOTATION_ONLY_STATUSES = {"reference_answer_wrong"}
FIXTURE_STATUSES = APPROVED_GRAPH_STATUSES | {
    REFERENCE_GRAPH_STATUS,
    MANUAL_REFERENCE_STATUS,
}
FIXTURE_KINDS = {"accepted_both", "approved_graph", "reference_graph", "manual_reference"}
FIXTURE_RECORD_KEYS = {
    "case_id",
    "kind",
    "structure_file",
    "row_index",
    "total_charge",
    "total_radical_electrons",
    "spin_multiplicity",
    "reference_smiles",
    "approved_smiles",
    "accepted_smiles",
    "source",
    "review_status",
}
_FIXTURE_LOCK = threading.Lock()


def resolve_xyz_path(stored_path: str, xyz_dir: Path | None) -> Path:
    path = Path(stored_path)
    if path.exists() or xyz_dir is None:
        return path
    remapped = xyz_dir / path.name
    if remapped.exists():
        return remapped
    return path


def _raw_payload(case: Mapping[str, Any]) -> dict[str, Any]:
    metadata_json = case.get("metadata_json")
    if isinstance(metadata_json, str) and metadata_json:
        try:
            payload = json.loads(metadata_json)
        except json.JSONDecodeError:
            return {}
        if isinstance(payload, dict):
            return payload
    return {}


def _optional_int(value: Any) -> int | None:
    if value in (None, ""):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _molecule_electronic_state(mol: Chem.Mol) -> tuple[int, int]:
    total_charge = 0
    total_radical_electrons = 0
    for atom in mol.GetAtoms():  # pyright: ignore[reportCallIssue]
        total_charge += int(atom.GetFormalCharge())
        total_radical_electrons += get_atom_unpaired_electrons(atom)
    return total_charge, total_radical_electrons


def case_electronic_state(case: Mapping[str, Any]) -> tuple[int, int, int]:
    """Return total charge, radical electrons, and spin multiplicity for one case."""

    raw = _raw_payload(case)
    total_charge = _optional_int(case.get("total_charge"))
    if total_charge is None:
        total_charge = _optional_int(case.get("charge")) or 0
    spin_multiplicity = _optional_int(case.get("spin_multiplicity"))
    if spin_multiplicity is None:
        spin_multiplicity = _optional_int(raw.get("spin_multiplicity_used"))
    total_radical_electrons = _optional_int(case.get("total_radical_electrons"))
    if total_radical_electrons is None:
        total_radical_electrons = _optional_int(raw.get("total_radical_electrons_used"))

    if total_radical_electrons is None and spin_multiplicity is not None:
        total_radical_electrons = max(0, spin_multiplicity - 1)
    if total_radical_electrons is None:
        smiles = str(case.get("reference_smiles") or "").strip()
        reference = Chem.MolFromSmiles(smiles, sanitize=False) if smiles else None
        if reference is not None:
            _, total_radical_electrons = _molecule_electronic_state(reference)
    if total_radical_electrons is None:
        total_radical_electrons = 0
    if spin_multiplicity is None:
        spin_multiplicity = total_radical_electrons + 1
    return total_charge, total_radical_electrons, max(1, spin_multiplicity)


def reconstruct_case_mol(
    case: Mapping[str, Any],
    *,
    xyz_dir: Path | None,
) -> Chem.Mol:
    xyz_path = resolve_xyz_path(str(case.get("xyz_path") or ""), xyz_dir)
    if not xyz_path.exists():
        raise FileNotFoundError(xyz_path)
    total_charge, _, spin_multiplicity = case_electronic_state(case)
    return xyz_to_rdmol(
        xyz_path.read_text(encoding="utf-8"),
        total_charge=total_charge,
        spin_multiplicity=spin_multiplicity,
        backend="cpp",
        make_dative_bonds=True,
        make_stereochemistry=True,
    )


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temp_path = Path(handle.name)
        handle.write(text)
    temp_path.replace(path)


def _atomic_copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=target.parent,
        prefix=f".{target.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temp_path = Path(handle.name)
    try:
        shutil.copyfile(source, temp_path)
        temp_path.replace(target)
    finally:
        temp_path.unlink(missing_ok=True)


def _sdf_text(mol: Chem.Mol, properties: Mapping[str, Any]) -> str:
    block = Chem.MolToMolBlock(mol, forceV3000=True)
    if not block.endswith("\n"):
        block += "\n"
    fields = []
    atom_property_list = f"atom.iprop.{METAL_UNPAIRED_ELECTRONS_PROP}"
    if mol.HasProp(atom_property_list):
        fields.append(f">  <{atom_property_list}>\n{mol.GetProp(atom_property_list)}\n")
    for key, value in properties.items():
        fields.append(f">  <{key}>\n{value}\n")
    return block + "\n".join(fields) + "$$$$\n"


def _load_manifest(fixtures_dir: Path) -> dict[str, Any]:
    path = fixtures_dir / MANIFEST_NAME
    if not path.exists():
        return {"schema_version": 1, "fixtures": []}
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("fixtures"), list):
        raise ValueError(f"invalid review fixture manifest: {path}")
    invalid_kinds = {
        str(record.get("kind") or "")
        for record in payload["fixtures"]
        if str(record.get("kind") or "") not in FIXTURE_KINDS
    }
    if invalid_kinds:
        raise ValueError(
            "review fixture manifest contains non-confirmed answers: "
            + ", ".join(sorted(invalid_kinds))
        )
    return payload


def load_fixture_records(fixtures_dir: Path = DEFAULT_FIXTURES_DIR) -> dict[str, dict[str, Any]]:
    """Return the managed fixture records indexed by case id."""

    manifest = _load_manifest(fixtures_dir)
    records: dict[str, dict[str, Any]] = {}
    for record in manifest["fixtures"]:
        case_id = str(record.get("case_id") or "").strip()
        if case_id:
            records[case_id] = dict(record)
    return records


def _write_manifest(
    fixtures_dir: Path,
    records: list[dict[str, Any]],
    *,
    existing: Mapping[str, Any] | None = None,
) -> None:
    payload = {key: value for key, value in (existing or {}).items() if key != "fixtures"}
    payload["schema_version"] = 1
    payload["fixtures"] = sorted(
        [
            {key: value for key, value in record.items() if key in FIXTURE_RECORD_KEYS}
            for record in records
        ],
        key=lambda item: (item["case_id"], item["kind"]),
    )
    _atomic_write_text(
        fixtures_dir / MANIFEST_NAME,
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )


def _remove_case_files(fixtures_dir: Path, case_id: str) -> None:
    (fixtures_dir / "accepted_both" / f"{case_id}.sdf").unlink(missing_ok=True)
    (fixtures_dir / "approved_graph" / f"{case_id}.sdf").unlink(missing_ok=True)
    (fixtures_dir / "reference_graph" / f"{case_id}.xyz").unlink(missing_ok=True)
    (fixtures_dir / "manual_reference" / f"{case_id}.xyz").unlink(missing_ok=True)


def sync_review_fixture(
    case: Mapping[str, Any],
    review: Mapping[str, Any],
    *,
    fixtures_dir: Path = DEFAULT_FIXTURES_DIR,
    xyz_dir: Path | None = None,
) -> dict[str, Any] | None:
    """Create, replace, or remove the fixture selected by one review decision."""

    case_id = str(case.get("case_id") or "").strip()
    if not case_id:
        raise ValueError("case_id_is_required")
    status = str(review.get("status") or "").strip()

    with _FIXTURE_LOCK:
        manifest = _load_manifest(fixtures_dir)
        existing_record = next(
            (record for record in manifest["fixtures"] if record.get("case_id") == case_id),
            None,
        )
        records = [record for record in manifest["fixtures"] if record.get("case_id") != case_id]
        if status in ANNOTATION_ONLY_STATUSES:
            return existing_record
        if status not in FIXTURE_STATUSES:
            _remove_case_files(fixtures_dir, case_id)
            _write_manifest(fixtures_dir, records, existing=manifest)
            return None

        row_index = _optional_int(case.get("row_index")) or 0
        reference_smiles = str(case.get("reference_smiles") or "").strip()
        xyz_path = resolve_xyz_path(str(case.get("xyz_path") or ""), xyz_dir)
        if not xyz_path.exists():
            raise FileNotFoundError(xyz_path)

        if status == REFERENCE_GRAPH_STATUS:
            if not reference_smiles:
                raise ValueError("reference_graph_smiles_is_required")
            total_charge, total_radical_electrons, spin_multiplicity = case_electronic_state(case)
            relative_path = Path("reference_graph") / f"{case_id}.xyz"
            source = "reference_graph"
            approved_smiles = reference_smiles
            kind = "reference_graph"
        elif status == MANUAL_REFERENCE_STATUS:
            approved_smiles = str(review.get("corrected_smiles") or "").strip()
            if not approved_smiles:
                raise ValueError("corrected_smiles_required_for_manual_reference_fixture")
            if Chem.MolFromSmiles(approved_smiles, sanitize=False) is None:
                raise ValueError("corrected_smiles_is_invalid")
            total_charge, total_radical_electrons, spin_multiplicity = case_electronic_state(case)
            relative_path = Path("manual_reference") / f"{case_id}.xyz"
            source = "corrected_smiles"
            kind = "manual_reference"
        else:
            if status == ACCEPT_BOTH_STATUS and not reference_smiles:
                raise ValueError("reference_graph_smiles_is_required")
            mol = reconstruct_case_mol(case, xyz_dir=xyz_dir)
            source = "molgr_reconstruction"
            # SDF-level electronic constraints describe the original case, not
            # metal-local unpaired electrons carried by the reconstructed graph.
            total_charge, total_radical_electrons, spin_multiplicity = case_electronic_state(case)
            try:
                approved_smiles = Chem.MolToSmiles(
                    mol, canonical=True, isomericSmiles=True
                )
            except Chem.KekulizeException:
                approved_smiles = Chem.MolToSmiles(
                    mol, canonical=True, isomericSmiles=True, kekuleSmiles=False
                )
            if status == ACCEPT_BOTH_STATUS:
                relative_path = Path("accepted_both") / f"{case_id}.sdf"
                source = "candidate_or_reference_graph"
                kind = "accepted_both"
            else:
                relative_path = Path("approved_graph") / f"{case_id}.sdf"
                kind = "approved_graph"
            properties = {
                "CASE_ID": case_id,
                "CANDIDATE_ID": case_id,
                "TOTAL_CHARGE": total_charge,
                "TOTAL_RADICAL_ELECTRONS": total_radical_electrons,
                "SPIN_MULTIPLICITY": spin_multiplicity,
                "REVIEW_STATUS": status,
                "APPROVED_SMILES": approved_smiles,
            }

        accepted_smiles = list(dict.fromkeys([approved_smiles, reference_smiles]))
        if status not in {REFERENCE_GRAPH_STATUS, ACCEPT_BOTH_STATUS}:
            accepted_smiles = [approved_smiles]

        record = {
            "case_id": case_id,
            "kind": kind,
            "structure_file": relative_path.as_posix(),
            "row_index": row_index,
            "total_charge": total_charge,
            "total_radical_electrons": total_radical_electrons,
            "spin_multiplicity": spin_multiplicity,
            "reference_smiles": reference_smiles,
            "approved_smiles": approved_smiles,
            "accepted_smiles": accepted_smiles,
            "source": source,
            "review_status": status,
        }
        _remove_case_files(fixtures_dir, case_id)
        if kind in {"accepted_both", "approved_graph"}:
            _atomic_write_text(fixtures_dir / relative_path, _sdf_text(mol, properties))
        else:
            _atomic_copy(xyz_path, fixtures_dir / relative_path)
        records.append(record)
        _write_manifest(fixtures_dir, records, existing=manifest)
        return record


def prune_unreviewed_fixtures(fixtures_dir: Path, reviewed_case_ids: set[str]) -> None:
    """Remove manifest entries that no longer have a corresponding review row."""

    with _FIXTURE_LOCK:
        manifest = _load_manifest(fixtures_dir)
        stale_ids = {
            str(record.get("case_id") or "")
            for record in manifest["fixtures"]
            if record.get("case_id") not in reviewed_case_ids
        }
        for case_id in stale_ids:
            _remove_case_files(fixtures_dir, case_id)
        records = [
            record for record in manifest["fixtures"] if record.get("case_id") in reviewed_case_ids
        ]
        _write_manifest(fixtures_dir, records, existing=manifest)


def remove_review_fixtures(fixtures_dir: Path, case_ids: set[str]) -> None:
    """Remove fixture records and files only for the requested review cases."""

    if not case_ids:
        return
    with _FIXTURE_LOCK:
        manifest = _load_manifest(fixtures_dir)
        for case_id in case_ids:
            _remove_case_files(fixtures_dir, case_id)
        records = [
            record for record in manifest["fixtures"] if record.get("case_id") not in case_ids
        ]
        _write_manifest(fixtures_dir, records, existing=manifest)


__all__ = [
    "ACCEPT_BOTH_STATUS",
    "APPROVED_GRAPH_STATUSES",
    "DEFAULT_FIXTURES_DIR",
    "FIXTURE_STATUSES",
    "MANUAL_REFERENCE_STATUS",
    "REFERENCE_GRAPH_STATUS",
    "case_electronic_state",
    "load_fixture_records",
    "prune_unreviewed_fixtures",
    "remove_review_fixtures",
    "reconstruct_case_mol",
    "resolve_xyz_path",
    "sync_review_fixture",
]
