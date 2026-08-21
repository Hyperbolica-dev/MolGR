#!/usr/bin/env python3
"""Local molecule graph review server."""
# ruff: noqa: E402, I001

from __future__ import annotations

import argparse
import csv
import html
import json
import mimetypes
import os
import sqlite3
import sys
import threading
import traceback
import uuid
from datetime import datetime, timezone
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from collections import Counter
from collections.abc import Mapping
from typing import Any
from urllib.parse import parse_qs, quote, unquote, urlparse
from urllib.request import Request, urlopen


APP_DIR = Path(__file__).resolve().parent
REPO_ROOT = next(
    (
        candidate
        for candidate in (APP_DIR.parent, Path(__file__).resolve().parents[2])
        if (candidate / "scripts" / "reconstruction_trace.py").is_file()
    ),
    Path(__file__).resolve().parents[2],
)
if __package__ in (None, ""):
    sys.path.insert(0, str(APP_DIR))
    sys.path.insert(0, str(REPO_ROOT))
    sys.path.insert(0, str(REPO_ROOT / "src"))

from rdkit import Chem, RDLogger

from fixture_builder import (
    case_electronic_state,
    load_fixture_records,
    reconstruct_case_mol,
    resolve_xyz_path,
    restore_review_fixture_snapshot,
    snapshot_review_fixture,
    sync_review_fixture,
)
from molgr.fallback.utils.consts import NON_METAL_DICT
from molgr.utils.equivalence import check_equivalence
from project_runtime import validate_project_runtime
from reference_diagnostics import classify_reference_problem, comparison_skip_reasons
from scripts.reconstruction_trace import TraceInputCase, render_trace_report
from triage_mapping import atom_h_count, is_metal, map_candidate_reference_xyz, parse_xyz_atoms


RDLogger.DisableLog("rdApp.*")  # type: ignore[arg-type]

STATIC_DIR = APP_DIR / "static"
DEFAULT_DB = REPO_ROOT / ".local" / "molgr_review" / "review.sqlite"
DEFAULT_XYZ_DIR = os.environ.get("MOLGR_XYZ_DIR", "")
KETCHER_BASE_URL = "https://lifescience.opensource.epam.com"
ALLOWED_STATUSES = {
    "accept_both",
    "accept_candidate",
    "accept_reference",
    "manual_reference",
    "reference_answer_wrong",
    "needs_followup",
    "skip",
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db",
        type=Path,
        default=Path(os.environ.get("MOLGR_REVIEW_DB", DEFAULT_DB)),
        help="SQLite database path.",
    )
    parser.add_argument(
        "--host",
        default=os.environ.get("MOLGR_REVIEW_HOST", "127.0.0.1"),
        help="Bind host.",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=int(os.environ.get("MOLGR_REVIEW_PORT", "8765")),
        help="Bind port.",
    )
    parser.add_argument(
        "--xyz-dir",
        type=Path,
        default=Path(DEFAULT_XYZ_DIR) if DEFAULT_XYZ_DIR else None,
        help="Optional directory used to resolve XYZ files by basename.",
    )
    parser.add_argument(
        "--fixtures-dir",
        type=Path,
        required=True,
        help="Directory updated immediately from fixture-producing review decisions.",
    )
    parser.add_argument(
        "--triage-csv",
        type=Path,
        default=None,
        help="Optional read-only triage CSV used for queue filtering and reviewer evidence.",
    )
    parser.add_argument(
        "--family-qa-csv",
        type=Path,
        default=None,
        help="Optional frozen representation family QA representative queue CSV.",
    )
    parser.add_argument(
        "--family-qa-manifest",
        type=Path,
        default=None,
        help="Pending family QA manifest updated by family decisions (never the review DB).",
    )
    return parser.parse_args()


def _connect(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def _row_dict(
    row: sqlite3.Row | Mapping[str, Any] | None,
    *,
    fixture: dict[str, Any] | None = None,
    triage: dict[str, Any] | None = None,
) -> dict[str, Any] | None:
    if row is None:
        return None
    payload = dict(row)
    payload["fixture"] = fixture
    payload["triage"] = triage
    if triage:
        payload["triage_bucket"] = triage.get("triage_bucket", "")
    metadata_json = payload.pop("metadata_json", None)
    if isinstance(metadata_json, str) and metadata_json:
        try:
            raw_payload = json.loads(metadata_json)
        except json.JSONDecodeError:
            raw_payload = None
        if isinstance(raw_payload, dict):
            for key, value in raw_payload.items():
                if key not in payload:
                    payload[key] = value

    diagnostic_group, diagnostic_reason = classify_reference_problem(
        reference_smiles=payload.get("reference_smiles"),
        skip_reasons=comparison_skip_reasons(payload.get("error")),
        formula_status=payload.get("reference_formula_check_status"),
    )
    payload["reference_diagnostic_group"] = diagnostic_group
    payload["reference_diagnostic_reason"] = diagnostic_reason
    if triage:
        payload["triage"] = {
            **triage,
            "reference_diagnostic_group": diagnostic_group,
            "reference_diagnostic_reason": diagnostic_reason,
        }

    runtime_label = f"py{sys.version_info.major}{sys.version_info.minor}"
    snapshot_smiles = ""
    snapshot_status = ""
    snapshot_found = False
    for method_id in ("candidate_cpp", "molgr_cpp"):
        smiles_key = f"{runtime_label}_{method_id}_smiles"
        status_key = f"{runtime_label}_{method_id}_status"
        if smiles_key in payload or status_key in payload:
            snapshot_smiles = str(payload.get(smiles_key) or "")
            snapshot_status = str(payload.get(status_key) or "")
            snapshot_found = True
            break
    payload["candidate_snapshot_runtime"] = runtime_label
    if snapshot_found:
        payload["candidate_snapshot_smiles"] = snapshot_smiles
        payload["candidate_snapshot_status"] = snapshot_status
    else:
        payload["candidate_snapshot_smiles"] = str(payload.get("candidate_smiles") or "")
        payload["candidate_snapshot_status"] = str(payload.get("candidate_status") or "")
    candidate_smiles = str(payload.get("candidate_smiles") or "").strip()
    reference_smiles = str(payload.get("reference_smiles") or "").strip()
    candidate_organic_smiles = str(payload.get("candidate_organic_smiles") or "").strip()
    reference_organic_smiles = str(payload.get("reference_organic_smiles") or "").strip()
    available_render_kinds = ["candidate"]
    if reference_smiles:
        available_render_kinds.append("reference")
    if candidate_organic_smiles and candidate_organic_smiles != candidate_smiles:
        available_render_kinds.append("candidate_organic")
    if reference_organic_smiles and reference_organic_smiles != reference_smiles:
        available_render_kinds.append("reference_organic")
    payload["available_render_kinds"] = available_render_kinds
    return payload


def _json_response(
    handler: BaseHTTPRequestHandler,
    payload: Any,
    *,
    status: HTTPStatus = HTTPStatus.OK,
) -> None:
    body = json.dumps(payload, ensure_ascii=False, sort_keys=True).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def _text_response(
    handler: BaseHTTPRequestHandler,
    text: str,
    *,
    status: HTTPStatus = HTTPStatus.OK,
    content_type: str = "text/plain; charset=utf-8",
) -> None:
    body = text.encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", content_type)
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def _trace_error_page(case_id: str, error: Exception) -> str:
    escaped_case = html.escape(case_id)
    escaped_error = html.escape(f"{type(error).__name__}: {error}")
    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><title>Molecule Trace - {escaped_case}</title>
<style>body{{font-family:system-ui,sans-serif;margin:2rem;color:#172033}}code{{white-space:pre-wrap}}</style>
</head><body><h1>Trace 生成失败: {escaped_case}</h1><code>{escaped_error}</code></body></html>"""


def _svg_fragment_from_image(image: Any) -> str:
    data = getattr(image, "data", None)
    if isinstance(data, str):
        return data
    repr_svg = getattr(image, "_repr_svg_", None)
    if callable(repr_svg):
        rendered = repr_svg()
        if isinstance(rendered, str):
            return rendered
    if isinstance(image, str):
        return image
    return str(image)


def _safe_smiles(mol: Chem.Mol) -> str:
    clone = Chem.RemoveHs(Chem.Mol(mol), sanitize=False)
    try:
        return Chem.MolToSmiles(clone, canonical=True, isomericSmiles=True)
    except Chem.KekulizeException:
        rw_mol = Chem.RWMol(clone)
        for atom in rw_mol.GetAtoms():  # pyright: ignore[reportCallIssue]
            atom.SetIsAromatic(False)
        for bond in rw_mol.GetBonds():  # pyright: ignore[reportCallIssue]
            bond.SetIsAromatic(False)
        fallback = rw_mol.GetMol()
        fallback.UpdatePropertyCache(strict=False)
        return Chem.MolToSmiles(fallback, canonical=True, isomericSmiles=True)


def _benchmark_candidate_mol(mol: Chem.Mol) -> Chem.Mol:
    """Apply the same successful-result postprocessing as the tmQMg benchmark."""

    # Candidate aromatic flags may be inconsistent after reconstruction edits;
    # RemoveHs must not trigger a fresh sanitization/Kekule pass here.
    return Chem.RemoveHs(Chem.Mol(mol), sanitize=False)


REVIEW_2D_RENDERER_VERSION = "review_2d_v2"


def _render_mol_from_smiles(smiles: str) -> Chem.Mol | None:
    if not smiles.strip():
        return None
    return Chem.MolFromSmiles(smiles)


def _mol_from_smiles_without_sanitize(smiles: str) -> Chem.Mol | None:
    if not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is not None:
        mol.UpdatePropertyCache(strict=False)
    return mol


def _atom_h_count(atom: Chem.Atom, *, implicit: bool) -> int | None:
    try:
        return int(atom.GetNumImplicitHs() if implicit else atom.GetNumExplicitHs())
    except (RuntimeError, ValueError):
        return None


def _review_graph_payload(mol: Chem.Mol, *, kind: str, smiles: str) -> dict[str, Any]:
    graph = Chem.Mol(mol)
    graph.UpdatePropertyCache(strict=False)
    atoms: list[dict[str, Any]] = []
    metal_atoms: list[dict[str, Any]] = []
    for atom in graph.GetAtoms():  # pyright: ignore[reportCallIssue]
        atomic_num = int(atom.GetAtomicNum())
        is_metal = atomic_num not in NON_METAL_DICT
        atom_payload = {
            "index": int(atom.GetIdx()),
            "element": atom.GetSymbol(),
            "formal_charge": int(atom.GetFormalCharge()),
            "radical_electrons": int(atom.GetNumRadicalElectrons()),
            "explicit_h": _atom_h_count(atom, implicit=False),
            "implicit_h": _atom_h_count(atom, implicit=True),
            "is_metal": is_metal,
            "neighbours": [
                {"index": int(neighbour.GetIdx()), "element": neighbour.GetSymbol()}
                for neighbour in atom.GetNeighbors()  # pyright: ignore[reportCallIssue]
            ],
        }
        atoms.append(atom_payload)
        if is_metal:
            metal_atoms.append(
                {
                    "index": atom_payload["index"],
                    "element": atom_payload["element"],
                    "formal_charge": atom_payload["formal_charge"],
                }
            )

    dative_types = {
        Chem.BondType.DATIVE,
        Chem.BondType.DATIVEL,
        Chem.BondType.DATIVER,
        Chem.BondType.DATIVEONE,
    }
    bond_names = {
        Chem.BondType.SINGLE: "single",
        Chem.BondType.DOUBLE: "double",
        Chem.BondType.TRIPLE: "triple",
        Chem.BondType.AROMATIC: "aromatic",
    }
    bonds: list[dict[str, Any]] = []
    for bond in graph.GetBonds():  # pyright: ignore[reportCallIssue]
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        bond_type = bond.GetBondType()
        kind_name = "dative" if bond_type in dative_types else bond_names.get(bond_type)
        if kind_name is None:
            kind_name = str(bond_type).lower()
        bonds.append(
            {
                "index": int(bond.GetIdx()),
                "begin_atom": int(begin.GetIdx()),
                "begin_element": begin.GetSymbol(),
                "end_atom": int(end.GetIdx()),
                "end_element": end.GetSymbol(),
                "type": kind_name,
                "directional": bond_type in dative_types,
            }
        )

    explicit_h_count = sum(atom.GetAtomicNum() == 1 for atom in graph.GetAtoms())
    explicit_h_count += sum(
        value
        for atom in graph.GetAtoms()
        if atom.GetAtomicNum() != 1
        for value in [_atom_h_count(atom, implicit=False)]
        if value is not None
    )
    return {
        "kind": kind,
        "smiles": smiles,
        "summary": {
            "total_formal_charge": sum(atom.GetFormalCharge() for atom in graph.GetAtoms()),
            "atom_count": graph.GetNumAtoms(),
            "explicit_h_count": explicit_h_count,
            "total_radical_electrons": sum(
                atom.GetNumRadicalElectrons() for atom in graph.GetAtoms()
            ),
            "metals": metal_atoms,
        },
        "atoms": atoms,
        "bonds": bonds,
    }


def _live_candidate_comparison(mol: Chem.Mol, snapshot_smiles: str) -> dict[str, Any]:
    live_smiles = _safe_smiles(mol)
    comparison: dict[str, Any] = {
        "candidate_snapshot_smiles": snapshot_smiles,
        "live_candidate_smiles_exact_match": (
            live_smiles == snapshot_smiles if snapshot_smiles else None
        ),
        "live_matches_candidate_snapshot": None,
        "live_candidate_equivalence_method": "",
        "live_candidate_equivalence_reason": "",
    }
    if not snapshot_smiles:
        comparison["live_candidate_equivalence_reason"] = "candidate_snapshot_smiles_unavailable"
        return comparison

    snapshot = _mol_from_smiles_without_sanitize(snapshot_smiles)
    if snapshot is None:
        comparison["live_candidate_equivalence_reason"] = "candidate_snapshot_smiles_invalid"
        return comparison

    equivalent, info = check_equivalence(mol, snapshot, use_chirality=False)
    comparison["live_matches_candidate_snapshot"] = equivalent
    comparison["live_candidate_equivalence_method"] = info.method.value if info.method else ""
    comparison["live_candidate_equivalence_reason"] = info.reason
    return comparison


def _prepare_review_2d_mol(mol: Chem.Mol, *, show_hydrogens: bool = False) -> Chem.Mol:
    """Prepare a display-only copy without changing the source graph."""

    copied = Chem.Mol(mol)
    if show_hydrogens:
        # Candidate reconstruction and the verified reference-XYZ display graph
        # already contain every source hydrogen as a real atom.  Render that
        # graph unchanged: AddHs would create atoms with no source-XYZ identity.
        copied.UpdatePropertyCache(strict=False)
        return copied
    return Chem.RemoveHs(copied, sanitize=False)


def _render_mol_svg(
    mol: Chem.Mol,
    *,
    legend: str,
    atom_notes: dict[int, str] | None = None,
    show_hydrogens: bool = False,
) -> str:
    from rdkit_dof import MolToDofImage

    annotated = Chem.Mol(mol)
    for atom_index, note in (atom_notes or {}).items():
        if 0 <= atom_index < annotated.GetNumAtoms():
            atom = annotated.GetAtomWithIdx(atom_index)
            atom.SetProp("atomNote", note)
            atom.SetBoolProp("_review_disputed_atom", True)
    render_mol = _prepare_review_2d_mol(annotated, show_hydrogens=show_hydrogens)
    highlight_atoms = [
        atom.GetIdx() for atom in render_mol.GetAtoms() if atom.HasProp("_review_disputed_atom")
    ]
    image = MolToDofImage(
        render_mol,
        size=(520, 360),
        legend=legend,
        use_svg=True,
        return_image=False,
        highlightAtoms=highlight_atoms,
        highlightColor=(0.10, 0.55, 0.64, 0.55),
    )
    return _svg_fragment_from_image(image)


def _render_cache_kind(kind: str, mode: str = "skeleton") -> str:
    return f"{REVIEW_2D_RENDERER_VERSION}:{kind}:{mode}"


def _triage_atom_notes(
    case: Mapping[str, Any],
    triage: dict[str, str] | None,
    candidate: Chem.Mol,
    reference: Chem.Mol | None,
    xyz_dir: Path | None,
) -> tuple[dict[int, str], dict[int, str]]:
    if not triage:
        return {}, {}
    candidate_notes: dict[int, list[str]] = {}
    reference_xyz_notes: dict[int, list[str]] = {}

    def add_note(target: dict[int, list[str]], index: Any, text: str) -> None:
        try:
            atom_index = int(index)
        except (TypeError, ValueError):
            return
        target.setdefault(atom_index, []).append(text)

    def atom_label(mol: Chem.Mol, index: Any) -> str:
        try:
            atom_index = int(index)
        except (TypeError, ValueError):
            return ""
        if not 0 <= atom_index < mol.GetNumAtoms():
            return ""
        return f"{mol.GetAtomWithIdx(atom_index).GetSymbol()} · #{atom_index}"

    for edge in _json_array(triage.get("metal_coordination_diff")):
        for index in edge.get("candidate_atoms", []):
            add_note(candidate_notes, index, atom_label(candidate, index))
            add_note(reference_xyz_notes, index, atom_label(candidate, index))
    for assignment in _json_array(triage.get("hydrogen_assignment_diff")):
        hydrogen = assignment.get("h_atom")
        add_note(candidate_notes, hydrogen, f"H · #{hydrogen}")
        candidate_center = assignment.get("candidate_center")
        reference_center = assignment.get("reference_center")
        add_note(candidate_notes, candidate_center, atom_label(candidate, candidate_center))
        add_note(reference_xyz_notes, hydrogen, f"H · #{hydrogen}")
        add_note(reference_xyz_notes, reference_center, atom_label(candidate, reference_center))
    candidate_result = {
        index: " · ".join(dict.fromkeys(notes)) for index, notes in candidate_notes.items()
    }
    if reference is None or not candidate_result:
        return candidate_result, {}
    try:
        xyz_path = resolve_xyz_path(str(case.get("xyz_path") or ""), xyz_dir)
        xyz_atoms = parse_xyz_atoms(xyz_path.read_text(encoding="utf-8"))
        mapping = map_candidate_reference_xyz(candidate, reference, xyz_atoms)
        if mapping.confidence not in {"exact", "unique_graph_mapping"}:
            return candidate_result, {}
        xyz_to_reference = {
            candidate_index: reference_index
            for reference_index, candidate_index in mapping.reference_to_candidate.items()
        }
        reference_result = {
            xyz_to_reference[xyz_index]: note
            for xyz_index, notes in reference_xyz_notes.items()
            if xyz_index in xyz_to_reference
            for note in [" · ".join(dict.fromkeys(notes))]
        }
        return candidate_result, reference_result
    except Exception:  # noqa: BLE001
        return candidate_result, {}


def _json_array(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, list):
        return [item for item in value if isinstance(item, dict)]
    if not value:
        return []
    try:
        parsed = json.loads(str(value))
    except (TypeError, ValueError, json.JSONDecodeError):
        return []
    return [item for item in parsed if isinstance(item, dict)] if isinstance(parsed, list) else []


def _triage_source_atom_notes(
    triage: dict[str, str] | None,
    source_mol: Chem.Mol,
    *,
    reference_assignment: bool,
) -> dict[int, str]:
    """Labels keyed by authoritative source XYZ atom index."""

    if not triage:
        return {}
    notes: dict[int, list[str]] = {}

    def add(index: Any, label: str) -> None:
        try:
            atom_index = int(index)
        except (TypeError, ValueError):
            return
        if not 0 <= atom_index < source_mol.GetNumAtoms():
            return
        notes.setdefault(atom_index, []).append(label)

    def label(index: Any) -> str:
        try:
            atom_index = int(index)
        except (TypeError, ValueError):
            return ""
        if not 0 <= atom_index < source_mol.GetNumAtoms():
            return ""
        return f"{source_mol.GetAtomWithIdx(atom_index).GetSymbol()} · #{atom_index}"

    for edge in _json_array(triage.get("metal_coordination_diff")):
        for index in edge.get("candidate_atoms", []):
            add(index, label(index))
    for assignment in _json_array(triage.get("hydrogen_assignment_diff")):
        hydrogen = assignment.get("h_atom")
        center_key = "reference_center" if reference_assignment else "candidate_center"
        center = assignment.get(center_key)
        add(hydrogen, f"H · #{hydrogen}")
        add(center, label(center))
    return {index: " · ".join(dict.fromkeys(values)) for index, values in notes.items()}


def _reference_xyz_mol(
    case: Mapping[str, Any],
    triage: dict[str, str] | None,
    candidate: Chem.Mol,
    reference: Chem.Mol,
    xyz_dir: Path | None,
    mapping_result: Any | None = None,
) -> tuple[Chem.Mol | None, str, str]:
    """Place reference connectivity on source XYZ using a valid mapping representative."""

    xyz_path = resolve_xyz_path(str(case.get("xyz_path") or ""), xyz_dir)
    xyz_atoms = parse_xyz_atoms(xyz_path.read_text(encoding="utf-8"))
    mapping = mapping_result or map_candidate_reference_xyz(candidate, reference, xyz_atoms)
    if mapping.confidence not in {"exact", "unique_graph_mapping", "ambiguous"}:
        return None, mapping.confidence, "atom_correspondence_not_reliable"
    if not mapping.reference_to_candidate:
        return None, mapping.confidence, "atom_correspondence_not_reliable"
    if candidate.GetNumAtoms() != len(xyz_atoms) or candidate.GetNumConformers() == 0:
        return None, mapping.confidence, "candidate_xyz_not_complete"

    reference_to_xyz = dict(mapping.reference_to_candidate)
    reference_copy = Chem.Mol(reference)
    reference_copy.UpdatePropertyCache(strict=False)
    candidate_copy = Chem.Mol(candidate)
    candidate_copy.UpdatePropertyCache(strict=False)
    editable = Chem.RWMol()
    for xyz_index, (symbol, _) in enumerate(xyz_atoms):
        source_atom = candidate_copy.GetAtomWithIdx(xyz_index)
        atom = Chem.Atom(symbol)
        atom.SetFormalCharge(source_atom.GetFormalCharge())
        atom.SetNumRadicalElectrons(source_atom.GetNumRadicalElectrons())
        atom.SetNoImplicit(True)
        editable.AddAtom(atom)

    for reference_atom in reference_copy.GetAtoms():
        if reference_atom.GetAtomicNum() == 1:
            continue
        xyz_index = reference_to_xyz.get(reference_atom.GetIdx())
        if xyz_index is None:
            return None, mapping.confidence, "reference_mapping_incomplete"
        target = editable.GetAtomWithIdx(xyz_index)
        target.SetFormalCharge(reference_atom.GetFormalCharge())
        target.SetNumRadicalElectrons(reference_atom.GetNumRadicalElectrons())
        target.SetIsAromatic(reference_atom.GetIsAromatic())

    for bond in reference_copy.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if (
            reference_copy.GetAtomWithIdx(begin).GetAtomicNum() == 1
            or reference_copy.GetAtomWithIdx(end).GetAtomicNum() == 1
        ):
            continue
        mapped_begin = reference_to_xyz.get(begin)
        mapped_end = reference_to_xyz.get(end)
        if mapped_begin is None or mapped_end is None:
            return None, mapping.confidence, "reference_bond_mapping_incomplete"
        editable.AddBond(mapped_begin, mapped_end, bond.GetBondType())
        editable.GetBondBetweenAtoms(mapped_begin, mapped_end).SetIsAromatic(bond.GetIsAromatic())

    hydrogen_centers: dict[int, int] = {}
    for atom in candidate_copy.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        heavy_neighbours = [
            neighbour.GetIdx() for neighbour in atom.GetNeighbors() if neighbour.GetAtomicNum() != 1
        ]
        if len(heavy_neighbours) == 1:
            hydrogen_centers[atom.GetIdx()] = heavy_neighbours[0]
    for assignment in _json_array((triage or {}).get("hydrogen_assignment_diff")):
        try:
            hydrogen_centers[int(assignment["h_atom"])] = int(assignment["reference_center"])
        except (KeyError, TypeError, ValueError):
            return None, mapping.confidence, "invalid_hydrogen_assignment"

    expected_h_counts = {
        reference_to_xyz[atom.GetIdx()]: atom_h_count(atom)
        for atom in reference_copy.GetAtoms()
        if atom.GetAtomicNum() != 1 and atom.GetIdx() in reference_to_xyz
    }
    actual_h_counts = Counter(hydrogen_centers.values())
    if any(actual_h_counts.get(index, 0) != count for index, count in expected_h_counts.items()):
        return None, mapping.confidence, "reference_hydrogen_mapping_incomplete"
    for hydrogen, center in hydrogen_centers.items():
        editable.AddBond(hydrogen, center, Chem.BondType.SINGLE)

    result = editable.GetMol()
    conformer = Chem.Conformer(len(xyz_atoms))
    source_conformer = candidate_copy.GetConformer()
    for index in range(len(xyz_atoms)):
        conformer.SetAtomPosition(index, source_conformer.GetAtomPosition(index))
    result.AddConformer(conformer)
    result.UpdatePropertyCache(strict=False)
    return result, mapping.confidence, ""


def _mapped_coordination_comparison(
    candidate: Chem.Mol,
    reference: Chem.Mol,
    reference_to_xyz: Mapping[int, int],
) -> dict[str, Any]:
    """Describe mapped metal-donor edges without changing either molecular graph."""

    xyz_to_reference = {xyz: reference_index for reference_index, xyz in reference_to_xyz.items()}

    def coordination_edges(
        mol: Chem.Mol, atom_map: Mapping[int, int] | None = None
    ) -> dict[tuple[int, int], tuple[int, int]]:
        result: dict[tuple[int, int], tuple[int, int]] = {}
        for bond in mol.GetBonds():
            begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            begin_metal = is_metal(mol.GetAtomWithIdx(begin))
            end_metal = is_metal(mol.GetAtomWithIdx(end))
            if begin_metal == end_metal:
                continue
            metal, donor = (begin, end) if begin_metal else (end, begin)
            mapped_metal = atom_map.get(metal) if atom_map is not None else metal
            mapped_donor = atom_map.get(donor) if atom_map is not None else donor
            if mapped_metal is None or mapped_donor is None:
                continue
            result[(mapped_metal, mapped_donor)] = (metal, donor)
        return result

    candidate_edges = coordination_edges(candidate)
    reference_edges = coordination_edges(reference, reference_to_xyz)
    conformer = candidate.GetConformer() if candidate.GetNumConformers() else None

    def mapped_ligand_group(donor_xyz: int) -> list[dict[str, Any]]:
        donor = candidate.GetAtomWithIdx(donor_xyz)
        group = {donor_xyz}
        for center in donor.GetNeighbors():
            if is_metal(center) or center.GetAtomicNum() == 1:
                continue
            group.update(
                neighbour.GetIdx()
                for neighbour in center.GetNeighbors()
                if neighbour.GetAtomicNum() == donor.GetAtomicNum() and not is_metal(neighbour)
            )
        return [
            {
                "element": candidate.GetAtomWithIdx(xyz_index).GetSymbol(),
                "candidate_xyz_index": xyz_index,
                "reference_atom_index": xyz_to_reference.get(xyz_index),
                "role": "donor" if xyz_index == donor_xyz else "mapped_equivalent",
            }
            for xyz_index in sorted(group)
            if xyz_index in xyz_to_reference
        ]

    edges = []
    for metal_xyz, donor_xyz in sorted(set(candidate_edges) | set(reference_edges)):
        candidate_source = candidate_edges.get((metal_xyz, donor_xyz))
        reference_source = reference_edges.get((metal_xyz, donor_xyz))
        presence = (
            "common"
            if candidate_source is not None and reference_source is not None
            else "candidate_only"
            if candidate_source is not None
            else "reference_only"
        )
        distance = None
        if conformer is not None:
            distance = float(
                conformer.GetAtomPosition(metal_xyz).Distance(
                    conformer.GetAtomPosition(donor_xyz)
                )
            )
        group = mapped_ligand_group(donor_xyz)
        edges.append(
            {
                "presence": presence,
                "metal_element": candidate.GetAtomWithIdx(metal_xyz).GetSymbol(),
                "donor_element": candidate.GetAtomWithIdx(donor_xyz).GetSymbol(),
                "candidate_metal_xyz_index": metal_xyz,
                "candidate_donor_xyz_index": donor_xyz,
                "reference_metal_atom_index": (
                    reference_source[0] if reference_source is not None else xyz_to_reference.get(metal_xyz)
                ),
                "reference_donor_atom_index": (
                    reference_source[1] if reference_source is not None else xyz_to_reference.get(donor_xyz)
                ),
                "distance": distance,
                "mapped_ligand_group": group if len(group) > 1 else [],
                "interpretation": (
                    "mapped_donor_preserved" if presence == "common" else "coordination_edge_missing"
                ),
            }
        )
    return {
        "candidate_to_xyz": {str(index): index for index in range(candidate.GetNumAtoms())},
        "reference_to_xyz": {str(key): value for key, value in reference_to_xyz.items()},
        "coordination_edges": edges,
    }


def _mapping_ambiguity_details(mapping: Any) -> tuple[str, str]:
    if getattr(mapping, "confidence", "") != "ambiguous":
        return "", ""
    if getattr(mapping, "timeout", False):
        return "mapping_timeout", "mapping_timeout_before_unique_correspondence"
    if getattr(mapping, "enumeration_truncated", False):
        return "mapping_enumeration_truncated", "mapping_enumeration_truncated_before_unique_correspondence"
    if getattr(mapping, "mapping_signature_count", 0) > 1:
        return "multiple_valid_mappings", "multiple_equally_valid_atom_mappings"
    if getattr(mapping, "equal_best_mapping_count", 0) > 1:
        return "symmetry_equivalent_atoms", "symmetry_equivalent_atoms_prevent_unique_correspondence"
    return "ambiguous_mapping", "unique_atom_correspondence_not_established"


def _mapping_ambiguity_locations(mapping: Any) -> dict[str, Any]:
    signatures = getattr(mapping, "decision_relevant_signatures", ())
    alternatives = []
    affected_atoms: set[int] = set()
    for alternative_index, signature in enumerate(signatures, start=1):
        differences = []
        for difference in signature:
            if not difference:
                continue
            kind = str(difference[0])
            if kind in {"metal_bond", "organic_bond"}:
                atoms = [int(index) for index in difference[1]]
                affected_atoms.update(atoms)
                differences.append(
                    {
                        "kind": kind,
                        "xyz_atoms": atoms,
                        "candidate_present": bool(difference[2]),
                        "reference_present": bool(difference[3]),
                        "candidate_bond": str(difference[4]),
                        "reference_bond": str(difference[5]),
                    }
                )
            elif kind == "hydrogen_assignment":
                hydrogen_atoms = [int(index) for index in difference[1]]
                candidate_center = difference[2]
                reference_centers = [int(index) for index in difference[3]]
                affected_atoms.update(hydrogen_atoms)
                affected_atoms.update(reference_centers)
                if candidate_center is not None:
                    affected_atoms.add(int(candidate_center))
                differences.append(
                    {
                        "kind": kind,
                        "hydrogen_xyz_atoms": hydrogen_atoms,
                        "candidate_center_xyz": candidate_center,
                        "reference_center_xyz": reference_centers,
                    }
                )
        alternatives.append({"alternative": alternative_index, "differences": differences})
    return {
        "affected_xyz_atoms": sorted(affected_atoms),
        "alternatives": alternatives,
        "location_proven": (
            not getattr(mapping, "enumeration_truncated", False)
            and not getattr(mapping, "timeout", False)
            and len(signatures) > 1
        ),
    }


def _mol_to_sdf_block(mol: Chem.Mol) -> str:
    # Candidate graphs may carry aromatic flags that are inconsistent after
    # charge/bond edits. Preserve the explicit graph without forcing Kekule.
    block = Chem.MolToMolBlock(mol, kekulize=False, forceV3000=True)
    if not block.endswith("\n"):
        block += "\n"
    return block + "$$$$\n"


def _mol_from_sdf_block(sdf: str) -> Chem.Mol:
    molblock = sdf.partition("$$$$")[0]
    mol = Chem.MolFromMolBlock(
        molblock,
        sanitize=False,
        removeHs=False,
        strictParsing=False,
    )
    if mol is None:
        raise ValueError("invalid_sdf")
    mol.UpdatePropertyCache(strict=False)
    return mol


def _dof_size(payload: dict[str, Any], key: str, default: tuple[int, int]) -> tuple[int, int]:
    raw = payload.get(key)
    if raw is None:
        return default
    if not isinstance(raw, list) or len(raw) != 2:
        raise ValueError(f"invalid_{key}")
    width, height = int(raw[0]), int(raw[1])
    if not (1 <= width <= 2000 and 1 <= height <= 2000):
        raise ValueError(f"invalid_{key}")
    return width, height


def _render_deferred_dof(payload: dict[str, Any]) -> str:
    render_type = str(payload.get("render_type") or "")
    legends = payload.get("legends") or []
    if not isinstance(legends, list) or not all(isinstance(item, str) for item in legends):
        raise ValueError("invalid_legends")

    if render_type == "single":
        sdf = payload.get("sdf")
        if not isinstance(sdf, str) or not sdf:
            raise ValueError("missing_sdf")
        from rdkit_dof import MolToDofImage

        image = MolToDofImage(
            _mol_from_sdf_block(sdf),
            size=_dof_size(payload, "size", (360, 300)),
            legend=legends[0] if legends else "",
            use_svg=True,
            return_image=False,
        )
    else:
        sdfs = payload.get("sdfs")
        if (
            not isinstance(sdfs, list)
            or not sdfs
            or not all(isinstance(item, str) for item in sdfs)
        ):
            raise ValueError("missing_sdfs")
        mols = [_mol_from_sdf_block(sdf) for sdf in sdfs]
        if legends and len(legends) != len(mols):
            raise ValueError("legend_count_mismatch")
        if render_type == "grid":
            from rdkit_dof import MolsToGridDofImage

            mols_per_row = int(payload.get("mols_per_row") or 3)
            if not 1 <= mols_per_row <= 12:
                raise ValueError("invalid_mols_per_row")
            image = MolsToGridDofImage(
                mols,
                molsPerRow=mols_per_row,
                subImgSize=_dof_size(payload, "sub_image_size", (320, 260)),
                legends=legends,
                use_svg=True,
                return_image=False,
            )
        elif render_type == "animation":
            from rdkit_dof import MolsToDofSvgAnimation

            duration = int(payload.get("duration") or 650)
            if not 50 <= duration <= 10_000:
                raise ValueError("invalid_duration")
            image = MolsToDofSvgAnimation(
                mols,
                size=_dof_size(payload, "size", (360, 300)),
                legends=legends,
                duration=duration,
                loop=0,
                return_image=False,
            )
        else:
            raise ValueError("invalid_render_type")
    svg = _svg_fragment_from_image(image)
    svg_start = svg.find("<svg")
    return svg[svg_start:] if svg_start >= 0 else svg


def _case_query(select_extra: str = "") -> str:
    return f"""
        SELECT
            c.case_id, c.row_index, c.source, c.category, c.xyz_path,
            c.total_charge, c.total_radical_electrons, c.spin_multiplicity,
            c.reference_smiles, c.candidate_smiles,
            c.candidate_organic_smiles, c.reference_organic_smiles,
            c.candidate_status, c.metadata_json,
            r.status AS review_status, r.corrected_smiles, r.corrected_molblock,
            r.notes, r.reviewer, r.updated_at
            {select_extra}
        FROM cases c
        LEFT JOIN reviews r ON r.case_id = c.case_id
    """


class ReviewServer(ThreadingHTTPServer):
    db_path: Path
    xyz_dir: Path | None
    fixtures_dir: Path
    runtime_info: dict[str, Any]
    triage_records: dict[str, dict[str, str]]
    review_histories: dict[str, list[dict[str, Any]]]
    review_mutation_lock: threading.RLock
    family_qa: dict[str, Any] | None
    family_qa_manifest_path: Path | None

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.review_histories = {}
        self.review_mutation_lock = threading.RLock()


def load_triage_records(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            row["case_id"]: row
            for row in csv.DictReader(handle)
            if str(row.get("case_id") or "").strip()
        }


def load_family_qa(path: Path | None, manifest_path: Path | None) -> dict[str, Any] | None:
    if path is None:
        if manifest_path is not None:
            raise ValueError("--family-qa-manifest requires --family-qa-csv")
        return None
    if not path.is_file():
        raise FileNotFoundError(path)
    if manifest_path is None:
        raise ValueError("--family-qa-manifest is required with --family-qa-csv")
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    with path.open(encoding="utf-8", newline="") as handle:
        representatives = [dict(row) for row in csv.DictReader(handle)]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict) or not isinstance(manifest.get("families"), list):
        raise ValueError("invalid family QA manifest")
    manifest_families = {
        str(family.get("family_id") or ""): family
        for family in manifest["families"]
        if isinstance(family, dict)
    }
    families: dict[str, dict[str, Any]] = {}
    for row in representatives:
        family_id = str(row.get("family_id") or "").strip()
        case_id = str(row.get("case_id") or "").strip()
        if not family_id or not case_id or family_id not in manifest_families:
            raise ValueError(f"invalid frozen family QA row: family={family_id!r}, case={case_id!r}")
        family = families.setdefault(
            family_id,
            {"metadata": dict(row), "representatives": []},
        )
        family["representatives"].append(row)
    for family_id, family in families.items():
        frozen = manifest_families[family_id]
        expected = {str(value) for value in frozen.get("representatives", [])}
        actual = {str(row["case_id"]) for row in family["representatives"]}
        if expected != actual:
            raise ValueError(f"frozen representative mismatch for {family_id}")
    return {"csv": path, "families": families}


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _family_manifest(server: ReviewServer) -> dict[str, Any]:
    path = server.family_qa_manifest_path
    if path is None:
        raise ValueError("family_qa_not_configured")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("families"), list):
        raise ValueError("invalid family QA manifest")
    return payload


def _public_family_qa(server: ReviewServer) -> dict[str, Any]:
    qa = server.family_qa
    if qa is None:
        return {"enabled": False, "families": [], "progress": {}}
    manifest = _family_manifest(server)
    manifest_by_id = {
        str(item.get("family_id") or ""): item
        for item in manifest["families"]
        if isinstance(item, dict)
    }
    families = []
    reviewed = approved_cases = 0
    total_cases = 0
    for family_id, source in qa["families"].items():
        pending = manifest_by_id[family_id]
        decision = str(pending.get("qa_decision") or "")
        family_size = int(source["metadata"].get("family_size") or pending.get("size") or 0)
        total_cases += family_size
        if decision:
            reviewed += 1
            if decision.startswith("approve_"):
                approved_cases += family_size
        marks = pending.get("representative_marks") or {}
        representatives = []
        for row in source["representatives"]:
            case_id = str(row["case_id"])
            representatives.append({**row, "qa_mark": str(marks.get(case_id) or "")})
        families.append(
            {
                **source["metadata"],
                "family_id": family_id,
                "family_size": family_size,
                "representatives": representatives,
                "decision": decision,
                "proposed_status": pending.get("proposed_status"),
                "proposed_reason": pending.get("proposed_reason"),
            }
        )
    return {
        "enabled": True,
        "approval_required": True,
        "approved": False,
        "families": families,
        "progress": {
            "reviewed_families": reviewed,
            "total_families": len(families),
            "approved_cases": approved_cases,
            "total_cases": total_cases,
        },
    }


REVIEW_COLUMNS = (
    "status",
    "corrected_smiles",
    "corrected_molblock",
    "notes",
    "reviewer",
    "updated_at",
)


def _review_state(row: sqlite3.Row | Mapping[str, Any] | None) -> dict[str, str] | None:
    if row is None:
        return None
    return {key: str(row[key] or "") for key in REVIEW_COLUMNS}


def _public_review_mutation(mutation: Mapping[str, Any]) -> dict[str, Any]:
    if mutation.get("mutation_type") == "family_qa":
        return {
            "mutation_id": mutation["mutation_id"],
            "mutation_type": "family_qa",
            "case_id": mutation.get("case_id", ""),
            "family_id": mutation["family_id"],
            "status": mutation["label"],
            "reviewer": "family QA",
            "notes": "",
            "timestamp": mutation["timestamp"],
            "undone": bool(mutation.get("undone")),
        }
    return {
        "mutation_id": mutation["mutation_id"],
        "case_id": mutation["case_id"],
        "status": mutation["after_review"]["status"],
        "reviewer": mutation["after_review"]["reviewer"],
        "notes": mutation["after_review"]["notes"],
        "timestamp": mutation["timestamp"],
        "undone": bool(mutation.get("undone")),
    }


class ReviewHandler(BaseHTTPRequestHandler):
    server: ReviewServer

    def log_message(self, format: str, *args: object) -> None:
        sys.stderr.write("[molgr-review] " + format % args + "\n")

    def do_GET(self) -> None:  # noqa: N802
        try:
            self._handle_get()
        except Exception as exc:  # noqa: BLE001
            traceback.print_exc()
            _json_response(
                self,
                {"error": f"{type(exc).__name__}: {exc}"},
                status=HTTPStatus.INTERNAL_SERVER_ERROR,
            )

    def do_POST(self) -> None:  # noqa: N802
        try:
            self._handle_post()
        except Exception as exc:  # noqa: BLE001
            traceback.print_exc()
            _json_response(
                self,
                {"error": f"{type(exc).__name__}: {exc}"},
                status=HTTPStatus.INTERNAL_SERVER_ERROR,
            )

    def _handle_get(self) -> None:
        parsed = urlparse(self.path)
        path = parsed.path
        query = parse_qs(parsed.query)
        if path in ("", "/"):
            self._serve_static("index.html")
        elif path.startswith("/static/"):
            self._serve_static(path[len("/static/") :])
        elif path.startswith("/KetcherDemoSA/"):
            self._proxy_ketcher(path, parsed.query)
        elif path.startswith("/ketcher/"):
            self._proxy_ketcher("/KetcherDemoSA/" + path[len("/ketcher/") :], parsed.query)
        elif path.startswith("/trace/"):
            self._trace_page(unquote(path[len("/trace/") :]))
        elif path == "/api/stats":
            self._api_stats()
        elif path == "/api/review-reasons":
            self._api_review_reasons()
        elif path == "/api/review-history":
            self._api_review_history()
        elif path == "/api/family-qa":
            self._api_family_qa()
        elif path == "/api/cases":
            self._api_cases(query)
        elif path.startswith("/api/cases/") and path.endswith("/render"):
            case_id = unquote(path.split("/")[3])
            self._api_render(case_id, query)
        elif path.startswith("/api/cases/") and path.endswith("/xyz"):
            case_id = unquote(path.split("/")[3])
            self._api_xyz(case_id)
        elif path.startswith("/api/cases/") and path.endswith("/candidate-sdf"):
            case_id = unquote(path.split("/")[3])
            self._api_candidate_sdf(case_id, query)
        elif path.startswith("/api/cases/") and path.endswith("/reference-xyz"):
            case_id = unquote(path.split("/")[3])
            self._api_reference_xyz(case_id)
        elif path.startswith("/api/cases/") and path.endswith("/graph"):
            case_id = unquote(path.split("/")[3])
            self._api_graph(case_id, query)
        elif path.startswith("/api/cases/"):
            case_id = unquote(path.split("/")[3])
            self._api_case(case_id)
        else:
            _json_response(self, {"error": "not_found"}, status=HTTPStatus.NOT_FOUND)

    def _trace_page(self, case_id: str) -> None:
        if not case_id:
            _text_response(self, "case_id is required", status=HTTPStatus.BAD_REQUEST)
            return
        with _connect(self.server.db_path) as conn:
            row = conn.execute("SELECT * FROM cases WHERE case_id = ?", (case_id,)).fetchone()
        if row is None:
            _text_response(self, "case_not_found", status=HTTPStatus.NOT_FOUND)
            return

        case = dict(row)
        xyz_path = resolve_xyz_path(str(case.get("xyz_path") or ""), self.server.xyz_dir)
        try:
            if not xyz_path.is_file():
                raise FileNotFoundError(xyz_path)
            total_charge, total_radicals, _ = case_electronic_state(case)
            fixture = load_fixture_records(self.server.fixtures_dir).get(case_id) or {}
            input_case = TraceInputCase(
                id=case_id,
                xyz_block=xyz_path.read_text(encoding="utf-8"),
                total_charge=total_charge,
                total_radical_electrons=total_radicals,
                xyz_path=xyz_path,
                xyz_source="review_page",
                fixture_kind=str(fixture.get("kind") or ""),
                fixture_structure_file=str(fixture.get("structure_file") or ""),
                expected_smiles=str(fixture.get("approved_smiles") or ""),
                expected_smiles_options=tuple(
                    str(value).strip()
                    for value in fixture.get("accepted_smiles", [])
                    if str(value).strip()
                ),
                reference_smiles=str(case.get("reference_smiles") or ""),
            )

            report = render_trace_report(
                [input_case],
                score_all_candidates=False,
                dof_max_images=None,
                defer_dof_images=True,
            )
        except Exception as exc:  # noqa: BLE001
            _text_response(
                self,
                _trace_error_page(case_id, exc),
                status=HTTPStatus.INTERNAL_SERVER_ERROR,
                content_type="text/html; charset=utf-8",
            )
            return
        _text_response(self, report, content_type="text/html; charset=utf-8")

    def _handle_post(self) -> None:
        parsed = urlparse(self.path)
        path = parsed.path
        if path == "/api/render-dof":
            self._api_render_dof()
        elif path == "/api/review-undo":
            self._api_review_undo()
        elif path == "/api/family-qa":
            self._api_family_qa_mutation()
        elif path.startswith("/api/cases/") and path.endswith("/review"):
            case_id = unquote(path.split("/")[3])
            self._api_review(case_id)
        else:
            _json_response(self, {"error": "not_found"}, status=HTTPStatus.NOT_FOUND)

    def _api_render_dof(self) -> None:
        content_length = int(self.headers.get("Content-Length") or 0)
        if not 0 < content_length <= 16 * 1024 * 1024:
            _json_response(self, {"error": "invalid_content_length"}, status=HTTPStatus.BAD_REQUEST)
            return
        try:
            payload = json.loads(self.rfile.read(content_length))
            if not isinstance(payload, dict):
                raise ValueError("payload_must_be_an_object")
            svg = _render_deferred_dof(payload)
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            _json_response(self, {"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)
            return
        _text_response(self, svg, content_type="image/svg+xml; charset=utf-8")

    def _serve_static(self, name: str) -> None:
        requested = (STATIC_DIR / name).resolve()
        if STATIC_DIR.resolve() not in requested.parents and requested != STATIC_DIR.resolve():
            _json_response(self, {"error": "invalid_static_path"}, status=HTTPStatus.BAD_REQUEST)
            return
        if not requested.exists() or not requested.is_file():
            _json_response(self, {"error": "static_not_found"}, status=HTTPStatus.NOT_FOUND)
            return
        body = requested.read_bytes()
        content_type = mimetypes.guess_type(str(requested))[0] or "application/octet-stream"
        if requested.suffix in {".html", ".js", ".css"}:
            content_type += "; charset=utf-8"
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _proxy_ketcher(self, path: str, query: str) -> None:
        if ".." in path:
            _json_response(self, {"error": "invalid_ketcher_path"}, status=HTTPStatus.BAD_REQUEST)
            return
        encoded_path = quote(path, safe="/._-")
        url = KETCHER_BASE_URL + encoded_path
        if query:
            url += "?" + query
        request = Request(
            url,
            headers={
                "User-Agent": "Molecule graph review client",
                "Accept": self.headers.get("Accept", "*/*"),
            },
        )
        with urlopen(request, timeout=30) as response:
            body = response.read()
            content_type = response.headers.get("Content-Type") or mimetypes.guess_type(path)[0]
            content_type = content_type or "application/octet-stream"
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "public, max-age=3600")
            self.end_headers()
            self.wfile.write(body)

    def _read_json_body(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0"))
        body = self.rfile.read(length)
        if not body:
            return {}
        value = json.loads(body.decode("utf-8"))
        if not isinstance(value, dict):
            raise ValueError("JSON request body must be an object")
        return value

    def _api_stats(self) -> None:
        with _connect(self.server.db_path) as conn:
            category_rows = conn.execute(
                "SELECT category, COUNT(*) AS count FROM cases GROUP BY category"
            ).fetchall()
            status_rows = conn.execute(
                """
                SELECT COALESCE(r.status, 'unreviewed') AS status, COUNT(*) AS count
                FROM cases c
                LEFT JOIN reviews r ON r.case_id = c.case_id
                GROUP BY COALESCE(r.status, 'unreviewed')
                """
            ).fetchall()
            metadata_rows = conn.execute("SELECT key, value FROM metadata").fetchall()
        _json_response(
            self,
            {
                "categories": {row["category"]: row["count"] for row in category_rows},
                "review_statuses": {row["status"]: row["count"] for row in status_rows},
                "metadata": {row["key"]: row["value"] for row in metadata_rows},
                "runtime": getattr(self.server, "runtime_info", {}),
                "storage": {
                    "review_db": str(self.server.db_path.resolve()),
                    "fixtures_dir": str(self.server.fixtures_dir.resolve()),
                },
                "triage_buckets": dict(
                    sorted(
                        Counter(
                            row.get("triage_bucket", "")
                            for row in getattr(self.server, "triage_records", {}).values()
                            if row.get("triage_bucket")
                        ).items()
                    )
                ),
            },
        )

    def _api_review_reasons(self) -> None:
        with _connect(self.server.db_path) as conn:
            rows = conn.execute(
                """
                SELECT reviewer, COUNT(*) AS count
                FROM reviews
                WHERE reviewer IS NOT NULL AND TRIM(reviewer) != ''
                GROUP BY reviewer
                ORDER BY COUNT(*) DESC, reviewer
                """
            ).fetchall()
        _json_response(
            self,
            {"items": [{"reviewer": row["reviewer"], "count": row["count"]} for row in rows]},
        )

    def _review_session_id(self) -> str:
        session_id = self.headers.get("X-MolGR-Review-Session", "").strip()
        return session_id if 0 < len(session_id) <= 128 else ""

    def _api_review_history(self) -> None:
        session_id = self._review_session_id()
        if not session_id:
            _json_response(self, {"error": "review_session_required"}, status=HTTPStatus.BAD_REQUEST)
            return
        with self.server.review_mutation_lock:
            history = self.server.review_histories.get(session_id, [])
            _json_response(
                self,
                {"items": [_public_review_mutation(item) for item in reversed(history)]},
            )

    def _api_family_qa(self) -> None:
        with self.server.review_mutation_lock:
            payload = _public_family_qa(self.server)
        if payload.get("enabled"):
            case_ids = [
                str(rep["case_id"])
                for family in payload["families"]
                for rep in family["representatives"]
            ]
            case_rows: dict[str, dict[str, Any]] = {}
            if case_ids:
                placeholders = ",".join("?" for _ in case_ids)
                with _connect(self.server.db_path) as conn:
                    rows = conn.execute(
                        _case_query() + f" WHERE c.case_id IN ({placeholders})", case_ids
                    ).fetchall()
                fixtures = load_fixture_records(self.server.fixtures_dir)
                case_rows = {
                    str(row["case_id"]): _row_dict(
                        row,
                        fixture=fixtures.get(str(row["case_id"])),
                        triage=self.server.triage_records.get(str(row["case_id"])),
                    )
                    for row in rows
                }
            missing = []
            for family in payload["families"]:
                for rep in family["representatives"]:
                    case_id = str(rep["case_id"])
                    rep["case"] = case_rows.get(case_id)
                    rep["triage"] = self.server.triage_records.get(case_id)
                    if rep["case"] is None or rep["triage"] is None:
                        missing.append(case_id)
            payload["missing_join_cases"] = sorted(set(missing))
        _json_response(self, payload)

    def _api_family_qa_mutation(self) -> None:
        if self.server.family_qa is None or self.server.family_qa_manifest_path is None:
            _json_response(self, {"error": "family_qa_not_configured"}, status=HTTPStatus.NOT_FOUND)
            return
        session_id = self._review_session_id()
        if not session_id:
            _json_response(self, {"error": "review_session_required"}, status=HTTPStatus.BAD_REQUEST)
            return
        payload = self._read_json_body()
        family_id = str(payload.get("family_id") or "").strip()
        action = str(payload.get("action") or "").strip()
        case_id = str(payload.get("case_id") or "").strip()
        value = str(payload.get("value") or "").strip()
        source = self.server.family_qa["families"].get(family_id)
        if source is None:
            _json_response(self, {"error": "family_not_in_frozen_queue"}, status=HTTPStatus.NOT_FOUND)
            return
        allowed_decisions = {"approve_resonance", "approve_redox", "reject_split"}
        allowed_marks = {"matches_family", "outlier_blocker", ""}
        if action == "decision" and value not in allowed_decisions:
            _json_response(self, {"error": "invalid_family_decision"}, status=HTTPStatus.BAD_REQUEST)
            return
        frozen_reps = {str(row["case_id"]) for row in source["representatives"]}
        if action == "representative_mark" and (case_id not in frozen_reps or value not in allowed_marks):
            _json_response(self, {"error": "invalid_representative_mark"}, status=HTTPStatus.BAD_REQUEST)
            return
        if action not in {"decision", "representative_mark"}:
            _json_response(self, {"error": "invalid_family_action"}, status=HTTPStatus.BAD_REQUEST)
            return
        with self.server.review_mutation_lock:
            before = _family_manifest(self.server)
            after = json.loads(json.dumps(before))
            family = next(
                (item for item in after["families"] if str(item.get("family_id") or "") == family_id),
                None,
            )
            if family is None:
                _json_response(self, {"error": "family_not_in_manifest"}, status=HTTPStatus.CONFLICT)
                return
            now = datetime.now(timezone.utc).isoformat()
            if action == "decision":
                family["qa_decision"] = value
                family["qa_passed"] = value.startswith("approve_")
                family["proposed_status"] = "accept_both" if value.startswith("approve_") else None
                family["proposed_reason"] = {
                    "approve_resonance": "resonance-representation",
                    "approve_redox": "redox-representation",
                }.get(value)
                label = value
            else:
                marks = family.setdefault("representative_marks", {})
                if value:
                    marks[case_id] = value
                else:
                    marks.pop(case_id, None)
                label = value or "clear representative mark"
            family["qa_updated_at"] = now
            after["approval_required"] = True
            after["approved"] = False
            _write_json_atomic(self.server.family_qa_manifest_path, after)
            mutation = {
                "mutation_id": uuid.uuid4().hex,
                "mutation_type": "family_qa",
                "family_id": family_id,
                "case_id": case_id,
                "label": label,
                "before_manifest": before,
                "after_manifest": after,
                "timestamp": now,
                "undone": False,
            }
            history = self.server.review_histories.setdefault(session_id, [])
            history.append(mutation)
            del history[:-20]
        _json_response(self, {"ok": True, "mutation": _public_review_mutation(mutation), "family_qa": _public_family_qa(self.server)})

    def _api_review_undo(self) -> None:
        session_id = self._review_session_id()
        if not session_id:
            _json_response(self, {"error": "review_session_required"}, status=HTTPStatus.BAD_REQUEST)
            return
        payload = self._read_json_body()
        requested_id = str(payload.get("mutation_id") or "").strip()
        with self.server.review_mutation_lock:
            history = self.server.review_histories.get(session_id, [])
            mutation = next((item for item in reversed(history) if not item.get("undone")), None)
            if mutation is None:
                _json_response(self, {"error": "no_review_mutation_to_undo"}, status=HTTPStatus.CONFLICT)
                return
            if requested_id and requested_id != mutation["mutation_id"]:
                _json_response(self, {"error": "undo_must_target_latest_mutation"}, status=HTTPStatus.CONFLICT)
                return
            if mutation.get("mutation_type") == "family_qa":
                current = _family_manifest(self.server)
                if current != mutation["after_manifest"]:
                    _json_response(self, {"error": "family_manifest_changed_since_session_action"}, status=HTTPStatus.CONFLICT)
                    return
                _write_json_atomic(self.server.family_qa_manifest_path, mutation["before_manifest"])
                mutation["undone"] = True
                _json_response(
                    self,
                    {
                        "ok": True,
                        "case_id": mutation.get("case_id") or "",
                        "family_id": mutation["family_id"],
                        "mutation": _public_review_mutation(mutation),
                    },
                )
                return
            case_id = str(mutation["case_id"])
            with _connect(self.server.db_path) as conn:
                current = conn.execute(
                    "SELECT * FROM reviews WHERE case_id = ?", (case_id,)
                ).fetchone()
                if _review_state(current) != mutation["after_review"]:
                    _json_response(
                        self,
                        {"error": "review_changed_since_session_save"},
                        status=HTTPStatus.CONFLICT,
                    )
                    return
                current_fixture = snapshot_review_fixture(self.server.fixtures_dir, case_id)
                if current_fixture != mutation["after_fixture"]:
                    _json_response(
                        self,
                        {"error": "fixture_changed_since_session_save"},
                        status=HTTPStatus.CONFLICT,
                    )
                    return
                try:
                    restore_review_fixture_snapshot(
                        self.server.fixtures_dir, case_id, mutation["before_fixture"]
                    )
                    before = mutation["before_review"]
                    if before is None:
                        conn.execute("DELETE FROM reviews WHERE case_id = ?", (case_id,))
                    else:
                        conn.execute(
                            """
                            INSERT INTO reviews(
                                case_id, status, corrected_smiles, corrected_molblock,
                                notes, reviewer, updated_at
                            ) VALUES(?, ?, ?, ?, ?, ?, ?)
                            ON CONFLICT(case_id) DO UPDATE SET
                                status=excluded.status,
                                corrected_smiles=excluded.corrected_smiles,
                                corrected_molblock=excluded.corrected_molblock,
                                notes=excluded.notes,
                                reviewer=excluded.reviewer,
                                updated_at=excluded.updated_at
                            """,
                            (case_id, *(before[key] for key in REVIEW_COLUMNS)),
                        )
                    conn.commit()
                except Exception:
                    restore_review_fixture_snapshot(
                        self.server.fixtures_dir, case_id, current_fixture
                    )
                    raise
            mutation["undone"] = True
            _json_response(
                self,
                {
                    "ok": True,
                    "case_id": case_id,
                    "restored_review": mutation["before_review"],
                    "mutation": _public_review_mutation(mutation),
                },
            )

    def _api_cases(self, query: dict[str, list[str]]) -> None:
        category = query.get("category", [""])[0]
        status = query.get("status", [""])[0]
        reviewer = query.get("reviewer", [""])[0].strip()
        search = query.get("q", [""])[0].strip()
        triage_bucket = query.get("triage_bucket", [""])[0]
        limit = max(1, min(int(query.get("limit", ["100"])[0]), 500))
        offset = max(0, int(query.get("offset", ["0"])[0]))

        where: list[str] = []
        params: list[Any] = []
        if category:
            where.append("c.category = ?")
            params.append(category)
        if status:
            if status == "unreviewed":
                where.append("r.status IS NULL")
            else:
                where.append("r.status = ?")
                params.append(status)
        if reviewer:
            where.append("TRIM(r.reviewer) = ?")
            params.append(reviewer)
        if search:
            where.append("(c.case_id LIKE ? OR CAST(c.row_index AS TEXT) = ?)")
            params.extend((f"%{search}%", search))
        where_sql = (" WHERE " + " AND ".join(where)) if where else ""

        triage_records = getattr(self.server, "triage_records", {})
        with _connect(self.server.db_path) as conn:
            count_rows = conn.execute(
                """
                SELECT c.case_id
                FROM cases c
                LEFT JOIN reviews r ON r.case_id = c.case_id
                """
                + where_sql,
                params,
            ).fetchall()
            triage_bucket_counts = dict(
                sorted(
                    Counter(
                        triage_records.get(str(row["case_id"]), {}).get("triage_bucket", "")
                        for row in count_rows
                        if triage_records.get(str(row["case_id"]), {}).get("triage_bucket")
                    ).items()
                )
            )
            if triage_bucket:
                matching = conn.execute(
                    _case_query() + where_sql + " ORDER BY c.row_index",
                    params,
                ).fetchall()
                matching = [
                    row
                    for row in matching
                    if triage_records.get(str(row["case_id"]), {}).get("triage_bucket")
                    == triage_bucket
                ]
                total = len(matching)
                rows = matching[offset : offset + limit]
            else:
                rows = conn.execute(
                    _case_query() + where_sql + " ORDER BY c.row_index LIMIT ? OFFSET ?",
                    (*params, limit, offset),
                ).fetchall()
                total = conn.execute(
                    """
                    SELECT COUNT(*) AS count
                    FROM cases c
                    LEFT JOIN reviews r ON r.case_id = c.case_id
                    """
                    + where_sql,
                    params,
                ).fetchone()["count"]
        fixture_records = load_fixture_records(self.server.fixtures_dir)
        _json_response(
            self,
            {
                "total": total,
                "triage_bucket_counts": triage_bucket_counts,
                "items": [
                    _row_dict(
                        row,
                        fixture=fixture_records.get(str(row["case_id"])),
                        triage=triage_records.get(str(row["case_id"])),
                    )
                    for row in rows
                ],
            },
        )

    def _api_case(self, case_id: str) -> None:
        with _connect(self.server.db_path) as conn:
            row = conn.execute(_case_query() + " WHERE c.case_id = ?", (case_id,)).fetchone()
        fixture_records = load_fixture_records(self.server.fixtures_dir)
        payload = _row_dict(
            row,
            fixture=fixture_records.get(case_id),
            triage=getattr(self.server, "triage_records", {}).get(case_id),
        )
        if payload is None:
            _json_response(self, {"error": "case_not_found"}, status=HTTPStatus.NOT_FOUND)
            return
        _json_response(self, payload)

    def _api_xyz(self, case_id: str) -> None:
        with _connect(self.server.db_path) as conn:
            row = conn.execute(
                "SELECT xyz_path FROM cases WHERE case_id = ?", (case_id,)
            ).fetchone()
        if row is None:
            _json_response(self, {"error": "case_not_found"}, status=HTTPStatus.NOT_FOUND)
            return
        path = resolve_xyz_path(row["xyz_path"], self.server.xyz_dir)
        if not path.exists():
            _json_response(
                self,
                {"error": "xyz_not_found", "xyz_path": str(path)},
                status=HTTPStatus.NOT_FOUND,
            )
            return
        _text_response(self, path.read_text(encoding="utf-8"), content_type="chemical/x-xyz")

    def _api_candidate_sdf(self, case_id: str, query: dict[str, list[str]]) -> None:
        render_mode = query.get("mode", ["skeleton"])[0]
        if render_mode not in {"skeleton", "hydrogen"}:
            _json_response(self, {"error": "invalid_render_mode"}, status=HTTPStatus.BAD_REQUEST)
            return
        with _connect(self.server.db_path) as conn:
            row = conn.execute(
                """
                SELECT *
                FROM cases
                WHERE case_id = ?
                """,
                (case_id,),
            ).fetchone()
        if row is None:
            _json_response(self, {"error": "case_not_found"}, status=HTTPStatus.NOT_FOUND)
            return

        case = _row_dict(row)
        assert case is not None
        snapshot_smiles = str(case["candidate_snapshot_smiles"] or "")

        try:
            reconstructed_mol = reconstruct_case_mol(dict(row), xyz_dir=self.server.xyz_dir)
            candidate_mol = _benchmark_candidate_mol(reconstructed_mol)
            sdf = _mol_to_sdf_block(reconstructed_mol)
        except Exception as exc:  # noqa: BLE001
            _json_response(
                self,
                {
                    "sdf": "",
                    "svg": "",
                    "smiles": "",
                    "available": False,
                    "source": "live_reconstruction",
                    "live_candidate_status": "failed",
                    "live_candidate_smiles": "",
                    "live_candidate_smiles_exact_match": None,
                    "candidate_snapshot_smiles": snapshot_smiles,
                    "live_candidate_equivalence_reason": "live_reconstruction_failed",
                    "error": f"{type(exc).__name__}: {exc}",
                    "render_error": "",
                },
            )
            return

        smiles = _safe_smiles(candidate_mol)
        svg = ""
        render_error = ""
        try:
            reference_mol = _render_mol_from_smiles(str(row["reference_smiles"] or ""))
            candidate_notes, _ = _triage_atom_notes(
                dict(row),
                getattr(self.server, "triage_records", {}).get(case_id),
                reconstructed_mol,
                reference_mol,
                self.server.xyz_dir,
            )
            svg = _render_mol_svg(
                reconstructed_mol,
                legend=f"{case_id} candidate current",
                atom_notes=candidate_notes,
                show_hydrogens=render_mode == "hydrogen",
            )
        except Exception as exc:  # noqa: BLE001
            render_error = f"{type(exc).__name__}: {exc}"
        _json_response(
            self,
            {
                "sdf": sdf,
                "svg": svg,
                "smiles": smiles,
                "available": True,
                "source": "live_reconstruction",
                "live_candidate_status": "ok",
                "live_candidate_smiles": smiles,
                "error": "",
                "render_error": render_error,
                **_live_candidate_comparison(candidate_mol, snapshot_smiles),
            },
        )

    def _api_reference_xyz(self, case_id: str) -> None:
        with _connect(self.server.db_path) as conn:
            row = conn.execute("SELECT * FROM cases WHERE case_id = ?", (case_id,)).fetchone()
        if row is None:
            _json_response(self, {"error": "case_not_found"}, status=HTTPStatus.NOT_FOUND)
            return
        reference = _render_mol_from_smiles(str(row["reference_smiles"] or ""))
        if reference is None:
            _json_response(
                self,
                {
                    "available": False,
                    "reference_xyz_status": "unavailable",
                    "mapping_confidence": "failed",
                    "reference_atom_to_xyz": {},
                    "failure_code": "reference_missing_or_invalid",
                    "error": "reference_missing_or_invalid",
                    "sdf": "",
                },
            )
            return
        try:
            candidate = reconstruct_case_mol(dict(row), xyz_dir=self.server.xyz_dir)
            xyz_atoms = parse_xyz_atoms(
                resolve_xyz_path(str(row["xyz_path"] or ""), self.server.xyz_dir).read_text(
                    encoding="utf-8"
                )
            )
            mapping = map_candidate_reference_xyz(candidate, reference, xyz_atoms)
            ambiguity_type, ambiguity_reason = _mapping_ambiguity_details(mapping)
            reference_xyz, confidence, error = _reference_xyz_mol(
                dict(row),
                getattr(self.server, "triage_records", {}).get(case_id),
                candidate,
                reference,
                self.server.xyz_dir,
                mapping,
            )
            representative_mapping = reference_xyz is not None and confidence == "ambiguous"
            _json_response(
                self,
                {
                    "available": reference_xyz is not None,
                    "sdf": _mol_to_sdf_block(reference_xyz) if reference_xyz is not None else "",
                    "reference_xyz_status": "available" if reference_xyz is not None else "unavailable",
                    "mapping_confidence": confidence,
                    "mapping_is_representative": representative_mapping,
                    "mapping_selection_method": (
                        "existing_valid_deterministic_representative"
                        if representative_mapping
                        else "unique_mapping"
                        if reference_xyz is not None
                        else ""
                    ),
                    "mapping_ambiguity_type": ambiguity_type,
                    "mapping_ambiguity_reason": ambiguity_reason,
                    "mapping_enumeration_truncated": mapping.enumeration_truncated,
                    "mapping_timeout": mapping.timeout,
                    "mapping_signature_count": mapping.mapping_signature_count,
                    "mapping_equal_best_count": mapping.equal_best_mapping_count,
                    "mapping_ambiguity_locations": _mapping_ambiguity_locations(mapping),
                    "reference_atom_to_xyz": (
                        {
                            str(reference_index): candidate_index
                            for reference_index, candidate_index in mapping.reference_to_candidate.items()
                        }
                        if reference_xyz is not None
                        else {}
                    ),
                    "coordinate_source": "source_xyz",
                    "connectivity_source": "reference_graph",
                    "mapped_comparison": (
                        _mapped_coordination_comparison(
                            candidate, reference, mapping.reference_to_candidate
                        )
                        if reference_xyz is not None
                        else {}
                    ),
                    "failure_code": error,
                    "error": error,
                },
            )
        except Exception as exc:  # noqa: BLE001
            failure = f"{type(exc).__name__}: {exc}"
            _json_response(
                self,
                {
                    "available": False,
                    "reference_xyz_status": "unavailable",
                    "mapping_confidence": "failed",
                    "reference_atom_to_xyz": {},
                    "failure_code": "mapping_failed",
                    "sdf": "",
                    "error": failure,
                },
            )

    def _api_graph(self, case_id: str, query: dict[str, list[str]]) -> None:
        kind = query.get("kind", ["candidate"])[0]
        if kind not in {"candidate", "reference"}:
            _json_response(self, {"error": "invalid_graph_kind"}, status=HTTPStatus.BAD_REQUEST)
            return
        with _connect(self.server.db_path) as conn:
            case = conn.execute("SELECT * FROM cases WHERE case_id = ?", (case_id,)).fetchone()
        if case is None:
            _json_response(self, {"error": "case_not_found"}, status=HTTPStatus.NOT_FOUND)
            return
        if kind == "candidate":
            mol = reconstruct_case_mol(dict(case), xyz_dir=self.server.xyz_dir)
            smiles = _safe_smiles(mol)
        else:
            source_smiles = str(case["reference_smiles"] or "")
            mol = _render_mol_from_smiles(source_smiles)
            if mol is None:
                _json_response(
                    self,
                    {"error": "reference_smiles_missing_or_invalid"},
                    status=HTTPStatus.UNPROCESSABLE_ENTITY,
                )
                return
            smiles = source_smiles
        _json_response(self, _review_graph_payload(mol, kind=kind, smiles=smiles))

    def _api_render(self, case_id: str, query: dict[str, list[str]]) -> None:
        kind = query.get("kind", ["candidate"])[0]
        localize = query.get("localize", [""])[0] == "1"
        render_mode = query.get("mode", ["skeleton"])[0]
        if render_mode not in {"skeleton", "hydrogen"}:
            _json_response(self, {"error": "invalid_render_mode"}, status=HTTPStatus.BAD_REQUEST)
            return
        if kind not in {"candidate", "reference", "candidate_organic", "reference_organic"}:
            _json_response(self, {"error": "invalid_render_kind"}, status=HTTPStatus.BAD_REQUEST)
            return

        with _connect(self.server.db_path) as conn:
            cache_kind = _render_cache_kind(kind, render_mode)
            if kind != "candidate" and not localize:
                cached = conn.execute(
                    "SELECT svg, smiles, error FROM render_cache WHERE case_id = ? AND kind = ?",
                    (case_id, cache_kind),
                ).fetchone()
                if cached is not None and not cached["error"]:
                    _json_response(self, {"kind": kind, **dict(cached)})
                    return

            case = conn.execute("SELECT * FROM cases WHERE case_id = ?", (case_id,)).fetchone()
            if case is None:
                _json_response(self, {"error": "case_not_found"}, status=HTTPStatus.NOT_FOUND)
                return

            svg = ""
            smiles = ""
            error = ""
            try:
                if kind == "candidate":
                    mol = reconstruct_case_mol(dict(case), xyz_dir=self.server.xyz_dir)
                    smiles = _safe_smiles(mol)
                    svg = _render_mol_svg(mol, legend=f"{case_id} candidate")
                elif kind == "reference":
                    mol = _render_mol_from_smiles(case["reference_smiles"] or "")
                    if mol is None:
                        raise ValueError("reference_smiles_missing_or_invalid")
                    smiles = _safe_smiles(mol)
                    atom_notes = None
                    if localize:
                        candidate = reconstruct_case_mol(dict(case), xyz_dir=self.server.xyz_dir)
                        triage = getattr(self.server, "triage_records", {}).get(case_id)
                        _, atom_notes = _triage_atom_notes(
                            dict(case),
                            triage,
                            candidate,
                            mol,
                            self.server.xyz_dir,
                        )
                        if render_mode == "hydrogen":
                            reference_xyz, _, _ = _reference_xyz_mol(
                                dict(case), triage, candidate, mol, self.server.xyz_dir
                            )
                            if reference_xyz is not None:
                                mol = reference_xyz
                                atom_notes = _triage_source_atom_notes(
                                    triage, candidate, reference_assignment=True
                                )
                    svg = _render_mol_svg(
                        mol,
                        legend=f"{case_id} Reference",
                        atom_notes=atom_notes,
                        show_hydrogens=render_mode == "hydrogen",
                    )
                elif kind == "candidate_organic":
                    mol = _render_mol_from_smiles(case["candidate_organic_smiles"] or "")
                    if mol is None:
                        raise ValueError("candidate_organic_smiles_missing_or_invalid")
                    smiles = _safe_smiles(mol)
                    svg = _render_mol_svg(mol, legend=f"{case_id} candidate organic")
                else:
                    mol = _render_mol_from_smiles(case["reference_organic_smiles"] or "")
                    if mol is None:
                        raise ValueError("reference_organic_smiles_missing_or_invalid")
                    smiles = _safe_smiles(mol)
                    svg = _render_mol_svg(mol, legend=f"{case_id} Reference organic")
            except Exception as exc:  # noqa: BLE001
                error = f"{type(exc).__name__}: {exc}"

            if kind != "candidate" and not localize:
                conn.execute(
                    """
                    INSERT OR REPLACE INTO render_cache(
                        case_id, kind, svg, smiles, error, generated_at
                    )
                    VALUES(?, ?, ?, ?, ?, ?)
                    """,
                    (
                        case_id,
                        cache_kind,
                        svg,
                        smiles,
                        error,
                        datetime.now(timezone.utc).isoformat(),
                    ),
                )
                conn.commit()
        _json_response(self, {"kind": kind, "svg": svg, "smiles": smiles, "error": error})

    def _api_review(self, case_id: str) -> None:
        payload = self._read_json_body()
        status = str(payload.get("status", "")).strip()
        if status not in ALLOWED_STATUSES:
            _json_response(
                self,
                {"error": "invalid_status", "allowed": sorted(ALLOWED_STATUSES)},
                status=HTTPStatus.BAD_REQUEST,
            )
            return

        corrected_smiles = str(payload.get("corrected_smiles", "")).strip()
        corrected_molblock = str(payload.get("corrected_molblock", ""))
        if status == "manual_reference" and not corrected_smiles:
            _json_response(
                self,
                {"error": "corrected_smiles_required_for_manual_reference"},
                status=HTTPStatus.BAD_REQUEST,
            )
            return

        updated_at = datetime.now(timezone.utc).isoformat()
        review = {
            "status": status,
            "corrected_smiles": corrected_smiles,
            "corrected_molblock": corrected_molblock,
            "notes": str(payload.get("notes", "")).strip(),
            "reviewer": str(payload.get("reviewer", "")).strip(),
            "updated_at": updated_at,
        }
        session_id = self._review_session_id()
        mutation = None
        with self.server.review_mutation_lock:
            with _connect(self.server.db_path) as conn:
                case = conn.execute("SELECT * FROM cases WHERE case_id = ?", (case_id,)).fetchone()
                if case is None:
                    _json_response(self, {"error": "case_not_found"}, status=HTTPStatus.NOT_FOUND)
                    return
                before_review = _review_state(
                    conn.execute("SELECT * FROM reviews WHERE case_id = ?", (case_id,)).fetchone()
                )
                before_fixture = snapshot_review_fixture(self.server.fixtures_dir, case_id)
                try:
                    fixture = sync_review_fixture(
                        dict(case),
                        review,
                        fixtures_dir=self.server.fixtures_dir,
                        xyz_dir=self.server.xyz_dir,
                    )
                    conn.execute(
                        """
                        INSERT INTO reviews(
                            case_id, status, corrected_smiles, corrected_molblock,
                            notes, reviewer, updated_at
                        )
                        VALUES(?, ?, ?, ?, ?, ?, ?)
                        ON CONFLICT(case_id) DO UPDATE SET
                            status = excluded.status,
                            corrected_smiles = excluded.corrected_smiles,
                            corrected_molblock = excluded.corrected_molblock,
                            notes = excluded.notes,
                            reviewer = excluded.reviewer,
                            updated_at = excluded.updated_at
                        """,
                        (
                            case_id,
                            status,
                            corrected_smiles,
                            corrected_molblock,
                            review["notes"],
                            review["reviewer"],
                            updated_at,
                        ),
                    )
                    conn.commit()
                    after_fixture = snapshot_review_fixture(self.server.fixtures_dir, case_id)
                except (FileNotFoundError, ValueError) as exc:
                    restore_review_fixture_snapshot(
                        self.server.fixtures_dir, case_id, before_fixture
                    )
                    _json_response(self, {"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)
                    return
                except Exception:
                    restore_review_fixture_snapshot(
                        self.server.fixtures_dir, case_id, before_fixture
                    )
                    raise
            if session_id:
                mutation = {
                    "mutation_id": uuid.uuid4().hex,
                    "case_id": case_id,
                    "before_review": before_review,
                    "before_fixture": before_fixture,
                    "after_fixture": after_fixture,
                    "after_review": dict(review),
                    "timestamp": updated_at,
                    "undone": False,
                }
                history = self.server.review_histories.setdefault(session_id, [])
                history.append(mutation)
                del history[:-20]
        _json_response(
            self,
            {
                "ok": True,
                "fixture": fixture,
                "mutation": _public_review_mutation(mutation) if mutation else None,
            },
        )


def main() -> None:
    args = _parse_args()
    if not args.db.exists():
        raise SystemExit(
            f"Database does not exist: {args.db}\n"
            "Create it first with: uv run python tools/molgr_review/import_cases.py --input <queue.csv> --db <path>"
        )
    try:
        runtime_info = validate_project_runtime(REPO_ROOT)
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    server = ReviewServer((args.host, args.port), ReviewHandler)
    server.db_path = args.db
    server.xyz_dir = args.xyz_dir
    server.fixtures_dir = args.fixtures_dir
    server.triage_records = load_triage_records(args.triage_csv)
    server.family_qa = load_family_qa(args.family_qa_csv, args.family_qa_manifest)
    server.family_qa_manifest_path = args.family_qa_manifest
    server.review_histories = {}
    server.review_mutation_lock = threading.RLock()
    server.runtime_info = runtime_info
    server.fixtures_dir.mkdir(parents=True, exist_ok=True)
    print(f"Molecule review server: http://{args.host}:{args.port}")
    print(f"database: {args.db}")
    if args.xyz_dir is not None:
        print(f"xyz directory: {args.xyz_dir}")
    print(f"review fixtures: {args.fixtures_dir}")
    if args.family_qa_csv is not None:
        print(f"family QA queue: {args.family_qa_csv}")
        print(f"family QA pending manifest: {args.family_qa_manifest}")
    print(f"molgr source: {runtime_info['molgr_source']}")
    print(f"C++ extension: {runtime_info['cpp_extension']}")
    print(
        "checkout: "
        f"{runtime_info['git_revision'] or 'unknown'}"
        f"{' (dirty)' if runtime_info['git_dirty'] else ''}"
    )
    server.serve_forever()


if __name__ == "__main__":
    main()
