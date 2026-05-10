#!/usr/bin/env python3
"""Trace tmQMg reconstruction by dataset id.

This is a tmQMg input adapter around scripts/reconstruction_trace.py. The
generic script owns the reconstruction trajectory tracing and report rendering;
this wrapper only resolves CSV rows, XYZ paths, charges, and radical counts.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any, Optional, Sequence


if __package__ in (None, ""):
    project_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(project_root))
    sys.path.insert(0, str(project_root / "src"))

from rdkit import Chem, RDLogger

from molgr.config import MolGRConfig
from scripts.reconstruction_trace import (
    TraceInputCase,
    _jsonable,
    _render_markdown_report,
    _resolve_output_format,
    split_repeated_values,
    trace_reconstruction_case,
)


RDLogger.DisableLog("rdApp.*")  # type: ignore[arg-type]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csv",
        type=Path,
        required=True,
        help="tmQMg metadata CSV path. Rows are selected by the 'id' column.",
    )
    parser.add_argument(
        "--xyz-dir",
        type=Path,
        required=True,
        help="Directory containing tmQMg XYZ files named <id>.xyz.",
    )
    parser.add_argument(
        "--id",
        dest="ids",
        action="append",
        required=True,
        help="tmQMg id to trace. May be repeated or comma-separated.",
    )
    parser.add_argument(
        "--spin-source",
        choices=("closed_shell", "reference_smiles"),
        default="closed_shell",
        help=(
            "How to choose total radical electrons when --total-radical-electrons is not set. "
            "Default matches the closed-shell tmQMg regression workflow."
        ),
    )
    parser.add_argument(
        "--total-radical-electrons",
        type=int,
        default=None,
        help="Explicit radical-electron count to use for every requested id.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Optional output path. Defaults to stdout.",
    )
    parser.add_argument(
        "--format",
        choices=("auto", "markdown", "json"),
        default="auto",
        help=(
            "Output format. In auto mode, .json writes JSON; every other output path and stdout "
            "write a Markdown report."
        ),
    )
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="JSON indentation when --format json is used. Use 0 for compact JSON.",
    )
    parser.add_argument(
        "--no-score-all-candidates",
        action="store_true",
        help=(
            "Do not compute analysis metal-candidate metrics for candidates filtered out by "
            "discordance. Production metadata is still reported."
        ),
    )
    args = parser.parse_args()
    if args.total_radical_electrons is not None and args.total_radical_electrons < 0:
        parser.error("--total-radical-electrons must be >= 0")
    return args


def _load_rows_by_id(path: Path) -> dict[str, tuple[int, dict[str, str]]]:
    rows: dict[str, tuple[int, dict[str, str]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "id" not in reader.fieldnames:
            raise ValueError(f"tmQMg CSV must contain an 'id' column: {path}")
        for row_index, row in enumerate(reader, start=1):
            case_id = row.get("id", "").strip()
            if case_id and case_id not in rows:
                rows[case_id] = (row_index, row)
    return rows


def _parse_int_field(row: dict[str, str], field_name: str) -> int:
    raw_value = row.get(field_name, "").strip()
    if raw_value == "":
        raise ValueError(
            f"Missing required integer field {field_name!r} for id={row.get('id', '')}"
        )
    return int(raw_value)


def _reference_smiles(row: dict[str, str]) -> str:
    for field_name in ("smiles", "SMILES", "canonical_smiles", "CanonicalSMILES"):
        value = row.get(field_name, "").strip()
        if value:
            return value
    raise ValueError(f"No reference SMILES column found for id={row.get('id', '')}")


def _radicals_from_reference_smiles(row: dict[str, str]) -> int:
    smiles = _reference_smiles(row)
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        raise ValueError(f"Could not parse reference SMILES for id={row.get('id', '')}: {smiles}")
    return int(sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()))


def _resolve_total_radicals(
    row: dict[str, str],
    *,
    spin_source: str,
    explicit_total_radicals: Optional[int],
) -> tuple[int, str]:
    if explicit_total_radicals is not None:
        return explicit_total_radicals, "explicit"
    if spin_source == "closed_shell":
        return 0, "closed_shell"
    if spin_source == "reference_smiles":
        return _radicals_from_reference_smiles(row), "reference_smiles"
    raise ValueError(f"Unsupported spin source: {spin_source!r}")


def _trace_tmqmg_case(
    row_index: int,
    row: dict[str, str],
    *,
    xyz_dir: Path,
    spin_source: str,
    explicit_total_radicals: Optional[int],
    score_all_candidates: bool,
    config: MolGRConfig | None = None,
) -> dict[str, Any]:
    started = time.perf_counter()
    case_id = row.get("id", "").strip()
    xyz_path = xyz_dir / f"{case_id}.xyz"
    total_charge = _parse_int_field(row, "charge")
    total_radicals, resolved_spin_source = _resolve_total_radicals(
        row,
        spin_source=spin_source,
        explicit_total_radicals=explicit_total_radicals,
    )
    input_case = TraceInputCase(
        id=case_id,
        xyz_block=xyz_path.read_text(encoding="utf-8"),
        total_charge=total_charge,
        total_radical_electrons=total_radicals,
        xyz_path=xyz_path,
        xyz_source="tmqmg_xyz_path",
    )
    trace = trace_reconstruction_case(
        input_case,
        score_all_candidates=score_all_candidates,
        config=config,
    )
    trace.update(
        {
            "row_index": row_index,
            "reference_smiles": row.get("smiles", ""),
            "spin_source": resolved_spin_source,
            "elapsed_seconds": time.perf_counter() - started,
        }
    )
    return trace


def _write_output(args: argparse.Namespace, output: dict[str, Any]) -> None:
    output_format = _resolve_output_format(args)
    if output_format == "json":
        output_text = json.dumps(
            _jsonable(output),
            ensure_ascii=False,
            indent=None if args.indent == 0 else args.indent,
            allow_nan=False,
        )
    else:
        output_text = _render_markdown_report(_jsonable(output))

    if args.out is None:
        print(output_text, end="" if output_text.endswith("\n") else "\n")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(
            output_text if output_text.endswith("\n") else output_text + "\n",
            encoding="utf-8",
        )


def _trace_requested_rows(
    case_ids: Sequence[str],
    rows_by_id: dict[str, tuple[int, dict[str, str]]],
    *,
    args: argparse.Namespace,
    config: MolGRConfig | None,
) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for case_id in case_ids:
        started = time.perf_counter()
        row_entry = rows_by_id.get(case_id)
        if row_entry is None:
            cases.append(
                {
                    "id": case_id,
                    "status": "missing_csv_row",
                    "trace_kind": "unknown",
                    "elapsed_seconds": time.perf_counter() - started,
                }
            )
            continue

        row_index, row = row_entry
        try:
            cases.append(
                _trace_tmqmg_case(
                    row_index,
                    row,
                    xyz_dir=args.xyz_dir,
                    spin_source=args.spin_source,
                    explicit_total_radicals=args.total_radical_electrons,
                    score_all_candidates=not args.no_score_all_candidates,
                    config=config,
                )
            )
        except Exception as exc:
            cases.append(
                {
                    "row_index": row_index,
                    "id": case_id,
                    "status": "error",
                    "trace_kind": "unknown",
                    "error": f"{type(exc).__name__}: {exc}",
                    "elapsed_seconds": time.perf_counter() - started,
                }
            )
    return cases


def main() -> int:
    args = _parse_args()
    case_ids = split_repeated_values(args.ids)
    rows_by_id = _load_rows_by_id(args.csv)
    config: MolGRConfig | None = None
    cases = _trace_requested_rows(case_ids, rows_by_id, args=args, config=config)

    output = {
        "input": {
            "source": "tmQMg",
            "csv": str(args.csv),
            "xyz_dir": str(args.xyz_dir),
            "ids": case_ids,
            "spin_source": args.spin_source,
            "total_radical_electrons": args.total_radical_electrons,
        },
        "case_count": len(cases),
        "cases": cases,
    }
    _write_output(args, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
