#!/usr/bin/env python3
# pyright: reportCallIssue=false
"""Benchmark tmQMg reconstruction speed for Python and C++ backends.

The script walks tmQMg rows in CSV order, resolves each XYZ file, and runs one
reconstruction call with backend="python" and one with backend="cpp" for each
row. It uses a fresh MolGR config object so each run is isolated from global
configuration mutations.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import Counter
from pathlib import Path
from statistics import median
from typing import Any, Literal


if __package__ in (None, ""):
    project_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(project_root / "src"))
    sys.path.insert(0, str(project_root))

from rdkit import Chem, RDLogger

from molgr.config import MolGRConfig
from molgr.interface import xyz_to_rdmol


RDLogger.DisableLog("rdApp.*")  # type: ignore[arg-type]

BackendName = Literal["python", "cpp"]

RESULT_FIELDNAMES = (
    "row_index",
    "id",
    "xyz_path",
    "charge",
    "reference_smiles_input",
    "spin_source",
    "total_radical_electrons_used",
    "spin_multiplicity_used",
    "python_status",
    "python_elapsed_ms",
    "python_error",
    "cpp_status",
    "cpp_elapsed_ms",
    "cpp_error",
    "speedup_ratio_python_over_cpp",
    "row_elapsed_ms",
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csv",
        type=Path,
        required=True,
        help="tmQMg metadata CSV path.",
    )
    parser.add_argument(
        "--xyz-dir",
        type=Path,
        default=None,
        help=(
            "Optional directory containing <id>.xyz files. If omitted, use the "
            "xyz_path column from the CSV."
        ),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional maximum number of rows to benchmark after filtering.",
    )
    parser.add_argument(
        "--start-row",
        type=int,
        default=1,
        help="1-based CSV row index to start from. Default: 1.",
    )
    parser.add_argument(
        "--end-row",
        type=int,
        default=None,
        help="Optional 1-based CSV row index to stop at, inclusive.",
    )
    parser.add_argument(
        "--spin-source",
        choices=("closed_shell", "reference_smiles"),
        default="closed_shell",
        help=(
            "How to determine total radical electrons. Default keeps the tmQMg "
            "closed-shell benchmark path."
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output CSV path.",
    )
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=None,
        help="Optional summary JSON path. Defaults to <out>.summary.json.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=10,
        help="Print one stderr progress line every N rows. Use 0 to silence.",
    )
    args = parser.parse_args()
    if args.start_row < 1:
        parser.error("--start-row must be >= 1")
    if args.end_row is not None and args.end_row < args.start_row:
        parser.error("--end-row must be >= --start-row")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be >= 1")
    if args.progress_every < 0:
        parser.error("--progress-every must be >= 0")
    return args


def _summary_path_from_output(output_path: Path) -> Path:
    if output_path.suffix:
        return output_path.with_suffix(f"{output_path.suffix}.summary.json")
    return output_path.with_name(f"{output_path.name}.summary.json")


def _load_rows(csv_path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with csv_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError(f"tmQMg CSV has no header row: {csv_path}")
        return list(reader.fieldnames), list(reader)


def _validate_columns(
    fieldnames: list[str],
    *,
    xyz_dir: Path | None,
    spin_source: str,
) -> None:
    fieldname_set = set(fieldnames)
    missing = [name for name in ("id", "charge") if name not in fieldname_set]
    if missing:
        raise ValueError(
            "tmQMg CSV is missing required columns: " + ", ".join(sorted(missing))
        )
    if xyz_dir is None and "xyz_path" not in fieldname_set:
        raise ValueError(
            "tmQMg CSV must contain an 'xyz_path' column when --xyz-dir is omitted."
        )
    if spin_source == "reference_smiles" and "smiles" not in fieldname_set:
        raise ValueError(
            "tmQMg CSV must contain a 'smiles' column when --spin-source reference_smiles is used."
        )


def _selected_rows(
    rows: list[dict[str, str]],
    *,
    start_row: int,
    end_row: int | None,
    limit: int | None,
) -> list[tuple[int, dict[str, str]]]:
    selected: list[tuple[int, dict[str, str]]] = []
    for row_index, row in enumerate(rows, start=1):
        if row_index < start_row:
            continue
        if end_row is not None and row_index > end_row:
            break
        selected.append((row_index, row))
        if limit is not None and len(selected) >= limit:
            break
    return selected


def _resolve_xyz_path(
    row: dict[str, str],
    *,
    xyz_dir: Path | None,
) -> Path:
    stored_xyz_path = row.get("xyz_path", "").strip()
    if xyz_dir is not None:
        if stored_xyz_path:
            return xyz_dir / Path(stored_xyz_path).name
        case_id = row.get("id", "").strip()
        if case_id:
            return xyz_dir / f"{case_id}.xyz"
        raise ValueError("Missing id for xyz-dir lookup.")
    if stored_xyz_path:
        return Path(stored_xyz_path)
    raise ValueError("tmQMg row has no xyz_path and no --xyz-dir was provided.")


def _reference_total_radical_electrons(row: dict[str, str]) -> int:
    reference_smiles = row.get("smiles", "").strip()
    if not reference_smiles:
        raise ValueError("missing_reference_smiles")

    reference_mol = Chem.MolFromSmiles(reference_smiles, sanitize=False)
    if reference_mol is None:
        raise ValueError("reference_parse_failed")
    return sum(atom.GetNumRadicalElectrons() for atom in reference_mol.GetAtoms())


def _resolve_total_radical_electrons(
    row: dict[str, str],
    *,
    spin_source: str,
) -> int:
    if spin_source == "closed_shell":
        return 0
    if spin_source == "reference_smiles":
        return _reference_total_radical_electrons(row)
    raise ValueError(f"Unsupported spin source: {spin_source!r}")


def _time_backend(
    xyz_block: str,
    *,
    total_charge: int,
    spin_multiplicity: int,
    backend: BackendName,
    config,
) -> dict[str, Any]:
    started_at = time.perf_counter()
    try:
        xyz_to_rdmol(
            xyz_block,
            total_charge,
            spin_multiplicity,
            backend=backend,
            config=config,
        )
    except Exception as exc:  # noqa: BLE001
        elapsed_ms = (time.perf_counter() - started_at) * 1000.0
        return {
            "status": f"failed:{type(exc).__name__}",
            "elapsed_ms": round(elapsed_ms, 6),
            "error": str(exc),
        }

    elapsed_ms = (time.perf_counter() - started_at) * 1000.0
    return {
        "status": "ok",
        "elapsed_ms": round(elapsed_ms, 6),
        "error": "",
    }


def _empty_result(
    row_index: int,
    row: dict[str, str],
    xyz_path: Path | None,
    spin_source: str,
) -> dict[str, Any]:
    return {
        "row_index": row_index,
        "id": row.get("id", "").strip(),
        "xyz_path": str(xyz_path) if xyz_path is not None else "",
        "charge": row.get("charge", "").strip(),
        "reference_smiles_input": row.get("smiles", "").strip(),
        "spin_source": spin_source,
        "total_radical_electrons_used": "",
        "spin_multiplicity_used": "",
        "python_status": "",
        "python_elapsed_ms": "",
        "python_error": "",
        "cpp_status": "",
        "cpp_elapsed_ms": "",
        "cpp_error": "",
        "speedup_ratio_python_over_cpp": "",
        "row_elapsed_ms": "",
    }


def _finalize_skipped(
    result: dict[str, Any],
    *,
    reason: str,
    error: str,
    started_at: float,
) -> dict[str, Any]:
    status = f"skipped:{reason}"
    result["python_status"] = status
    result["python_error"] = error
    result["cpp_status"] = status
    result["cpp_error"] = error
    result["row_elapsed_ms"] = round((time.perf_counter() - started_at) * 1000.0, 6)
    return result


def _benchmark_row(
    row_index: int,
    row: dict[str, str],
    *,
    xyz_dir: Path | None,
    spin_source: str,
    config,
) -> dict[str, Any]:
    started_at = time.perf_counter()
    try:
        xyz_path = _resolve_xyz_path(row, xyz_dir=xyz_dir)
    except Exception as exc:  # noqa: BLE001
        result = _empty_result(row_index, row, None, spin_source)
        return _finalize_skipped(
            result,
            reason="xyz_path_unresolved",
            error=str(exc),
            started_at=started_at,
        )

    result = _empty_result(row_index, row, xyz_path, spin_source)
    if not xyz_path.exists():
        return _finalize_skipped(
            result,
            reason="xyz_missing",
            error=f"Missing XYZ file: {xyz_path}",
            started_at=started_at,
        )

    charge_raw = row.get("charge", "").strip()
    if not charge_raw:
        return _finalize_skipped(
            result,
            reason="missing_charge",
            error=f"Missing charge for id={row.get('id', '').strip() or '<unknown>'}",
            started_at=started_at,
        )
    try:
        total_charge = int(charge_raw)
    except ValueError:
        return _finalize_skipped(
            result,
            reason="invalid_charge",
            error=f"Invalid charge for id={row.get('id', '').strip() or '<unknown>'}: {charge_raw!r}",
            started_at=started_at,
        )

    try:
        total_radical_electrons = _resolve_total_radical_electrons(
            row,
            spin_source=spin_source,
        )
    except Exception as exc:  # noqa: BLE001
        return _finalize_skipped(
            result,
            reason=str(exc),
            error=f"{str(exc)} for id={row.get('id', '').strip() or '<unknown>'}",
            started_at=started_at,
        )

    result["total_radical_electrons_used"] = total_radical_electrons
    result["spin_multiplicity_used"] = total_radical_electrons + 1

    xyz_block = xyz_path.read_text(encoding="utf-8")
    python_result = _time_backend(
        xyz_block,
        total_charge=total_charge,
        spin_multiplicity=total_radical_electrons + 1,
        backend="python",
        config=config,
    )
    cpp_result = _time_backend(
        xyz_block,
        total_charge=total_charge,
        spin_multiplicity=total_radical_electrons + 1,
        backend="cpp",
        config=config,
    )

    result["python_status"] = python_result["status"]
    result["python_elapsed_ms"] = python_result["elapsed_ms"]
    result["python_error"] = python_result["error"]
    result["cpp_status"] = cpp_result["status"]
    result["cpp_elapsed_ms"] = cpp_result["elapsed_ms"]
    result["cpp_error"] = cpp_result["error"]
    if (
        result["python_status"] == "ok"
        and result["cpp_status"] == "ok"
        and float(result["cpp_elapsed_ms"]) > 0.0
    ):
        result["speedup_ratio_python_over_cpp"] = round(
            float(result["python_elapsed_ms"]) / float(result["cpp_elapsed_ms"]),
            6,
        )
    result["row_elapsed_ms"] = round((time.perf_counter() - started_at) * 1000.0, 6)
    return result


def _timing_stats(values: list[float]) -> dict[str, float]:
    if not values:
        return {
            "count": 0,
            "total_ms": 0.0,
            "mean_ms": 0.0,
            "median_ms": 0.0,
            "min_ms": 0.0,
            "max_ms": 0.0,
        }
    return {
        "count": len(values),
        "total_ms": round(sum(values), 6),
        "mean_ms": round(sum(values) / len(values), 6),
        "median_ms": round(median(values), 6),
        "min_ms": round(min(values), 6),
        "max_ms": round(max(values), 6),
    }


def _status_counts(results: list[dict[str, Any]], field_name: str) -> dict[str, int]:
    counts = Counter(str(row[field_name]) for row in results if row.get(field_name))
    return dict(sorted(counts.items()))


def _summarize(
    results: list[dict[str, Any]],
    *,
    csv_path: Path,
    xyz_dir: Path | None,
    spin_source: str,
    config,
    wall_seconds: float,
) -> dict[str, Any]:
    python_successes = [
        float(row["python_elapsed_ms"])
        for row in results
        if row.get("python_status") == "ok"
    ]
    cpp_successes = [
        float(row["cpp_elapsed_ms"])
        for row in results
        if row.get("cpp_status") == "ok"
    ]
    paired_speedups = [
        float(row["speedup_ratio_python_over_cpp"])
        for row in results
        if row.get("speedup_ratio_python_over_cpp") != ""
    ]

    summary = {
        "source_csv": str(csv_path),
        "xyz_dir": str(xyz_dir) if xyz_dir is not None else None,
        "spin_source": spin_source,
        "record_count": len(results),
        "wall_seconds": round(wall_seconds, 6),
        "config": {
            "cpp_backend": {
                "max_threads": config.cpp_backend.max_threads,
                "enable_target_bucket_parallelism": config.cpp_backend.enable_target_bucket_parallelism,
                "enable_candidate_scoring_parallelism": config.cpp_backend.enable_candidate_scoring_parallelism,
                "enable_uff_atom_typing_cache": config.cpp_backend.enable_uff_atom_typing_cache,
                "target_bucket_parallel_threshold": config.cpp_backend.target_bucket_parallel_threshold,
                "candidate_score_parallel_threshold": config.cpp_backend.candidate_score_parallel_threshold,
            }
        },
        "python": {
            **_timing_stats(python_successes),
            "ok_count": len(python_successes),
            "failed_count": sum(
                1 for row in results if str(row.get("python_status", "")).startswith("failed:")
            ),
            "skipped_count": sum(
                1 for row in results if str(row.get("python_status", "")).startswith("skipped:")
            ),
        },
        "cpp": {
            **_timing_stats(cpp_successes),
            "ok_count": len(cpp_successes),
            "failed_count": sum(
                1 for row in results if str(row.get("cpp_status", "")).startswith("failed:")
            ),
            "skipped_count": sum(
                1 for row in results if str(row.get("cpp_status", "")).startswith("skipped:")
            ),
        },
        "paired": {
            "ok_count": len(paired_speedups),
            "mean_speedup_ratio_python_over_cpp": round(sum(paired_speedups) / len(paired_speedups), 6)
            if paired_speedups
            else 0.0,
            "median_speedup_ratio_python_over_cpp": round(median(paired_speedups), 6)
            if paired_speedups
            else 0.0,
        },
        "status_counts": {
            "python": _status_counts(results, "python_status"),
            "cpp": _status_counts(results, "cpp_status"),
        },
    }
    return summary


def _write_results_csv(path: Path, results: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=RESULT_FIELDNAMES)
        writer.writeheader()
        writer.writerows(results)


def _write_summary_json(path: Path, summary: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def run(
    *,
    csv_path: Path,
    out: Path,
    xyz_dir: Path | None = None,
    limit: int | None = None,
    start_row: int = 1,
    end_row: int | None = None,
    spin_source: str = "closed_shell",
    progress_every: int = 10,
    summary_out: Path | None = None,
) -> list[dict[str, Any]]:
    fieldnames, rows = _load_rows(csv_path)
    _validate_columns(fieldnames, xyz_dir=xyz_dir, spin_source=spin_source)
    selected_rows = _selected_rows(
        rows,
        start_row=start_row,
        end_row=end_row,
        limit=limit,
    )

    config = MolGRConfig()
    started_at = time.perf_counter()
    results: list[dict[str, Any]] = []

    for processed, (row_index, row) in enumerate(selected_rows, start=1):
        result = _benchmark_row(
            row_index,
            row,
            xyz_dir=xyz_dir,
            spin_source=spin_source,
            config=config,
        )
        results.append(result)
        if progress_every and processed % progress_every == 0:
            print(
                "[tmqmg-speed] processed "
                f"{processed}/{len(selected_rows)} rows; latest id={result['id']} "
                f"python_ms={result['python_elapsed_ms']} cpp_ms={result['cpp_elapsed_ms']}",
                file=sys.stderr,
            )

    _write_results_csv(out, results)
    total_wall_seconds = time.perf_counter() - started_at
    summary = _summarize(
        results,
        csv_path=csv_path,
        xyz_dir=xyz_dir,
        spin_source=spin_source,
        config=config,
        wall_seconds=total_wall_seconds,
    )
    _write_summary_json(summary_out or _summary_path_from_output(out), summary)
    print(
        f"Wrote {len(results)} tmQMg speed rows to {out} and summary to "
        f"{summary_out or _summary_path_from_output(out)}",
        file=sys.stderr,
    )
    return results


def main() -> int:
    args = _parse_args()
    run(
        csv_path=args.csv,
        xyz_dir=args.xyz_dir,
        limit=args.limit,
        start_row=args.start_row,
        end_row=args.end_row,
        spin_source=args.spin_source,
        out=args.out,
        summary_out=args.summary_out,
        progress_every=args.progress_every,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
