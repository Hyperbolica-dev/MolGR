from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from benchmarks.smiles_xyz_benchmark.schema import BenchmarkResult


SUMMARY_COLUMNS: tuple[str, ...] = (
    "method_id",
    "count",
    "success_count",
    "fail_count",
    "skip_count",
    "avg_ms_total",
    "p50_ms_total",
    "p95_ms_total",
)

RESULT_COLUMNS: tuple[str, ...] = (
    "case_idx",
    "method_id",
    "input_smiles",
    "ground_truth_smiles",
    "status",
    "error",
    "predicted_smiles",
    "equivalent",
    "equivalence_method",
    "timing_ms_total",
)

TIMING_PRIORITY_COLUMNS: tuple[str, ...] = (
    "method_ms",
    "equivalence_ms",
)


def _percentile(values: list[float], percentile: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    index = (len(ordered) - 1) * percentile
    low = int(index)
    high = min(low + 1, len(ordered) - 1)
    weight = index - low
    return ordered[low] * (1.0 - weight) + ordered[high] * weight


def summarize_results(results: list[BenchmarkResult]) -> list[dict[str, str | int | float]]:
    by_method: dict[str, list[BenchmarkResult]] = {}
    for result in results:
        by_method.setdefault(result.method_id, []).append(result)

    summary_rows: list[dict[str, str | int | float]] = []
    for method_id in sorted(by_method):
        rows = by_method[method_id]
        timings = [row.timing_ms_total for row in rows]
        count = len(rows)
        success_count = sum(1 for row in rows if row.equivalent is True)
        fail_count = sum(1 for row in rows if row.status == "error" or row.equivalent is False)
        skip_count = sum(1 for row in rows if row.status == "skipped")

        avg_ms = sum(timings) / count if count else 0.0
        summary_rows.append(
            {
                "method_id": method_id,
                "count": count,
                "success_count": success_count,
                "fail_count": fail_count,
                "skip_count": skip_count,
                "avg_ms_total": avg_ms,
                "p50_ms_total": _percentile(timings, 0.50),
                "p95_ms_total": _percentile(timings, 0.95),
            }
        )
    return summary_rows


def write_results_csv(path: Path, results: list[BenchmarkResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    extra_timing_columns: set[str] = set()
    for result in results:
        breakdown = result.timing_ms_breakdown or {}
        extra_timing_columns.update(key for key in breakdown if key not in TIMING_PRIORITY_COLUMNS)

    timing_columns = list(TIMING_PRIORITY_COLUMNS) + sorted(extra_timing_columns)
    fieldnames = [*RESULT_COLUMNS, *timing_columns, "timing_ms_breakdown_json"]
    rows: list[dict[str, str | int | float | bool | None]] = []

    try:
        for result in results:
            breakdown = result.timing_ms_breakdown or {}
            row: dict[str, str | int | float | bool | None] = {
                "case_idx": result.case_idx,
                "method_id": result.method_id,
                "input_smiles": result.input_smiles,
                "ground_truth_smiles": result.ground_truth_smiles,
                "status": result.status,
                "error": result.error,
                "predicted_smiles": result.predicted_smiles,
                "equivalent": result.equivalent,
                "equivalence_method": result.equivalence_method,
                "timing_ms_total": result.timing_ms_total,
                "timing_ms_breakdown_json": json.dumps(breakdown, ensure_ascii=True),
            }
            for column in timing_columns:
                row[column] = breakdown.get(column)
            rows.append(row)

        pd.DataFrame(rows, columns=pd.Index(fieldnames)).fillna("").to_csv(path, index=False)
    except OSError as exc:
        raise RuntimeError(f"unable to write benchmark results CSV: {path}") from exc


def write_summary_csv(path: Path, summary_rows: list[dict[str, str | int | float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summary_rows, columns=pd.Index(SUMMARY_COLUMNS)).to_csv(path, index=False)
