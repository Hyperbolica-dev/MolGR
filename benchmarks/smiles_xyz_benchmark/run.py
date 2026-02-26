from __future__ import annotations

import argparse
import sys
import time
from importlib import import_module
from pathlib import Path


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from benchmarks.smiles_xyz_benchmark.io import (
    summarize_results,
    write_results_csv,
    write_summary_csv,
)
from benchmarks.smiles_xyz_benchmark.methods import get_method_registry
from benchmarks.smiles_xyz_benchmark.schema import BenchmarkResult
from scripts.molgr_cases_smiles_csv import load_smiles_csv_cases


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run SMILES XYZ benchmark skeleton.")
    parser.add_argument("--input", type=Path, required=True, help="Input CSV path.")
    parser.add_argument("--limit", type=int, default=None, help="Optional case limit.")
    parser.add_argument("--out", type=Path, required=True, help="Output directory.")
    return parser.parse_args()


def _run_case_method(case: dict, method_id: str, method_runner) -> BenchmarkResult:
    check_equivalence = import_module("molgr.utils.equivalence").check_equivalence

    started = time.perf_counter()
    output = method_runner(case)
    method_elapsed_ms = (time.perf_counter() - started) * 1000.0

    breakdown = dict(output.timing_ms_breakdown or {})
    breakdown.setdefault("method_ms", method_elapsed_ms)
    breakdown.setdefault("equivalence_ms", 0.0)

    status = output.status
    error = output.error
    equivalent = output.equivalent
    equivalence_method = output.equivalence_method

    if case.get("provider_error"):
        total_elapsed_ms = (time.perf_counter() - started) * 1000.0
        return BenchmarkResult(
            case_idx=int(case["case_idx"]),
            method_id=method_id,
            input_smiles=str(case["input_smiles"]),
            ground_truth_smiles=case.get("ground_truth_smiles"),
            status="skipped",
            error=str(case["provider_error"]),
            predicted_smiles=None,
            equivalent=None,
            equivalence_method=None,
            timing_ms_total=total_elapsed_ms,
            timing_ms_breakdown=breakdown,
        )

    ground_truth_rdmol = case.get("ground_truth_rdmol")
    if ground_truth_rdmol is not None and output.rdkit_mol is not None:
        equivalence_started = time.perf_counter()
        try:
            is_equivalent, info = check_equivalence(
                ground_truth_rdmol,
                output.rdkit_mol,
                use_chirality=True,
                max_resonance=100,
            )
            equivalent = is_equivalent
            equivalence_method = info.method.value if info.method is not None else None
        except Exception as exc:  # noqa: BLE001
            status = "error"
            equivalence_error = f"equivalence check failed: {exc}"
            error = f"{error}; {equivalence_error}" if error else equivalence_error
            equivalent = None
            equivalence_method = None
        finally:
            breakdown["equivalence_ms"] = (time.perf_counter() - equivalence_started) * 1000.0

    total_elapsed_ms = (time.perf_counter() - started) * 1000.0

    return BenchmarkResult(
        case_idx=int(case["case_idx"]),
        method_id=method_id,
        input_smiles=str(case["input_smiles"]),
        ground_truth_smiles=case.get("ground_truth_smiles"),
        status=status,
        error=error,
        predicted_smiles=output.predicted_smiles,
        equivalent=equivalent,
        equivalence_method=equivalence_method,
        timing_ms_total=total_elapsed_ms,
        timing_ms_breakdown=breakdown,
    )


def run(input_path: Path, out_dir: Path, limit: int | None = None) -> list[BenchmarkResult]:
    cases = load_smiles_csv_cases(input_path=input_path, limit=limit)
    methods = get_method_registry()

    results: list[BenchmarkResult] = []
    for case in cases:
        for method in methods:
            results.append(_run_case_method(case, method.method_id, method.run))

    out_dir.mkdir(parents=True, exist_ok=True)
    write_results_csv(out_dir / "results.csv", results)
    write_summary_csv(out_dir / "summary.csv", summarize_results(results))
    return results


def main() -> int:
    args = _parse_args()
    run(input_path=args.input, out_dir=args.out, limit=args.limit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
