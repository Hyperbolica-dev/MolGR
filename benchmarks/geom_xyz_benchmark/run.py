from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
import time
from pathlib import Path
from typing import Any


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from rdkit import Chem

from benchmarks._timeout import CaseTimeoutError, case_timeout
from benchmarks.geom_xyz_benchmark.adapter import load_cases
from benchmarks.smiles_xyz_benchmark.methods.molgr_cpp import MolGRCppMethod
from molgr.utils.equivalence import evaluate_equivalence


def _args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the molecule-level GEOM READY benchmark.")
    parser.add_argument("--input", type=Path, required=True, help="Derived GEOM JSONL fixture.")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--seed", type=int, default=20260825)
    parser.add_argument("--case-timeout-seconds", type=float, default=10.0)
    return parser.parse_args()


def _canonical(mol: Chem.Mol, *, stereo: bool) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=stereo)


def _size_stratum(heavy_atoms: int) -> str:
    if heavy_atoms <= 15:
        return "01_15"
    if heavy_atoms <= 25:
        return "16_25"
    if heavy_atoms <= 35:
        return "26_35"
    if heavy_atoms <= 50:
        return "36_50"
    return "51_plus"


def _percentile(values: list[float], percentile: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    position = (len(ordered) - 1) * percentile
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def run(
    input_path: Path, out_dir: Path, *, limit: int | None, seed: int, timeout: float
) -> list[dict[str, Any]]:
    method = MolGRCppMethod()
    rows: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    for case in load_cases(input_path, limit=limit, seed=seed):
        started = time.perf_counter()
        row = {key: case.get(key) for key in ("case_idx", "case_id", "molecule_id", "conformer_id")}
        source_metadata = case.get("source_metadata", {})
        reference = case.get("ground_truth_rdmol")
        row.update(
            reference_smiles=case.get("ground_truth_smiles"),
            heavy_atom_count=reference.GetNumHeavyAtoms() if reference is not None else None,
            source_conformer_count=source_metadata.get("source_conformer_count"),
            relativeenergy_kcal_mol=source_metadata.get("relativeenergy_kcal_mol"),
            status="loader_failure" if case.get("provider_error") else "pending",
            error=case.get("provider_error"),
            evaluator_decision=None,
            evaluator_relation=None,
            graph_equivalent=None,
            exact_smiles=None,
            stereo_equivalent=None,
            charge_consistent=None,
            radical_consistent=None,
            predicted_smiles=None,
        )
        if not case.get("provider_error"):
            try:
                with case_timeout(timeout or None, f"molgr_cpp case {case['case_id']}"):
                    output = method.run(case)
                row["status"], row["error"] = output.status, output.error
                row["predicted_smiles"] = output.predicted_smiles
                if output.rdkit_mol is not None:
                    reference = case["ground_truth_rdmol"]
                    primary = evaluate_equivalence(output.rdkit_mol, reference, use_chirality=False)
                    stereo = evaluate_equivalence(output.rdkit_mol, reference, use_chirality=True)
                    row["evaluator_decision"] = primary.decision.value
                    row["evaluator_relation"] = primary.relation.value
                    row["graph_equivalent"] = primary.equivalent
                    row["stereo_equivalent"] = stereo.equivalent
                    row["exact_smiles"] = _canonical(reference, stereo=False) == _canonical(
                        output.rdkit_mol, stereo=False
                    )
                    row["charge_consistent"] = primary.checks.formal_charge.passed
                    row["radical_consistent"] = primary.checks.radical_electrons.passed
            except CaseTimeoutError as exc:
                row["status"], row["error"] = "timeout", str(exc)
            except Exception as exc:  # noqa: BLE001
                row["status"], row["error"] = "exception", f"{type(exc).__name__}: {exc}"
        row["runtime_ms"] = (time.perf_counter() - started) * 1000
        rows.append(row)
        if row["status"] != "ok" or row["graph_equivalent"] is not True:
            failures.append(row)

    out_dir.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0]) if rows else []
    for name, payload in (("results.csv", rows), ("failures.csv", failures)):
        with (out_dir / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(payload)
    successful = [row for row in rows if row["status"] == "ok"]
    decisions = ("equivalent", "not_equivalent", "inconclusive")
    relations = ("normalized_graph_identity", "resonance_equivalence")
    strata: dict[str, dict[str, int]] = {}
    for name in ("01_15", "16_25", "26_35", "36_50", "51_plus"):
        members = [row for row in rows if _size_stratum(int(row["heavy_atom_count"])) == name]
        strata[name] = {
            "n": len(members),
            "reconstruction_failure": sum(row["status"] != "ok" for row in members),
            "equivalent": sum(row["evaluator_decision"] == "equivalent" for row in members),
            "not_equivalent": sum(row["evaluator_decision"] == "not_equivalent" for row in members),
            "inconclusive": sum(row["evaluator_decision"] == "inconclusive" for row in members),
        }
    runtimes = [float(row["runtime_ms"]) for row in rows]
    summary = {
        "protocol": "geom-molecule-one-conformer-v1",
        "seed": seed,
        "case_count": len(rows),
        "heavy_atom_strata": strata,
        "reconstruction_success": len(successful),
        "reconstruction_failure": len(rows) - len(successful),
        **{
            f"decision_{decision}": sum(row["evaluator_decision"] == decision for row in rows)
            for decision in decisions
        },
        **{
            f"relation_{relation}": sum(row["evaluator_relation"] == relation for row in rows)
            for relation in relations
        },
        "graph_equivalent": sum(row["graph_equivalent"] is True for row in rows),
        "exact_smiles": sum(row["exact_smiles"] is True for row in rows),
        "exact_smiles_mismatch": sum(row["exact_smiles"] is False for row in rows),
        "stereo_equivalent": sum(row["stereo_equivalent"] is True for row in rows),
        "stereo_not_equivalent": sum(row["stereo_equivalent"] is False for row in rows),
        "charge_consistent": sum(row["charge_consistent"] is True for row in rows),
        "charge_inconsistent": sum(row["charge_consistent"] is False for row in rows),
        "radical_consistent": sum(row["radical_consistent"] is True for row in rows),
        "radical_inconsistent": sum(row["radical_consistent"] is False for row in rows),
        "timeout": sum(row["status"] == "timeout" for row in rows),
        "exception": sum(row["status"] == "exception" for row in rows),
        "failures": len(failures),
        "runtime_ms_total": sum(float(row["runtime_ms"]) for row in rows),
        "runtime_ms": {
            "p50": statistics.median(runtimes) if runtimes else math.nan,
            "p95": _percentile(runtimes, 0.95),
            "p99": _percentile(runtimes, 0.99),
            "max": max(runtimes, default=math.nan),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return rows


def main() -> int:
    args = _args()
    run(args.input, args.out, limit=args.limit, seed=args.seed, timeout=args.case_timeout_seconds)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
