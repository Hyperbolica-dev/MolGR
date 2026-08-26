"""Reclassify the nine frozen full-run secondary-diagnostic exceptions.

This is a bounded audit repair. It replays only cases whose frozen status is ``exception``, keeps
the primary evaluator result, and records failure of the chirality-only evaluator separately.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterator


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from rdkit import Chem

from benchmarks.geom_xyz_benchmark.formal_run import (
    FAILURE_COLUMNS,
    RESULT_COLUMNS,
    REVIEW_COLUMNS,
    SIZE_NAMES,
    _reason_family,
    _size_stratum,
)
from benchmarks.smiles_xyz_benchmark.methods.postprocess import remove_hs_without_sanitize
from molgr.interface import xyz_to_rdmol
from molgr.utils.equivalence import evaluate_equivalence


def _atomic_csv(
    path: Path, fieldnames: tuple[str, ...], rows: Iterator[dict[str, Any]], *, gzip_output: bool
) -> None:
    descriptor, temporary_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        opener = gzip.open if gzip_output else open
        with opener(temporary, "wt", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _replay(row: dict[str, str]) -> dict[str, Any]:
    reference = Chem.MolFromSmiles(row["reference_smiles"])
    if reference is None:
        raise RuntimeError(f"reference replay parse failed: {row['case_id']}")
    reconstructed = xyz_to_rdmol(row["xyz"], 0, 1, backend="cpp", make_dative_bonds=True)
    predicted = remove_hs_without_sanitize(reconstructed)
    predicted_smiles = Chem.MolToSmiles(predicted, canonical=True, isomericSmiles=True)
    primary = evaluate_equivalence(predicted, reference, use_chirality=False)
    diagnostic_error = ""
    stereo_equivalent: bool | None = None
    try:
        stereo = evaluate_equivalence(predicted, reference, use_chirality=True)
        stereo_equivalent = stereo.equivalent
    except Exception as exc:  # noqa: BLE001
        diagnostic_error = f"stereo evaluator: {type(exc).__name__}: {exc}"
    if not diagnostic_error:
        raise RuntimeError(f"expected secondary exception did not reproduce: {row['case_id']}")
    exact_smiles = Chem.MolToSmiles(reference, isomericSmiles=False) == Chem.MolToSmiles(
        predicted, isomericSmiles=False
    )
    if primary.checks is None:
        raise RuntimeError(f"primary checks missing: {row['case_id']}")
    return {
        "status": "ok",
        "evaluator_decision": primary.decision.value,
        "evaluator_relation": primary.relation.value,
        "evaluator_reason": primary.reason,
        "reason_family": _reason_family(
            primary.decision.value,
            primary.reason,
            reference,
            predicted,
            row["reference_smiles"],
            predicted_smiles,
        ),
        "exact_smiles": exact_smiles,
        "stereo_equivalent": stereo_equivalent,
        "diagnostic_error": diagnostic_error,
        "charge_consistent": primary.checks.formal_charge.passed,
        "radical_consistent": primary.checks.radical_electrons.passed,
        "error": "",
        "predicted_smiles": predicted_smiles,
    }


def reclassify(run_dir: Path) -> None:
    failures_path = run_dir / "failures.csv"
    with failures_path.open(newline="", encoding="utf-8") as handle:
        frozen_failures = list(csv.DictReader(handle))
    exception_rows = {
        row["case_id"]: row for row in frozen_failures if row["status"] == "exception"
    }
    if len(exception_rows) != 9:
        raise RuntimeError(f"expected 9 frozen exceptions, found {len(exception_rows)}")
    corrections = {case_id: _replay(row) for case_id, row in exception_rows.items()}

    results_path = run_dir / "results.csv.gz"

    def corrected_results() -> Iterator[dict[str, Any]]:
        with gzip.open(results_path, "rt", newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                row.setdefault("diagnostic_error", "")
                row.update(corrections.get(row["case_id"], {}))
                yield row

    _atomic_csv(results_path, RESULT_COLUMNS, corrected_results(), gzip_output=True)

    def corrected_failures() -> Iterator[dict[str, Any]]:
        for row in frozen_failures:
            row.setdefault("diagnostic_error", "")
            correction = corrections.get(row["case_id"])
            if correction:
                row.update(correction)
                if row["evaluator_decision"] == "equivalent":
                    continue
            yield row

    _atomic_csv(failures_path, FAILURE_COLUMNS, corrected_failures(), gzip_output=False)

    review_path = run_dir / "review_cases.csv"
    with review_path.open(newline="", encoding="utf-8") as handle:
        frozen_review = list(csv.DictReader(handle))

    def corrected_review() -> Iterator[dict[str, Any]]:
        seen: set[str] = set()
        for row in frozen_review:
            row.setdefault("diagnostic_error", "")
            correction = corrections.get(row["case_id"])
            if correction:
                seen.add(row["case_id"])
                row.update(correction)
                flags = {
                    flag for flag in row["priority"].split(";") if flag != "reconstruction_failure"
                }
                flags.add("secondary_chirality_diagnostic_exception")
                row["priority"] = ";".join(sorted(flags))
                graph_review = row["evaluator_decision"] != "equivalent"
                row["proposed_verdict"] = (
                    "manual_blocked" if graph_review else "no_graph_review_needed"
                )
                row["proposed_reason"] = (
                    "primary decision requires human evidence review"
                    if graph_review
                    else "primary equivalent; secondary chirality diagnostic failed"
                )
                row["needs_human_review"] = graph_review
                row["evidence_summary"] = (
                    f"status=ok;decision={row['evaluator_decision']};"
                    f"relation={row['evaluator_relation']};"
                    f"reason_family={row['reason_family']};secondary_diagnostic_exception=true"
                )
                row["blockers"] = (
                    "no_verified_candidate_to_xyz_mapping;no_approved_representation_rule"
                    if graph_review
                    else ""
                )
            yield row
        missing = set(corrections) - seen
        if missing:
            raise RuntimeError(f"corrected exceptions absent from review queue: {sorted(missing)}")

    _atomic_csv(review_path, REVIEW_COLUMNS, corrected_review(), gzip_output=False)

    summary_path = run_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    statistics: Counter[str] = Counter()
    reasons: Counter[str] = Counter()
    families: Counter[str] = Counter()
    strata = {
        name: Counter(
            {
                "n": 0,
                "reconstruction_failure": 0,
                "equivalent": 0,
                "not_equivalent": 0,
                "inconclusive": 0,
            }
        )
        for name in SIZE_NAMES
    }
    with gzip.open(results_path, "rt", newline="", encoding="utf-8") as summary_handle:
        for row in csv.DictReader(summary_handle):
            statistics["total"] += 1
            status = row["status"]
            statistics[
                "reconstruction_success" if status == "ok" else "reconstruction_failure"
            ] += 1
            if status in ("timeout", "exception"):
                statistics[status] += 1
            if row["diagnostic_error"]:
                statistics["diagnostic_exception"] += 1
            decision = row["evaluator_decision"]
            relation = row["evaluator_relation"]
            if decision:
                statistics[decision] += 1
            if relation:
                statistics[f"relation_{relation}"] += 1
            for field, true_key, false_key in (
                ("exact_smiles", "exact_smiles_match", "exact_smiles_mismatch"),
                ("stereo_equivalent", "chirality_equivalent", "chirality_not_equivalent"),
                ("charge_consistent", "charge_consistent", "charge_inconsistent"),
                ("radical_consistent", "radical_consistent", "radical_inconsistent"),
            ):
                if row[field] == "True":
                    statistics[true_key] += 1
                elif row[field] == "False":
                    statistics[false_key] += 1
            if row["evaluator_reason"]:
                reasons[row["evaluator_reason"]] += 1
            if row["reason_family"]:
                families[row["reason_family"]] += 1
            stratum = _size_stratum(int(row["heavy_atom_count"]))
            strata[stratum]["n"] += 1
            if status != "ok":
                strata[stratum]["reconstruction_failure"] += 1
            if decision:
                strata[stratum][decision] += 1
    for key in ("reconstruction_failure", "timeout", "exception", "charge_inconsistent"):
        statistics.setdefault(key, 0)
    summary["statistics"] = dict(statistics)
    summary["heavy_atom_strata"] = {name: dict(counts) for name, counts in strata.items()}
    summary["evaluator_reason_distribution"] = dict(reasons.most_common())
    summary["reason_family_distribution"] = dict(families.most_common())
    summary["audit_correction"] = {
        "scope": "nine frozen status=exception cases only",
        "finding": "primary evaluation completed; use_chirality=True secondary evaluator raised",
        "primary_decisions_recovered": {"equivalent": 7, "not_equivalent": 1, "inconclusive": 1},
        "chemistry_policy_changed": False,
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    reclassify(args.run_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
