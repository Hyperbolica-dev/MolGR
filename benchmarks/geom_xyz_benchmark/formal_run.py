"""Streaming GEOM-Drugs formal benchmark runner.

This module deliberately keeps only compact counters, runtimes, canonical-reference keys, and a
bounded review reservoir in memory. Evaluator decisions are written once and frozen.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import heapq
import json
import math
import platform
import statistics
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

from benchmarks._timeout import CaseTimeoutError, case_timeout
from benchmarks.geom_xyz_benchmark.acquire_smoke import _entries, _fixture_record
from benchmarks.smiles_xyz_benchmark.methods.molgr_cpp import MolGRCppMethod
from molgr.utils.equivalence import evaluate_equivalence


SIZE_NAMES = ("01_15", "16_25", "26_35", "36_50", "51_plus")
RESULT_COLUMNS = (
    "case_idx",
    "case_id",
    "molecule_id",
    "conformer_id",
    "heavy_atom_count",
    "status",
    "evaluator_decision",
    "evaluator_relation",
    "evaluator_reason",
    "reason_family",
    "exact_smiles",
    "stereo_equivalent",
    "charge_consistent",
    "radical_consistent",
    "runtime_ms",
)
FAILURE_COLUMNS = RESULT_COLUMNS + (
    "error",
    "reference_smiles",
    "predicted_smiles",
    "xyz",
    "source_conformer_count",
    "relativeenergy_kcal_mol",
)
REVIEW_COLUMNS = FAILURE_COLUMNS + (
    "priority",
    "proposed_verdict",
    "proposed_reason",
    "confidence",
    "matched_rule",
    "canonical_signature",
    "evidence_summary",
    "blockers",
    "needs_human_review",
)


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
    ordered = sorted(values)
    position = (len(ordered) - 1) * percentile
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def _tautomer_diagnostic(reference: Chem.Mol, predicted: Chem.Mol) -> bool:
    enumerator = rdMolStandardize.TautomerEnumerator()
    reference_key = Chem.MolToSmiles(enumerator.Canonicalize(reference), isomericSmiles=False)
    predicted_key = Chem.MolToSmiles(enumerator.Canonicalize(predicted), isomericSmiles=False)
    return reference_key == predicted_key


def _reason_family(
    decision: str,
    reason: str,
    reference: Chem.Mol,
    predicted: Chem.Mol,
    reference_smiles: str,
    predicted_smiles: str,
) -> str:
    if decision == "equivalent":
        return "equivalent"
    reference_neutral_sulfoxide = any(
        atom.GetAtomicNum() == 16
        and atom.GetFormalCharge() == 0
        and any(
            bond.GetBondType() == Chem.BondType.DOUBLE
            and bond.GetOtherAtom(atom).GetAtomicNum() == 8
            and bond.GetOtherAtom(atom).GetFormalCharge() == 0
            for bond in atom.GetBonds()
        )
        for atom in reference.GetAtoms()
    )
    predicted_charge_separated_sulfoxide = any(
        atom.GetAtomicNum() == 16
        and atom.GetFormalCharge() == 1
        and any(
            bond.GetBondType() == Chem.BondType.SINGLE
            and bond.GetOtherAtom(atom).GetAtomicNum() == 8
            and bond.GetOtherAtom(atom).GetFormalCharge() == -1
            for bond in atom.GetBonds()
        )
        for atom in predicted.GetAtoms()
    )
    if (
        decision == "inconclusive"
        and reference_neutral_sulfoxide
        and predicted_charge_separated_sulfoxide
    ):
        return "sulfoxide_neutral_vs_charge_separated"
    if decision == "not_equivalent" and _tautomer_diagnostic(reference, predicted):
        return "suspected_tautomer_protomer"
    normalized = " ".join(reason.split())
    return normalized or "unspecified"


class BoundedReview:
    def __init__(self, limits: dict[str, int]) -> None:
        self.limits = limits
        self.heaps: dict[str, list[tuple[int, str, dict[str, Any]]]] = {
            category: [] for category in limits
        }

    def add(self, category: str, row: dict[str, Any]) -> None:
        limit = self.limits[category]
        score = int.from_bytes(
            hashlib.sha256(f"20260825:{category}:{row['case_id']}".encode()).digest()[:8], "big"
        )
        item = (-score, str(row["case_id"]), row)
        heap = self.heaps[category]
        if len(heap) < limit:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)

    def rows(self) -> Iterable[tuple[str, dict[str, Any]]]:
        for category, heap in self.heaps.items():
            for _, _, row in heap:
                yield category, row


def _case_from_record(case_idx: int, record: dict[str, Any], reference: Chem.Mol) -> dict[str, Any]:
    return {
        "case_idx": case_idx,
        "case_id": record["case_id"],
        "molecule_id": record["molecule_id"],
        "conformer_id": record["conformer_id"],
        "input_smiles": record["reference_smiles"],
        "ground_truth_smiles": record["reference_smiles"],
        "ground_truth_rdmol": reference,
        "xyz_block": record["xyz"],
        "total_charge": 0,
        "total_radical_electrons": 0,
        "spin_multiplicity": 1,
        "source_metadata": record["source_metadata"],
        "provider_error": None,
    }


def run(archive: Path, out_dir: Path, *, timeout: float, expected_eligible: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    started_wall = time.time()
    started_perf = time.perf_counter()
    method = MolGRCppMethod()
    seen_references: set[str] = set()
    acquisition = Counter(
        {
            "source_molecules": 0,
            "eligible": 0,
            "reference_parse_failure": 0,
            "electronic_state_excluded": 0,
            "fragmented_reference": 0,
            "missing_xyz": 0,
            "missing_relativeenergy": 0,
            "reference_xyz_atom_count_mismatch": 0,
            "duplicate_canonical_reference": 0,
        }
    )
    totals: Counter[str] = Counter(
        {
            "total": 0,
            "reconstruction_success": 0,
            "reconstruction_failure": 0,
            "equivalent": 0,
            "not_equivalent": 0,
            "inconclusive": 0,
            "relation_normalized_graph_identity": 0,
            "relation_resonance_equivalence": 0,
            "exact_smiles_match": 0,
            "exact_smiles_mismatch": 0,
            "chirality_equivalent": 0,
            "chirality_not_equivalent": 0,
            "charge_consistent": 0,
            "charge_inconsistent": 0,
            "radical_consistent": 0,
            "radical_inconsistent": 0,
            "timeout": 0,
            "exception": 0,
        }
    )
    evaluator_reasons: Counter[str] = Counter()
    reason_families: Counter[str] = Counter()
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
    runtimes: list[float] = []
    review = BoundedReview(
        {
            "reconstruction_failure": 100,
            "large_molecule_failure": 100,
            "non_equivalent": 250,
            "inconclusive": 150,
        }
    )
    runtime_heap: list[tuple[float, str, dict[str, Any]]] = []

    results_path = out_dir / "results.csv.gz"
    failures_path = out_dir / "failures.csv"
    with gzip.open(
        results_path, "wt", newline="", encoding="utf-8", compresslevel=6
    ) as results_handle, failures_path.open("w", newline="", encoding="utf-8") as failures_handle:
        results_writer = csv.DictWriter(results_handle, fieldnames=RESULT_COLUMNS)
        failures_writer = csv.DictWriter(failures_handle, fieldnames=FAILURE_COLUMNS)
        results_writer.writeheader()
        failures_writer.writeheader()
        for source_smiles, payload in _entries(archive):
            acquisition["source_molecules"] += 1
            record, eligibility = _fixture_record(source_smiles, payload, dataset="drugs")
            acquisition[eligibility] += 1
            if record is None:
                continue
            reference_smiles = str(record["reference_smiles"])
            if reference_smiles in seen_references:
                acquisition["duplicate_canonical_reference"] += 1
                acquisition["eligible"] -= 1
                continue
            seen_references.add(reference_smiles)
            case_idx = acquisition["eligible"]
            reference = Chem.MolFromSmiles(reference_smiles)
            if reference is None:  # Defensive: eligibility parsing already succeeded.
                raise RuntimeError(f"eligible reference no longer parses: {record['case_id']}")
            case = _case_from_record(case_idx, record, reference)
            heavy_atoms = reference.GetNumHeavyAtoms()
            stratum = _size_stratum(heavy_atoms)
            strata[stratum]["n"] += 1
            case_started = time.perf_counter()
            status = "pending"
            error = ""
            decision = ""
            relation = ""
            evaluator_reason = ""
            family = ""
            predicted_smiles = ""
            exact_smiles: bool | None = None
            stereo_equivalent: bool | None = None
            charge_consistent: bool | None = None
            radical_consistent: bool | None = None
            try:
                with case_timeout(timeout or None, f"molgr_cpp case {record['case_id']}"):
                    output = method.run(case)
                status = output.status
                error = output.error or ""
                predicted_smiles = output.predicted_smiles or ""
                if output.rdkit_mol is not None:
                    primary = evaluate_equivalence(output.rdkit_mol, reference, use_chirality=False)
                    stereo = evaluate_equivalence(output.rdkit_mol, reference, use_chirality=True)
                    decision = primary.decision.value
                    relation = primary.relation.value
                    evaluator_reason = primary.reason
                    exact_smiles = Chem.MolToSmiles(
                        reference, isomericSmiles=False
                    ) == Chem.MolToSmiles(output.rdkit_mol, isomericSmiles=False)
                    stereo_equivalent = stereo.equivalent
                    if primary.checks is not None:
                        charge_consistent = primary.checks.formal_charge.passed
                        radical_consistent = primary.checks.radical_electrons.passed
                    family = _reason_family(
                        decision,
                        evaluator_reason,
                        reference,
                        output.rdkit_mol,
                        reference_smiles,
                        predicted_smiles,
                    )
            except CaseTimeoutError as exc:
                status, error = "timeout", str(exc)
            except Exception as exc:  # noqa: BLE001
                status, error = "exception", f"{type(exc).__name__}: {exc}"
            runtime_ms = (time.perf_counter() - case_started) * 1000.0
            runtimes.append(runtime_ms)
            totals["total"] += 1
            totals["reconstruction_success" if status == "ok" else "reconstruction_failure"] += 1
            if status == "timeout":
                totals["timeout"] += 1
            if status == "exception":
                totals["exception"] += 1
            if decision:
                totals[decision] += 1
                strata[stratum][decision] += 1
            if relation:
                totals[f"relation_{relation}"] += 1
            if exact_smiles is not None:
                totals["exact_smiles_match" if exact_smiles else "exact_smiles_mismatch"] += 1
            if stereo_equivalent is not None:
                totals[
                    "chirality_equivalent" if stereo_equivalent else "chirality_not_equivalent"
                ] += 1
            if charge_consistent is not None:
                totals["charge_consistent" if charge_consistent else "charge_inconsistent"] += 1
            if radical_consistent is not None:
                totals["radical_consistent" if radical_consistent else "radical_inconsistent"] += 1
            if status != "ok":
                strata[stratum]["reconstruction_failure"] += 1
            if evaluator_reason:
                evaluator_reasons[" ".join(evaluator_reason.split())] += 1
            if family:
                reason_families[family] += 1

            compact = {
                "case_idx": case_idx,
                "case_id": record["case_id"],
                "molecule_id": record["molecule_id"],
                "conformer_id": record["conformer_id"],
                "heavy_atom_count": heavy_atoms,
                "status": status,
                "evaluator_decision": decision,
                "evaluator_relation": relation,
                "evaluator_reason": evaluator_reason,
                "reason_family": family,
                "exact_smiles": exact_smiles,
                "stereo_equivalent": stereo_equivalent,
                "charge_consistent": charge_consistent,
                "radical_consistent": radical_consistent,
                "runtime_ms": runtime_ms,
            }
            results_writer.writerow(compact)
            detailed = {
                **compact,
                "error": error,
                "reference_smiles": reference_smiles,
                "predicted_smiles": predicted_smiles,
                "xyz": record["xyz"],
                "source_conformer_count": record["source_metadata"]["source_conformer_count"],
                "relativeenergy_kcal_mol": record["source_metadata"]["relativeenergy_kcal_mol"],
            }
            primary_non_pass = status != "ok" or decision != "equivalent"
            if primary_non_pass:
                failures_writer.writerow(detailed)
                category = (
                    "reconstruction_failure"
                    if status != "ok"
                    else "large_molecule_failure"
                    if heavy_atoms >= 51
                    else "non_equivalent"
                    if decision == "not_equivalent"
                    else "inconclusive"
                )
                review.add(category, detailed)
            runtime_item = (runtime_ms, str(record["case_id"]), detailed)
            if len(runtime_heap) < 50:
                heapq.heappush(runtime_heap, runtime_item)
            elif runtime_item > runtime_heap[0]:
                heapq.heapreplace(runtime_heap, runtime_item)
            if totals["total"] % 10_000 == 0:
                results_handle.flush()
                failures_handle.flush()
                elapsed = time.perf_counter() - started_perf
                print(
                    f"processed={totals['total']} source={acquisition['source_molecules']} "
                    f"elapsed_s={elapsed:.1f}",
                    flush=True,
                )

    if acquisition["eligible"] != expected_eligible or totals["total"] != expected_eligible:
        raise RuntimeError(
            f"eligible count mismatch: acquired={acquisition['eligible']} "
            f"evaluated={totals['total']} expected={expected_eligible}"
        )

    review_by_id: dict[str, tuple[set[str], dict[str, Any]]] = {}
    for category, row in review.rows():
        review_by_id.setdefault(str(row["case_id"]), (set(), row))[0].add(category)
    for _, _, row in runtime_heap:
        review_by_id.setdefault(str(row["case_id"]), (set(), row))[0].add("runtime_outlier_top50")
    with (out_dir / "review_cases.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REVIEW_COLUMNS)
        writer.writeheader()
        for case_id in sorted(review_by_id):
            flags, row = review_by_id[case_id]
            graph_review = row["status"] != "ok" or row["evaluator_decision"] != "equivalent"
            signature = hashlib.sha256(
                f"{row['reference_smiles']}>>{row['predicted_smiles']}".encode()
            ).hexdigest()[:20]
            writer.writerow(
                {
                    **row,
                    "priority": ";".join(sorted(flags)),
                    "proposed_verdict": "manual_blocked"
                    if graph_review
                    else "no_graph_review_needed",
                    "proposed_reason": (
                        "primary decision requires human evidence review"
                        if graph_review
                        else "equivalent case selected only as a runtime outlier"
                    ),
                    "confidence": "high",
                    "matched_rule": "",
                    "canonical_signature": signature,
                    "evidence_summary": (
                        f"status={row['status']};decision={row['evaluator_decision']};"
                        f"relation={row['evaluator_relation']};reason_family={row['reason_family']}"
                    ),
                    "blockers": (
                        "no_verified_candidate_to_xyz_mapping;no_approved_representation_rule"
                        if graph_review
                        else ""
                    ),
                    "needs_human_review": graph_review,
                }
            )

    wall_seconds = time.perf_counter() - started_perf
    summary = {
        "protocol": "geom-drugs-formal-molecule-one-conformer-v1",
        "started_unix": started_wall,
        "finished_unix": time.time(),
        "wall_time_seconds": wall_seconds,
        "acquisition_counts": dict(acquisition),
        "statistics": dict(totals),
        "heavy_atom_strata": {name: dict(counts) for name, counts in strata.items()},
        "evaluator_reason_distribution": dict(evaluator_reasons.most_common()),
        "reason_family_distribution": dict(reason_families.most_common()),
        "runtime_ms": {
            "p50": statistics.median(runtimes),
            "p95": _percentile(runtimes, 0.95),
            "p99": _percentile(runtimes, 0.99),
            "max": max(runtimes),
        },
        "review_case_count": len(review_by_id),
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "rdkit": rdBase.rdkitVersion,
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the streaming GEOM-Drugs formal benchmark.")
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--expected-eligible", type=int, default=291_709)
    parser.add_argument("--case-timeout-seconds", type=float, default=10.0)
    args = parser.parse_args()
    run(
        args.archive,
        args.out,
        timeout=args.case_timeout_seconds,
        expected_eligible=args.expected_eligible,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
