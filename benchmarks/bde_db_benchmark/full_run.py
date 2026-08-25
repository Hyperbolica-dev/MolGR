from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import heapq
import importlib.metadata
import json
import platform
import shlex
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from contextlib import suppress
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterator


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from rdkit import rdBase

import molgr
from benchmarks._timeout import case_timeout
from benchmarks.bde_db_benchmark.adapter import BDECase, iter_bde_cases
from benchmarks.bde_db_benchmark.run import BDEResult, _run_case
from benchmarks.smiles_xyz_benchmark.methods.molgr_cpp import MolGRCppMethod


COMPACT_FIELDS = [
    "case_id",
    "source_record_index",
    "stratum",
    "reference_smiles",
    "predicted_smiles",
    "total_charge",
    "spin_multiplicity",
    "radical_site",
    "reconstruction_success",
    "status",
    "failure_kind",
    "error",
    "equivalent",
    "evaluator_decision",
    "evaluator_relation",
    "evaluator_reason",
    "equivalence_method",
    "evaluator_inconclusive",
    "bounded_search_attempted",
    "bounded_search_limit",
    "bounded_search_limit_reached",
    "bounded_search_exhaustive",
    "bounded_search_candidate_count",
    "bounded_search_reference_count",
    "exact_smiles_match",
    "charge_consistent",
    "radical_electron_consistent",
    "atom_order_preserved",
    "atom_identity_guard_reason",
    "formal_radical_atom_index_match",
    "runtime_ms",
]
HEAVY_FIELDS = list(BDEResult.__dataclass_fields__)
REVIEW_LIMIT_PER_CATEGORY = 200
RUNTIME_REVIEW_LIMIT = 100


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_sha() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()


def _git_clean() -> bool:
    return not subprocess.check_output(["git", "status", "--short"], text=True).strip()


def _percentile(values: list[float], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * fraction
    low = int(position)
    high = min(low + 1, len(ordered) - 1)
    weight = position - low
    return ordered[low] * (1.0 - weight) + ordered[high] * weight


class Aggregate:
    def __init__(self) -> None:
        self.counts: Counter[str] = Counter()
        self.decisions: Counter[str] = Counter()
        self.relations: Counter[str] = Counter()
        self.methods: Counter[str] = Counter()
        self.failure_kinds: Counter[str] = Counter()
        self.runtimes: list[float] = []

    def add(self, result: BDEResult) -> None:
        self.counts["total"] += 1
        self.counts["reconstruction_success"] += result.reconstruction_success
        self.counts["reconstruction_failure"] += not result.reconstruction_success
        self.counts["charge_agreement"] += result.charge_consistent is True
        self.counts["formal_radical_electron_agreement"] += (
            result.radical_electron_consistent is True
        )
        self.counts["formal_radical_atom_index_match"] += (
            result.formal_radical_atom_index_match is True
        )
        self.counts["formal_radical_atom_index_mismatch"] += (
            result.formal_radical_atom_index_match is False
        )
        self.counts["formal_radical_atom_index_unknown"] += (
            result.formal_radical_atom_index_match is None
        )
        self.counts["exact_smiles_match"] += result.exact_smiles_match is True
        self.counts["timeouts"] += result.failure_kind in {"timeout", "evaluator_timeout"}
        self.counts["exceptions"] += result.failure_kind in {
            "exception",
            "evaluator_exception",
        }
        self.counts["method_errors"] += result.failure_kind == "method_error"
        if result.evaluator_decision:
            self.decisions[result.evaluator_decision] += 1
        if result.evaluator_relation:
            self.relations[result.evaluator_relation] += 1
        if result.equivalence_method:
            self.methods[result.equivalence_method] += 1
        if result.failure_kind:
            self.failure_kinds[result.failure_kind] += 1
        self.runtimes.append(result.runtime_ms)

    def to_dict(self) -> dict[str, Any]:
        runtime_sum = sum(self.runtimes)
        return {
            **dict(self.counts),
            "equivalent": self.decisions["equivalent"],
            "not_equivalent": self.decisions["not_equivalent"],
            "inconclusive": self.decisions["inconclusive"],
            "normalized_graph_identity": self.relations["normalized_graph_identity"],
            "resonance_equivalence": self.relations["resonance_equivalence"],
            "decisions": dict(sorted(self.decisions.items())),
            "relations": dict(sorted(self.relations.items())),
            "methods": dict(sorted(self.methods.items())),
            "failure_kinds": dict(sorted(self.failure_kinds.items())),
            "runtime_mean_ms": runtime_sum / len(self.runtimes) if self.runtimes else 0.0,
            "runtime_p50_ms": _percentile(self.runtimes, 0.50),
            "runtime_p95_ms": _percentile(self.runtimes, 0.95),
            "runtime_p99_ms": _percentile(self.runtimes, 0.99),
            "runtime_max_ms": max(self.runtimes, default=0.0),
        }


class ReviewCollector:
    def __init__(self, seed: int) -> None:
        self.seed = seed
        self.critical: dict[int, tuple[str, BDEResult]] = {}
        self.samples: dict[str, list[tuple[int, int, BDEResult]]] = defaultdict(list)
        self.runtime: list[tuple[float, int, BDEResult]] = []

    def _score(self, category: str, result: BDEResult) -> int:
        payload = f"{self.seed}:{category}:{result.source_record_index}:{result.case_id}".encode()
        return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")

    def _sample(self, category: str, result: BDEResult) -> None:
        score = self._score(category, result)
        heap = self.samples[category]
        item = (-score, result.source_record_index, result)
        if len(heap) < REVIEW_LIMIT_PER_CATEGORY:
            heapq.heappush(heap, item)
        elif score < -heap[0][0]:
            heapq.heapreplace(heap, item)

    def add(self, result: BDEResult) -> None:
        critical_reasons = []
        if result.evaluator_decision == "not_equivalent":
            critical_reasons.append("non_equivalent")
        if result.evaluator_decision == "inconclusive":
            critical_reasons.append("inconclusive")
        if result.status != "ok":
            critical_reasons.append(result.failure_kind or "failure")
        if result.charge_consistent is False or result.radical_electron_consistent is False:
            critical_reasons.append("charge_or_radical_electron_mismatch")
        if critical_reasons:
            self.critical[result.source_record_index] = (";".join(critical_reasons), result)
        if result.evaluator_relation == "resonance_equivalence":
            self._sample("resonance_equivalent", result)
        if result.formal_radical_atom_index_match is False:
            self._sample("formal_radical_atom_index_mismatch", result)
        if result.formal_radical_atom_index_match is None and result.reconstruction_success:
            self._sample("formal_radical_atom_index_unknown", result)
        item = (result.runtime_ms, result.source_record_index, result)
        if len(self.runtime) < RUNTIME_REVIEW_LIMIT:
            heapq.heappush(self.runtime, item)
        elif result.runtime_ms > self.runtime[0][0]:
            heapq.heapreplace(self.runtime, item)

    def rows(self) -> list[dict[str, Any]]:
        collected: dict[int, tuple[int, set[str], BDEResult]] = {}

        def add(priority: int, reason: str, result: BDEResult) -> None:
            current = collected.get(result.source_record_index)
            if current is None:
                collected[result.source_record_index] = (priority, {reason}, result)
            else:
                collected[result.source_record_index] = (
                    min(priority, current[0]),
                    current[1] | {reason},
                    result,
                )

        for reason, result in self.critical.values():
            add(1, reason, result)
        priorities = {
            "resonance_equivalent": 4,
            "formal_radical_atom_index_mismatch": 5,
            "formal_radical_atom_index_unknown": 6,
        }
        for category, heap in self.samples.items():
            for _, _, result in heap:
                add(priorities[category], category, result)
        for _, _, result in self.runtime:
            add(8, "runtime_top_100", result)
        ordered = sorted(
            collected.values(),
            key=lambda item: (item[0], -item[2].runtime_ms, item[2].source_record_index),
        )
        return [
            {
                "review_priority": priority,
                "review_reasons": ";".join(sorted(reasons)),
                **asdict(result),
            }
            for priority, reasons, result in ordered
        ]


def _environment() -> dict[str, str]:
    packages = {}
    for name in ("rdkit", "openbabel", "openbabel-wheel"):
        with suppress(importlib.metadata.PackageNotFoundError):
            packages[name] = importlib.metadata.version(name)
    return {
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "molgr": molgr.__version__,
        "rdkit_runtime": rdBase.rdkitVersion,
        **packages,
    }


def _write_csv(path: Path, rows: Iterator[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _copy_file(source: Path, destination: Path) -> None:
    with source.open("rb") as input_stream, destination.open("wb") as output_stream:
        shutil.copyfileobj(input_stream, output_stream, length=1024 * 1024)


def _paper_bundle(
    run_dir: Path,
    input_path: Path,
    summary: dict[str, Any],
    command: str,
) -> Path:
    bundle = run_dir / "bde_db_paper_export"
    bundle.mkdir(parents=True, exist_ok=True)
    for name in ("results.csv.gz", "failures.csv", "review_cases.csv"):
        _copy_file(run_dir / name, bundle / name)
    (bundle / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    provenance = f"""# BDE-db full benchmark provenance

- Dataset: `{input_path.name}`
- Dataset SHA-256: `{summary['dataset']['sha256']}`
- Dataset records: {summary['dataset']['record_count']}
- MolGR git SHA: `{summary['molgr_git_sha']}`
- Evaluator: evaluator v1, `evaluate_equivalence(Candidate, Reference)`
- Candidate: fresh internal `MolGRCppMethod` RDKit molecule
- Reference: official SDF graph with benchmark-equivalent hydrogen removal
- Primary policy: `use_chirality=False`, `max_resonance=100`
- Exact command: `{command}`

The compact case-level results were copied directly from the frozen full run. The paper
export does not rerun reconstruction or equivalence evaluation. XYZ, bond dumps, and full
source metadata are restricted to failure and bounded manual-review records.
"""
    (bundle / "provenance.md").write_text(provenance, encoding="utf-8")
    files = {}
    for path in sorted(bundle.iterdir()):
        if path.name == "manifest.json":
            continue
        files[path.name] = {"sha256": _sha256(path), "size": path.stat().st_size}
    manifest = {
        "bundle_schema": "bde-db-paper-export-v1",
        "created_utc": summary["completed_utc"],
        "molgr_git_sha": summary["molgr_git_sha"],
        "dataset_sha256": summary["dataset"]["sha256"],
        "files": files,
    }
    (bundle / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    return bundle


def run_full(
    input_path: Path,
    out_dir: Path,
    *,
    expected_records: int = 289_639,
    seed: int = 0,
    timeout_seconds: float | None = 5.0,
    command_record: str | None = None,
) -> dict[str, Any]:
    if not _git_clean():
        raise RuntimeError("full benchmark requires a clean Git worktree")
    out_dir.mkdir(parents=True, exist_ok=False)
    git_sha = _git_sha()
    command = command_record or shlex.join(sys.argv)
    dataset_sha256 = _sha256(input_path)
    aggregates: dict[str, Aggregate] = defaultdict(Aggregate)
    review = ReviewCollector(seed)
    loader_failures: list[dict[str, str | int]] = []
    failures: list[dict[str, Any]] = []
    method = MolGRCppMethod()
    case_iter = iter_bde_cases(input_path)
    buffered: list[tuple[BDECase | None, dict[str, str | int] | None]] = []
    warmup_case = None
    for item in case_iter:
        buffered.append(item)
        if item[0] is not None and item[0].stratum == "closed_shell":
            warmup_case = item[0]
            break
    warmup_started = time.perf_counter()
    warmup_status = "skipped"
    warmup_error = None
    if warmup_case is not None:
        try:
            with case_timeout(timeout_seconds, f"warm-up case {warmup_case.case_id}"):
                warmup_output = method.run(warmup_case.to_method_case())
            warmup_status = warmup_output.status
            warmup_error = warmup_output.error
        except Exception as exc:
            warmup_status = "error"
            warmup_error = f"{type(exc).__name__}: {exc}"
    warmup_ms = (time.perf_counter() - warmup_started) * 1000.0
    started = time.time()
    processed = 0

    def items() -> Iterator[tuple[BDECase | None, dict[str, str | int] | None]]:
        yield from buffered
        yield from case_iter

    with gzip.open(out_dir / "results.csv.gz", "wt", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=COMPACT_FIELDS)
        writer.writeheader()
        for case, loader_failure in items():
            if loader_failure is not None:
                loader_failures.append(loader_failure)
                continue
            assert case is not None
            result = _run_case(case, method, timeout_seconds)
            row = asdict(result)
            writer.writerow({field: row[field] for field in COMPACT_FIELDS})
            processed += 1
            aggregates["overall"].add(result)
            aggregates[result.stratum].add(result)
            review.add(result)
            if result.status != "ok":
                failures.append(row)
            if processed % 1000 == 0:
                stream.flush()
                elapsed = time.time() - started
                print(
                    f"processed={processed} elapsed_s={elapsed:.1f} "
                    f"rate={processed / elapsed:.1f}/s",
                    flush=True,
                )
    if processed + len(loader_failures) != expected_records:
        raise RuntimeError(
            f"record count mismatch: processed={processed}, loader_failures={len(loader_failures)}, "
            f"expected={expected_records}"
        )
    failure_rows = failures + [
        {
            **dict.fromkeys(HEAVY_FIELDS, ""),
            "source_record_index": failure["record_index"],
            "status": "loader_error",
            "failure_kind": "loader_error",
            "error": failure["error"],
        }
        for failure in loader_failures
    ]
    _write_csv(out_dir / "failures.csv", iter(failure_rows), HEAVY_FIELDS)
    review_rows = review.rows()
    _write_csv(
        out_dir / "review_cases.csv",
        iter(review_rows),
        ["review_priority", "review_reasons", *HEAVY_FIELDS],
    )
    completed = time.time()
    summary = {
        "schema": "bde-db-full-run-v1",
        "molgr_git_sha": git_sha,
        "git_worktree_clean_at_start": True,
        "dataset": {
            "filename": input_path.name,
            "path": str(input_path.resolve()),
            "sha256": dataset_sha256,
            "record_count": expected_records,
        },
        "protocol": {
            "method": "molgr_cpp",
            "evaluator": "v1",
            "candidate": "fresh internal MolGR RDKit molecule",
            "reference": "official SDF graph",
            "use_chirality": False,
            "max_resonance": 100,
            "seed": seed,
            "case_timeout_seconds": timeout_seconds,
            "compact_results": True,
            "review_limit_per_diagnostic_category": REVIEW_LIMIT_PER_CATEGORY,
            "runtime_review_limit": RUNTIME_REVIEW_LIMIT,
        },
        "command": command,
        "environment": _environment(),
        "warmup": {
            "case_id": warmup_case.case_id if warmup_case else None,
            "status": warmup_status,
            "error": warmup_error,
            "runtime_ms": warmup_ms,
        },
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(started)),
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(completed)),
        "wall_time_seconds": completed - started,
        "loader_failure_count": len(loader_failures),
        "review_case_count": len(review_rows),
        "metrics": {name: aggregate.to_dict() for name, aggregate in aggregates.items()},
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    bundle = _paper_bundle(out_dir, input_path, summary, command)
    summary["paper_export_bundle"] = str(bundle.resolve())
    return summary


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the streaming BDE-db full benchmark.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--expected-records", type=int, default=289_639)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--case-timeout-seconds", type=float, default=5.0)
    parser.add_argument("--command-record")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    summary = run_full(
        args.input,
        args.out,
        expected_records=args.expected_records,
        seed=args.seed,
        timeout_seconds=args.case_timeout_seconds or None,
        command_record=args.command_record,
    )
    print(json.dumps(summary["metrics"]["overall"], indent=2, sort_keys=True))
    print(f"paper_export_bundle={summary['paper_export_bundle']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
