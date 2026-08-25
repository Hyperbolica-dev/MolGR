#!/usr/bin/env python3
"""Resumable, internal-molecule tmQMg evaluator-v1 classification refresh.

The runner never persists an RDKit molecule or reparses serialized Candidate
SMILES for Candidate-to-Reference classification.  Each native reconstruction
chunk is evaluated immediately, reduced to lightweight JSONL records, and
released before the next chunk.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sqlite3
import statistics
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rdkit import Chem, RDLogger  # noqa: E402

from benchmarks.tmqmg_xyz_benchmark.methods import get_method_registry  # noqa: E402
from benchmarks.tmqmg_xyz_benchmark.run import _build_case  # noqa: E402
from molgr.utils.equivalence import EquivalenceResult, evaluate_equivalence  # noqa: E402


DEFAULT_CSV = ROOT / ".local" / "tmQMg" / "data" / "tmQMg_properties_and_targets.csv"
DEFAULT_XYZ_DIR = ROOT / ".local" / "tmQMg" / "data" / "xyz"
DEFAULT_OUT_DIR = ROOT / ".local" / "tmqmg_evaluator_v1_refresh"
DEFAULT_REVIEW_DB = ROOT / ".local" / "molgr_review" / "review.sqlite"
EVALUATOR_PATH = ROOT / "src" / "molgr" / "utils" / "equivalence.py"
PHASE_CAPS = (100, 500, 1000)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--xyz-dir", type=Path, default=DEFAULT_XYZ_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--review-db", type=Path, default=DEFAULT_REVIEW_DB)
    parser.add_argument("--chunk-size", type=int, default=750)
    parser.add_argument("--max-workers", type=int, default=8)
    parser.add_argument(
        "--cap-1000-threshold",
        type=int,
        default=5000,
        help="Run cap 1000 only when no more than this many cap-500 cases remain.",
    )
    resume = parser.add_mutually_exclusive_group()
    resume.add_argument("--resume", dest="resume", action="store_true")
    resume.add_argument("--no-resume", dest="resume", action="store_false")
    parser.set_defaults(resume=True)
    parser.add_argument(
        "--stop-after-chunks",
        type=int,
        default=None,
        help="Testing aid: stop after writing this many new chunks across all phases.",
    )
    args = parser.parse_args()
    if args.chunk_size < 1:
        parser.error("--chunk-size must be positive")
    if args.max_workers < 1:
        parser.error("--max-workers must be positive")
    return args


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _fingerprint(args: argparse.Namespace, total: int) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "evaluator_policy": "v1",
        "evaluator_sha256": _sha256(EVALUATOR_PATH),
        "csv_path": str(args.csv.resolve()),
        "csv_sha256": _sha256(args.csv),
        "xyz_dir": str(args.xyz_dir.resolve()),
        "total_rows": total,
        "chunk_size": args.chunk_size,
        "spin_policy": "tmqmg_closed_shell",
        "use_chirality": False,
        "phase_caps": list(PHASE_CAPS),
        "cap_1000_threshold": args.cap_1000_threshold,
    }


def _atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def _atomic_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def _load_rows(path: Path) -> list[tuple[int, dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return [(index, row) for index, row in enumerate(csv.DictReader(stream), start=1)]


def _chunks(values: Sequence[Any], size: int) -> Iterator[Sequence[Any]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def _chunk_path(out_dir: Path, cap: int, ordinal: int, items: Sequence[Any]) -> Path:
    first = int(items[0][0])
    last = int(items[-1][0])
    return out_dir / f"cap_{cap}" / f"chunk_{ordinal:05d}_{first:06d}_{last:06d}.jsonl"


def _read_jsonl(path: Path) -> Iterator[dict[str, Any]]:
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.strip():
                yield json.loads(line)


def _all_phase_records(out_dir: Path, cap: int) -> Iterator[dict[str, Any]]:
    for path in sorted((out_dir / f"cap_{cap}").glob("chunk_*.jsonl")):
        yield from _read_jsonl(path)


def _percentile(values: Sequence[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    return ordered[max(0, min(len(ordered) - 1, int(fraction * len(ordered)) - 1))]


def _result_fields(result: EquivalenceResult, elapsed_ms: float) -> dict[str, Any]:
    bounded = result.bounded_search
    return {
        "decision": result.decision.value,
        "relation": result.relation.value,
        "method": result.method.value if result.method else None,
        "reason": result.reason,
        "contradictions": list(result.contradictions),
        "bounded_attempted": bool(bounded and bounded.attempted),
        "bounded_limit": bounded.limit if bounded else None,
        "bounded_limit_reached": bool(bounded and bounded.limit_reached),
        "bounded_exhaustive": bounded.exhaustive if bounded else None,
        "equivalence_ms": round(elapsed_ms, 6),
    }


def _is_bounded_only(record: dict[str, Any]) -> bool:
    return bool(
        record.get("decision") == "inconclusive"
        and record.get("bounded_limit_reached")
        and not record.get("contradictions")
        and str(record.get("reason", "")).startswith("Inconclusive: bounded resonance search")
    )


def _roundtrip(internal: Chem.Mol, predicted_smiles: str) -> dict[str, Any]:
    try:
        reparsed = Chem.MolFromSmiles(predicted_smiles)
    except Exception as exc:  # noqa: BLE001
        return {
            "roundtrip_status": "reparse_failed",
            "roundtrip_decision": None,
            "roundtrip_reason": f"{type(exc).__name__}: {exc}",
        }
    if reparsed is None:
        return {
            "roundtrip_status": "reparse_failed",
            "roundtrip_decision": None,
            "roundtrip_reason": "predicted SMILES could not be reparsed",
        }
    started = time.perf_counter()
    result = evaluate_equivalence(internal, reparsed, use_chirality=False, max_resonance=100)
    return {
        "roundtrip_status": "preserved" if result.equivalent else "roundtrip_changed",
        "roundtrip_decision": result.decision.value,
        "roundtrip_reason": result.reason,
        "roundtrip_ms": round((time.perf_counter() - started) * 1000.0, 6),
    }


def _classify_chunk(
    items: Sequence[tuple[int, dict[str, str]]],
    *,
    cap: int,
    xyz_dir: Path,
    max_workers: int,
) -> list[dict[str, Any]]:
    cases = [_build_case(index, row, xyz_dir=xyz_dir) for index, row in items]
    method = {item.method_id: item for item in get_method_registry(("molgr_cpp",))}["molgr_cpp"]
    chunk_started = time.perf_counter()
    outputs = method.run_batch(cases, max_workers=max_workers)
    reconstruction_wall_ms = (time.perf_counter() - chunk_started) * 1000.0
    per_case_reconstruction_ms = reconstruction_wall_ms / len(items)
    records: list[dict[str, Any]] = []
    for (row_index, row), case in zip(items, cases):
        case_id = row.get("id", "").strip()
        output = outputs.get(row_index)
        record: dict[str, Any] = {
            "row_index": row_index,
            "case_id": case_id,
            "cap": cap,
            "charge": row.get("charge", ""),
            "reconstruction_ms_share": round(per_case_reconstruction_ms, 6),
            "reference_status": "ok",
            "reference_error": case.get("reference_error"),
            "reconstruction_status": "ok",
            "reconstruction_error": None,
        }
        reference_smiles = row.get("smiles", "").strip()
        if not reference_smiles:
            record["reference_status"] = "unavailable"
        reference = Chem.MolFromSmiles(reference_smiles) if reference_smiles else None
        if reference_smiles and reference is None:
            record["reference_status"] = "unparseable"
        if output is None or output.status != "ok" or output.rdkit_mol is None:
            record["reconstruction_status"] = "failed"
            record["reconstruction_error"] = (
                output.error if output is not None else "native batch returned no result"
            )
            records.append(record)
            continue
        internal = output.rdkit_mol
        if output.predicted_smiles:
            record.update(_roundtrip(internal, output.predicted_smiles))
        else:
            record.update(
                {
                    "roundtrip_status": "reparse_failed",
                    "roundtrip_decision": None,
                    "roundtrip_reason": "reconstruction produced no serialized SMILES",
                }
            )
        if reference is not None:
            started = time.perf_counter()
            result = evaluate_equivalence(
                internal,
                reference,
                use_chirality=False,
                max_resonance=cap,
            )
            record.update(_result_fields(result, (time.perf_counter() - started) * 1000.0))
        records.append(record)
        # The only strong reference to the internal molecule is released when
        # this loop advances and the chunk-local output dictionary is dropped.
    return records


def _run_phase(
    items: Sequence[tuple[int, dict[str, str]]],
    *,
    cap: int,
    args: argparse.Namespace,
    new_chunk_budget: list[int | None],
) -> bool:
    phase_dir = args.out_dir / f"cap_{cap}"
    phase_dir.mkdir(parents=True, exist_ok=True)
    total_chunks = (len(items) + args.chunk_size - 1) // args.chunk_size
    phase_started = time.perf_counter()
    completed_cases = 0
    for ordinal, chunk in enumerate(_chunks(items, args.chunk_size)):
        path = _chunk_path(args.out_dir, cap, ordinal, chunk)
        if args.resume and path.exists():
            existing = list(_read_jsonl(path))
            if len(existing) != len(chunk):
                raise RuntimeError(f"Incomplete checkpoint: {path}")
            completed_cases += len(existing)
            continue
        if new_chunk_budget[0] is not None and new_chunk_budget[0] <= 0:
            return False
        started = time.perf_counter()
        records = _classify_chunk(
            chunk,
            cap=cap,
            xyz_dir=args.xyz_dir,
            max_workers=args.max_workers,
        )
        chunk_wall_ms = (time.perf_counter() - started) * 1000.0
        chunk_wall_ms_share = chunk_wall_ms / len(records)
        for record in records:
            record["chunk_wall_ms_share"] = round(chunk_wall_ms_share, 6)
        _atomic_jsonl(path, records)
        completed_cases += len(records)
        elapsed = time.perf_counter() - phase_started
        rate = completed_cases / elapsed if elapsed else 0.0
        remaining = len(items) - completed_cases
        eta = remaining / rate if rate else 0.0
        decisions = Counter(record.get("decision", "not_evaluated") for record in records)
        print(
            f"[cap={cap}] chunk {ordinal + 1}/{total_chunks} "
            f"cases={completed_cases}/{len(items)} chunk_s={time.perf_counter() - started:.1f} "
            f"rate={rate:.2f}/s eta_s={eta:.0f} decisions={dict(decisions)}",
            file=sys.stderr,
            flush=True,
        )
        if new_chunk_budget[0] is not None:
            new_chunk_budget[0] -= 1
    return True


def _latest_records(out_dir: Path) -> dict[str, dict[str, Any]]:
    latest: dict[str, dict[str, Any]] = {}
    for cap in PHASE_CAPS:
        for record in _all_phase_records(out_dir, cap):
            latest[str(record["case_id"])] = record
    return latest


def _effective_records(out_dir: Path) -> list[dict[str, Any]]:
    phase = {
        cap: {str(record["case_id"]): record for record in _all_phase_records(out_dir, cap)}
        for cap in PHASE_CAPS
    }
    effective: list[dict[str, Any]] = []
    for case_id, initial in phase[100].items():
        if case_id in phase[1000]:
            record = dict(phase[1000][case_id])
        elif case_id in phase[500]:
            record = dict(phase[500][case_id])
        else:
            record = dict(initial)
        record["initial_decision_at_cap_100"] = initial.get("decision")
        record["retried"] = case_id in phase[500]
        record["resolved_at_cap"] = (
            int(record["cap"])
            if record.get("decision") and not _is_bounded_only(record)
            else None
        )
        effective.append(record)
    return sorted(effective, key=lambda item: int(item["row_index"]))


def _migration_report(
    review_db: Path,
    latest: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    if not review_db.exists():
        return {"status": "review_db_missing", "path": str(review_db)}
    connection = sqlite3.connect(f"file:{review_db.resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        rows = connection.execute(
            """
            SELECT c.case_id, c.category, c.metadata_json,
                   r.status AS review_status, r.notes, r.reviewer
            FROM cases c
            LEFT JOIN reviews r ON r.case_id = c.case_id
            WHERE c.category = 'graph_not_equivalent'
            """
        ).fetchall()
    finally:
        connection.close()
    transitions: Counter[str] = Counter()
    reviewed: Counter[tuple[str, str, str]] = Counter()
    unavailable = 0
    for row in rows:
        record = latest.get(str(row["case_id"]))
        if record is None or not record.get("decision"):
            unavailable += 1
            continue
        decision = str(record["decision"])
        transitions[decision] += 1
        if row["review_status"] is not None:
            # Historical automation commonly stores its reason/rule label in
            # reviewer, while human prose remains in notes.  Preserve both
            # dimensions without rewriting or interpreting either field.
            reviewed[(str(row["review_status"]), str(row["reviewer"] or ""), decision)] += 1
    return {
        "status": "ok",
        "review_db": str(review_db.resolve()),
        "old_graph_not_equivalent_total": len(rows),
        "classified": sum(transitions.values()),
        "unavailable": unavailable,
        "transitions": {
            "old_graph_not_equivalent_to_equivalent": transitions["equivalent"],
            "old_graph_not_equivalent_to_inconclusive": transitions["inconclusive"],
            "old_graph_not_equivalent_to_still_not_equivalent": transitions["not_equivalent"],
        },
        "reviewed_changes_by_existing_verdict_reason": [
            {
                "existing_verdict": verdict,
                "existing_reason_or_reviewer": reviewer,
                "new_decision": decision,
                "count": count,
            }
            for (verdict, reviewer, decision), count in sorted(reviewed.items())
        ],
    }


def _aggregate(args: argparse.Namespace, total: int) -> dict[str, Any]:
    phase_records = {cap: list(_all_phase_records(args.out_dir, cap)) for cap in PHASE_CAPS}
    latest = _latest_records(args.out_dir)
    reconstruction = Counter(record.get("reconstruction_status", "missing") for record in phase_records[100])
    references = Counter(record.get("reference_status", "missing") for record in phase_records[100])
    roundtrip = Counter(record.get("roundtrip_status", "missing") for record in phase_records[100])
    decisions = Counter(record.get("decision", "not_evaluated") for record in latest.values())
    relations = Counter(record.get("relation") or "none" for record in latest.values() if record.get("decision"))
    methods = Counter(record.get("method") or "none" for record in latest.values() if record.get("decision"))
    equivalence_ms = [
        float(record["equivalence_ms"])
        for record in latest.values()
        if record.get("equivalence_ms") is not None
    ]
    phase_wall_seconds = {
        str(cap): sum(float(record.get("chunk_wall_ms_share", 0.0)) for record in records)
        / 1000.0
        for cap, records in phase_records.items()
    }
    bounded_by_cap = {
        str(cap): sum(_is_bounded_only(record) for record in phase_records[cap])
        for cap in PHASE_CAPS
    }
    resolved_at_cap = Counter()
    cap100 = {str(record["case_id"]): record for record in phase_records[100]}
    cap500 = {str(record["case_id"]): record for record in phase_records[500]}
    cap1000 = {str(record["case_id"]): record for record in phase_records[1000]}
    for case_id in cap500:
        if not _is_bounded_only(cap500[case_id]):
            resolved_at_cap[500] += 1
        elif case_id in cap1000 and not _is_bounded_only(cap1000[case_id]):
            resolved_at_cap[1000] += 1
        else:
            resolved_at_cap["unresolved"] += 1
    summary = {
        "evaluator_policy": "v1",
        "total": total,
        "checkpointed_at_cap_100": len(cap100),
        "reconstruction": dict(reconstruction),
        "reference": dict(references),
        "decisions": dict(decisions),
        "relation_counts": dict(relations),
        "method_counts": dict(methods),
        "bounded_resonance_inconclusive_by_cap": bounded_by_cap,
        "attempted_by_cap": {str(cap): len(records) for cap, records in phase_records.items()},
        "resolved_at_cap": {str(key): value for key, value in resolved_at_cap.items()},
        "serialization_roundtrip": dict(roundtrip),
        "equivalence_runtime_ms": {
            "count": len(equivalence_ms),
            "median": statistics.median(equivalence_ms) if equivalence_ms else None,
            "p95": _percentile(equivalence_ms, 0.95),
            "mean": statistics.mean(equivalence_ms) if equivalence_ms else None,
            "total": sum(equivalence_ms),
        },
        "checkpoint_runtime_seconds_by_cap": phase_wall_seconds,
        "checkpoint_runtime_seconds_total": sum(phase_wall_seconds.values()),
        "first_pass_throughput_cases_per_second": (
            len(phase_records[100]) / phase_wall_seconds["100"]
            if phase_wall_seconds["100"]
            else None
        ),
        "overall_checkpoint_throughput_cases_per_second": (
            sum(len(records) for records in phase_records.values()) / sum(phase_wall_seconds.values())
            if sum(phase_wall_seconds.values())
            else None
        ),
        "migration": _migration_report(args.review_db, latest),
    }
    return summary


def main() -> int:
    args = _parse_args()
    RDLogger.DisableLog("rdApp.*")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = _load_rows(args.csv)
    fingerprint = _fingerprint(args, len(rows))
    metadata_path = args.out_dir / "run_metadata.json"
    if metadata_path.exists():
        existing = json.loads(metadata_path.read_text(encoding="utf-8"))
        if existing != fingerprint:
            raise RuntimeError(
                "Checkpoint fingerprint differs from this run. Use a new --out-dir; "
                "evaluator policy and dataset are frozen within one refresh."
            )
    else:
        _atomic_json(metadata_path, fingerprint)
    budget: list[int | None] = [args.stop_after_chunks]
    run_started = time.perf_counter()
    if not _run_phase(rows, cap=100, args=args, new_chunk_budget=budget):
        print("Stopped after requested chunk budget; resume with the same command.", file=sys.stderr)
        return 75
    cap100 = list(_all_phase_records(args.out_dir, 100))
    retry_500_ids = {str(record["case_id"]) for record in cap100 if _is_bounded_only(record)}
    retry_500 = [item for item in rows if item[1].get("id", "").strip() in retry_500_ids]
    if retry_500 and not _run_phase(retry_500, cap=500, args=args, new_chunk_budget=budget):
        print("Stopped after requested chunk budget; resume with the same command.", file=sys.stderr)
        return 75
    cap500 = list(_all_phase_records(args.out_dir, 500))
    retry_1000_ids = {str(record["case_id"]) for record in cap500 if _is_bounded_only(record)}
    if retry_1000_ids and len(retry_1000_ids) <= args.cap_1000_threshold:
        retry_1000 = [item for item in rows if item[1].get("id", "").strip() in retry_1000_ids]
        if not _run_phase(retry_1000, cap=1000, args=args, new_chunk_budget=budget):
            print("Stopped after requested chunk budget; resume with the same command.", file=sys.stderr)
            return 75
    summary = _aggregate(args, len(rows))
    summary["wall_seconds_this_invocation"] = round(time.perf_counter() - run_started, 6)
    _atomic_jsonl(args.out_dir / "effective_results.jsonl", _effective_records(args.out_dir))
    _atomic_json(args.out_dir / "summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
