#!/usr/bin/env python
"""Stress MolGR native-batch and spawn-process concurrency on XYZ data.

The native-batch mode keeps one process and gives its worker budget to the C++
batch scheduler. The spawn mode deliberately avoids importing MolGR in the
parent; every worker imports the native extension after spawn and directly
calls ``xyz_to_rdmol`` with the default per-molecule parallelism. Run one mode
and worker count per invocation to keep measurements in a fresh process tree.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gc
import json
import math
import multiprocessing as mp
import os
import platform
import time
from pathlib import Path
from queue import Empty
from typing import Any, Tuple


Item = Tuple[str, str, int, int]


def _read_cpu_quota() -> int | None:
    """Return the integer CPU budget imposed by a Linux cgroup, if present."""
    quota_paths = (
        Path("/sys/fs/cgroup/cpu.max"),
        Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us"),
    )
    period_paths = (
        Path("/sys/fs/cgroup/cpu.max"),
        Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us"),
    )
    try:
        v2_quota = quota_paths[0].read_text(encoding="ascii").split()
        if len(v2_quota) == 2:
            if v2_quota[0] == "max":
                return None
            quota, period = int(v2_quota[0]), int(v2_quota[1])
            if quota > 0 and period > 0:
                return max(1, math.floor(quota / period))
    except (OSError, ValueError, IndexError):
        pass
    try:
        quota = int(quota_paths[1].read_text(encoding="ascii").strip())
        period = int(period_paths[1].read_text(encoding="ascii").strip())
        if quota > 0 and period > 0:
            return max(1, math.floor(quota / period))
    except (OSError, ValueError):
        return None
    return None


def available_cpu_count() -> int:
    """Estimate CPUs available to this process, including affinity/cgroup limits."""
    counts: list[int] = []
    process_cpu_count = getattr(os, "process_cpu_count", None)
    if process_cpu_count is not None:
        count = process_cpu_count()
        if count is not None and count > 0:
            counts.append(int(count))
    if hasattr(os, "sched_getaffinity"):
        with contextlib.suppress(OSError):
            counts.append(len(os.sched_getaffinity(0)))
    system_count = os.cpu_count()
    if system_count is not None and system_count > 0:
        counts.append(int(system_count))
    quota_count = _read_cpu_quota()
    if quota_count is not None:
        counts.append(quota_count)
    return max(1, min(counts)) if counts else 1


def stress_worker_counts(*, max_workers: int = 8) -> tuple[int, ...]:
    """Return bounded powers of two plus the runner's exact available CPU budget."""
    if max_workers < 1:
        raise ValueError("max_workers must be positive")
    available = min(max_workers, available_cpu_count())
    counts = [1]
    candidate = 2
    while candidate <= available:
        counts.append(candidate)
        candidate *= 2
    if available not in counts:
        counts.append(available)
    return tuple(counts)


def _load_items(
    csv_path: Path,
    xyz_dir: Path,
    limit: int,
    *,
    repeat_inputs: bool,
) -> list[Item]:
    with csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    available: list[Item] = []
    for row in rows:
        case_id = (row.get("id") or row.get("case_id") or "").strip()
        if not case_id:
            continue
        xyz_file = (row.get("xyz_file") or f"{case_id}.xyz").strip()
        path = xyz_dir / xyz_file
        if not path.is_file():
            continue
        if row.get("spin_multiplicity"):
            spin_multiplicity = int(row["spin_multiplicity"])
        elif row.get("n_electrons"):
            spin_multiplicity = 1 if int(row["n_electrons"]) % 2 == 0 else 2
        else:
            spin_multiplicity = 1
        available.append((case_id, str(path), int(row.get("charge") or 0), spin_multiplicity))
        if len(available) >= limit:
            break
    if len(available) >= limit:
        return available[:limit]
    if repeat_inputs and available:
        return [available[index % len(available)] for index in range(limit)]
    raise SystemExit(
        f"requested {limit} items but found only {len(available)} CSV/XYZ pairs; "
        "pass --repeat-inputs to cycle the available fixtures"
    )


def _rss_peak_kib() -> int | None:
    try:
        import resource

        return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    except (ImportError, OSError):
        return None


def _thread_count() -> int | None:
    try:
        return len(os.listdir("/proc/self/task"))
    except OSError:
        return None


def _worker(worker_id: int, items: list[Item], result_queue: Any) -> None:
    # Native imports must happen in the spawned child, after process creation.
    from molgr import MolGRConfig
    from molgr.interface import xyz_to_rdmol

    config = MolGRConfig()
    started = time.perf_counter()
    success = 0
    failures = 0
    errors: list[dict[str, str]] = []
    atom_count = 0
    peak_thread_count = _thread_count()
    for case_id, xyz_path, total_charge, spin_multiplicity in items:
        try:
            molecule = xyz_to_rdmol(
                Path(xyz_path).read_text(encoding="utf-8"),
                total_charge=total_charge,
                spin_multiplicity=spin_multiplicity,
                backend="cpp",
                config=config,
            )
            success += 1
            atom_count += molecule.GetNumAtoms()
            current_thread_count = _thread_count()
            if current_thread_count is not None:
                peak_thread_count = max(peak_thread_count or 0, current_thread_count)
        except Exception as exc:  # noqa: BLE001
            failures += 1
            if len(errors) < 10:
                error_type = type(exc).__name__
                errors.append({"id": case_id, "error": f"{error_type}: {exc}"})
    result_queue.put(
        {
            "worker_id": worker_id,
            "pid": os.getpid(),
            "assigned": len(items),
            "success": success,
            "failures": failures,
            "atom_count": atom_count,
            "elapsed_seconds": time.perf_counter() - started,
            "peak_rss_kib": _rss_peak_kib(),
            "peak_native_thread_count": peak_thread_count,
            "errors": errors,
            "config": {
                "max_threads": config.cpp_backend.max_threads,
                "enable_target_bucket_parallelism": config.cpp_backend.enable_target_bucket_parallelism,
                "target_bucket_parallel_max_threads": config.cpp_backend.target_bucket_parallel_max_threads,
                "enable_candidate_scoring_parallelism": config.cpp_backend.enable_candidate_scoring_parallelism,
            },
        }
    )


def _partition(items: list[Item], worker_count: int) -> list[list[Item]]:
    return [items[index::worker_count] for index in range(worker_count)]


def _run_native_batch(
    items: list[Item], worker_count: int, timeout_seconds: float
) -> dict[str, Any]:
    del timeout_seconds
    from molgr import MolGRConfig, ReconstructionBatchRequest, iter_xyz_to_rdmol_batch

    requests = [
        ReconstructionBatchRequest(
            xyz_block=Path(xyz_path).read_text(encoding="utf-8"),
            total_charge=total_charge,
            spin_multiplicity=spin_multiplicity,
        )
        for _, xyz_path, total_charge, spin_multiplicity in items
    ]
    config = MolGRConfig()
    effective_max_threads = config.cpp_backend.max_threads
    effective_target_parallelism = config.cpp_backend.enable_target_bucket_parallelism
    effective_target_max_threads = config.cpp_backend.target_bucket_parallel_max_threads
    effective_candidate_parallelism = config.cpp_backend.enable_candidate_scoring_parallelism
    if worker_count > 1:
        effective_max_threads = 1
        effective_target_parallelism = False
        effective_target_max_threads = 1
        effective_candidate_parallelism = False
    started = time.perf_counter()
    success = 0
    failures = 0
    atom_count = 0
    errors: list[dict[str, str]] = []
    peak_thread_count = _thread_count()
    try:
        results = iter_xyz_to_rdmol_batch(
            requests,
            backend="cpp",
            max_workers=worker_count,
            queue_size=16,
            ordered=False,
            config=config,
        )
        for result in results:
            if result.molecule is None:
                failures += 1
                if len(errors) < 10 and result.status is not None:
                    errors.append(
                        {
                            "error": f"{result.status.code}: {result.status.message}",
                        }
                    )
            else:
                success += 1
                atom_count += result.molecule.GetNumAtoms()
            current_thread_count = _thread_count()
            if current_thread_count is not None:
                peak_thread_count = max(peak_thread_count or 0, current_thread_count)
    finally:
        gc.collect()
    elapsed = time.perf_counter() - started
    return {
        "mode": "native_batch",
        "start_method": "native_threads",
        "process_count": 1,
        "native_batch_worker_count": worker_count,
        "worker_count": worker_count,
        "requested": len(items),
        "reported": success + failures,
        "success": success,
        "failures": failures,
        "unreported": len(items) - success - failures,
        "timed_out": False,
        "terminated_pids": [],
        "crashed_pids": [],
        "nonzero_exit_pids": [],
        "wall_seconds": elapsed,
        "items_per_second": success / elapsed if elapsed else 0.0,
        "peak_rss_kib_by_worker": {"native": _rss_peak_kib()},
        "peak_native_thread_count": peak_thread_count,
        "workers": [
            {
                "worker_id": 0,
                "pid": os.getpid(),
                "assigned": len(items),
                "success": success,
                "failures": failures,
                "atom_count": atom_count,
                "peak_rss_kib": _rss_peak_kib(),
                "peak_native_thread_count": peak_thread_count,
                "errors": errors,
                "config": {
                    "max_threads": effective_max_threads,
                    "enable_target_bucket_parallelism": effective_target_parallelism,
                    "target_bucket_parallel_max_threads": effective_target_max_threads,
                    "enable_candidate_scoring_parallelism": effective_candidate_parallelism,
                },
            }
        ],
    }


def _run(items: list[Item], worker_count: int, timeout_seconds: float) -> dict[str, Any]:
    context = mp.get_context("spawn")
    result_queue = context.Queue()
    chunks = _partition(items, worker_count)
    processes = [
        context.Process(target=_worker, args=(worker_id, chunk, result_queue))
        for worker_id, chunk in enumerate(chunks)
    ]
    started = time.perf_counter()
    for process in processes:
        process.start()
    deadline = started + timeout_seconds
    results: dict[int, dict[str, Any]] = {}
    crashed_pids: list[int] = []
    while len(results) < worker_count and time.perf_counter() < deadline:
        for process in processes:
            if process.is_alive():
                continue
            if process.exitcode not in (None, 0):
                crashed_pids.append(process.pid or -1)
                break
        if crashed_pids:
            break
        try:
            result = result_queue.get(timeout=min(0.5, max(0.01, deadline - time.perf_counter())))
        except Empty:
            continue
        results[int(result["worker_id"])] = result
    timed_out = len(results) < worker_count and not crashed_pids and time.perf_counter() >= deadline
    abort_remaining = bool(crashed_pids)
    for process in processes:
        if abort_remaining or timed_out:
            process.join(timeout=0.1)
        else:
            process.join(timeout=max(0.0, deadline - time.perf_counter()))
    terminated: list[int] = []
    for process in processes:
        if process.is_alive():
            terminated.append(process.pid or -1)
            process.terminate()
    for process in processes:
        process.join(timeout=5)
    elapsed = time.perf_counter() - started
    completed = sum(int(result["assigned"]) for result in results.values())
    success = sum(int(result["success"]) for result in results.values())
    failures = sum(int(result["failures"]) for result in results.values())
    exited_nonzero = [
        process.pid or -1 for process in processes if process.exitcode not in (0, None)
    ]
    result_queue.close()
    return {
        "mode": "spawn_xyz_to_rdmol",
        "start_method": "spawn",
        "process_count": worker_count,
        "native_batch_worker_count": 0,
        "worker_count": worker_count,
        "requested": len(items),
        "reported": completed,
        "success": success,
        "failures": failures,
        "unreported": len(items) - completed,
        "timed_out": timed_out,
        "terminated_pids": terminated,
        "crashed_pids": crashed_pids,
        "nonzero_exit_pids": exited_nonzero,
        "wall_seconds": elapsed,
        "items_per_second": success / elapsed if elapsed else 0.0,
        "peak_rss_kib_by_worker": {
            str(worker_id): result["peak_rss_kib"] for worker_id, result in sorted(results.items())
        },
        "workers": [results[key] for key in sorted(results)],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=("native_batch", "spawn_xyz_to_rdmol"),
        default="spawn_xyz_to_rdmol",
    )
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--xyz-dir", type=Path, required=True)
    parser.add_argument("--items", type=int, default=1000)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    parser.add_argument("--repeat-inputs", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.items < 1 or args.workers < 1 or args.timeout_seconds <= 0:
        parser.error("--items, --workers, and --timeout-seconds must be positive")
    items = _load_items(
        args.csv,
        args.xyz_dir,
        args.items,
        repeat_inputs=args.repeat_inputs,
    )
    if args.mode == "native_batch":
        summary = _run_native_batch(items, args.workers, args.timeout_seconds)
    else:
        summary = _run(items, args.workers, args.timeout_seconds)
    summary.update(
        {
            "csv": str(args.csv),
            "xyz_dir": str(args.xyz_dir),
            "items": args.items,
            "repeat_inputs": args.repeat_inputs,
            "platform": platform.platform(),
            "python": platform.python_version(),
            "cpu_count": os.cpu_count(),
            "available_cpu_count": available_cpu_count(),
            "rss_unit": "KiB on Linux/macOS; bytes on Windows",
        }
    )
    payload = json.dumps(summary, ensure_ascii=True, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")
    print(payload)
    if summary["timed_out"] or summary["nonzero_exit_pids"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
