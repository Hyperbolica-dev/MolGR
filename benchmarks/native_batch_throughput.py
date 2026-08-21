#!/usr/bin/env python
"""Measure native reconstruction throughput in a fresh Python process.

Each invocation measures one execution mode so the caller can compare worker
counts without sharing Open Babel/RDKit state between configurations.
"""

from __future__ import annotations

import argparse
import csv
import gc
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable, Iterable, cast

from molgr import _core as core
from molgr.batch import ReconstructionBatchRequest, iter_xyz_to_rdmol_batch
from molgr.diagnostics import ReconstructionError
from molgr.interface import xyz_to_rdmol


def _load_requests(fixture_dir: Path, item_count: int) -> list[ReconstructionBatchRequest]:
    paths = sorted(fixture_dir.glob("*.xyz"))
    if not paths:
        raise SystemExit(f"no XYZ fixtures found in {fixture_dir}")
    fixtures = [ReconstructionBatchRequest(path.read_text()) for path in paths]
    return [fixtures[index % len(fixtures)] for index in range(item_count)]


def _load_tmqmg_requests(
    xyz_dir: Path,
    properties_csv: Path,
    item_count: int,
) -> list[ReconstructionBatchRequest]:
    with properties_csv.open(newline="") as stream:
        properties = {row["id"]: row for row in csv.DictReader(stream)}
    paths = sorted(xyz_dir.glob("*.xyz"))[:item_count]
    if len(paths) < item_count:
        raise SystemExit(f"only found {len(paths)} XYZ files in {xyz_dir}")
    requests: list[ReconstructionBatchRequest] = []
    for path in paths:
        row = properties.get(path.stem)
        if row is None:
            raise SystemExit(f"missing CSV row for {path.stem}")
        total_charge = int(row["charge"])
        electron_count = int(row["n_electrons"])
        spin_multiplicity = 1 if electron_count % 2 == 0 else 2
        requests.append(
            ReconstructionBatchRequest(
                path.read_text(),
                total_charge=total_charge,
                spin_multiplicity=spin_multiplicity,
            )
        )
    return requests


def _consume_cpp_native(
    requests: Iterable[ReconstructionBatchRequest], workers: int
) -> tuple[int, int]:
    native_requests = (
        (request.xyz_block, request.total_charge, request.spin_multiplicity - 1)
        for request in requests
    )
    iterator = core.pipeline.reconstruct_with_metals.batch_xyz2omol(
        native_requests,
        max_workers=workers,
        queue_size=16,
        ordered=False,
    )
    success = 0
    atom_count = 0
    try:
        for result in iterator:
            molecule_data = cast(Any, result["molecule_data"])
            if molecule_data is not None:
                success += 1
                atom_count += len(molecule_data.atoms)
    finally:
        iterator.close()
    return success, atom_count


def _consume_cpp_batch(
    requests: Iterable[ReconstructionBatchRequest], workers: int
) -> tuple[int, int]:
    success = 0
    atom_count = 0
    max_workers = None if workers == 0 else workers
    for _, molecule, _ in iter_xyz_to_rdmol_batch(requests, backend="cpp", max_workers=max_workers):
        if molecule is not None:
            success += 1
            atom_count += molecule.GetNumAtoms()
    return success, atom_count


def _consume_python_batch(requests: Iterable[ReconstructionBatchRequest]) -> tuple[int, int]:
    success = 0
    atom_count = 0
    for _, molecule, _ in iter_xyz_to_rdmol_batch(requests, backend="python"):
        if molecule is not None:
            success += 1
            atom_count += molecule.GetNumAtoms()
    return success, atom_count


def _consume_single_cpp(requests: Iterable[ReconstructionBatchRequest]) -> tuple[int, int]:
    success = 0
    atom_count = 0
    for request in requests:
        try:
            molecule = xyz_to_rdmol(
                request.xyz_block,
                request.total_charge,
                request.spin_multiplicity,
                backend="cpp",
            )
        except ReconstructionError:
            continue
        if molecule is not None:
            success += 1
            atom_count += molecule.GetNumAtoms()
    return success, atom_count


def _git_revision() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=("cpp_native_batch", "cpp_batch", "cpp_single", "python_batch"),
        required=True,
    )
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--items", type=int, default=60)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--fixture-dir", type=Path, default=Path("tests/data/tmqmg/reconstruction"))
    parser.add_argument("--tmqmg-xyz-dir", type=Path)
    parser.add_argument("--tmqmg-csv", type=Path)
    args = parser.parse_args()
    if args.items < 1 or args.repeats < 1 or args.warmup < 0:
        parser.error("items and repeats must be positive; warmup must be non-negative")
    if args.mode in ("cpp_native_batch", "cpp_batch") and args.workers < 0:
        parser.error("workers must be non-negative for C++ batch modes (0 means auto)")
    if (args.tmqmg_xyz_dir is None) != (args.tmqmg_csv is None):
        parser.error("--tmqmg-xyz-dir and --tmqmg-csv must be supplied together")

    if args.tmqmg_xyz_dir is not None and args.tmqmg_csv is not None:
        requests = _load_tmqmg_requests(args.tmqmg_xyz_dir, args.tmqmg_csv, args.items)
        fixture_source = str(args.tmqmg_xyz_dir)
    else:
        requests = _load_requests(args.fixture_dir, args.items)
        fixture_source = str(args.fixture_dir)
    mode_fn: Callable[..., tuple[int, int]]
    if args.mode == "cpp_native_batch":
        mode_fn = _consume_cpp_native
    elif args.mode == "cpp_batch":
        mode_fn = _consume_cpp_batch
    elif args.mode == "cpp_single":
        mode_fn = _consume_single_cpp
    else:
        mode_fn = _consume_python_batch

    for _ in range(args.warmup):
        warmup_requests = requests[: min(6, len(requests))]
        if args.mode in ("cpp_native_batch", "cpp_batch"):
            mode_fn(warmup_requests, args.workers)
        else:
            mode_fn(warmup_requests)
    gc.collect()

    samples: list[float] = []
    successes: list[int] = []
    atom_counts: list[int] = []
    for _ in range(args.repeats):
        started = time.perf_counter()
        if args.mode in ("cpp_native_batch", "cpp_batch"):
            success, atom_count = mode_fn(requests, args.workers)
        else:
            success, atom_count = mode_fn(requests)
        elapsed = time.perf_counter() - started
        samples.append(elapsed)
        successes.append(success)
        atom_counts.append(atom_count)

    median_seconds = statistics.median(samples)
    result: dict[str, Any] = {
        "mode": args.mode,
        "workers": args.workers if args.mode in ("cpp_native_batch", "cpp_batch") else 1,
        "items": args.items,
        "repeats": args.repeats,
        "warmup": args.warmup,
        "samples_seconds": samples,
        "median_seconds": median_seconds,
        "mean_seconds": statistics.mean(samples),
        "min_seconds": min(samples),
        "throughput_items_per_second": args.items / median_seconds,
        "success_counts": successes,
        "atom_count_checksums": atom_counts,
        "all_success": all(value == args.items for value in successes),
        "consistent_atom_count": len(set(atom_counts)) == 1,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "git_revision": _git_revision(),
        "fixture_dir": fixture_source,
    }
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
