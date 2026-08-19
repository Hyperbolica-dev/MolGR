from __future__ import annotations

import json
import multiprocessing as mp
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pytest


pytest.importorskip("openbabel")
pytest.importorskip("rdkit")


_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from benchmarks.tmqmg_concurrency_stress import available_cpu_count  # noqa: E402


_WATER_XYZ = """3
water
O 0.0000 0.0000 0.0000
H 0.7586 0.0000 0.5043
H -0.7586 0.0000 0.5043
"""

_TMQMG_ABACAL_XYZ = (
    Path(__file__).parent / "data" / "tmqmg" / "reconstruction" / "ABACAL.xyz"
).read_text(encoding="utf-8")


def _spawn_stress_process_count() -> int:
    return min(4, max(1, available_cpu_count()))


def _spawn_stress_native_worker_count(process_count: int) -> int:
    if sys.platform == "win32":
        # Windows defaults to serial per-molecule native execution. Keep the
        # process-level stress test within the documented safety boundary.
        return 1
    return max(1, min(2, available_cpu_count() // process_count))


def _spawn_native_batch_worker(rounds_and_size: tuple[int, int, int]) -> dict[str, int]:
    from molgr import ReconstructionBatchRequest, iter_xyz_to_rdmol_batch

    rounds, batch_size, native_worker_count = rounds_and_size
    completed = 0
    for _ in range(rounds):
        requests = [ReconstructionBatchRequest(_WATER_XYZ) for _ in range(batch_size)]
        results = list(
            iter_xyz_to_rdmol_batch(
                requests,
                backend="cpp",
                max_workers=native_worker_count,
                queue_size=2,
                ordered=False,
            )
        )
        assert len(results) == batch_size
        assert all(result.molecule is not None for result in results)
        completed += len(results)
    return {"pid": os.getpid(), "completed": completed}


def _run_spawn_stress() -> None:
    process_count = _spawn_stress_process_count()
    native_worker_count = _spawn_stress_native_worker_count(process_count)
    rounds = 5
    batch_size = 24
    context = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=process_count, mp_context=context) as executor:
        summaries = list(
            executor.map(
                _spawn_native_batch_worker,
                [(rounds, batch_size, native_worker_count)] * process_count,
            )
        )
    assert len({summary["pid"] for summary in summaries}) == process_count
    assert sum(summary["completed"] for summary in summaries) == (
        process_count * rounds * batch_size
    )
    print(
        json.dumps(
            {
                "start_method": "spawn",
                "process_count": process_count,
                "native_batch_worker_count": native_worker_count,
                "available_cpu_count": available_cpu_count(),
                "workers": summaries,
            },
            sort_keys=True,
        )
    )


def _spawn_public_xyz_worker(rounds_and_size: tuple[int, int]) -> dict[str, object]:
    from molgr import MolGRConfig
    from molgr.interface import xyz_to_rdmol

    rounds, batch_size = rounds_and_size
    config = MolGRConfig()
    completed = 0
    peak_thread_count = 0
    for _ in range(rounds):
        for _ in range(batch_size):
            molecule = xyz_to_rdmol(_TMQMG_ABACAL_XYZ, backend="cpp", config=config)
            assert molecule.GetNumAtoms() == 27
            completed += 1
            if sys.platform.startswith("linux"):
                peak_thread_count = max(peak_thread_count, len(os.listdir("/proc/self/task")))
    return {
        "pid": os.getpid(),
        "completed": completed,
        "peak_thread_count": peak_thread_count or None,
        "config": {
            "max_threads": config.cpp_backend.max_threads,
            "enable_target_bucket_parallelism": config.cpp_backend.enable_target_bucket_parallelism,
            "target_bucket_parallel_max_threads": config.cpp_backend.target_bucket_parallel_max_threads,
        },
    }


def _run_public_xyz_stress() -> None:
    process_count = _spawn_stress_process_count()
    rounds = 4
    batch_size = 20
    context = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=process_count, mp_context=context) as executor:
        summaries = list(
            executor.map(
                _spawn_public_xyz_worker,
                [(rounds, batch_size)] * process_count,
            )
        )
    assert len({summary["pid"] for summary in summaries}) == process_count
    assert sum(summary["completed"] for summary in summaries) == (
        process_count * rounds * batch_size
    )
    print(
        json.dumps(
            {
                "api": "xyz_to_rdmol",
                "process_count": process_count,
                "available_cpu_count": available_cpu_count(),
                "workers": summaries,
            },
            sort_keys=True,
        )
    )


def test_spawn_processes_can_run_native_batches_under_stress() -> None:
    completed = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "--spawn-stress"],
        check=False,
        capture_output=True,
        text=True,
        timeout=180,
    )

    assert completed.returncode == 0, completed.stderr
    summary = json.loads(completed.stdout.strip().splitlines()[-1])
    assert summary["start_method"] == "spawn"
    process_count = _spawn_stress_process_count()
    assert summary["available_cpu_count"] == available_cpu_count()
    assert summary["process_count"] == process_count
    assert summary["native_batch_worker_count"] == _spawn_stress_native_worker_count(process_count)
    assert len(summary["workers"]) == process_count
    assert sum(worker["completed"] for worker in summary["workers"]) == process_count * 120


def test_spawn_processes_can_call_public_xyz_to_rdmol_under_stress() -> None:
    completed = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "--public-xyz-stress"],
        check=False,
        capture_output=True,
        text=True,
        timeout=180,
    )

    assert completed.returncode == 0, completed.stderr
    summary = json.loads(completed.stdout.strip().splitlines()[-1])
    assert summary["api"] == "xyz_to_rdmol"
    process_count = _spawn_stress_process_count()
    assert summary["available_cpu_count"] == available_cpu_count()
    assert summary["process_count"] == process_count
    assert len(summary["workers"]) == process_count
    assert sum(worker["completed"] for worker in summary["workers"]) == process_count * 80
    assert all(
        worker["config"]["enable_target_bucket_parallelism"] for worker in summary["workers"]
    )
    if sys.platform.startswith("linux"):
        assert all(worker["config"]["max_threads"] is None for worker in summary["workers"])
        assert all(worker["peak_thread_count"] > 1 for worker in summary["workers"])


@pytest.mark.skipif(not hasattr(os, "fork"), reason="POSIX fork is unavailable")
def test_fork_after_import_is_rejected_repeatedly_after_native_pool_warmup() -> None:
    script = f"""\
import gc
import os
import sys
from molgr import _core

xyz = {_WATER_XYZ!r}
requests = [(xyz, 0, 0)] * 32
warmup = _core.pipeline.reconstruct_with_metals.batch_xyz2omol(
    requests, max_workers=4, queue_size=2, ordered=False
)
assert len(list(warmup)) == len(requests)
warmup.close()

for index in range(32):
    inherited = _core.pipeline.reconstruct_with_metals.batch_xyz2omol(
        requests, max_workers=4, queue_size=1, ordered=False
    )
    assert next(inherited)["molecule_data"] is not None
    pid = os.fork()
    if pid == 0:
        try:
            _core.pipeline.reconstruct_with_metals.xyz2omol(xyz, 0, 0)
        except RuntimeError as exc:
            assert "forked child after MolGR/Open Babel was initialized" in str(exc)
            del inherited
            gc.collect()
            sys.exit(0)
        sys.exit(3)
    _, status = os.waitpid(pid, 0)
    inherited.close()
    assert os.WIFEXITED(status), status
    assert os.WEXITSTATUS(status) == 0, status

print("fork-rejected", 32)
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().endswith("fork-rejected 32")


@pytest.mark.skipif(not hasattr(os, "fork"), reason="POSIX fork is unavailable")
def test_fork_before_import_remains_supported() -> None:
    script = f"""\
import os
import sys

pid = os.fork()
if pid == 0:
    from molgr.interface import xyz_to_rdmol
    mol = xyz_to_rdmol({_WATER_XYZ!r})
    sys.exit(0 if mol.GetNumAtoms() == 3 else 4)

_, status = os.waitpid(pid, 0)
assert os.WIFEXITED(status), status
assert os.WEXITSTATUS(status) == 0, status
print("fork-before-import-ok")
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().endswith("fork-before-import-ok")


if __name__ == "__main__":
    if sys.argv[1:] == ["--spawn-stress"]:
        _run_spawn_stress()
    elif sys.argv[1:] == ["--public-xyz-stress"]:
        _run_public_xyz_stress()
