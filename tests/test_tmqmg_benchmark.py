from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from benchmarks.tmqmg_xyz_benchmark.io import summarize_results
from benchmarks.tmqmg_xyz_benchmark.run import (
    TmqmgBenchmarkInput,
    _build_case,
    _build_cpp_backend_config_payload,
    _run_case_method,
    _run_method_cases_worker,
    _run_method_subprocesses,
    _select_methods,
    _selected_rows,
    get_method_registry,
    run,
)
from benchmarks.tmqmg_xyz_benchmark.schema import BenchmarkResult


def test_tmqmg_subset_selection_respects_ids_and_row_bounds() -> None:
    rows = [
        {"id": "A", "smiles": "C", "charge": "0"},
        {"id": "B", "smiles": "CC", "charge": "0"},
        {"id": "C", "smiles": "CCC", "charge": "0"},
        {"id": "D", "smiles": "CCCC", "charge": "0"},
    ]

    selected = _selected_rows(
        rows,
        ids={"B", "D"},
        start_row=2,
        end_row=4,
        limit=1,
    )

    assert [item.row["id"] for item in selected] == ["B"]
    assert [item.row_index for item in selected] == [2]


def test_tmqmg_run_writes_results_for_selected_subset(tmp_path: Path, monkeypatch) -> None:
    csv_path = tmp_path / "tmqmg.csv"
    xyz_dir = tmp_path / "xyz"
    out_dir = tmp_path / "out"
    xyz_dir.mkdir()

    csv_path.write_text(
        "id,smiles,charge\n"
        "A,C,0\n"
        "B,CC,0\n"
        "C,CCC,0\n",
        encoding="utf-8",
    )
    (xyz_dir / "B.xyz").write_text("1\nB\nB 0 0 0\n", encoding="utf-8")

    def _fake_methods(method_ids=None):
        class _Method:
            method_id = "fake"

            def run(self, case):
                return _Output(
                    status="ok",
                    error=None,
                    predicted_smiles=case["input_smiles"],
                    rdkit_mol=None,
                    equivalent=True,
                    equivalence_method="ideal",
                    timing_ms_breakdown={"method_ms": 1.0},
                )

        return [_Method()]

    class _Output:
        def __init__(
            self,
            *,
            status: str,
            error: str | None,
            predicted_smiles: str | None,
            rdkit_mol,
            equivalent: bool | None,
            equivalence_method: str | None,
            timing_ms_breakdown: dict[str, float],
        ) -> None:
            self.status = status
            self.error = error
            self.predicted_smiles = predicted_smiles
            self.rdkit_mol = rdkit_mol
            self.equivalent = equivalent
            self.equivalence_method = equivalence_method
            self.timing_ms_breakdown = timing_ms_breakdown

    monkeypatch.setattr("benchmarks.tmqmg_xyz_benchmark.run.get_method_registry", _fake_methods)
    monkeypatch.setattr(
        "benchmarks.tmqmg_xyz_benchmark.run._run_method_subprocess",
        lambda **kwargs: [
            type(
                "_Result",
                (),
                {
                    "case_idx": 2,
                    "case_id": "B",
                    "method_id": "fake",
                    "input_smiles": "CC",
                    "ground_truth_smiles": "CC",
                    "status": "ok",
                    "error": None,
                    "predicted_smiles": "CC",
                    "equivalent": True,
                    "equivalence_method": "ideal",
                    "comparison_skipped": False,
                    "comparison_skip_reason": None,
                    "timing_ms_total": 1.0,
                    "timing_ms_breakdown": {"method_ms": 1.0},
                },
            )
        ],
    )

    results = run(
        csv_path=csv_path,
        xyz_dir=xyz_dir,
        out_dir=out_dir,
        limit=1,
        start_row=2,
        end_row=3,
        ids=["B"],
    )

    assert len(results) == 1
    assert (out_dir / "results.csv").exists()
    assert (out_dir / "summary.csv").exists()

    with (out_dir / "results.csv").open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert [row["id"] for row in rows] == ["B"]
    assert [row["input_smiles"] for row in rows] == ["CC"]


def test_tmqmg_run_filters_selected_methods(tmp_path: Path, monkeypatch) -> None:
    csv_path = tmp_path / "tmqmg.csv"
    xyz_dir = tmp_path / "xyz"
    out_dir = tmp_path / "out"
    xyz_dir.mkdir()

    csv_path.write_text("id,smiles,charge\nA,C,0\n", encoding="utf-8")
    (xyz_dir / "A.xyz").write_text("1\nA\nC 0 0 0\n", encoding="utf-8")

    class _Method:
        def __init__(self, method_id: str) -> None:
            self.method_id = method_id

    monkeypatch.setattr(
        "benchmarks.tmqmg_xyz_benchmark.run.get_method_registry",
        lambda method_ids=None: [
            _Method(method_id)
            for method_id in (method_ids or ("skip_me", "molgr_fallback", "molgr_cpp"))
        ],
    )

    called_method_ids: list[str] = []

    def _fake_run_method_subprocess(**kwargs):
        method_id = kwargs["method_id"]
        called_method_ids.append(method_id)
        return [
            BenchmarkResult(
                case_idx=1,
                method_id=method_id,
                input_smiles="C",
                ground_truth_smiles="C",
                status="ok",
                error=None,
                predicted_smiles="C",
                equivalent=True,
                equivalence_method="ideal",
                timing_ms_total=1.0,
                timing_ms_breakdown={"method_ms": 1.0},
                case_id="A",
            )
        ]

    monkeypatch.setattr(
        "benchmarks.tmqmg_xyz_benchmark.run._run_method_subprocess",
        _fake_run_method_subprocess,
    )

    results = run(
        csv_path=csv_path,
        xyz_dir=xyz_dir,
        out_dir=out_dir,
        method_ids=["molgr_cpp", "molgr_fallback"],
    )

    assert called_method_ids == ["molgr_cpp", "molgr_fallback"]
    assert [result.method_id for result in results] == ["molgr_cpp", "molgr_fallback"]


def test_tmqmg_method_filter_rejects_unknown_method_id() -> None:
    class _Method:
        method_id = "known"

    try:
        _select_methods([_Method()], ["missing"])
    except ValueError as exc:
        assert "missing" in str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_tmqmg_method_registry_only_imports_selected_method_modules(monkeypatch) -> None:
    imported_modules: list[str] = []

    def fake_import_module(module_name: str):
        imported_modules.append(module_name)

        class _Method:
            method_id = "molgr_cpp"

        return type("_Module", (), {"MolGRCppMethod": _Method})

    monkeypatch.setattr(
        "benchmarks.smiles_xyz_benchmark.methods.import_module",
        fake_import_module,
    )

    methods = get_method_registry(("molgr_cpp",))

    assert [method.method_id for method in methods] == ["molgr_cpp"]
    assert imported_modules == ["benchmarks.smiles_xyz_benchmark.methods.molgr_cpp"]


def test_cpp_acceleration_preset_enables_available_cpp_backend_switches() -> None:
    cpp_backend_config = _build_cpp_backend_config_payload(use_all_accelerations=True)

    assert cpp_backend_config["enable_target_bucket_parallelism"] is True
    assert cpp_backend_config["enable_candidate_scoring_parallelism"] is False
    assert cpp_backend_config["enable_uff_atom_typing_cache"] is True
    assert cpp_backend_config["enable_target_bucket_score_bundle_preheat"] is True
    assert cpp_backend_config["target_bucket_parallel_threshold"] == 1
    assert cpp_backend_config["target_bucket_parallel_max_threads"] is None
    assert cpp_backend_config["candidate_score_parallel_threshold"] == 32


def test_cpp_uff_atom_typing_cache_can_be_enabled_with_default_preset() -> None:
    cpp_backend_config = _build_cpp_backend_config_payload(
        use_all_accelerations=False,
        enable_uff_atom_typing_cache=True,
    )

    assert cpp_backend_config["enable_target_bucket_parallelism"] is True
    assert cpp_backend_config["enable_candidate_scoring_parallelism"] is False
    assert cpp_backend_config["enable_uff_atom_typing_cache"] is True
    assert cpp_backend_config["target_bucket_parallel_threshold"] == 1
    assert cpp_backend_config["target_bucket_parallel_max_threads"] is None


def test_tmqmg_worker_applies_cpp_backend_config_before_running_methods(
    tmp_path: Path,
    monkeypatch,
) -> None:
    seen_configs = []
    xyz_dir = tmp_path / "xyz"
    xyz_dir.mkdir()
    (xyz_dir / "A.xyz").write_text("2\nA\nC 0 0 0\nC 0 0 1\n", encoding="utf-8")

    class _Output:
        def __init__(self, config):
            self.status = "ok"
            self.error = None
            self.predicted_smiles = "CC"
            self.rdkit_mol = None
            self.equivalent = True
            self.equivalence_method = "ideal"
            self.timing_ms_breakdown = {"method_ms": 1.0}
            self._config = config

    class _Method:
        method_id = "fake"

        def run(self, case):
            from molgr.config import CONFIG

            config = CONFIG
            seen_configs.append(config.cpp_backend)
            return _Output(config)

    monkeypatch.setattr(
        "benchmarks.tmqmg_xyz_benchmark.run.get_method_registry",
        lambda method_ids=None: [_Method()],
    )

    payload = {
        "method_id": "fake",
        "cases": [{"row_index": 1, "row": {"id": "A", "smiles": "CC", "charge": "0"}}],
        "xyz_dir": str(xyz_dir),
        "case_timeout_seconds": 1.0,
        "cpp_backend_config": {
            "max_threads": None,
            "enable_target_bucket_parallelism": True,
            "enable_candidate_scoring_parallelism": True,
            "enable_uff_atom_typing_cache": False,
            "enable_target_bucket_score_bundle_preheat": True,
            "target_bucket_parallel_threshold": 4,
            "target_bucket_parallel_max_threads": 1,
            "candidate_score_parallel_threshold": 1,
        },
    }

    results = _run_method_cases_worker(payload)

    assert len(results) == 1
    assert seen_configs
    assert seen_configs[0].enable_target_bucket_parallelism is True
    assert seen_configs[0].enable_candidate_scoring_parallelism is True
    assert seen_configs[0].enable_uff_atom_typing_cache is False
    assert seen_configs[0].enable_target_bucket_score_bundle_preheat is True
    assert seen_configs[0].target_bucket_parallel_threshold == 4
    assert seen_configs[0].target_bucket_parallel_max_threads == 1


def test_tmqmg_method_subprocesses_split_cases_across_process_workers(
    monkeypatch,
) -> None:
    class _FakeProcess:
        def __init__(self, *, env):
            self.returncode = 0
            self._out_path = Path(env["MOLGR_TMQMG_SUBPROCESS_OUT"])
            payload = json.loads(
                Path(env["MOLGR_TMQMG_SUBPROCESS_PAYLOAD_PATH"]).read_text(
                    encoding="utf-8"
                )
            )
            with self._out_path.open("w", encoding="utf-8") as fh:
                for item in payload["cases"]:
                    row_index = int(item["row_index"])
                    fh.write(
                        json.dumps(
                            BenchmarkResult(
                                case_idx=row_index,
                                method_id=payload["method_id"],
                                input_smiles="C",
                                ground_truth_smiles="C",
                                status="ok",
                                error=None,
                                predicted_smiles="C",
                                equivalent=True,
                                equivalence_method="ideal",
                                timing_ms_total=float(row_index),
                                timing_ms_breakdown={"method_ms": float(row_index)},
                                case_id=item["row"]["id"],
                            ).to_dict(),
                            ensure_ascii=True,
                        )
                        + "\n"
                    )

        def communicate(self):
            return "", ""

        def poll(self):
            return self.returncode

        def kill(self):
            self.returncode = -9

        def wait(self):
            return self.returncode

    popen_calls = []

    def fake_popen(*args, **kwargs):
        popen_calls.append((args, kwargs))
        return _FakeProcess(env=kwargs["env"])

    monkeypatch.setattr(subprocess, "Popen", fake_popen)

    cases = [
        TmqmgBenchmarkInput(row_index=idx, row={"id": f"C{idx}"})
        for idx in range(1, 6)
    ]
    results = _run_method_subprocesses(
        method_id="molgr_cpp",
        cases=cases,
        xyz_dir=Path("xyz"),
        case_timeout_seconds=1.0,
        cpp_backend_config=_build_cpp_backend_config_payload(use_all_accelerations=True),
        process_workers=3,
    )

    assert len(popen_calls) == 3
    assert [result.case_idx for result in results] == [1, 2, 3, 4, 5]
    assert [result.case_id for result in results] == ["C1", "C2", "C3", "C4", "C5"]


def test_build_case_marks_provider_error_for_missing_xyz(tmp_path: Path) -> None:
    case = _build_case(
        2,
        {"id": "B", "smiles": "CC", "charge": "0"},
        xyz_dir=tmp_path / "xyz",
    )

    assert case["provider_error"] is not None
    assert case["xyz_block"] is None


def test_build_case_marks_reference_error_for_reference_xyz_element_mismatch(
    tmp_path: Path,
) -> None:
    xyz_dir = tmp_path / "xyz"
    xyz_dir.mkdir()
    (xyz_dir / "A.xyz").write_text(
        "2\nA\nC 0 0 0\nH 0 0 1\n",
        encoding="utf-8",
    )

    case = _build_case(
        1,
        {"id": "A", "smiles": "C", "charge": "0"},
        xyz_dir=xyz_dir,
    )

    assert case["provider_error"] is None
    assert case["xyz_block"] is not None
    assert case["reference_error"] is not None
    assert "Reference SMILES element counts differ from XYZ" in case["reference_error"]


def test_summary_count_excludes_uncomparable_reference_errors_but_keeps_timing() -> None:
    summary = summarize_results(
        [
            BenchmarkResult(
                case_idx=1,
                method_id="fake",
                input_smiles="C",
                ground_truth_smiles=None,
                status="ok",
                error=None,
                predicted_smiles="C",
                equivalent=None,
                equivalence_method=None,
                timing_ms_total=4.0,
                timing_ms_breakdown={"method_ms": 4.0},
                case_id="A",
                comparison_skipped=True,
                comparison_skip_reason="reference error",
            ),
            BenchmarkResult(
                case_idx=2,
                method_id="fake",
                input_smiles="C",
                ground_truth_smiles="C",
                status="ok",
                error=None,
                predicted_smiles="C",
                equivalent=True,
                equivalence_method="ideal",
                timing_ms_total=2.0,
                timing_ms_breakdown={"method_ms": 2.0},
                case_id="B",
            ),
            BenchmarkResult(
                case_idx=3,
                method_id="fake",
                input_smiles="C",
                ground_truth_smiles="C",
                status="error",
                error="method failed",
                predicted_smiles=None,
                equivalent=None,
                equivalence_method=None,
                timing_ms_total=3.0,
                timing_ms_breakdown={"method_ms": 3.0},
                case_id="C",
            ),
        ]
    )

    assert summary[0]["count"] == 2
    assert summary[0]["success_count"] == 1
    assert summary[0]["fail_count"] == 1
    assert summary[0]["skip_count"] == 0
    assert summary[0]["comparison_skip_count"] == 1
    assert summary[0]["avg_ms_total"] == 3.0


def test_reference_element_mismatch_skips_only_comparison_not_method(
    tmp_path: Path,
) -> None:
    xyz_dir = tmp_path / "xyz"
    xyz_dir.mkdir()
    (xyz_dir / "A.xyz").write_text(
        "2\nA\nC 0 0 0\nH 0 0 1\n",
        encoding="utf-8",
    )
    case = _build_case(
        1,
        {"id": "A", "smiles": "C", "charge": "0"},
        xyz_dir=xyz_dir,
    )
    calls = {"count": 0}

    class _Output:
        status = "ok"
        error = None
        predicted_smiles = "C"
        rdkit_mol = None
        equivalent = True
        equivalence_method = "method-native"
        timing_ms_breakdown = {"method_ms": 1.0}

    def _runner(_case):
        calls["count"] += 1
        return _Output()

    result = _run_case_method(
        case,
        "fake",
        _runner,
        case_timeout_seconds=1.0,
    )

    assert calls["count"] == 1
    assert result.status == "ok"
    assert result.error is None
    assert result.predicted_smiles == "C"
    assert result.case_id == "A"
    assert result.equivalent is None
    assert result.equivalence_method is None
    assert result.comparison_skipped is True
    assert "Reference SMILES element counts differ from XYZ" in str(
        result.comparison_skip_reason
    )


def test_tmqmg_cpp_all_accelerations_worker_survives_target_bucket_parallelism(
    tmp_path: Path,
) -> None:
    pytest.importorskip("molgr._core.pipeline")
    csv_path = Path("/mnt/e/download/tmQMg_properties_and_targets.csv")
    xyz_dir = Path("/mnt/e/download/tmQMg_xyz/xyz")
    if not csv_path.exists() or not xyz_dir.exists():
        pytest.skip("tmQMg benchmark dataset is not available")

    with csv_path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    case_row = next((row for row in rows if row.get("id") == "ABAJAP"), None)
    if case_row is None or not (xyz_dir / "ABAJAP.xyz").exists():
        pytest.skip("tmQMg ABAJAP crash-regression case is not available")

    payload_path = tmp_path / "payload.json"
    out_jsonl = tmp_path / "results.jsonl"
    repeats = 8
    payload_path.write_text(
        json.dumps(
            {
                "method_id": "molgr_cpp",
                "cases": [
                    {"row_index": 5, "row": case_row}
                    for _ in range(repeats)
                ],
                "xyz_dir": str(xyz_dir),
                "case_timeout_seconds": 2.0,
                "cpp_backend_config": _build_cpp_backend_config_payload(
                    use_all_accelerations=True
                ),
            },
            ensure_ascii=True,
        ),
        encoding="utf-8",
    )
    out_jsonl.write_text("", encoding="utf-8")

    env = os.environ.copy()
    env["MOLGR_TMQMG_SUBPROCESS_PAYLOAD_PATH"] = str(payload_path)
    env["MOLGR_TMQMG_SUBPROCESS_OUT"] = str(out_jsonl)
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "benchmarks.tmqmg_xyz_benchmark.run",
            "--subprocess-worker",
        ],
        cwd=Path(__file__).resolve().parents[1],
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert completed.returncode == 0, completed.stderr
    result_rows = [
        json.loads(line)
        for line in out_jsonl.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert len(result_rows) == repeats
    assert {row["method_id"] for row in result_rows} == {"molgr_cpp"}
    assert {row["case_id"] for row in result_rows} == {"ABAJAP"}
