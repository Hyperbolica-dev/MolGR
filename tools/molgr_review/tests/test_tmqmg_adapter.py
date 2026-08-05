from __future__ import annotations

import csv
import sys
from dataclasses import replace
from pathlib import Path


APP_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_DIR))

from adapters.tmqmg import (  # noqa: E402
    BENCHMARK_RDKIT_REQUIREMENT,
    BENCHMARK_RDKIT_RUNTIME_VERSION,
    REFERENCE_FORMULA_MISMATCH_PREFIX,
    BenchmarkRow,
    PythonRun,
    _benchmark_build_dependencies,
    _build_review_rows,
    _case_issue,
    _merge_partial_review_rows,
    _python_results_equivalent,
    _reference_formula_mismatch_fields,
    _validate_benchmark_environment_versions,
)


def _row(case_id: str, row_index: int, candidate_smiles: str) -> dict[str, str]:
    return {
        "case_id": case_id,
        "row_index": str(row_index),
        "candidate_smiles": candidate_smiles,
    }


def _benchmark_row(label: str) -> BenchmarkRow:
    return BenchmarkRow(
        label=label,
        method_id="molgr_cpp",
        case_idx=1,
        case_id="FORMULA",
        input_smiles="C",
        ground_truth_smiles="C",
        status="ok",
        error="",
        predicted_smiles="[CH3]",
        equivalent="",
        equivalence_method="",
        comparison_skipped="True",
        comparison_skip_reason=(
            REFERENCE_FORMULA_MISMATCH_PREFIX + " reference=C:1,H:4; xyz=C:1,H:3"
        ),
        timing_ms_total=1.0,
    )


def _write_xyz(path: Path) -> None:
    path.write_text(
        "4\nformula mismatch\nC 0 0 0\nH 1 0 0\nH 0 1 0\nH 0 0 1\n",
        encoding="utf-8",
    )


def _write_result(path: Path, row: BenchmarkRow) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=(
                "case_idx",
                "id",
                "method_id",
                "input_smiles",
                "ground_truth_smiles",
                "status",
                "error",
                "predicted_smiles",
                "equivalent",
                "equivalence_method",
                "comparison_skipped",
                "comparison_skip_reason",
                "timing_ms_total",
            ),
        )
        writer.writeheader()
        writer.writerow(
            {
                "case_idx": row.case_idx,
                "id": row.case_id,
                "method_id": row.method_id,
                "input_smiles": row.input_smiles,
                "ground_truth_smiles": row.ground_truth_smiles,
                "status": row.status,
                "error": row.error,
                "predicted_smiles": row.predicted_smiles,
                "equivalent": row.equivalent,
                "equivalence_method": row.equivalence_method,
                "comparison_skipped": row.comparison_skipped,
                "comparison_skip_reason": row.comparison_skip_reason,
                "timing_ms_total": row.timing_ms_total,
            }
        )


def test_partial_refresh_replaces_only_cases_in_scope() -> None:
    existing = [
        _row("KEEP", 1, "C"),
        _row("REFRESH", 2, "N"),
        _row("REMOVE", 3, "O"),
    ]
    refreshed = [_row("REFRESH", 2, "[N+]")]

    merged = _merge_partial_review_rows(
        existing,
        refreshed,
        refreshed_case_ids={"REFRESH", "REMOVE"},
    )

    assert [(row["case_id"], row["candidate_smiles"]) for row in merged] == [
        ("KEEP", "C"),
        ("REFRESH", "[N+]"),
    ]


def test_partial_refresh_adds_new_review_case_in_dataset_order() -> None:
    merged = _merge_partial_review_rows(
        [_row("LATE", 20, "C")],
        [_row("EARLY", 10, "N")],
        refreshed_case_ids={"EARLY"},
    )

    assert [row["case_id"] for row in merged] == ["EARLY", "LATE"]


def test_cross_openbabel_benchmark_pins_shared_rdkit() -> None:
    dependencies = _benchmark_build_dependencies()

    assert BENCHMARK_RDKIT_REQUIREMENT in dependencies
    assert "openbabel-wheel==3.1.1.22; python_version < '3.9' and sys_platform == 'win32'" in dependencies
    assert "openbabel-wheel; python_version < '3.10' and (sys_platform != 'win32' or python_version >= '3.9')" in dependencies
    assert "openbabel>=3.2.0; python_version >= '3.10'" in dependencies


def test_cross_openbabel_benchmark_rejects_rdkit_drift(tmp_path: Path) -> None:
    run = PythonRun("py310", "python3.10", tmp_path / "py310")

    try:
        _validate_benchmark_environment_versions(
            run,
            {
                "python": "3.10.18",
                "openbabel": "3.2.1",
                "rdkit": "2026.03.4",
            },
        )
    except SystemExit as exc:
        assert BENCHMARK_RDKIT_RUNTIME_VERSION in str(exc)
    else:
        raise AssertionError("RDKit version drift was not rejected")


def test_reference_formula_mismatch_has_dedicated_category() -> None:
    rows = {
        ("py38", "molgr_cpp"): _benchmark_row("py38"),
        ("py310", "molgr_cpp"): _benchmark_row("py310"),
    }

    category, details = _case_issue(
        rows,
        labels=("py38", "py310"),
        method_ids=("molgr_cpp",),
    )

    assert category == "reference_formula_mismatch"
    assert len(details["reference_formula_mismatches"]) == 2


def test_reference_formula_mismatch_fields_recompute_hydrogen_counts(tmp_path: Path) -> None:
    xyz_path = tmp_path / "FORMULA.xyz"
    _write_xyz(xyz_path)

    fields = _reference_formula_mismatch_fields("C", xyz_path)

    assert fields["xyz_formula"] == "C:1,H:3"
    assert fields["reference_formula_with_h"] == "C:1,H:4"
    assert fields["reference_formula_mismatch_detail"] == "H:xyz=3,ref=4"
    assert fields["reference_answer_wrong"] == "True"
    assert fields["reference_answer_status"] == "formula_mismatch"


def test_review_row_marks_reference_formula_mismatch(tmp_path: Path) -> None:
    xyz_dir = tmp_path / "xyz"
    xyz_dir.mkdir()
    _write_xyz(xyz_dir / "FORMULA.xyz")
    result_paths = {
        "py38": tmp_path / "py38.csv",
        "py310": tmp_path / "py310.csv",
    }
    for label, path in result_paths.items():
        _write_result(path, _benchmark_row(label))

    rows, summary = _build_review_rows(
        metadata={(1, "FORMULA"): {"smiles": "C", "n_atoms": "4"}},
        result_paths=result_paths,
        xyz_dir=xyz_dir,
        method_ids=("molgr_cpp",),
        python_version_comparison="graph",
    )

    assert summary["category_counts"] == {"reference_formula_mismatch": 1}
    assert len(rows) == 1
    assert rows[0]["category"] == "reference_formula_mismatch"
    assert rows[0]["reference_formula_check_status"] == "formula_mismatch"
    assert rows[0]["reference_answer_wrong"] == "True"
    assert rows[0]["reference_formula_mismatch_detail"] == "H:xyz=3,ref=4"


def test_review_row_marks_both_boron_cluster_answers_not_assessable(tmp_path: Path) -> None:
    xyz_dir = tmp_path / "xyz"
    xyz_dir.mkdir()
    (xyz_dir / "ADOCOL.xyz").write_text(
        "4\nboron cluster\nB 0 0 0\nB 1 0 0\nB 0 1 0\nB 0 0 1\n",
        encoding="utf-8",
    )
    result_paths = {
        "py38": tmp_path / "py38.csv",
        "py310": tmp_path / "py310.csv",
    }
    for label, path in result_paths.items():
        benchmark_row = replace(
            _benchmark_row(label),
            case_id="ADOCOL",
            input_smiles="B.B.B.B",
            ground_truth_smiles="B.B.B.B",
            predicted_smiles="B.B.B.B",
        )
        _write_result(path, benchmark_row)

    rows, summary = _build_review_rows(
        metadata={(1, "ADOCOL"): {"smiles": "B.B.B.B", "n_atoms": "4"}},
        result_paths=result_paths,
        xyz_dir=xyz_dir,
        method_ids=("molgr_cpp",),
        python_version_comparison="graph",
    )

    assert summary["category_counts"] == {"no_clear_evidence_boron_cluster": 1}
    assert len(rows) == 1
    assert rows[0]["category"] == "no_clear_evidence_boron_cluster"
    assert rows[0]["reference_answer_wrong"] == "False"
    assert rows[0]["reference_answer_status"] == "not_assessable"
    assert rows[0]["tmqmg_answer_assessment"] == "not_assessable"
    assert rows[0]["molgr_answer_assessment"] == "not_assessable"
    assert rows[0]["equivalent"] == ""
    assert rows[0]["effective_equivalent"] == ""
    assert "3-center-2-electron" in rows[0]["accuracy_assessment_reason"]


def test_python_versions_compare_candidate_graphs_not_reference_verdicts() -> None:
    py38 = replace(
        _benchmark_row("py38"),
        predicted_smiles="C[N+](C)(C)C",
        equivalent="",
        comparison_skipped="True",
        comparison_skip_reason="equivalence timed out",
    )
    py310 = replace(
        _benchmark_row("py310"),
        predicted_smiles="C[N+](C)(C)C",
        equivalent="True",
        comparison_skipped="False",
        comparison_skip_reason="",
    )

    equivalent, reason = _python_results_equivalent(py38, py310)

    assert equivalent is True
    assert reason == "identical_smiles"
