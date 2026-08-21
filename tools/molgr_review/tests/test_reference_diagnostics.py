import sys
from pathlib import Path


APP_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_DIR))

from tools.molgr_review.reference_diagnostics import (  # noqa: E402
    classify_reference_problem,
    comparison_skip_reasons,
    reference_formula_conserved,
    reference_metal_charges,
    site_metal_valences,
)
from tools.molgr_review.server import _row_dict  # noqa: E402


_XYZ_ZN = "\n".join(
    [
        "3",
        "",
        "Zn 0.0 0.0 0.0",
        "O 0.0 0.0 1.9",
        "O 0.0 0.0 -1.9",
    ]
)


def test_formula_mismatch_is_kept_separate_from_metal_valence() -> None:
    conserved, detail = reference_formula_conserved("[O-2].[O-2].[O-2].[Zn+2]", _XYZ_ZN)
    assert not conserved
    assert "O:xyz=2,ref=3" in detail


def test_metal_valence_out_of_range_uses_site_enumeration() -> None:
    assert site_metal_valences(_XYZ_ZN) == {"Zn": [0, 1, 2, 4]}
    assert reference_metal_charges("[O-2].[O-2].[Zn+8]") == [("Zn", 8)]


def test_metal_valence_within_enumeration_is_not_flagged() -> None:
    valences = site_metal_valences(_XYZ_ZN)
    assert valences is not None
    assert reference_metal_charges("[O-2].[O-2].[Zn+2]") == [("Zn", 2)]
    assert 2 in valences["Zn"]


def test_reference_problem_diagnostics_are_not_formula_statuses() -> None:
    assert classify_reference_problem(reference_smiles="", skip_reasons=[])[0] == (
        "missing_reference"
    )
    assert (
        classify_reference_problem(
            reference_smiles="C",
            skip_reasons=["molgr_cpp equivalence 42 timed out after 1.000s"] * 2,
        )[0]
        == "equivalence_timeout"
    )
    assert (
        classify_reference_problem(
            reference_smiles="C",
            skip_reasons=["equivalence check failed: predicted_smiles could not be reparsed"] * 2,
        )[0]
        == "candidate_reparse_failure"
    )
    assert (
        classify_reference_problem(
            reference_smiles="C", skip_reasons=[], formula_status="formula_mismatch"
        )[0]
        == "formula_mismatch"
    )


def test_comparison_skip_reasons_reads_retained_benchmark_rows() -> None:
    error = {
        "rows": {
            "py38/molgr_cpp": {
                "comparison_skipped": "True",
                "comparison_skip_reason": "molgr_cpp equivalence 42 timed out after 1.000s",
            },
            "py310/molgr_cpp": {
                "comparison_skipped": "False",
                "comparison_skip_reason": "",
            },
        }
    }

    assert comparison_skip_reasons(error) == ["molgr_cpp equivalence 42 timed out after 1.000s"]


def test_legacy_case_payload_gets_runtime_diagnostic_without_rewriting_metadata() -> None:
    row = {
        "case_id": "TIMEOUT",
        "reference_smiles": "C",
        "metadata_json": '{"reference_formula_check_status":"comparison_skipped",'
        '"error":"{\\"rows\\":{\\"py38/molgr_cpp\\":{\\"comparison_skipped\\":'
        '\\"True\\",\\"comparison_skip_reason\\":\\"molgr_cpp equivalence 1 timed out '
        'after 1.000s\\"}}}"}',
    }

    payload = _row_dict(row)

    assert payload is not None
    assert payload["reference_formula_check_status"] == "comparison_skipped"
    assert payload["reference_diagnostic_group"] == "equivalence_timeout"
