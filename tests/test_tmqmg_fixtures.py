from __future__ import annotations

import csv
from pathlib import Path

import pytest
from rdkit import Chem

from molgr.fallback.pipeline.reconstruct_with_metals import xyz2omol_state
from molgr.interface import xyz_to_rdmol
from scripts.tmqmg_reference_formula_check import (
    _read_xyz_element_counts,
    _smiles_element_counts_with_h,
)


_FIXTURE_ROOT = Path(__file__).parent / "data" / "tmqmg"


def _manifest(group: str) -> list[dict[str, str]]:
    with (_FIXTURE_ROOT / group / "manifest.csv").open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _canonical_smiles(xyz_block: str, charge: int, backend: str) -> str:
    mol = xyz_to_rdmol(
        xyz_block,
        charge,
        1,
        backend=backend,
        make_dative_bonds=True,
        make_stereochemistry=True,
    )
    return Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True)


@pytest.mark.parametrize("case", _manifest("reconstruction"), ids=lambda case: case["case_id"])
def test_tmqmg_reconstruction_fixture_exercises_expected_path(case: dict[str, str]) -> None:
    xyz_path = _FIXTURE_ROOT / "reconstruction" / case["xyz_file"]
    xyz_block = xyz_path.read_text(encoding="utf-8")
    state = xyz2omol_state(xyz_block, int(case["charge"]), 0)

    assert state is not None
    assert state.no_metal_state is not None
    no_metal_state = state.no_metal_state
    actions = str(no_metal_state.metadata.get("neighbor_radical_actions", ""))
    classification = case["classification"]
    if classification == "primary_no_neighbor_branch":
        assert not actions
        assert int(no_metal_state.metadata.get("recovery_tier", 0)) == 0
    elif classification == "neighbor_bond_order":
        assert "bond_order" in actions
    elif classification == "radical_charge_localization":
        assert "assign_negative_charges_from_radicals" in no_metal_state.phase_history
    elif classification == "neighbor_charge_separation":
        assert "charge_separation" in actions
    elif classification.startswith("recovery_deformed_pi"):
        assert no_metal_state.metadata["recovery_tier"] == 1
    elif classification.startswith("recovery_bond_break"):
        assert no_metal_state.metadata["recovery_tier"] == 2
    else:
        raise AssertionError(f"Unhandled fixture classification: {classification}")


@pytest.mark.parametrize("case", _manifest("reconstruction"), ids=lambda case: case["case_id"])
def test_tmqmg_reconstruction_fixture_cpp_python_parity(case: dict[str, str]) -> None:
    xyz_block = (_FIXTURE_ROOT / "reconstruction" / case["xyz_file"]).read_text(encoding="utf-8")
    charge = int(case["charge"])

    assert _canonical_smiles(xyz_block, charge, "cpp") == _canonical_smiles(
        xyz_block,
        charge,
        "python",
    )


def test_tmqmg_source_issue_fixtures_are_separate_and_documented() -> None:
    reconstruction_ids = {case["case_id"] for case in _manifest("reconstruction")}
    source_issues = _manifest("source_issues")

    assert source_issues
    assert reconstruction_ids.isdisjoint(case["case_id"] for case in source_issues)
    assert all(case["reason"].strip() for case in source_issues)
    assert all(
        (_FIXTURE_ROOT / "source_issues" / case["xyz_file"]).exists() for case in source_issues
    )


def test_tmqmg_missing_reference_graph_fixture_is_reproducible() -> None:
    missing_references = [
        case
        for case in _manifest("source_issues")
        if case["classification"] == "missing_reference_graph"
    ]

    assert [case["case_id"] for case in missing_references] == ["ABARAX"]
    assert missing_references[0]["reference_smiles"] == ""


def test_tmqmg_manually_accepted_reference_issues_are_explicit() -> None:
    accepted_ids = {
        case["case_id"]
        for case in _manifest("source_issues")
        if case["classification"] == "accepted_reference_issue"
    }

    assert accepted_ids == {"ABATEC", "ABEGOD", "ABETIK", "ABETOQ"}


@pytest.mark.parametrize(
    "case",
    [
        case
        for case in _manifest("source_issues")
        if case["classification"] == "reference_formula_mismatch"
    ],
    ids=lambda case: case["case_id"],
)
def test_tmqmg_reference_formula_mismatch_fixture_is_reproducible(case: dict[str, str]) -> None:
    xyz_path = _FIXTURE_ROOT / "source_issues" / case["xyz_file"]

    assert _read_xyz_element_counts(xyz_path) != _smiles_element_counts_with_h(
        case["reference_smiles"]
    )
