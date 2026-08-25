from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from benchmarks.bde_db_benchmark.adapter import (
    load_bde_cases,
    load_bde_cases_by_record_index,
)
from benchmarks.bde_db_benchmark.run import _atom_identity_mapping, _run_case
from benchmarks.smiles_xyz_benchmark.methods.base import BenchmarkMethod, MethodRunOutput
from benchmarks.smiles_xyz_benchmark.methods.postprocess import remove_hs_without_sanitize


def _write_sdf(path: Path, smiles: list[str]) -> None:
    plain_path = path.with_suffix("")
    writer = Chem.SDWriter(str(plain_path))
    for index, value in enumerate(smiles, start=1):
        mol = Chem.AddHs(Chem.MolFromSmiles(value))
        assert AllChem.EmbedMolecule(mol, randomSeed=index) == 0
        mol.SetProp("_Name", f"case-{index}")
        mol.SetProp("SMILES", value)
        writer.write(mol)
    writer.close()
    with plain_path.open("rb") as source, gzip.open(path, "wb") as destination:
        shutil.copyfileobj(source, destination)


def test_stratified_loader_is_deterministic(tmp_path: Path) -> None:
    input_path = tmp_path / "test.sdf.gz"
    _write_sdf(input_path, ["CC", "[CH3]", "[NH2]", "[OH]", "[O-]", "[H]", "CCC"])

    first, first_diagnostics = load_bde_cases(input_path, limit=5, seed=17)
    second, second_diagnostics = load_bde_cases(input_path, limit=5, seed=17)

    assert [case.case_id for case in first] == [case.case_id for case in second]
    assert {case.reference_smiles for case in first}.issuperset({"[H]", "[O-]"})
    assert first_diagnostics.selected_records == 5
    assert second_diagnostics.selected_records == 5


def test_loader_respects_inclusive_record_range(tmp_path: Path) -> None:
    input_path = tmp_path / "test.sdf.gz"
    _write_sdf(input_path, ["C", "CC", "CCC", "CCCC"])

    cases, diagnostics = load_bde_cases(input_path, limit=10, start=2, end=3, seed=0)

    assert [case.source_record_index for case in cases] == [2, 3]
    assert diagnostics.scanned_records == 4


def test_exact_record_replay_preserves_old_order_and_appends_required_smiles(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "test.sdf.gz"
    _write_sdf(input_path, ["C", "[H]", "[CH3]", "CC"])

    cases, diagnostics = load_bde_cases_by_record_index(
        input_path,
        [4, 1, 3],
        required_smiles=["[H]"],
    )

    assert [case.source_record_index for case in cases] == [4, 1, 3, 2]
    assert [case.case_idx for case in cases] == [1, 2, 3, 4]
    assert cases[-1].reference_smiles == "[H]"
    assert diagnostics.selected_records == 4


def test_loader_records_ambiguous_multiradical_multiplicity(tmp_path: Path) -> None:
    input_path = tmp_path / "test.sdf.gz"
    _write_sdf(input_path, ["[CH2][CH2]", "CC"])

    cases, diagnostics = load_bde_cases(input_path, limit=10, seed=0)

    assert [case.case_id for case in cases] == ["case-2"]
    assert len(diagnostics.failures) == 1
    assert "not uniquely supported" in str(diagnostics.failures[0]["error"])


def test_atom_identity_mapping_guards_index_comparison(tmp_path: Path) -> None:
    input_path = tmp_path / "test.sdf.gz"
    _write_sdf(input_path, ["[CH2]C"])
    cases, _ = load_bde_cases(input_path, limit=1, seed=0)
    reference = cases[0].reference_mol
    predicted = Chem.RemoveHs(reference, sanitize=False)

    preserved, mapping = _atom_identity_mapping(reference, predicted)

    assert preserved is True
    assert mapping == {0: 0, 1: 1}

    reordered = Chem.RenumberAtoms(predicted, [1, 0])
    preserved, mapping = _atom_identity_mapping(reference, reordered)

    assert preserved is False
    assert mapping == {}


class _ThrowingMethod(BenchmarkMethod):
    def run(self, case: dict) -> MethodRunOutput:
        del case
        raise RuntimeError("expected test exception")


class _IdentityMethod(BenchmarkMethod):
    def run(self, case: dict) -> MethodRunOutput:
        mol = remove_hs_without_sanitize(case["ground_truth_rdmol"])
        return MethodRunOutput(
            status="ok",
            predicted_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
            rdkit_mol=mol,
        )


def test_runner_records_ordinary_method_exception(tmp_path: Path) -> None:
    input_path = tmp_path / "test.sdf.gz"
    _write_sdf(input_path, ["[CH3]"])
    cases, _ = load_bde_cases(input_path, limit=1, seed=0)

    result = _run_case(cases[0], _ThrowingMethod("throwing"), None)

    assert result.reconstruction_success is False
    assert result.failure_kind == "exception"
    assert "expected test exception" in str(result.error)


def test_runner_checks_formal_radical_atom_identity(tmp_path: Path) -> None:
    input_path = tmp_path / "test.sdf.gz"
    _write_sdf(input_path, ["[CH2]C"])
    cases, _ = load_bde_cases(input_path, limit=1, seed=0)

    result = _run_case(cases[0], _IdentityMethod("identity"), None)

    assert result.atom_order_preserved is True
    assert result.formal_radical_atom_index_match is True
    assert result.evaluator_decision == "equivalent"
    assert result.evaluator_relation == "normalized_graph_identity"
    assert result.evaluator_inconclusive is False
    assert result.equivalent is True
    assert result.evaluator_reason


def test_runner_records_evaluator_exception_and_continues(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_path = tmp_path / "test.sdf.gz"
    _write_sdf(input_path, ["[CH3]"])
    cases, _ = load_bde_cases(input_path, limit=1, seed=0)

    def fail_evaluation(*args: object, **kwargs: object) -> None:
        del args, kwargs
        raise RuntimeError("expected evaluator exception")

    monkeypatch.setattr("benchmarks.bde_db_benchmark.run.evaluate_equivalence", fail_evaluation)
    result = _run_case(cases[0], _IdentityMethod("identity"), None)

    assert result.reconstruction_success is True
    assert result.status == "error"
    assert result.failure_kind == "evaluator_exception"
    assert result.equivalent is None
    assert result.evaluator_decision is None
    assert "expected evaluator exception" in str(result.error)
