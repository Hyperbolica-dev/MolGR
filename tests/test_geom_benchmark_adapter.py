from __future__ import annotations

import json
from pathlib import Path

from rdkit import Chem

from benchmarks.geom_xyz_benchmark.acquire_smoke import _drugs_stratum, _fixture_record
from benchmarks.geom_xyz_benchmark.adapter import load_cases, stable_molecule_id
from benchmarks.geom_xyz_benchmark.formal_run import (
    BoundedReview,
)
from benchmarks.geom_xyz_benchmark.formal_run import (
    _percentile as _formal_percentile,
)
from benchmarks.geom_xyz_benchmark.formal_run import (
    _size_stratum as _formal_size_stratum,
)
from benchmarks.geom_xyz_benchmark.run import _percentile, _size_stratum


def test_geom_adapter_is_deterministic_and_preserves_state(tmp_path: Path) -> None:
    path = tmp_path / "fixture.jsonl"
    rows = [
        {
            "case_id": f"c{i}",
            "molecule_id": f"m{i}",
            "conformer_id": "0",
            "xyz": "3\n\nO 0 0 0\nH 0 0 1\nH 0 1 0\n",
            "total_charge": 0,
            "spin_multiplicity": 1,
            "reference_smiles": "O",
            "source_metadata": {"source": "test"},
        }
        for i in range(4)
    ]
    path.write_text("\n".join(json.dumps(row) for row in rows), encoding="utf-8")
    first = load_cases(path, limit=2, seed=7)
    second = load_cases(path, limit=2, seed=7)
    assert [case["case_id"] for case in first] == [case["case_id"] for case in second]
    assert all(case["total_charge"] == 0 for case in first)
    assert all(case["spin_multiplicity"] == 1 for case in first)
    assert stable_molecule_id("O") == stable_molecule_id("O")


def test_geom_adapter_records_invalid_reference_without_repair(tmp_path: Path) -> None:
    path = tmp_path / "bad.jsonl"
    path.write_text(
        json.dumps(
            {
                "case_id": "bad",
                "reference_smiles": "not-smiles",
                "total_charge": 0,
                "spin_multiplicity": 1,
                "xyz": "",
            }
        ),
        encoding="utf-8",
    )
    case = load_cases(path)[0]
    assert case["ground_truth_rdmol"] is None
    assert "failed to parse" in case["provider_error"]


def test_drugs_fixture_selects_lowest_energy_not_source_order() -> None:
    xyz_high = [
        [6, 0.0, 0.0, 0.0],
        [1, 1.0, 0.0, 0.0],
        [1, 0.0, 1.0, 0.0],
        [1, 0.0, 0.0, 1.0],
        [1, -1.0, 0.0, 0.0],
    ]
    xyz_low = [[row[0], row[1] + 0.1, row[2], row[3]] for row in xyz_high]
    record, reason = _fixture_record(
        "C",
        {
            "conformers": [
                {"xyz": xyz_high, "relativeenergy": 2.0},
                {"xyz": xyz_low, "relativeenergy": 0.0},
            ]
        },
        dataset="drugs",
    )
    assert reason == "eligible"
    assert record is not None
    assert record["conformer_id"] == "1"
    assert record["source_metadata"]["relativeenergy_kcal_mol"] == 0.0
    assert record["source_metadata"]["dataset"] == "GEOM-DRUGS"


def test_drugs_size_strata_boundaries() -> None:
    expected = {
        15: "heavy_01_15",
        16: "heavy_16_25",
        26: "heavy_26_35",
        36: "heavy_36_50",
        51: "heavy_51_plus",
    }
    for heavy_atoms, stratum in expected.items():
        mol = Chem.MolFromSmiles("C" * heavy_atoms)
        assert mol is not None
        assert _drugs_stratum(mol) == stratum


def test_pilot_size_strata_and_percentile() -> None:
    assert [_size_stratum(value) for value in (15, 16, 25, 26, 35, 36, 50, 51)] == [
        "01_15",
        "16_25",
        "16_25",
        "26_35",
        "26_35",
        "36_50",
        "36_50",
        "51_plus",
    ]
    assert _percentile([1.0, 2.0, 3.0, 4.0], 0.5) == 2.5
    assert _formal_size_stratum(51) == "51_plus"
    assert _formal_percentile([1.0, 2.0, 3.0, 4.0], 0.95) == 3.8499999999999996


def test_formal_review_reservoir_is_bounded_and_order_independent() -> None:
    rows = [{"case_id": f"case-{index}"} for index in range(10)]
    first = BoundedReview({"non_equivalent": 3})
    second = BoundedReview({"non_equivalent": 3})
    for row in rows:
        first.add("non_equivalent", row)
    for row in reversed(rows):
        second.add("non_equivalent", row)
    assert sorted(row["case_id"] for _, row in first.rows()) == sorted(
        row["case_id"] for _, row in second.rows()
    )
