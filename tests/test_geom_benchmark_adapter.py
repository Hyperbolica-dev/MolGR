from __future__ import annotations

import json
from pathlib import Path

from benchmarks.geom_xyz_benchmark.adapter import load_cases, stable_molecule_id


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
