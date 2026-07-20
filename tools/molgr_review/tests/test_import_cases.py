from __future__ import annotations

import json
import sys
from pathlib import Path


APP_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_DIR))

from import_cases import _normalize_row  # noqa: E402


def test_normalize_generic_review_queue_row() -> None:
    row = _normalize_row(
        {
            "case_id": "WATER",
            "source": "smoke",
            "category": "graph_check",
            "xyz_path": "/tmp/WATER.xyz",
            "total_charge": "0",
            "total_radical_electrons": "0",
            "spin_multiplicity": "1",
            "reference_smiles": "O",
            "candidate_smiles": "O",
            "candidate_status": "ok",
        },
        1,
    )

    assert row["case_id"] == "WATER"
    assert row["source"] == "smoke"
    assert row["candidate_smiles"] == "O"
    assert row["spin_multiplicity"] == 1


def test_dataset_specific_columns_remain_metadata() -> None:
    row = _normalize_row(
        {
            "case_id": "ABEGOD",
            "source": "tmqmg",
            "category": "candidate_failed",
            "total_charge": "-1",
            "reference_smiles": "[Cl-]",
            "candidate_smiles": "[Cl-]",
            "candidate_status": "error",
            "review_category": "molgr_failed",
            "xyz_path": "/tmp/ABEGOD.xyz",
            "charge": "-1",
            "reference_smiles_input": "[Cl-]",
            "molgr_smiles_canonical": "[Cl-]",
            "molgr_status": "error",
        },
        7,
    )

    assert row["case_id"] == "ABEGOD"
    assert row["category"] == "candidate_failed"
    assert row["total_charge"] == -1
    assert row["reference_smiles"] == "[Cl-]"
    assert row["candidate_smiles"] == "[Cl-]"
    assert row["candidate_status"] == "error"
    assert json.loads(row["metadata_json"])["molgr_status"] == "error"
