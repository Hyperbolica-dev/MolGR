from __future__ import annotations

import hashlib
import json
import random
from pathlib import Path
from typing import Any

from rdkit import Chem


def _electronic_state(mol: Chem.Mol) -> tuple[int, int]:
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    radicals = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    return charge, radicals


def load_cases(
    path: Path, *, limit: int | None = None, seed: int = 20260825
) -> list[dict[str, Any]]:
    """Load an immutable, line-delimited fixture and deterministically select molecules.

    The fixture contains exactly one source conformer per molecule.  Sampling is by
    molecule, never by conformer, and does not repair invalid records.
    """
    records: list[dict[str, Any]] = []
    for line_no, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        try:
            record = json.loads(line)
            records.append(record)
        except Exception as exc:  # noqa: BLE001
            records.append({"loader_error": f"line {line_no}: {type(exc).__name__}: {exc}"})

    rng = random.Random(seed)
    rng.shuffle(records)
    if limit is not None:
        records = records[: max(0, limit)]

    cases: list[dict[str, Any]] = []
    for case_idx, record in enumerate(records, start=1):
        case: dict[str, Any] = {
            "case_idx": case_idx,
            "case_id": record.get("case_id"),
            "molecule_id": record.get("molecule_id"),
            "conformer_id": record.get("conformer_id"),
            "input_smiles": record.get("reference_smiles", ""),
            "ground_truth_smiles": record.get("reference_smiles"),
            "ground_truth_rdmol": None,
            "xyz_block": record.get("xyz"),
            "total_charge": record.get("total_charge"),
            "total_radical_electrons": record.get("spin_multiplicity", 1) - 1,
            "spin_multiplicity": record.get("spin_multiplicity"),
            "source_metadata": record.get("source_metadata", {}),
            "provider_error": record.get("loader_error"),
        }
        try:
            if case["provider_error"]:
                raise ValueError(case["provider_error"])
            mol = Chem.MolFromSmiles(str(case["ground_truth_smiles"]))
            if mol is None:
                raise ValueError("RDKit failed to parse reference_smiles")
            charge, radicals = _electronic_state(mol)
            if charge != case["total_charge"]:
                raise ValueError(
                    f"reference charge {charge} != fixture charge {case['total_charge']}"
                )
            if radicals != case["total_radical_electrons"]:
                raise ValueError(
                    f"reference radicals {radicals} != multiplicity-derived "
                    f"radicals {case['total_radical_electrons']}"
                )
            case["ground_truth_rdmol"] = mol
        except Exception as exc:  # noqa: BLE001
            case["provider_error"] = f"{type(exc).__name__}: {exc}"
        cases.append(case)
    return cases


def stable_molecule_id(smiles: str) -> str:
    return hashlib.sha256(smiles.encode("utf-8")).hexdigest()[:16]
