from __future__ import annotations

from pathlib import Path
from typing import Any

from rdkit import Chem, RDLogger


RDLogger.DisableLog("rdApp.*")  # type: ignore


_MOLFILE_SUFFIXES = {".mol", ".molfile", ".sdf"}


def _total_charge_and_radicals(mol: Chem.Mol) -> tuple[int, int]:
    charge = 0
    radical = 0
    for atom in mol.GetAtoms():  # pyright: ignore[reportCallIssue]
        charge += int(atom.GetFormalCharge())
        radical += int(atom.GetNumRadicalElectrons())
    return charge, radical


def _canonical_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _iter_molfile_paths(input_path: Path) -> list[Path]:
    if input_path.is_file():
        return [input_path]
    if not input_path.is_dir():
        raise FileNotFoundError(f"input path does not exist: {input_path}")

    paths = [
        path
        for path in input_path.rglob("*")
        if path.is_file() and path.suffix.lower() in _MOLFILE_SUFFIXES
    ]
    paths.sort(key=lambda path: str(path))
    return paths


def _build_case(molfile_path: Path, case_idx: int) -> dict[str, Any]:
    case: dict[str, Any] = {
        "case_idx": case_idx,
        "input_smiles": "",
        "xyz_block": None,
        "total_charge": None,
        "total_radical_electrons": None,
        "provider_error": None,
        "ground_truth_rdmol": None,
        "ground_truth_smiles": "",
        "ground_truth_path": str(molfile_path),
    }

    try:
        mol = Chem.MolFromMolFile(
            str(molfile_path),
            sanitize=True,
            removeHs=False,
            strictParsing=False,
        )
        if mol is None:
            raise ValueError("RDKit failed to parse molfile")
        if mol.GetNumConformers() == 0:
            raise ValueError("molfile has no conformer coordinates")

        total_charge, total_radical_electrons = _total_charge_and_radicals(mol)
        smiles = _canonical_smiles(mol)

        case["input_smiles"] = smiles
        case["xyz_block"] = Chem.MolToXYZBlock(mol)
        case["total_charge"] = total_charge
        case["total_radical_electrons"] = total_radical_electrons
        case["ground_truth_rdmol"] = mol
        case["ground_truth_smiles"] = smiles
    except Exception as exc:
        case["provider_error"] = f"{type(exc).__name__}: {exc}"

    return case


def load_molfile_cases(input_path: Path, limit: int | None = None) -> list[dict[str, Any]]:
    molfile_paths = _iter_molfile_paths(input_path)
    if limit is not None:
        molfile_paths = molfile_paths[: max(limit, 0)]

    cases: list[dict[str, Any]] = []
    for case_idx, molfile_path in enumerate(molfile_paths, start=1):
        cases.append(_build_case(molfile_path=molfile_path, case_idx=case_idx))
    return cases
