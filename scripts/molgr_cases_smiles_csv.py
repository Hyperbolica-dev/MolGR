"""
Author: TMJ
Date: 2026-02-24 11:44:45
LastEditors: TMJ
LastEditTime: 2026-02-24 23:33:23
Description: 请填写简介
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom


RDLogger.DisableLog("rdApp.*")  # type: ignore


def _total_charge_and_radicals(mol: Chem.Mol) -> tuple[int, int]:
    charge = 0
    radical = 0
    for atom in mol.GetAtoms():  # pyright: ignore[reportCallIssue]
        charge += int(atom.GetFormalCharge())
        radical += int(atom.GetNumRadicalElectrons())
    return charge, radical


def _build_case(smiles: str, case_idx: int) -> dict[str, Any]:
    case: dict[str, Any] = {
        "case_idx": case_idx,
        "input_smiles": smiles,
        "ground_truth_rdmol": None,
        "ground_truth_smiles": None,
        "xyz_block": None,
        "total_charge": None,
        "total_radical_electrons": None,
        "provider_error": None,
    }

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("RDKit failed to parse SMILES")
        ground_truth_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

        mol_h = Chem.AddHs(mol)
        embed_code = rdDistGeom.EmbedMolecule(mol_h)  # pyright: ignore[reportCallIssue]
        if int(embed_code) != 0:
            raise ValueError(f"RDKit EmbedMolecule failed: code={embed_code}")

        total_charge, total_radical_electrons = _total_charge_and_radicals(mol_h)
        xyz_block = Chem.MolToXYZBlock(mol_h)

        case["ground_truth_rdmol"] = mol
        case["ground_truth_smiles"] = ground_truth_smiles
        case["xyz_block"] = xyz_block
        case["total_charge"] = total_charge
        case["total_radical_electrons"] = total_radical_electrons
    except Exception as exc:
        case["provider_error"] = f"{type(exc).__name__}: {exc}"

    return case


def load_smiles_csv_cases(input_path: Path, limit: int | None = None) -> list[dict[str, Any]]:
    raw_lines = input_path.read_text(encoding="utf-8").splitlines()
    lines = [line.strip() for line in raw_lines if line.strip()]

    if not lines:
        return []

    header0 = lines[0].strip().lower()
    header1 = lines[1].strip().lower() if len(lines) > 1 else ""

    if header0 == "general" and header1 in {"canonicalsmiles", "canonicalsmi", "smiles"}:
        smiles_lines = lines[2:]
    elif header0 in {"canonicalsmiles", "canonicalsmi", "smiles"}:
        smiles_lines = lines[1:]
    else:
        smiles_lines = lines

    if limit is not None:
        smiles_lines = smiles_lines[: max(limit, 0)]

    cases: list[dict[str, Any]] = []
    for case_idx, smiles in enumerate(smiles_lines, start=1):
        cases.append(_build_case(smiles=smiles, case_idx=case_idx))
    return cases
