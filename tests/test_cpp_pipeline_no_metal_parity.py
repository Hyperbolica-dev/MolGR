# pyright: reportMissingImports=false

from __future__ import annotations

import csv
from pathlib import Path

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr import _core
from molgr.fallback.pipeline.reconstruct_without_metals import (
    xyz_to_omol_no_metal as py_xyz_to_omol_no_metal,
)
from molgr.fallback.utils.tools import typed_lru_cache
from molgr.utils.converter import mol_data_to_rdkit
from molgr.utils.equivalence import check_equivalence


@typed_lru_cache(maxsize=1024, typed=True)
def cpp_xyz_to_omol_no_metal(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
) -> _core.utils.MoleculeData | None:
    return _core.pipeline.reconstruct_without_metals.xyz_to_omol_no_metal(
        xyz_block,
        total_charge,
        total_radical_electrons,
    )


def _get_ptr(obmol: ob.OBMol) -> int:
    this = getattr(obmol, "this", None)
    if this is not None:
        return int(this)  # type: ignore[arg-type]
    return int(obmol)  # type: ignore[arg-type]


def _charge_and_radicals(mol: pybel.Molecule) -> tuple[int, int]:
    charge = 0
    radicals = 0
    for atom in mol.atoms:
        charge += atom.OBAtom.GetFormalCharge()
        radicals += atom.OBAtom.GetSpinMultiplicity()
    return charge, radicals


def _load_curated_smiles() -> list[str]:
    csv_path = Path(__file__).with_name("test_cases.csv")
    with csv_path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    curated_rows = [1, 2, 5, 10, 17, 26, 34, 42]
    return [rows[idx - 1]["smiles"] for idx in curated_rows]


@pytest.mark.parametrize("smiles", _load_curated_smiles())
def test_xyz_to_omol_no_metal_cpp_matches_fallback(smiles: str) -> None:
    seed = pybel.readstring("smi", smiles)
    xyz_block = str(seed.write("xyz"))
    total_charge, total_radical_electrons = _charge_and_radicals(seed)

    py_result = py_xyz_to_omol_no_metal(xyz_block, total_charge, total_radical_electrons)
    cpp_result = cpp_xyz_to_omol_no_metal(xyz_block, total_charge, total_radical_electrons)

    assert (py_result is None) == (cpp_result is None)
    if py_result is None or cpp_result is None:
        return

    py_mol_data = _core.utils.extract_molecule_data(_get_ptr(py_result.OBMol))

    py_rdmol = mol_data_to_rdkit(py_mol_data)
    cpp_rdmol = mol_data_to_rdkit(cpp_result)

    equivalent, info = check_equivalence(py_rdmol, cpp_rdmol)
    assert equivalent, info.reason


def test_cpp_xyz_to_omol_no_metal_cache_clear() -> None:
    xyz_block = """2
H2
H 0.0 0.0 0.0
H 0.0 0.0 0.74
"""

    first = cpp_xyz_to_omol_no_metal(xyz_block, 0, 0)
    second = cpp_xyz_to_omol_no_metal(xyz_block, 0, 0)
    assert first is second

    cpp_xyz_to_omol_no_metal.cache_clear()  # pyright: ignore[reportFunctionMemberAccess]

    third = cpp_xyz_to_omol_no_metal(xyz_block, 0, 0)
    assert third is not second
