# pyright: reportMissingImports=false

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr import _core
from molgr.fallback.pipeline.reconstruct_with_metals import (
    _build_metal_states as py_build_metal_states,
)
from molgr.fallback.pipeline.reconstruct_with_metals import (
    combine_metal_with_omol as py_combine_metal_with_omol,
)
from molgr.fallback.pipeline.reconstruct_with_metals import (
    get_possible_metal_radicals as py_get_possible_metal_radicals,
)
from molgr.fallback.pipeline.reconstruct_with_metals import (
    xyz2omol as py_xyz2omol,
)
from molgr.fallback.utils.tools import typed_lru_cache
from molgr.utils.converter import mol_data_to_rdkit
from molgr.utils.equivalence import check_equivalence


_with_metals: Any = _core.pipeline.reconstruct_with_metals


@typed_lru_cache(maxsize=1024, typed=True)
def cpp_xyz2omol(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
) -> _core.utils.MoleculeData | None:
    return _with_metals.xyz2omol(
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


def _load_parity_cases() -> list[object]:
    csv_path = Path(__file__).with_name("test_cases.csv")
    with csv_path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    curated_rows = [1, 2, 5, 10]
    cases: list[object] = []
    for idx in curated_rows:
        smiles = rows[idx - 1]["smiles"]
        seed = pybel.readstring("smi", smiles)
        xyz_block = str(seed.write("xyz"))
        total_charge, total_radical_electrons = _charge_and_radicals(seed)
        cases.append(
            pytest.param(
                xyz_block,
                total_charge,
                total_radical_electrons,
                id=f"curated-{idx}",
            )
        )

    synthetic_metal_xyz = """3
LiCO
Li 0.0 0.0 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
"""
    cases.append(pytest.param(synthetic_metal_xyz, 0, 0, id="synthetic-li-co"))
    return cases


@pytest.mark.parametrize(
    ("xyz_block", "total_charge", "total_radical_electrons"),
    _load_parity_cases(),
)
def test_xyz2omol_cpp_matches_fallback(
    xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
) -> None:
    py_result = py_xyz2omol(xyz_block, total_charge, total_radical_electrons)
    cpp_result = cpp_xyz2omol(xyz_block, total_charge, total_radical_electrons)

    assert (py_result is None) == (cpp_result is None)
    if py_result is None or cpp_result is None:
        return

    py_mol_data = _core.utils.extract_molecule_data(_get_ptr(py_result.OBMol))
    py_rdmol = mol_data_to_rdkit(py_mol_data)
    cpp_rdmol = mol_data_to_rdkit(cpp_result)

    equivalent, info = check_equivalence(py_rdmol, cpp_rdmol)
    assert equivalent, info.reason


def test_cpp_xyz2omol_cache_clear() -> None:
    xyz_block = """2
H2
H 0.0 0.0 0.0
H 0.0 0.0 0.74
"""

    first = cpp_xyz2omol(xyz_block, 0, 0)
    second = cpp_xyz2omol(xyz_block, 0, 0)
    assert first is not None
    assert first is second

    cpp_xyz2omol.cache_clear()  # pyright: ignore[reportFunctionMemberAccess]

    third = cpp_xyz2omol(xyz_block, 0, 0)
    assert third is not None
    assert third is not second


@pytest.mark.parametrize(
    ("metal", "valence"),
    [
        ("Fe", 2),
        ("Cu", 2),
        ("Li", 10),
        ("UUnobtainium", 1),
    ],
)
def test_get_possible_metal_radicals_cpp_matches_fallback(metal: str, valence: int) -> None:
    py_result = py_get_possible_metal_radicals(metal, valence)
    cpp_result = _with_metals.get_possible_metal_radicals(
        metal,
        valence,
    )
    assert py_result == cpp_result


def test_build_metal_states_ptr_cpp_matches_fallback() -> None:
    metal_seed = pybel.readstring(
        "xyz",
        """3
LiCO
Li 0.0 0.0 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    metal_atoms = [atom for atom in metal_seed.atoms if atom.OBAtom.IsMetal()]
    assert len(metal_atoms) == 1

    py_states = py_build_metal_states(metal_atoms[0].OBAtom)
    cpp_states = _with_metals.build_metal_states_ptr(
        _get_ptr(metal_seed.OBMol),
        metal_atoms[0].idx,
    )

    assert [
        (
            state.idx,
            state.symbol,
            state.element_idx,
            state.valence,
            state.radical_num,
            state.position_x,
            state.position_y,
            state.position_z,
        )
        for state in py_states
    ] == [
        (
            state.idx,
            state.symbol,
            state.element_idx,
            state.valence,
            state.radical_num,
            state.position_x,
            state.position_y,
            state.position_z,
        )
        for state in cpp_states
    ]


def test_combine_metal_with_omol_ptr_cpp_matches_fallback() -> None:
    metal_seed_py = pybel.readstring(
        "xyz",
        """3
LiCO
Li 0.0 0.0 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    metal_seed_cpp = pybel.readstring(
        "xyz",
        """3
LiCO
Li 0.0 0.0 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    py_state = py_build_metal_states(metal_seed_py.OBMol.GetAtom(1))[0]
    cpp_state = _with_metals.build_metal_states_ptr(
        _get_ptr(metal_seed_cpp.OBMol),
        1,
    )[0]

    organic_py = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    organic_cpp = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )

    py_combined = py_combine_metal_with_omol(organic_py, [py_state])
    _with_metals.combine_metal_with_omol_ptr(
        _get_ptr(organic_cpp.OBMol),
        [cpp_state],
    )

    assert [atom.atomicnum for atom in py_combined.atoms] == [
        atom.atomicnum for atom in organic_cpp.atoms
    ]
