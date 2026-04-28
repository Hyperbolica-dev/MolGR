"""
Author: TMJ
Date: 2026-02-27 14:45:40
LastEditors: TMJ
LastEditTime: 2026-02-27 22:34:34
Description: 请填写简介
"""
# pyright: reportMissingImports=false

import sys
from pathlib import Path

import pytest

from molgr.utils.converter import mol_data_to_rdkit


pytest.importorskip("rdkit")
pytest.importorskip("openbabel")


def test_mol_data_to_rdkit_sets_positions_by_atom_index() -> None:
    from openbabel import pybel

    from molgr import _core  # type: ignore
    from molgr.utils.converter import mol_data_to_rdkit

    xyz_block = """3
pos
C 0.0 0.0 0.0
O 1.5 -2.0 3.0
N -4.0 5.0 -6.0
"""

    omol = pybel.readstring("xyz", xyz_block)

    def _get_ptr(obmol) -> int:
        this = getattr(obmol, "this", None)
        if this is not None:
            return int(this)  # type: ignore[arg-type]
        return int(obmol)  # type: ignore[arg-type]

    mol_ptr = _get_ptr(omol.OBMol)
    md = _core.utils.extract_molecule_data(mol_ptr)

    mol = mol_data_to_rdkit(md, sanitize=True)
    assert mol.GetNumAtoms() == 3
    conf = mol.GetConformer()

    p0 = conf.GetAtomPosition(0)
    p1 = conf.GetAtomPosition(1)
    p2 = conf.GetAtomPosition(2)

    assert p0.x == pytest.approx(0.0)
    assert p0.y == pytest.approx(0.0)
    assert p0.z == pytest.approx(0.0)

    assert p1.x == pytest.approx(1.5)
    assert p1.y == pytest.approx(-2.0)
    assert p1.z == pytest.approx(3.0)

    assert p2.x == pytest.approx(-4.0)
    assert p2.y == pytest.approx(5.0)
    assert p2.z == pytest.approx(-6.0)


@pytest.mark.parametrize("case_idx", [1])
def test_mol_data_to_rdkit_matches_fallback_pybel_conversion_for_hard_cases(case_idx: int) -> None:
    from rdkit import Chem

    from molgr import _core  # type: ignore
    from molgr.fallback import xyz2omol
    from molgr.utils.converter import pybel_to_rdmol

    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    from scripts.molgr_cases_smiles_csv import load_smiles_csv_cases

    cases = load_smiles_csv_cases(Path("tests/test_cases.csv"))
    case = next(row for row in cases if int(row["case_idx"]) == case_idx)

    xyz_block = str(case["xyz_block"])
    total_charge = int(case["total_charge"])
    total_radical_electrons = int(case["total_radical_electrons"])

    fallback_omol = xyz2omol(
        xyz_block,
        total_charge=total_charge,
        total_radical_electrons=total_radical_electrons,
    )
    assert fallback_omol is not None

    cpp_mol_data = _core.pipeline.reconstruct_with_metals.xyz2omol(
        xyz_block,
        total_charge,
        total_radical_electrons,
    )
    assert cpp_mol_data is not None

    fallback_rdmol = Chem.RemoveHs(pybel_to_rdmol(fallback_omol))
    cpp_rdmol = Chem.RemoveHs(mol_data_to_rdkit(cpp_mol_data))

    fallback_smiles = Chem.MolToSmiles(fallback_rdmol, canonical=True, isomericSmiles=True)
    cpp_smiles = Chem.MolToSmiles(cpp_rdmol, canonical=True, isomericSmiles=True)
    assert cpp_smiles == fallback_smiles
