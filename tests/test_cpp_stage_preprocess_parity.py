# pyright: reportMissingImports=false

from typing import Any, Tuple

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr import _core  # type: ignore
from molgr.fallback.stages.preprocess import make_connections, pre_clean, validate_omol


_stages: Any = _core.stages


def _get_ptr(obmol) -> int:
    this = getattr(obmol, "this", None)
    if this is not None:
        return int(this)  # type: ignore[arg-type]
    return int(obmol)  # type: ignore[arg-type]


def _total_charge(mol: pybel.Molecule) -> int:
    return sum(atom.OBAtom.GetFormalCharge() for atom in mol.atoms)


def _total_radicals(mol: pybel.Molecule) -> Tuple[int, int]:
    raw = sum(atom.OBAtom.GetSpinMultiplicity() for atom in mol.atoms)
    singlet_compatible = sum(atom.OBAtom.GetSpinMultiplicity() % 2 for atom in mol.atoms)
    return raw, singlet_compatible


def _smiles_token(mol: pybel.Molecule) -> str:
    smi = mol.write("smi")
    assert smi is not None
    return smi.split()[0]


def _assert_stage_parity(py_mol: pybel.Molecule, cpp_mol: pybel.Molecule) -> None:
    assert _total_charge(py_mol) == _total_charge(cpp_mol)

    py_raw, py_singlet = _total_radicals(py_mol)
    cpp_raw, cpp_singlet = _total_radicals(cpp_mol)
    assert py_raw == cpp_raw and py_singlet == cpp_singlet

    assert _smiles_token(py_mol) == _smiles_token(cpp_mol)


def _clear_all_bonds(obmol: ob.OBMol) -> None:
    obmol.BeginModify()
    bonds = list(ob.OBMolBondIter(obmol))
    for bond in bonds:
        obmol.DeleteBond(bond)
    obmol.EndModify()


def test_make_connections_stage_parity_from_xyz() -> None:
    xyz_block = """2
no_bond
N 0.000 0.000 0.000
O 1.200 0.000 0.000
"""
    py_mol = pybel.readstring("xyz", xyz_block)
    cpp_mol = pybel.readstring("xyz", xyz_block)

    _clear_all_bonds(py_mol.OBMol)
    _clear_all_bonds(cpp_mol.OBMol)

    make_connections(py_mol, factor=1.4)
    _stages.preprocess.make_connections_ptr(_get_ptr(cpp_mol.OBMol), factor=1.4)

    _assert_stage_parity(py_mol, cpp_mol)


def test_pre_clean_stage_parity_from_xyz() -> None:
    xyz_block = """6
sif5_setup
Si 0.000 0.000 0.000
F 1.600 0.000 0.000
F -1.600 0.000 0.000
F 0.000 1.600 0.000
F 0.000 -1.600 0.000
F 0.000 0.000 1.600
"""
    py_mol = pybel.readstring("xyz", xyz_block)
    cpp_mol = pybel.readstring("xyz", xyz_block)

    for mol in (py_mol, cpp_mol):
        _clear_all_bonds(mol.OBMol)
        mol.OBMol.BeginModify()
        for idx in range(2, 7):
            mol.OBMol.AddBond(1, idx, 1)
        mol.OBMol.EndModify()

    pre_clean(py_mol)
    _stages.preprocess.pre_clean_ptr(_get_ptr(cpp_mol.OBMol))

    _assert_stage_parity(py_mol, cpp_mol)


def test_validate_omol_cpp_matches_python_spin_mod2_rule() -> None:
    xyz_block = """1
carbene
C 0.000 0.000 0.000
"""
    py_mol = pybel.readstring("xyz", xyz_block)
    cpp_mol = pybel.readstring("xyz", xyz_block)

    py_atom = py_mol.atoms[0].OBAtom
    cpp_atom = cpp_mol.atoms[0].OBAtom

    py_atom.SetFormalCharge(0)
    cpp_atom.SetFormalCharge(0)

    py_atom.SetSpinMultiplicity(2)
    cpp_atom.SetSpinMultiplicity(2)
    assert validate_omol(py_mol, 0, 0) is True
    assert _stages.preprocess.validate_omol_ptr(_get_ptr(cpp_mol.OBMol), 0, 0) is True

    py_atom.SetSpinMultiplicity(1)
    cpp_atom.SetSpinMultiplicity(1)
    assert validate_omol(py_mol, 0, 1) is True
    assert _stages.preprocess.validate_omol_ptr(_get_ptr(cpp_mol.OBMol), 0, 1) is True
