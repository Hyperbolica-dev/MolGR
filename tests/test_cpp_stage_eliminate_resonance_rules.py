# pyright: reportMissingImports=false

from typing import Any, Iterable, Tuple

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr import _core  # type: ignore
from molgr.fallback.stages.eliminate import (
    eliminate_1_3_dipole,
    eliminate_negative_charges,
    eliminate_positive_charges,
)


_stages: Any = _core.stages


def _get_ptr(obmol) -> int:
    this = getattr(obmol, "this", None)
    if this is not None:
        return int(this)  # type: ignore[arg-type]
    return int(obmol)  # type: ignore[arg-type]


def _smiles_token(mol: pybel.Molecule) -> str:
    smi = mol.write("smi")
    assert smi is not None
    return smi.split()[0]


def _assert_atom_state_parity(py_mol: pybel.Molecule, cpp_mol: pybel.Molecule) -> None:
    assert py_mol.OBMol.NumAtoms() == cpp_mol.OBMol.NumAtoms()
    for idx in range(1, py_mol.OBMol.NumAtoms() + 1):
        py_atom = py_mol.OBMol.GetAtom(idx)
        cpp_atom = cpp_mol.OBMol.GetAtom(idx)
        assert py_atom.GetFormalCharge() == cpp_atom.GetFormalCharge()
        assert py_atom.GetSpinMultiplicity() == cpp_atom.GetSpinMultiplicity()


def _assert_key_bond_orders(
    py_mol: pybel.Molecule, cpp_mol: pybel.Molecule, bonds: Iterable[Tuple[int, int]]
) -> None:
    for begin_idx, end_idx in bonds:
        py_bond = py_mol.OBMol.GetBond(begin_idx, end_idx)
        cpp_bond = cpp_mol.OBMol.GetBond(begin_idx, end_idx)
        assert py_bond is not None and cpp_bond is not None
        assert py_bond.GetBondOrder() == cpp_bond.GetBondOrder()


def _build_1_3_dipole_seed() -> ob.OBMol:
    obmol = ob.OBMol()
    obmol.BeginModify()

    atom1 = obmol.NewAtom()
    atom1.SetAtomicNum(6)
    atom1.SetFormalCharge(-1)
    atom1.SetSpinMultiplicity(0)

    atom2 = obmol.NewAtom()
    atom2.SetAtomicNum(7)
    atom2.SetFormalCharge(0)
    atom2.SetSpinMultiplicity(0)

    atom3 = obmol.NewAtom()
    atom3.SetAtomicNum(6)
    atom3.SetFormalCharge(0)
    atom3.SetSpinMultiplicity(1)

    obmol.AddBond(1, 2, 1)
    obmol.AddBond(2, 3, 2)
    obmol.EndModify()
    return obmol


def _build_positive_charges_seed() -> ob.OBMol:
    obmol = ob.OBMol()
    obmol.BeginModify()

    n1 = obmol.NewAtom()
    n1.SetAtomicNum(7)
    n1.SetFormalCharge(0)
    n1.SetSpinMultiplicity(0)

    n2 = obmol.NewAtom()
    n2.SetAtomicNum(7)
    n2.SetFormalCharge(0)
    n2.SetSpinMultiplicity(1)

    h = obmol.NewAtom()
    h.SetAtomicNum(1)
    h.SetFormalCharge(0)
    h.SetSpinMultiplicity(0)

    obmol.AddBond(1, 2, 2)
    obmol.AddBond(1, 3, 1)
    obmol.EndModify()
    return obmol


def _build_negative_charges_seed() -> ob.OBMol:
    obmol = ob.OBMol()
    obmol.BeginModify()

    oxygen = obmol.NewAtom()
    oxygen.SetAtomicNum(8)
    oxygen.SetFormalCharge(0)
    oxygen.SetSpinMultiplicity(1)

    carbon = obmol.NewAtom()
    carbon.SetAtomicNum(6)
    carbon.SetFormalCharge(0)
    carbon.SetSpinMultiplicity(1)

    h1 = obmol.NewAtom()
    h1.SetAtomicNum(1)
    h2 = obmol.NewAtom()
    h2.SetAtomicNum(1)
    h3 = obmol.NewAtom()
    h3.SetAtomicNum(1)

    obmol.AddBond(2, 3, 1)
    obmol.AddBond(2, 4, 1)
    obmol.AddBond(2, 5, 1)
    obmol.EndModify()
    return obmol


def test_eliminate_1_3_dipole_cpp_matches_python() -> None:
    py_mol = pybel.Molecule(_build_1_3_dipole_seed())
    cpp_mol = pybel.Molecule(_build_1_3_dipole_seed())

    py_mol, py_charge = eliminate_1_3_dipole(py_mol, 0)
    cpp_charge = _stages.eliminate.eliminate_1_3_dipole_ptr(_get_ptr(cpp_mol.OBMol), 0)

    assert py_charge == cpp_charge
    _assert_atom_state_parity(py_mol, cpp_mol)
    _assert_key_bond_orders(py_mol, cpp_mol, [(1, 2), (2, 3)])
    assert _smiles_token(py_mol) == _smiles_token(cpp_mol)


def test_eliminate_positive_charges_cpp_matches_python() -> None:
    py_mol = pybel.Molecule(_build_positive_charges_seed())
    cpp_mol = pybel.Molecule(_build_positive_charges_seed())

    py_mol, py_charge = eliminate_positive_charges(py_mol, 1)
    cpp_charge = _stages.eliminate.eliminate_positive_charges_ptr(_get_ptr(cpp_mol.OBMol), 1)

    assert py_charge == cpp_charge
    _assert_atom_state_parity(py_mol, cpp_mol)
    _assert_key_bond_orders(py_mol, cpp_mol, [(1, 2), (1, 3)])
    assert _smiles_token(py_mol) == _smiles_token(cpp_mol)


def test_eliminate_negative_charges_cpp_matches_python() -> None:
    py_mol = pybel.Molecule(_build_negative_charges_seed())
    cpp_mol = pybel.Molecule(_build_negative_charges_seed())

    py_mol, py_charge = eliminate_negative_charges(py_mol, -2)
    cpp_charge = _stages.eliminate.eliminate_negative_charges_ptr(_get_ptr(cpp_mol.OBMol), -2)

    assert py_charge == cpp_charge
    _assert_atom_state_parity(py_mol, cpp_mol)
    _assert_key_bond_orders(py_mol, cpp_mol, [(2, 3), (2, 4), (2, 5)])
    assert _smiles_token(py_mol) == _smiles_token(cpp_mol)
