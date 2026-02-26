# pyright: reportMissingImports=false

from typing import Any

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr import _core  # type: ignore
from molgr.fallback.stages.fresh import (
    assign_charge_radical_for_atom,
    assign_radical_dots,
    fresh_omol_charge_radical,
)


_stages: Any = _core.stages


def _get_ptr(obmol) -> int:
    this = getattr(obmol, "this", None)
    if this is not None:
        return int(this)  # type: ignore[arg-type]
    return int(obmol)  # type: ignore[arg-type]


def _build_seed_obmol() -> ob.OBMol:
    obmol = ob.OBMol()
    obmol.BeginModify()

    b = obmol.NewAtom()
    b.SetAtomicNum(5)
    b.SetFormalCharge(0)
    b.SetSpinMultiplicity(0)

    h1 = obmol.NewAtom()
    h1.SetAtomicNum(1)
    h2 = obmol.NewAtom()
    h2.SetAtomicNum(1)
    h3 = obmol.NewAtom()
    h3.SetAtomicNum(1)
    h4 = obmol.NewAtom()
    h4.SetAtomicNum(1)

    ne = obmol.NewAtom()
    ne.SetAtomicNum(10)
    ne.SetFormalCharge(7)
    ne.SetSpinMultiplicity(0)

    o = obmol.NewAtom()
    o.SetAtomicNum(8)
    o.SetFormalCharge(1)
    o.SetSpinMultiplicity(1)

    c = obmol.NewAtom()
    c.SetAtomicNum(6)
    c.SetFormalCharge(0)
    c.SetSpinMultiplicity(0)

    obmol.AddBond(1, 2, 1)
    obmol.AddBond(1, 3, 1)
    obmol.AddBond(1, 4, 1)
    obmol.AddBond(1, 5, 1)
    obmol.AddBond(8, 7, 2)

    obmol.EndModify()
    return obmol


def test_fresh_stage_cpp_matches_python_per_atom_and_totals() -> None:
    py_mol = pybel.Molecule(_build_seed_obmol())
    cpp_mol = pybel.Molecule(_build_seed_obmol())

    fresh_omol_charge_radical(py_mol)
    _stages.fresh.fresh_omol_charge_radical_ptr(_get_ptr(cpp_mol.OBMol))

    py_charge_total = 0
    cpp_charge_total = 0
    py_rad_total = 0
    cpp_rad_total = 0

    assert py_mol.OBMol.NumAtoms() == cpp_mol.OBMol.NumAtoms()
    for idx in range(1, py_mol.OBMol.NumAtoms() + 1):
        py_atom = py_mol.OBMol.GetAtom(idx)
        cpp_atom = cpp_mol.OBMol.GetAtom(idx)

        assert py_atom.GetAtomicNum() == cpp_atom.GetAtomicNum()
        assert py_atom.GetFormalCharge() == cpp_atom.GetFormalCharge()
        assert py_atom.GetSpinMultiplicity() == cpp_atom.GetSpinMultiplicity()

        py_charge_total += py_atom.GetFormalCharge()
        cpp_charge_total += cpp_atom.GetFormalCharge()
        py_rad_total += py_atom.GetSpinMultiplicity()
        cpp_rad_total += cpp_atom.GetSpinMultiplicity()

    assert py_charge_total == cpp_charge_total
    assert py_rad_total == cpp_rad_total


def test_assign_radical_dots_ptr_matches_fallback_per_atom() -> None:
    py_mol = pybel.Molecule(_build_seed_obmol())
    cpp_mol = pybel.Molecule(_build_seed_obmol())

    assert py_mol.OBMol.NumAtoms() == cpp_mol.OBMol.NumAtoms()
    for idx in range(1, py_mol.OBMol.NumAtoms() + 1):
        py_atom = py_mol.OBMol.GetAtom(idx)
        py_radical_dots = assign_radical_dots(py_atom)
        cpp_radical_dots = _stages.fresh.assign_radical_dots_ptr(_get_ptr(cpp_mol.OBMol), idx)
        assert py_radical_dots == cpp_radical_dots


def test_assign_charge_radical_for_atom_ptr_matches_fallback_per_atom() -> None:
    py_mol = pybel.Molecule(_build_seed_obmol())
    cpp_mol = pybel.Molecule(_build_seed_obmol())

    assert py_mol.OBMol.NumAtoms() == cpp_mol.OBMol.NumAtoms()
    for idx in range(1, py_mol.OBMol.NumAtoms() + 1):
        py_atom = py_mol.OBMol.GetAtom(idx)
        assign_charge_radical_for_atom(py_atom)

        _stages.fresh.assign_charge_radical_for_atom_ptr(_get_ptr(cpp_mol.OBMol), idx)
        cpp_atom = cpp_mol.OBMol.GetAtom(idx)

        assert py_atom.GetFormalCharge() == cpp_atom.GetFormalCharge()
        assert py_atom.GetSpinMultiplicity() == cpp_atom.GetSpinMultiplicity()
