"""
Author: TMJ
Date: 2026-02-25 00:41:17
LastEditors: TMJ
LastEditTime: 2026-02-25 13:45:25
Description: 请填写简介
"""
# pyright: reportMissingImports=false

from typing import Tuple

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.stages.preprocess import validate_omol
from molgr.fallback.utils.electrons import (
    get_lone_pair_count,
    get_unpaired_electron_count,
    has_unresolved_two_electron_center,
    set_lone_pair_count,
    set_unpaired_electron_count,
    set_unresolved_two_electron_center,
)


def test_validate_omol_distinguishes_lone_pair_from_two_unpaired_electrons() -> None:
    obmol = ob.OBMol()
    a = obmol.NewAtom()
    a.SetAtomicNum(6)
    a.SetFormalCharge(0)

    set_lone_pair_count(a, 1)
    mol = pybel.Molecule(obmol)
    assert validate_omol(mol, 0, 0) is True

    set_lone_pair_count(a, 0)
    set_unpaired_electron_count(a, 2)
    mol = pybel.Molecule(obmol)
    assert validate_omol(mol, 0, 0) is False
    assert validate_omol(mol, 0, 2) is True

    set_unpaired_electron_count(a, 1)
    mol = pybel.Molecule(obmol)
    assert validate_omol(mol, 0, 1) is True


@pytest.mark.parametrize(
    ("total_radical_electrons", "expected_occupancy"),
    [(0, (0, 1)), (2, (2, 0))],
)
def test_validate_omol_resolves_two_electron_center_for_python_and_cpp(
    total_radical_electrons: int,
    expected_occupancy: Tuple[int, int],
) -> None:
    def make_molecule() -> pybel.Molecule:
        obmol = ob.OBMol()
        atom = obmol.NewAtom()
        atom.SetAtomicNum(6)
        set_unresolved_two_electron_center(atom, True)
        return pybel.Molecule(obmol)

    mol = make_molecule()
    assert validate_omol(mol, 0, total_radical_electrons) is True
    atom = mol.OBMol.GetAtom(1)
    assert (get_unpaired_electron_count(atom), get_lone_pair_count(atom)) == expected_occupancy
    assert not has_unresolved_two_electron_center(atom)

    from molgr import _core  # type: ignore

    cpp_mol = make_molecule()
    mol_ptr = int(getattr(cpp_mol.OBMol, "this", cpp_mol.OBMol))
    assert (
        _core.dev.stages.preprocess.validate_omol_ptr(
            mol_ptr,
            0,
            total_radical_electrons,
        )
        is True
    )
    atom = cpp_mol.OBMol.GetAtom(1)
    assert (get_unpaired_electron_count(atom), get_lone_pair_count(atom)) == expected_occupancy
    assert not has_unresolved_two_electron_center(atom)


def test_validate_omol_leaves_unresolved_center_unchanged_when_budget_is_impossible() -> None:
    obmol = ob.OBMol()
    atom = obmol.NewAtom()
    atom.SetAtomicNum(6)
    set_unresolved_two_electron_center(atom, True)
    mol = pybel.Molecule(obmol)

    assert validate_omol(mol, 0, 1) is False
    assert has_unresolved_two_electron_center(atom)
    assert (get_unpaired_electron_count(atom), get_lone_pair_count(atom)) == (0, 0)

    from molgr import _core  # type: ignore

    mol_ptr = int(getattr(mol.OBMol, "this", mol.OBMol))
    assert _core.dev.stages.preprocess.validate_omol_ptr(mol_ptr, 0, 1) is False
    assert has_unresolved_two_electron_center(atom)
    assert (get_unpaired_electron_count(atom), get_lone_pair_count(atom)) == (0, 0)
