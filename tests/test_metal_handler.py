# pyright: reportMissingImports=false

import pytest
from openbabel import pybel

from molgr import _core  # type: ignore


_dev_with_metals = _core.dev.pipeline.reconstruct_with_metals


def get_ptr(mol):
    if isinstance(mol, pybel.Molecule):
        mol = mol.OBMol
    if hasattr(mol, "this"):
        return int(mol.this)
    return int(mol)


def test_build_metal_states_from_dev_api():
    xyz = """1
Pd
Pd 0.0 0.0 0.0
"""
    mol = pybel.readstring("xyz", xyz)

    states = _dev_with_metals.build_metal_states_ptr(get_ptr(mol), 1)

    assert len(states) > 0
    assert all(state.idx == 1 for state in states)
    assert all(state.symbol == "Pd" for state in states)
    assert all(isinstance(state.valence, int) for state in states)
    assert all(isinstance(state.radical_num, int) for state in states)


def test_combine_and_renumber_with_dev_reinsertion():
    xyz = """3
Renumber Test
C  1.0 0.0 0.0
Li 2.0 0.0 0.0
O  3.0 0.0 0.0
"""
    mol = pybel.readstring("xyz", xyz)
    ptr = get_ptr(mol)

    original_c_x = mol.atoms[0].coords[0]
    original_li_x = mol.atoms[1].coords[0]
    original_o_x = mol.atoms[2].coords[0]

    li_atom = mol.OBMol.GetAtom(2)
    metal = _dev_with_metals.MetalAtomPosition()
    metal.idx = 2
    metal.symbol = "Li"
    metal.element_idx = li_atom.GetAtomicNum()
    metal.valence = 1
    metal.radical_num = 0
    metal.position_x = li_atom.GetX()
    metal.position_y = li_atom.GetY()
    metal.position_z = li_atom.GetZ()

    mol.OBMol.BeginModify()
    mol.OBMol.DeleteAtom(li_atom)
    mol.OBMol.EndModify()
    assert mol.OBMol.NumAtoms() == 2

    _dev_with_metals.combine_metal_with_omol_ptr(ptr, [metal])

    assert mol.OBMol.NumAtoms() == 3

    atom1 = mol.atoms[0]
    atom2 = mol.atoms[1]
    atom3 = mol.atoms[2]

    assert atom1.atomicnum == 6
    assert atom2.atomicnum == 3
    assert atom3.atomicnum == 8

    assert atom1.coords[0] == pytest.approx(original_c_x)
    assert atom2.coords[0] == pytest.approx(original_li_x)
    assert atom3.coords[0] == pytest.approx(original_o_x)
