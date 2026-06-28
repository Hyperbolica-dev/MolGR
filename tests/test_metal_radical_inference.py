# pyright: reportMissingImports=false

from __future__ import annotations

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.utils.metal_radical_inference import infer_metal_radical_state
from molgr.fallback.utils.metals.preparation import prepare_metal_state


_RU_H2_P4_LOCAL_XYZ = """8
ADIYET RuH2P4 local environment
Ru         0.03520       -0.02870        0.47350
H          0.15310       -1.23300        1.59540
H          0.00770        0.98500        1.77670
P         -2.14440       -0.22840        1.01360
P          2.24650        0.14520        0.88020
P          0.12810       -1.64740       -1.13160
P         -0.23420        1.71230       -0.98260
H         -0.88470        2.79020        1.51270
"""


def _metal_atom_from_mol(mol: pybel.Molecule):
    for atom in ob.OBMolAtomIter(mol.OBMol):
        if atom.IsMetal():
            return atom
    raise AssertionError("metal atom not found")


def test_metal_radical_inference_counts_direct_hydrides_but_not_outer_hydrogen() -> None:
    mol = pybel.readstring("xyz", _RU_H2_P4_LOCAL_XYZ)
    metal_atom = _metal_atom_from_mol(mol)

    state = infer_metal_radical_state(metal_atom, 2)

    assert state.coordination_number == 6
    assert state.geometry == "octahedral_like"
    assert 0 in state.radical_counts


def test_prepare_metal_state_includes_low_spin_ru_ii_when_hydrides_are_direct() -> None:
    state = prepare_metal_state(_RU_H2_P4_LOCAL_XYZ, total_charge=0, total_radical_electrons=0)
    metal_states = state.available_valence_radical_states[0]

    assert ("Ru", 2, 0) in {
        (metal_state.symbol, metal_state.valence, metal_state.radical_num)
        for metal_state in metal_states
    }
