from __future__ import annotations

# pyright: reportMissingImports=false
from typing import cast

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.utils.metal_radical_inference import (
    infer_metal_radical_counts,
    infer_metal_radical_state,
)
from molgr.fallback.utils.metals.preparation import _build_metal_states


def test_infer_metal_radical_state_square_planar_co_i_keeps_closed_shell() -> None:
    xyz = """5
CoN4
Co 0.0000 0.0000 0.0000
N  1.9000 0.0000 0.0000
N -1.9000 0.0000 0.0000
N  0.0000 1.9000 0.0000
N  0.0000 -1.9000 0.0000
"""
    mol = pybel.readstring("xyz", xyz)
    metal = cast(ob.OBAtom, mol.OBMol.GetAtom(1))

    result = infer_metal_radical_state(metal, 1)

    assert result.geometry == "square_planar"
    assert result.field_strength == "strong"
    assert result.effective_d_electrons == 8
    assert result.radical_counts == (0,)


def test_infer_metal_radical_state_tetrahedral_co_i_prefers_open_shell() -> None:
    xyz = """5
CoCl4
Co 0.0000 0.0000 0.0000
Cl 1.1500 1.1500 1.1500
Cl -1.1500 -1.1500 1.1500
Cl -1.1500 1.1500 -1.1500
Cl 1.1500 -1.1500 -1.1500
"""
    mol = pybel.readstring("xyz", xyz)
    metal = cast(ob.OBAtom, mol.OBMol.GetAtom(1))

    result = infer_metal_radical_state(metal, 1)

    assert result.geometry == "tetrahedral"
    assert result.field_strength == "weak"
    assert result.effective_d_electrons == 8
    assert result.radical_counts == (2,)


def test_infer_metal_radical_state_octahedral_pt_iv_keeps_strong_and_weak_field_states() -> None:
    xyz = """7
PtO2Cl4
Pt 0.0000 0.0000 0.0000
O  2.1000 0.0000 0.0000
O -2.1000 0.0000 0.0000
Cl 0.0000 2.3000 0.0000
Cl 0.0000 -2.3000 0.0000
Cl 0.0000 0.0000 2.3500
Cl 0.0000 0.0000 -2.3500
"""
    mol = pybel.readstring("xyz", xyz)
    metal = cast(ob.OBAtom, mol.OBMol.GetAtom(1))

    result = infer_metal_radical_state(metal, 4)

    assert result.geometry == "octahedral_like"
    assert result.field_strength == "weak"
    assert result.effective_d_electrons == 6
    assert result.radical_counts == (0, 4)


def test_build_metal_states_uses_environment_sensitive_radical_inference() -> None:
    xyz = """5
CoN4
Co 0.0000 0.0000 0.0000
N  1.9000 0.0000 0.0000
N -1.9000 0.0000 0.0000
N  0.0000 1.9000 0.0000
N  0.0000 -1.9000 0.0000
"""
    mol = pybel.readstring("xyz", xyz)
    metal = cast(ob.OBAtom, mol.OBMol.GetAtom(1))

    states = _build_metal_states(metal)
    co_i_radicals = {state.radical_num for state in states if state.valence == 1}

    assert co_i_radicals == {0}
    assert infer_metal_radical_counts(metal, 1) == (0,)
