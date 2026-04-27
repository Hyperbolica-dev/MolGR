from __future__ import annotations

# pyright: reportMissingImports=false
from typing import cast

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.pipeline.reconstruct_with_metals import xyz2omol
from molgr.fallback.utils.metals.preparation import _build_metal_states


def test_build_metal_states_unknown_symbol_defaults() -> None:
    xyz = """1
U
U 0.0 0.0 0.0
"""
    mol = pybel.readstring("xyz", xyz)
    obatom = cast(ob.OBAtom, mol.atoms[0].OBAtom)

    states = _build_metal_states(obatom)

    assert states
    assert any(state.valence == 0 and state.radical_num == 0 for state in states)


def test_xyz2omol_unknown_symbol_does_not_raise() -> None:
    xyz = """1
U
U 0.0 0.0 0.0
"""

    xyz2omol(xyz, total_charge=0, total_radical_electrons=0)
