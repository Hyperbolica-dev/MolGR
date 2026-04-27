# pyright: reportMissingImports=false

import pytest


pytest.importorskip("openbabel")

from openbabel import pybel

from molgr.fallback.utils.dataclasses import MetalAtomPosition
from molgr.fallback.utils.metals.preparation import combine_metal_with_omol


def test_combine_metal_with_omol_inserts_li_into_original_slot() -> None:
    organic_xyz = """2

C 1.0 0.0 0.0
O 3.0 0.0 0.0
"""
    organic = pybel.readstring("xyz", organic_xyz)
    li = MetalAtomPosition(
        idx=2,
        symbol="Li",
        element_idx=3,
        valence=1,
        radical_num=0,
        position_x=2.0,
        position_y=0.0,
        position_z=0.0,
    )

    combined = combine_metal_with_omol(organic, [li])

    assert [atom.atomicnum for atom in combined.atoms] == [6, 3, 8]
