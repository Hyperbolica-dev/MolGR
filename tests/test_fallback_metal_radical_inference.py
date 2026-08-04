from __future__ import annotations

# pyright: reportMissingImports=false
from typing import cast

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.utils.consts import NON_METAL_DICT
from molgr.fallback.utils.metal_radical_inference import (
    _DONOR_FIELD_STRENGTH,
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


def test_infer_metal_radical_state_tetrahedral_co_i_keeps_both_spin_branches() -> None:
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
    assert result.radical_counts == (2, 0)


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
    assert result.field_strength == "ambiguous"
    assert result.effective_d_electrons == 6
    assert result.radical_counts == (4, 0)


def test_weak_field_octahedral_pt_iv_keeps_low_spin_fallback_for_python_and_cpp() -> None:
    xyz = """7
PtBr6
Pt 0.0000 0.0000 0.0000
Br 2.3000 0.0000 0.0000
Br -2.3000 0.0000 0.0000
Br 0.0000 2.3000 0.0000
Br 0.0000 -2.3000 0.0000
Br 0.0000 0.0000 2.3000
Br 0.0000 0.0000 -2.3000
"""
    mol = pybel.readstring("xyz", xyz)
    metal = cast(ob.OBAtom, mol.OBMol.GetAtom(1))

    result = infer_metal_radical_state(metal, 4)

    assert result.geometry == "octahedral_like"
    assert result.field_strength == "weak"
    assert result.radical_counts == (4, 0)

    from molgr import _core

    cpp_states = _core.dev.pipeline.reconstruct_with_metals.build_metal_states_ptr(
        int(mol.OBMol.this),
        1,
    )
    assert {state.radical_num for state in cpp_states if state.valence == 4} == {0, 4}


def test_tetrahedral_field_score_in_margin_keeps_both_spin_branches() -> None:
    xyz = """5
NiP2Cl2
Ni 0.0000 0.0000 0.0000
P  1.30 1.30 1.30
P -1.30 -1.30 1.30
Cl -1.30 1.30 -1.30
Cl 1.30 -1.30 -1.30
"""
    mol = pybel.readstring("xyz", xyz)
    metal = cast(ob.OBAtom, mol.OBMol.GetAtom(1))
    result = infer_metal_radical_state(metal, 2)
    assert result.geometry == "tetrahedral"
    assert result.field_strength == "ambiguous"
    assert result.radical_counts == (2, 0)


def test_all_metals_keep_low_spin_branch_in_weak_tetrahedral_field() -> None:
    def tetrahedral_chloride(symbol: str) -> pybel.Molecule:
        return pybel.readstring(
            "xyz",
            f"""5
{symbol}Cl4
{symbol} 0.0000 0.0000 0.0000
Cl 1.1500 1.1500 1.1500
Cl -1.1500 -1.1500 1.1500
Cl -1.1500 1.1500 -1.1500
Cl 1.1500 -1.1500 -1.1500
""",
        )

    nickel = tetrahedral_chloride("Ni")
    copper = tetrahedral_chloride("Cu")
    nickel_result = infer_metal_radical_state(cast(ob.OBAtom, nickel.OBMol.GetAtom(1)), 2)
    copper_result = infer_metal_radical_state(cast(ob.OBAtom, copper.OBMol.GetAtom(1)), 3)

    assert nickel_result.geometry == copper_result.geometry == "tetrahedral"
    assert nickel_result.field_strength == copper_result.field_strength == "weak"
    assert nickel_result.radical_counts == (2, 0)
    assert copper_result.radical_counts == (2, 0)

    from molgr import _core

    cpp_states = _core.dev.pipeline.reconstruct_with_metals.build_metal_states_ptr(
        int(copper.OBMol.this),
        1,
    )
    assert {state.radical_num for state in cpp_states if state.valence == 3} == {0, 2}


def test_donor_field_strength_covers_every_supported_non_metal() -> None:
    assert set(_DONOR_FIELD_STRENGTH) == set(NON_METAL_DICT)


def test_tetrahedral_selenium_coordination_is_strong_field_for_python_and_cpp() -> None:
    mol = pybel.readstring(
        "xyz",
        """5
CuSe4
Cu 0.0000 0.0000 0.0000
Se 1.3000 1.3000 1.3000
Se -1.3000 -1.3000 1.3000
Se -1.3000 1.3000 -1.3000
Se 1.3000 -1.3000 -1.3000
""",
    )
    metal = cast(ob.OBAtom, mol.OBMol.GetAtom(1))
    result = infer_metal_radical_state(metal, 3)

    assert result.geometry == "tetrahedral"
    assert result.field_strength == "strong"
    assert result.radical_counts == (0,)

    from molgr import _core

    cpp_states = _core.dev.pipeline.reconstruct_with_metals.build_metal_states_ptr(
        int(mol.OBMol.this),
        1,
    )
    assert {state.radical_num for state in cpp_states if state.valence == 3} == {0}


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
