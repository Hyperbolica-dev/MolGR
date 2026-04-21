# pyright: reportMissingImports=false

from __future__ import annotations

import csv
from pathlib import Path

import pytest


pytest.importorskip("openbabel")

from openbabel import pybel

from molgr.fallback.pipeline import reconstruct_without_metals as no_metal_module
from molgr.fallback.pipeline.reconstruct_without_metals import (
    xyz_to_omol_no_metal,
    xyz_to_omol_no_metal_state,
)
from molgr.fallback.state import ReconstructionState
from molgr.fallback.utils.tools import typed_lru_cache
from molgr.utils.converter import pybel_to_rdmol
from molgr.utils.equivalence import check_equivalence


@typed_lru_cache(maxsize=128, typed=True)
def _seed_case(smiles: str) -> tuple[str, int, int]:
    seed = pybel.readstring("smi", smiles)
    xyz_block = str(seed.write("xyz"))
    charge = 0
    radicals = 0
    for atom in seed.atoms:
        charge += atom.OBAtom.GetFormalCharge()
        radicals += atom.OBAtom.GetSpinMultiplicity()
    return xyz_block, charge, radicals


def _load_curated_smiles() -> list[str]:
    csv_path = Path(__file__).with_name("test_cases.csv")
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    curated_rows = [1, 2, 5, 10, 17, 26]
    return [rows[idx - 1]["smiles"] for idx in curated_rows]


@pytest.mark.parametrize("smiles", _load_curated_smiles())
def test_fallback_no_metal_reconstructs_curated_cases(smiles: str) -> None:
    xyz_block, total_charge, total_radical_electrons = _seed_case(smiles)

    result = xyz_to_omol_no_metal(xyz_block, total_charge, total_radical_electrons)

    assert result is not None

    expected = pybel_to_rdmol(pybel.readstring("smi", smiles))
    equivalent, info = check_equivalence(expected, pybel_to_rdmol(result))
    assert equivalent, info.reason


def test_fallback_no_metal_exposes_staged_history() -> None:
    xyz_block = """2
H2
H 0.0 0.0 0.0
H 0.0 0.0 0.74
"""
    state = xyz_to_omol_no_metal_state(xyz_block, 0, 0)

    assert state is not None
    assert state.phase_history[:6] == (
        "read_xyz",
        "make_connections",
        "pre_clean",
        "fresh_omol_charge_radical_initial",
        "initialize_charge_budget",
        "eliminate_NNN_negative",
    )
    assert "break_one_bond" in state.phase_history
    assert state.phase_history[-1] in {"clean_resonances", "select_best_resonance_candidate"}


def test_run_linear_pipeline_passes_current_charge_into_break_stages(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    xyz_block = """2
H2
H 0.0 0.0 0.0
H 0.0 0.0 0.74
"""
    state = ReconstructionState(
        omol=pybel.readstring("xyz", xyz_block),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=3,
        phase_history=("read_xyz",),
    )
    recorded: dict[str, tuple[int, int]] = {}

    monkeypatch.setattr(no_metal_module, "make_connections", lambda omol: (omol, False))
    monkeypatch.setattr(no_metal_module, "pre_clean", lambda omol: (omol, False))
    monkeypatch.setattr(no_metal_module, "fresh_omol_charge_radical", lambda omol: (omol, False))
    monkeypatch.setattr(
        no_metal_module,
        "eliminate_NNN",
        lambda omol, given_charge, positive: (omol, 7 if not positive else given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_module,
        "eliminate_high_positive_charge_atoms",
        lambda omol, given_charge: (omol, given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_module,
        "eliminate_CN_in_doubt",
        lambda omol, given_charge: (omol, given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_module,
        "eliminate_carboxyl",
        lambda omol, given_charge: (omol, given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_module,
        "eliminate_carbene_neighbor_heteroatom",
        lambda omol, given_charge: (omol, given_charge, False),
    )
    monkeypatch.setattr(no_metal_module, "clean_carbene_neighbor_unsaturated", lambda omol: (omol, False))
    monkeypatch.setattr(no_metal_module, "clean_neighbor_radicals", lambda omol: (omol, False))
    monkeypatch.setattr(
        no_metal_module,
        "eliminate_charge_spliting",
        lambda omol, given_charge: (omol, given_charge, False),
    )

    def record_break_deformed_ene(omol, given_charge, total_radical_electrons):
        recorded["break_deformed_ene"] = (given_charge, total_radical_electrons)
        return omol, False

    def record_break_one_bond(omol, given_charge, total_radical_electrons):
        recorded["break_one_bond"] = (given_charge, total_radical_electrons)
        return omol, given_charge, False

    monkeypatch.setattr(no_metal_module, "break_deformed_ene", record_break_deformed_ene)
    monkeypatch.setattr(no_metal_module, "break_one_bond", record_break_one_bond)

    next_state = no_metal_module._run_linear_pipeline(state)

    assert recorded["break_deformed_ene"] == (7, 3)
    assert recorded["break_one_bond"] == (7, 3)
    assert next_state.given_charge == 7


def test_fallback_no_metal_reuses_resonance_score_metadata(monkeypatch: pytest.MonkeyPatch) -> None:
    xyz_to_omol_no_metal_state.cache_clear()

    base_state = ReconstructionState(
        omol=object(),  # type: ignore[arg-type]
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz",),
    )
    weaker_candidate = ReconstructionState(
        omol=object(),  # type: ignore[arg-type]
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz", "validate_resonance_candidate"),
        metadata={
            "organic_core_score": 3.0,
            "post_reinsertion_base_key": "weak",
            "score": 5.0,
        },
    )
    stronger_candidate = ReconstructionState(
        omol=object(),  # type: ignore[arg-type]
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz", "validate_resonance_candidate"),
        metadata={
            "organic_core_score": 1.0,
            "post_reinsertion_base_key": "strong",
            "score": 2.0,
        },
    )

    monkeypatch.setattr(no_metal_module, "_seed_state", lambda *args, **kwargs: base_state)
    monkeypatch.setattr(no_metal_module, "_run_linear_pipeline", lambda state: state)
    monkeypatch.setattr(no_metal_module, "validate_omol", lambda *args, **kwargs: False)
    monkeypatch.setattr(
        no_metal_module,
        "_recover_resonance_candidates",
        lambda state, **kwargs: [weaker_candidate, stronger_candidate],
    )

    result = xyz_to_omol_no_metal_state("cached-resonance", 0, 0)

    assert result is not None
    assert result.metadata["score"] == 2.0
    assert result.metadata["organic_core_score"] == 1.0
    assert result.metadata["post_reinsertion_base_key"] == "strong"
    assert result.phase_history[-1] == "select_best_resonance_candidate"
