# pyright: reportMissingImports=false

from __future__ import annotations

import pytest


pytest.importorskip("openbabel")

from openbabel import pybel

from molgr.fallback.state import (
    MetalCandidateStateMachine,
    OmolStateMachine,
    ReconstructionState,
    make_metal_candidate_state,
)
from molgr.fallback.utils.dataclasses import MetalAtomPosition
from molgr.fallback.utils.force_field import (
    _build_score_key,
    build_force_field_score_key,
    force_field_evaluation_cache_clear,
    force_field_evaluation_cache_info,
    selection_force_field_energy,
)


def test_selection_force_field_cache_hits_for_equivalent_reloads() -> None:
    force_field_evaluation_cache_clear()
    xyz_block = """3
LiCO
Li 0.0 0.0 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
"""

    first = pybel.readstring("xyz", xyz_block)
    second = pybel.readstring("xyz", xyz_block)

    score_first = selection_force_field_energy(first)
    hits_after_first, misses_after_first, size_after_first = force_field_evaluation_cache_info()
    score_second = selection_force_field_energy(second)
    hits_after_second, misses_after_second, size_after_second = force_field_evaluation_cache_info()

    assert score_first == score_second
    assert (hits_after_first, misses_after_first, size_after_first) == (0, 1, 1)
    assert (hits_after_second, misses_after_second, size_after_second) == (1, 1, 1)


def test_force_field_cache_clear_resets_stats() -> None:
    xyz_block = """2
H2
H 0.0 0.0 0.0
H 0.0 0.0 0.74
"""

    force_field_evaluation_cache_clear()
    selection_force_field_energy(pybel.readstring("xyz", xyz_block))
    assert force_field_evaluation_cache_info() == (0, 1, 1)

    force_field_evaluation_cache_clear()
    assert force_field_evaluation_cache_info() == (0, 0, 0)


def test_score_cache_quantizes_coordinates_at_six_decimal_places() -> None:
    force_field_evaluation_cache_clear()
    first = pybel.readstring(
        "xyz",
        """2
CO
C 2.0000004 0.0 0.0
O 3.2000004 0.0 0.0
""",
    )
    second = pybel.readstring(
        "xyz",
        """2
CO
C 2.00000049 0.0 0.0
O 3.20000049 0.0 0.0
""",
    )

    score_first = selection_force_field_energy(first)
    assert force_field_evaluation_cache_info() == (0, 1, 1)

    score_second = selection_force_field_energy(second)
    assert score_first == score_second
    assert force_field_evaluation_cache_info() == (1, 1, 1)


def test_score_cache_distinguishes_different_geometries_with_same_topology() -> None:
    force_field_evaluation_cache_clear()
    first = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    second = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 4.2 0.0 0.0
""",
    )

    first_key = build_force_field_score_key(first)
    second_key = build_force_field_score_key(second)
    assert first_key != second_key

    selection_force_field_energy(first)
    assert force_field_evaluation_cache_info() == (0, 1, 1)

    selection_force_field_energy(second)
    assert force_field_evaluation_cache_info() == (0, 2, 2)


def test_reconstruction_state_caches_force_field_score_key_until_omol_changes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    force_field_evaluation_cache_clear()
    omol = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    state = ReconstructionState(
        omol=omol,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
    )
    score_key_calls = 0
    original_builder = _build_score_key

    def tracking_builder(current_omol):
        nonlocal score_key_calls
        score_key_calls += 1
        return original_builder(current_omol)

    monkeypatch.setattr("molgr.fallback.utils.force_field._build_score_key", tracking_builder)

    first_organic_score = state.organic_core_score()
    first_organic_key = state.force_field_score_key()
    second_organic_score = state.organic_core_score()
    second_organic_key = state.force_field_score_key()
    assert first_organic_score == second_organic_score
    assert first_organic_key == second_organic_key
    assert score_key_calls == 1

    first_full_score = state.full_score()
    second_full_score = state.full_score()
    assert first_full_score == second_full_score
    assert score_key_calls == 1

    machine = OmolStateMachine.from_reconstruction_state(state)
    machine.run_omol_stage("clone", lambda current_omol: (current_omol.clone, True))
    next_state = machine.freeze_like(state)

    cloned_organic_score = next_state.organic_core_score()
    cloned_organic_key = next_state.force_field_score_key()
    cloned_full_score = next_state.full_score()
    assert cloned_organic_score == first_organic_score
    assert cloned_organic_key == first_organic_key
    assert cloned_full_score == first_full_score
    assert score_key_calls == 2


def test_state_machine_preserves_score_key_cache_when_stage_reports_no_hit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    omol = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    state = ReconstructionState(
        omol=omol,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
    )
    score_key_calls = 0
    original_builder = _build_score_key

    def tracking_builder(current_omol):
        nonlocal score_key_calls
        score_key_calls += 1
        return original_builder(current_omol)

    def no_hit_stage(current_omol):
        return current_omol, False

    monkeypatch.setattr("molgr.fallback.utils.force_field._build_score_key", tracking_builder)

    state.organic_core_score()
    assert score_key_calls == 1

    machine = OmolStateMachine.from_reconstruction_state(state)
    machine.run_omol_stage("noop", no_hit_stage)
    next_state = machine.freeze_like(state)

    next_state.organic_core_score()
    assert score_key_calls == 1


def test_state_machine_invalidates_derived_metadata_on_hit() -> None:
    omol = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    state = ReconstructionState(
        omol=omol,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        metadata={
            "force_field_energy": 1.0,
            "force_field_requested": "auto",
            "force_field_resolved_force_field": "uff",
            "force_field_score_key": "cached-key",
            "organic_core_score": 1.0,
            "score": 1.0,
            "keep_me": "yes",
        },
    )

    machine = OmolStateMachine.from_reconstruction_state(state)
    machine.run_omol_stage("clone", lambda current_omol: (current_omol.clone, True))
    next_state = machine.freeze_like(state)

    assert next_state.omol_revision == 1
    assert "force_field_energy" not in next_state.metadata
    assert "force_field_requested" not in next_state.metadata
    assert "force_field_resolved_force_field" not in next_state.metadata
    assert "force_field_score_key" not in next_state.metadata
    assert "organic_core_score" not in next_state.metadata
    assert "score" not in next_state.metadata
    assert next_state.metadata["keep_me"] == "yes"


def test_metal_candidate_state_caches_combined_score_until_no_metal_state_changes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    no_metal_a = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    no_metal_b = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 3.4 0.0 0.0
""",
    )
    state_a = ReconstructionState(
        omol=no_metal_a,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        metadata={
            "organic_core_score": 1.0,
            "force_field_score_key": build_force_field_score_key(no_metal_a),
        },
    )
    state_b = ReconstructionState(
        omol=no_metal_b,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        metadata={
            "organic_core_score": 2.0,
            "force_field_score_key": build_force_field_score_key(no_metal_b),
        },
    )
    candidate = make_metal_candidate_state(
        (),
        (
            MetalAtomPosition(1, "Li", 3, 1, 0, 0.0, 0.0, 0.0),
            MetalAtomPosition(2, "Li", 3, 1, 0, 4.0, 0.0, 0.0),
        ),
        -2,
        0,
        combination_index=0,
    )
    calls = {"count": 0}

    def tracking_score(self, profile="full"):
        assert profile == "organic_core"
        calls["count"] += 1
        return float(self.metadata["organic_core_score"])

    monkeypatch.setattr(ReconstructionState, "score", tracking_score)

    candidate_machine = MetalCandidateStateMachine.from_candidate_state(candidate)
    candidate_machine.set_no_metal_state("attach_state_a", state_a)
    candidate = candidate_machine.freeze()

    assert candidate.combined_score() == 1.0
    assert candidate.combined_score() == 1.0
    assert calls["count"] == 1

    candidate_machine = MetalCandidateStateMachine.from_candidate_state(candidate)
    candidate_machine.set_no_metal_state("attach_state_b", state_b)
    next_candidate = candidate_machine.freeze()

    assert next_candidate.score is None
    assert next_candidate.combined_omol is None
    assert "score" not in next_candidate.metadata
    assert next_candidate.combined_score() == 2.0
    assert calls["count"] == 2
