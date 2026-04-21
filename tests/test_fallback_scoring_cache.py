# pyright: reportMissingImports=false

from __future__ import annotations

import pytest


pytest.importorskip("openbabel")

from openbabel import pybel

from molgr.fallback.pipeline.reconstruct_with_metals import combine_metal_with_omol
from molgr.fallback.state import (
    MetalCandidateStateMachine,
    OmolStateMachine,
    ReconstructionState,
    make_metal_candidate_state,
)
from molgr.fallback.utils import scoring as scoring_module
from molgr.fallback.utils.dataclasses import MetalAtomPosition
from molgr.fallback.utils.scoring import (
    build_post_reinsertion_base_components,
    build_post_reinsertion_base_key,
    combined_candidate_score_from_metal_states,
    omol_score,
    omol_score_cache_clear,
    omol_score_cache_info,
    omol_score_from_parts,
    organic_core_score,
    organic_core_score_cache_info,
    post_reinsertion_score_cache_info,
)


def test_fallback_omol_score_cache_hits_for_equivalent_reloads() -> None:
    omol_score_cache_clear()
    xyz_block = """3
LiCO
Li 0.0 0.0 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
"""

    first = pybel.readstring("xyz", xyz_block)
    second = pybel.readstring("xyz", xyz_block)

    score_first = omol_score(first)
    hits_after_first, misses_after_first, size_after_first = omol_score_cache_info()
    score_second = omol_score(second)
    hits_after_second, misses_after_second, size_after_second = omol_score_cache_info()

    assert score_first == score_second
    assert (hits_after_first, misses_after_first, size_after_first) == (0, 1, 1)
    assert hits_after_second == 1
    assert misses_after_second == 1
    assert size_after_second == 1


def test_fallback_omol_score_cache_clear_resets_stats() -> None:
    xyz_block = """2
H2
H 0.0 0.0 0.0
H 0.0 0.0 0.74
"""

    omol_score_cache_clear()
    omol_score(pybel.readstring("xyz", xyz_block))
    assert omol_score_cache_info() == (0, 1, 1)

    omol_score_cache_clear()
    assert omol_score_cache_info() == (0, 0, 0)


def test_omol_score_from_parts_reuses_no_metal_post_key() -> None:
    omol_score_cache_clear()
    omol = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )

    organic_score = organic_core_score(omol)
    post_base_key = build_post_reinsertion_base_key(omol)
    assert organic_core_score_cache_info() == (0, 1, 1)

    score_a = omol_score_from_parts(organic_score, post_base_key, omol)
    assert post_reinsertion_score_cache_info() == (0, 1, 1)

    score_b = omol_score_from_parts(organic_score, post_base_key, omol)
    assert score_a == score_b
    assert post_reinsertion_score_cache_info() == (1, 1, 1)
    assert score_a == omol_score(omol)


def test_combined_candidate_score_from_metal_states_matches_full_score_without_recombining() -> None:
    omol_score_cache_clear()
    no_metal = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    metal_states = (
        MetalAtomPosition(
            idx=1,
            symbol="Li",
            element_idx=3,
            valence=1,
            radical_num=0,
            position_x=0.0,
            position_y=0.0,
            position_z=0.0,
        ),
        MetalAtomPosition(
            idx=2,
            symbol="Li",
            element_idx=3,
            valence=1,
            radical_num=0,
            position_x=4.0,
            position_y=0.0,
            position_z=0.0,
        ),
    )
    combined = combine_metal_with_omol(no_metal, metal_states)
    organic_score = organic_core_score(no_metal)
    post_base_key = build_post_reinsertion_base_key(no_metal)
    base_symmetry_penalty, charged_atom_snapshots = build_post_reinsertion_base_components(no_metal)

    direct_score = combined_candidate_score_from_metal_states(
        organic_score,
        post_base_key,
        base_symmetry_penalty,
        charged_atom_snapshots,
        metal_states,
    )

    assert direct_score == omol_score(combined)
    assert post_reinsertion_score_cache_info() == (0, 1, 1)

    repeated_score = combined_candidate_score_from_metal_states(
        organic_score,
        post_base_key,
        base_symmetry_penalty,
        charged_atom_snapshots,
        metal_states,
    )
    assert repeated_score == direct_score
    assert post_reinsertion_score_cache_info() == (1, 1, 1)


def test_score_cache_quantizes_coordinates_at_six_decimal_places() -> None:
    omol_score_cache_clear()
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

    score_first = omol_score(first)
    assert omol_score_cache_info() == (0, 1, 1)

    score_second = omol_score(second)
    assert score_first == score_second
    assert omol_score_cache_info() == (1, 1, 1)


def test_reconstruction_state_caches_score_keys_until_omol_changes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    omol_score_cache_clear()
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
    original_builder = scoring_module._build_score_key

    def tracking_builder(current_omol):
        nonlocal score_key_calls
        score_key_calls += 1
        return original_builder(current_omol)

    monkeypatch.setattr(scoring_module, "_build_score_key", tracking_builder)

    first_organic_score = organic_core_score(state)
    first_organic_key = build_post_reinsertion_base_key(state)
    second_organic_score = organic_core_score(state)
    second_organic_key = build_post_reinsertion_base_key(state)
    assert first_organic_score == second_organic_score
    assert first_organic_key == second_organic_key
    assert score_key_calls == 1

    first_full_score = omol_score(state)
    second_full_score = omol_score(state)
    assert first_full_score == second_full_score
    assert score_key_calls == 1

    machine = OmolStateMachine.from_reconstruction_state(state)
    machine.run_omol_stage("clone", lambda current_omol: (current_omol.clone, True))
    next_state = machine.freeze_like(state)

    cloned_organic_score = organic_core_score(next_state)
    cloned_organic_key = build_post_reinsertion_base_key(next_state)
    cloned_full_score = omol_score(next_state)
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
    original_builder = scoring_module._build_score_key

    def tracking_builder(current_omol):
        nonlocal score_key_calls
        score_key_calls += 1
        return original_builder(current_omol)

    def no_hit_stage(current_omol):
        return current_omol, False

    monkeypatch.setattr(scoring_module, "_build_score_key", tracking_builder)

    organic_core_score(state)
    assert score_key_calls == 1

    machine = OmolStateMachine.from_reconstruction_state(state)
    machine.run_omol_stage("noop", no_hit_stage)
    next_state = machine.freeze_like(state)

    organic_core_score(next_state)
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
            "organic_core_score": 1.0,
            "post_reinsertion_base_key": "cached-key",
            "post_reinsertion_base_symmetry_penalty": 2.0,
            "post_reinsertion_charged_atom_snapshots": ((-1, 0.0, 0.0, 0.0),),
            "score": 3.0,
            "keep_me": "yes",
        },
    )

    machine = OmolStateMachine.from_reconstruction_state(state)
    machine.run_omol_stage("clone", lambda current_omol: (current_omol.clone, True))
    next_state = machine.freeze_like(state)

    assert next_state.omol_revision == 1
    assert "organic_core_score" not in next_state.metadata
    assert "post_reinsertion_base_key" not in next_state.metadata
    assert "post_reinsertion_base_symmetry_penalty" not in next_state.metadata
    assert "post_reinsertion_charged_atom_snapshots" not in next_state.metadata
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
    score_a = organic_core_score(no_metal_a)
    key_a = build_post_reinsertion_base_key(no_metal_a)
    symmetry_a, charged_atoms_a = build_post_reinsertion_base_components(no_metal_a)
    score_b = organic_core_score(no_metal_b)
    key_b = build_post_reinsertion_base_key(no_metal_b)
    symmetry_b, charged_atoms_b = build_post_reinsertion_base_components(no_metal_b)
    state_a = ReconstructionState(
        omol=no_metal_a,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        metadata={
            "organic_core_score": score_a,
            "post_reinsertion_base_key": key_a,
            "post_reinsertion_base_symmetry_penalty": symmetry_a,
            "post_reinsertion_charged_atom_snapshots": charged_atoms_a,
        },
    )
    state_b = ReconstructionState(
        omol=no_metal_b,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        metadata={
            "organic_core_score": score_b,
            "post_reinsertion_base_key": key_b,
            "post_reinsertion_base_symmetry_penalty": symmetry_b,
            "post_reinsertion_charged_atom_snapshots": charged_atoms_b,
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

    def fake_combined_score_from_metal_states(*args, **kwargs):
        calls["count"] += 1
        return float(calls["count"])

    monkeypatch.setattr(
        scoring_module,
        "combined_candidate_score_from_metal_states",
        fake_combined_score_from_metal_states,
    )

    candidate_machine = MetalCandidateStateMachine.from_candidate_state(candidate)
    candidate_machine.set_no_metal_state("attach_state_a", state_a)
    candidate = candidate_machine.freeze()

    assert candidate.combined_score() == 1.0
    assert candidate.combined_score() == 1.0
    assert calls["count"] == 1
    materialize_calls = {"count": 0}

    def tracked_combine(omol, metal_states):
        materialize_calls["count"] += 1
        return combine_metal_with_omol(omol, metal_states)

    first_combined = candidate.materialize_combined_omol(tracked_combine)
    second_combined = candidate.materialize_combined_omol(tracked_combine)
    assert first_combined is second_combined
    assert candidate.combined_omol is first_combined
    assert materialize_calls["count"] == 1

    candidate_machine = MetalCandidateStateMachine.from_candidate_state(candidate)
    candidate_machine.set_no_metal_state("attach_state_b", state_b)
    next_candidate = candidate_machine.freeze()

    assert next_candidate.score is None
    assert next_candidate.combined_omol is None
    assert "score" not in next_candidate.metadata
    assert next_candidate.combined_score() == 2.0
    assert calls["count"] == 2
