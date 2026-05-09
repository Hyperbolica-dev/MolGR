# pyright: reportMissingImports=false

from __future__ import annotations

from collections import deque
from dataclasses import replace

import pytest


pytest.importorskip("openbabel")

from openbabel import pybel

from molgr.config import make_default_config
from molgr.fallback.pipeline import resonance as resonance_module
from molgr.fallback.pipeline.resonance import (
    ResonanceTraversalContext,
    ResonanceTraversalMove,
    build_processed_resonance_key,
    build_resonance_state_key,
    get_radical_resonances,
    process_resonance,
)
from molgr.fallback.stages import clean as clean_module
from molgr.fallback.state import OmolStateMachine, ReconstructionState
from molgr.fallback.utils import resonance as resonance_utils_module
from molgr.fallback.utils.force_field import selection_force_field_energy
from molgr.fallback.utils.no_metals import resonance as no_metal_resonance_module


def _make_seed(smiles: str, radical_atom_indices: tuple[int, ...]) -> pybel.Molecule:
    mol = pybel.readstring("smi", smiles)
    for idx in radical_atom_indices:
        mol.OBMol.GetAtom(idx).SetSpinMultiplicity(1)
    return mol


def _naive_resonance_state_keys(seed: pybel.Molecule, max_depth: int = 2) -> set[tuple]:
    root_key, bond_index_map = resonance_utils_module._build_resonance_search_context(seed)
    seen = {root_key}
    frontier = deque([(seed, root_key, 0)])
    while frontier:
        current, current_key, depth = frontier.popleft()
        if depth >= max_depth:
            continue
        for move in resonance_utils_module._enumerate_one_step_resonance_moves(
            current,
            current_key,
            bond_index_map,
        ):
            next_resonance = resonance_utils_module._materialize_one_step_resonance(
                current, move.idxs
            )
            next_key = move.next_state_key
            if next_key in seen:
                continue
            seen.add(next_key)
            frontier.append((next_resonance, next_key, depth + 1))
    return seen


def test_build_resonance_state_key_matches_clone_and_changes_with_state() -> None:
    seed = _make_seed("C=CC=C", (2,))
    clone = seed.clone
    moved = get_radical_resonances(seed, max_depth=1)[1]

    assert build_resonance_state_key(seed) == build_resonance_state_key(clone)
    assert build_resonance_state_key(seed) != build_resonance_state_key(moved)


def test_build_processed_resonance_key_matches_clone_and_changes_with_state() -> None:
    seed = _make_seed("C=CC=C", (2,))
    clone = seed.clone
    moved = get_radical_resonances(seed, max_depth=1)[1]

    assert build_processed_resonance_key(seed) == build_processed_resonance_key(clone)
    assert build_processed_resonance_key(seed) != build_processed_resonance_key(moved)


def test_get_radical_resonances_avoids_smiles_serialization_for_dedup(monkeypatch) -> None:
    original_write = pybel.Molecule.write

    def patched_write(self, format: str = "smi", filename=None, overwrite=False, opt=None):
        if format == "smi":
            raise AssertionError("SMILES serialization should not be used for resonance dedup")
        return original_write(self, format=format, filename=filename, overwrite=overwrite, opt=opt)

    monkeypatch.setattr(pybel.Molecule, "write", patched_write)

    seed = _make_seed("C=CC=C", (2,))
    resonances = get_radical_resonances(seed)

    assert len(resonances) > 1
    assert len({build_resonance_state_key(mol) for mol in resonances}) == len(resonances)


def test_incremental_resonance_keys_match_rebuilt_keys() -> None:
    seed = _make_seed("C=CC=C", (2,))
    root_key, bond_index_map = resonance_utils_module._build_resonance_search_context(seed)

    for move in resonance_utils_module._enumerate_one_step_resonance_moves(
        seed,
        root_key,
        bond_index_map,
    ):
        moved = resonance_utils_module._materialize_one_step_resonance(seed, move.idxs)
        incremental_key = move.next_state_key
        assert incremental_key == build_resonance_state_key(moved)


def test_get_radical_resonances_builds_root_key_only_once(monkeypatch) -> None:
    original = resonance_module.build_resonance_state_key
    calls = []

    def wrapped(omol):
        calls.append(1)
        return original(omol)

    monkeypatch.setattr(resonance_module, "build_resonance_state_key", wrapped)
    seed = _make_seed("C=CC=C", (2,))
    resonances = get_radical_resonances(seed)

    assert len(resonances) > 1
    assert len(calls) == 0


def test_get_radical_resonances_matches_naive_bfs_on_independent_components() -> None:
    seed = _make_seed("C=CC=C.C=CC=C", (2, 6))

    resonances = get_radical_resonances(seed, max_depth=2)
    actual_keys = {build_resonance_state_key(mol) for mol in resonances}
    expected_keys = _naive_resonance_state_keys(seed, max_depth=2)

    assert actual_keys == expected_keys


def test_get_radical_resonances_accepts_traversal_policy_to_prune_directions() -> None:
    seed = _make_seed("C=CC=C.C=CC=C", (2, 6))
    calls = []

    def traversal_policy(
        context: ResonanceTraversalContext,
        moves: tuple[ResonanceTraversalMove, ...],
    ) -> list[ResonanceTraversalMove]:
        calls.append((context.depth, tuple(move.path for move in moves)))
        if context.depth == 0:
            return [move for move in moves if move.path == (2, 3, 4)]
        return []

    resonances = get_radical_resonances(
        seed,
        max_depth=2,
        traversal_policy=traversal_policy,
    )

    assert len(resonances) == 2
    assert calls[0] == (0, ((2, 3, 4), (6, 7, 8)))
    assert calls[1][0] == 1
    assert (6, 7, 8) in calls[1][1]


def test_clean_resonances_returns_false_when_no_rule_hits(monkeypatch) -> None:
    ordered_rules = (11, 0, 1, 2, 3, 4, 9, 5, 6, 7, 8, 9, 10, 12, 13)
    calls = []

    def make_stage(rule_id):
        def stage(omol):
            calls.append(rule_id)
            return omol, False

        return stage

    for rule_id in set(ordered_rules):
        monkeypatch.setattr(clean_module, f"clean_resonances_{rule_id}", make_stage(rule_id))

    omol = pybel.readstring("smi", "CC")
    result, hit = clean_module.clean_resonances(omol)

    assert result is omol
    assert hit is False
    assert calls == list(ordered_rules)


def test_clean_resonances_accumulates_stage_hits_without_global_signature(monkeypatch) -> None:
    ordered_rules = (11, 0, 1, 2, 3, 4, 9, 5, 6, 7, 8, 9, 10, 12, 13)
    calls = []

    def make_stage(rule_id):
        def stage(omol):
            calls.append(rule_id)
            return omol, rule_id in {3, 10}

        return stage

    for rule_id in set(ordered_rules):
        monkeypatch.setattr(clean_module, f"clean_resonances_{rule_id}", make_stage(rule_id))

    omol = pybel.readstring("smi", "CC")
    result, hit = clean_module.clean_resonances(omol)

    assert result is omol
    assert hit is True
    assert calls == list(ordered_rules)


def test_get_radical_resonances_limited_discrepancy_policy_prefers_low_discrepancy_paths(
    monkeypatch,
) -> None:
    seed = _make_seed("C=CC=C", (2,))
    state_key_a = (((6, 0, 0, False),), ((1, 2, 1, False),))
    state_key_b = (((6, 0, 0, False),), ((1, 2, 2, False),))
    state_key_c = (((6, 0, 0, False),), ((1, 2, 3, False),))
    state_key_a1 = (((6, 0, 0, False),), ((1, 2, 4, False),))
    state_key_b1 = (((6, 0, 0, False),), ((1, 2, 5, False),))
    move_a = resonance_utils_module._IndexedResonanceTraversalMove(  # type: ignore[attr-defined]
        idxs=(1, 2, 3),
        next_state_key=state_key_a,
    )
    move_b = resonance_utils_module._IndexedResonanceTraversalMove(  # type: ignore[attr-defined]
        idxs=(4, 5, 6),
        next_state_key=state_key_b,
    )
    move_c = resonance_utils_module._IndexedResonanceTraversalMove(  # type: ignore[attr-defined]
        idxs=(7, 8, 9),
        next_state_key=state_key_c,
    )
    child_move_a = resonance_utils_module._IndexedResonanceTraversalMove(  # type: ignore[attr-defined]
        idxs=(10, 11, 12),
        next_state_key=state_key_a1,
    )
    child_move_b = resonance_utils_module._IndexedResonanceTraversalMove(  # type: ignore[attr-defined]
        idxs=(20, 21, 22),
        next_state_key=state_key_b1,
    )
    policy = resonance_module.make_limited_discrepancy_force_field_traversal_policy(
        max_discrepancy=1
    )

    scores_by_state_and_move = {
        ("seed", (1, 2, 3)): 1.0,
        ("seed", (4, 5, 6)): 2.0,
        ("seed", (7, 8, 9)): 3.0,
        ("state-a", (10, 11, 12)): 1.0,
        ("state-b", (20, 21, 22)): 1.0,
    }

    def fake_force_field_score(omol, move_path, *, config=None):
        label = "seed" if omol is seed else omol
        return scores_by_state_and_move[(label, move_path)]

    def fake_enumerate(omol, state_key, bond_index_map):
        if omol is seed:
            return [move_c, move_b, move_a]
        if omol == "state-a":
            return [child_move_a]
        if omol == "state-b":
            return [child_move_b]
        return []

    def fake_materialize(omol, idxs):
        if omol is seed:
            return {
                (1, 2, 3): "state-a",
                (4, 5, 6): "state-b",
                (7, 8, 9): "state-c",
            }[idxs]
        if omol == "state-a":
            return "state-a1"
        if omol == "state-b":
            return "state-b1"
        return f"terminal-{idxs}"

    monkeypatch.setattr(
        resonance_utils_module,
        "_score_one_step_resonance_with_force_field",
        fake_force_field_score,
    )
    monkeypatch.setattr(
        resonance_utils_module, "_enumerate_one_step_resonance_moves", fake_enumerate
    )
    monkeypatch.setattr(resonance_utils_module, "_materialize_one_step_resonance", fake_materialize)

    resonances = get_radical_resonances(seed, max_depth=2, traversal_policy=policy)

    assert resonances[0] is seed
    assert resonances[1:] == ["state-a", "state-a1", "state-b", "state-b1"]


def test_force_field_limited_discrepancy_prunes_high_rank_resonance_branch(
    monkeypatch,
) -> None:
    seed = _make_seed("C=CC=C", (2,))
    state_key_a = (((6, 0, 0, False),), ((1, 2, 1, False),))
    state_key_b = (((6, 0, 0, False),), ((1, 2, 2, False),))
    state_key_c = (((6, 0, 0, False),), ((1, 2, 3, False),))
    move_a = resonance_utils_module._IndexedResonanceTraversalMove(  # type: ignore[attr-defined]
        idxs=(1, 2, 3),
        next_state_key=state_key_a,
    )
    move_b = resonance_utils_module._IndexedResonanceTraversalMove(  # type: ignore[attr-defined]
        idxs=(4, 5, 6),
        next_state_key=state_key_b,
    )
    move_c = resonance_utils_module._IndexedResonanceTraversalMove(  # type: ignore[attr-defined]
        idxs=(7, 8, 9),
        next_state_key=state_key_c,
    )
    policy = resonance_module.make_limited_discrepancy_force_field_traversal_policy(
        max_discrepancy=1
    )

    monkeypatch.setattr(
        resonance_utils_module,
        "_score_one_step_resonance_with_force_field",
        lambda omol, move_path, *, config=None: {
            (1, 2, 3): 1.0,
            (4, 5, 6): 2.0,
            (7, 8, 9): 3.0,
        }[move_path],
    )
    monkeypatch.setattr(
        resonance_utils_module,
        "_enumerate_one_step_resonance_moves",
        lambda omol, state_key, bond_index_map: [move_c, move_b, move_a] if omol is seed else [],
    )
    monkeypatch.setattr(
        resonance_utils_module,
        "_materialize_one_step_resonance",
        lambda omol, idxs: f"state-{idxs[0]}",
    )

    resonances = get_radical_resonances(seed, max_depth=1, traversal_policy=policy)

    assert resonances == [seed, "state-1", "state-4"]


def test_limited_discrepancy_search_prefers_lower_cost_duplicate_state(monkeypatch) -> None:
    seed = _make_seed("C=CC=C", (2,))
    shared_state_key = (((6, 0, 0, False),), ((1, 2, 1, False),))
    state_key_a = (((6, 0, 0, False),), ((1, 2, 2, False),))
    move_a = resonance_utils_module._IndexedResonanceTraversalMove(  # type: ignore[attr-defined]
        idxs=(1, 2, 3),
        next_state_key=state_key_a,
    )
    move_b = resonance_utils_module._IndexedResonanceTraversalMove(  # type: ignore[attr-defined]
        idxs=(4, 5, 6),
        next_state_key=shared_state_key,
    )
    child_move_a = resonance_utils_module._IndexedResonanceTraversalMove(  # type: ignore[attr-defined]
        idxs=(10, 11, 12),
        next_state_key=shared_state_key,
    )
    policy = resonance_module.make_limited_discrepancy_force_field_traversal_policy(
        max_discrepancy=1
    )

    scores_by_state_and_move = {
        ("seed", (1, 2, 3)): 1.0,
        ("seed", (4, 5, 6)): 2.0,
        ("state-a", (10, 11, 12)): 1.0,
    }

    def fake_force_field_score(omol, move_path, *, config=None):
        label = "seed" if omol is seed else omol
        return scores_by_state_and_move[(label, move_path)]

    def fake_enumerate(omol, state_key, bond_index_map):
        if omol is seed:
            return [move_b, move_a]
        if omol == "state-a":
            return [child_move_a]
        return []

    def fake_materialize(omol, idxs):
        if omol is seed and idxs == (1, 2, 3):
            return "state-a"
        return "shared-state"

    monkeypatch.setattr(
        resonance_utils_module,
        "_score_one_step_resonance_with_force_field",
        fake_force_field_score,
    )
    monkeypatch.setattr(
        resonance_utils_module, "_enumerate_one_step_resonance_moves", fake_enumerate
    )
    monkeypatch.setattr(resonance_utils_module, "_materialize_one_step_resonance", fake_materialize)

    resonances = get_radical_resonances(seed, max_depth=2, traversal_policy=policy)

    assert resonances[0] is seed
    assert resonances[1:] == ["state-a", "shared-state"]


def test_process_resonance_matches_equivalent_clone_states_without_cache() -> None:
    seed = _make_seed("C=CC=C", (2,))
    first_resonance, first_charge, _ = process_resonance(seed, 0)
    second_resonance, second_charge, _ = process_resonance(seed.clone, 0)

    assert first_charge == second_charge
    assert build_resonance_state_key(first_resonance) == build_resonance_state_key(second_resonance)


def test_resonance_move_score_cache_uses_config_in_key_and_execution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    resonance_utils_module.resonance_move_score_cache_clear()
    seed = _make_seed("C=CC=C", (2,))
    default_config = make_default_config()
    config_a = replace(
        default_config,
        force_field=replace(default_config.force_field, selection_force_field="uff"),
    )
    config_b = replace(
        default_config,
        force_field=replace(default_config.force_field, selection_force_field="auto"),
    )
    seen_configs = []

    def fake_selection_force_field_energy(omol, *, config=None):
        seen_configs.append(config)
        return 1.0

    monkeypatch.setattr(
        resonance_utils_module,
        "selection_force_field_energy",
        fake_selection_force_field_energy,
    )

    first = resonance_utils_module._score_one_step_resonance_with_force_field(
        seed,
        (2, 3, 4),
        config=config_a,
    )
    second = resonance_utils_module._score_one_step_resonance_with_force_field(
        seed,
        (2, 3, 4),
        config=config_a,
    )
    third = resonance_utils_module._score_one_step_resonance_with_force_field(
        seed,
        (2, 3, 4),
        config=config_b,
    )

    assert first == second == third == 1.0
    assert seen_configs == [config_a, config_b]


def test_omol_state_machine_caches_resonance_key_until_omol_changes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seed = _make_seed("C=CC=C", (2,))
    machine = OmolStateMachine(seed)
    calls = 0
    original_builder = resonance_module.build_resonance_state_key

    def tracking_builder(omol):
        nonlocal calls
        calls += 1
        return original_builder(omol)

    monkeypatch.setattr(resonance_module, "build_resonance_state_key", tracking_builder)

    first_key = machine.get_cached_omol_value(
        "resonance_state_key",
        resonance_module.build_resonance_state_key,
    )
    second_key = machine.get_cached_omol_value(
        "resonance_state_key",
        resonance_module.build_resonance_state_key,
    )
    assert first_key == second_key
    assert calls == 1

    machine.run_omol_stage("clone", lambda current_omol: (current_omol.clone, True))
    third_key = machine.get_cached_omol_value(
        "resonance_state_key",
        resonance_module.build_resonance_state_key,
    )
    assert third_key == first_key
    assert calls == 2


def test_recover_resonance_candidates_returns_valid_candidates() -> None:
    seed = _make_seed("C=CC=C", (2,))
    state = ReconstructionState(
        seed,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=1,
    )
    candidates = no_metal_resonance_module._recover_resonance_candidates(state)
    assert candidates
    scores = [selection_force_field_energy(candidate.omol) for candidate in candidates]
    assert min(scores) == sorted(scores)[0]


def test_recover_resonance_candidates_returns_candidates_without_shared_region_scoring() -> None:
    seed = _make_seed("C=CC=C", (2,))
    state = ReconstructionState(
        seed,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=1,
    )
    candidates = no_metal_resonance_module._recover_resonance_candidates(state)

    assert candidates


def test_recover_resonance_candidates_dedups_processed_states(monkeypatch) -> None:
    seed = _make_seed("C=CC=C", (2,))
    state = ReconstructionState(
        seed,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=1,
    )

    resonance_a = object()
    resonance_b = object()
    processed = object()
    validate_calls = 0
    organic_score_calls = 0
    force_field_score_calls = 0

    monkeypatch.setattr(
        no_metal_resonance_module,
        "get_radical_resonances",
        lambda omol: [resonance_a, resonance_b],
    )
    monkeypatch.setattr(
        resonance_utils_module,
        "process_resonance",
        lambda omol, given_charge: (processed, given_charge, True),
    )
    monkeypatch.setattr(
        resonance_utils_module,
        "build_processed_resonance_key",
        lambda omol: "same-processed-state",
    )

    def validate(*args, **kwargs):
        nonlocal validate_calls
        validate_calls += 1
        return True

    def organic_core_score(self):
        nonlocal organic_score_calls
        organic_score_calls += 1
        return 1.0

    def full_score(self):
        nonlocal force_field_score_calls
        force_field_score_calls += 1
        return 1.0

    monkeypatch.setattr(no_metal_resonance_module, "validate_omol", validate)
    monkeypatch.setattr(ReconstructionState, "organic_core_score", organic_core_score)
    monkeypatch.setattr(ReconstructionState, "full_score", full_score)

    candidates = no_metal_resonance_module._recover_resonance_candidates(state)

    assert len(candidates) == 1
    assert validate_calls == 1
    assert organic_score_calls == 1
    assert force_field_score_calls == 1


def test_recover_resonance_candidates_forwards_traversal_policy(monkeypatch) -> None:
    seed = _make_seed("C=CC=C", (2,))
    state = ReconstructionState(
        seed,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=1,
    )
    policy = object()
    recorded_policy = None

    def fake_get_radical_resonances(omol, *, traversal_policy=None, max_depth=2):
        nonlocal recorded_policy
        recorded_policy = traversal_policy
        return []

    monkeypatch.setattr(
        no_metal_resonance_module,
        "get_radical_resonances",
        fake_get_radical_resonances,
    )

    candidates = no_metal_resonance_module._recover_resonance_candidates(
        state,
        resonance_traversal_policy=policy,
    )

    assert candidates == []
    assert recorded_policy is policy


def test_resonance_pipeline_no_longer_exposes_incumbent_bound_helpers() -> None:
    from molgr.fallback.pipeline import reconstruct_without_metals as no_metal_module

    assert not hasattr(
        resonance_module,
        "estimate_remaining_resonance_score_improvement_upper_bound",
    )
    assert not hasattr(no_metal_module, "_recover_resonance_candidates_with_incumbent_bound")
    assert not hasattr(no_metal_module, "_recover_resonance_candidates")
