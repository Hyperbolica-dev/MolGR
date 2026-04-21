# pyright: reportMissingImports=false

from __future__ import annotations

import csv
from pathlib import Path

import pytest


pytest.importorskip("openbabel")

from openbabel import pybel

from molgr.fallback.pipeline import reconstruct_with_metals as with_metals_module
from molgr.fallback.pipeline.reconstruct_with_metals import (
    _group_candidates_by_target_dp,
    prepare_metal_state,
    xyz2omol,
)
from molgr.fallback.state import ReconstructionState, make_metal_candidate_state
from molgr.fallback.utils import scoring as scoring_module
from molgr.fallback.utils.dataclasses import MetalAtomPosition
from molgr.fallback.utils.scoring import (
    build_post_reinsertion_base_components,
    build_post_reinsertion_base_key,
    organic_core_score,
)


def _charge_and_radicals(mol: pybel.Molecule) -> tuple[int, int]:
    charge = 0
    radicals = 0
    for atom in mol.atoms:
        charge += atom.OBAtom.GetFormalCharge()
        radicals += atom.OBAtom.GetSpinMultiplicity()
    return charge, radicals


def _load_parity_cases() -> list[object]:
    csv_path = Path(__file__).with_name("test_cases.csv")
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    cases: list[object] = []
    for idx in [1, 2, 5, 10]:
        seed = pybel.readstring("smi", rows[idx - 1]["smiles"])
        xyz_block = str(seed.write("xyz"))
        total_charge, total_radical_electrons = _charge_and_radicals(seed)
        cases.append(
            pytest.param(
                xyz_block,
                total_charge,
                total_radical_electrons,
                id=f"curated-{idx}",
            )
        )

    cases.append(
        pytest.param(
            """3
LiCO
Li 0.0 0.0 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
            0,
            0,
            id="synthetic-li-co",
        )
    )
    return cases


@pytest.mark.parametrize(
    ("xyz_block", "total_charge", "total_radical_electrons"),
    _load_parity_cases(),
)
def test_fallback_with_metals_reconstructs_curated_cases(
    xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
) -> None:
    result = xyz2omol(xyz_block, total_charge, total_radical_electrons)
    assert result is not None


def test_fallback_with_metals_tracks_top_level_phases() -> None:
    state = prepare_metal_state(
        """3
LiCO
Li 0.0 0.0 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
        0,
        0,
    )

    assert state is not None
    assert state.phase_history[:4] == (
        "read_xyz",
        "build_metal_state_options",
        "remove_metal_atoms",
        "serialize_no_metal_xyz",
    )
    assert state.metadata["metal_atom_count"] == 1
def test_group_candidates_by_target_dp_caps_candidates_per_target() -> None:
    state_options = tuple(
        (
            MetalAtomPosition(idx, "Fe", 26, 0, 0, float(idx), 0.0, 0.0),
            MetalAtomPosition(idx, "Fe", 26, 2, 0, float(idx), 0.0, 0.0),
        )
        for idx in range(1, 5)
    )

    grouped = _group_candidates_by_target_dp(
        (),
        state_options,
        total_charge=0,
        total_radical_electrons=0,
        max_mixed_valence_spread=None,
        max_total_metal_radicals=0,
        max_assignments_per_target=1,
    )

    assert sorted(grouped) == [(-8, 0), (-6, 0), (-4, 0), (-2, 0), (0, 0)]
    assert sum(len(bucket) for bucket in grouped.values()) == 5
    assert all(len(bucket) == 1 for bucket in grouped.values())


def test_group_candidates_by_target_dp_uses_meet_in_the_middle_split(monkeypatch) -> None:
    state_options = tuple(
        (MetalAtomPosition(idx, "Fe", 26, 2, 0, float(idx), 0.0, 0.0),)
        for idx in range(1, 5)
    )
    calls: list[int] = []
    original_frontier = with_metals_module._enumerate_partial_assignment_frontier

    def tracking_frontier(*args, **kwargs):
        calls.append(len(args[0]))
        return original_frontier(*args, **kwargs)

    monkeypatch.setattr(
        with_metals_module,
        "_enumerate_partial_assignment_frontier",
        tracking_frontier,
    )

    _group_candidates_by_target_dp(
        (),
        state_options,
        total_charge=0,
        total_radical_electrons=0,
        max_mixed_valence_spread=3,
        max_total_metal_radicals=0,
        max_assignments_per_target=1,
    )

    assert calls == [2, 2]


def test_group_candidates_by_target_dp_preserves_cross_half_valence_spread_pruning() -> None:
    state_options = (
        (MetalAtomPosition(1, "Fe", 26, 2, 0, 1.0, 0.0, 0.0),),
        (MetalAtomPosition(2, "Fe", 26, 6, 0, 2.0, 0.0, 0.0),),
    )

    grouped = _group_candidates_by_target_dp(
        (),
        state_options,
        total_charge=0,
        total_radical_electrons=0,
        max_mixed_valence_spread=3,
        max_total_metal_radicals=0,
        max_assignments_per_target=1,
    )

    assert grouped == {}


def test_combine_partial_assignment_frontiers_only_checks_radical_compatible_buckets(
    monkeypatch,
) -> None:
    partial_assignment = with_metals_module._PartialMetalAssignment
    left_frontier = {
        (0, 1, ()): [
            partial_assignment(
                metal_states=(),
                total_metal_charge=0,
                total_metal_radicals=1,
                metal_assignment_rank=0.0,
                valence_bounds=(),
                order=0,
            )
        ]
    }
    right_frontier = {
        (1, 0, ()): [
            partial_assignment(
                metal_states=(),
                total_metal_charge=1,
                total_metal_radicals=0,
                metal_assignment_rank=0.0,
                valence_bounds=(),
                order=1,
            )
        ],
        (2, 3, ()): [
            partial_assignment(
                metal_states=(),
                total_metal_charge=2,
                total_metal_radicals=3,
                metal_assignment_rank=0.0,
                valence_bounds=(),
                order=2,
            )
        ],
    }
    calls = {"count": 0}

    def tracking_merge(left_bounds, right_bounds, max_mixed_valence_spread):
        calls["count"] += 1
        return ()

    monkeypatch.setattr(with_metals_module, "_merge_valence_bounds", tracking_merge)

    grouped = with_metals_module._combine_partial_assignment_frontiers(
        left_frontier,
        right_frontier,
        total_charge=0,
        total_radical_electrons=1,
        max_mixed_valence_spread=3,
        max_total_metal_radicals=None,
        max_assignments_per_target=1,
    )

    assert calls["count"] == 1
    assert sorted(grouped) == [(-1, 0)]


def test_iter_reachable_radical_buckets_uses_prefix_cutoff() -> None:
    prefix_reachability = with_metals_module._build_radical_prefix_reachability(
        {
            0: {1: []},
            2: {3: []},
            5: {7: []},
        }
    )

    reachable = with_metals_module._iter_reachable_radical_buckets(prefix_reachability, 3)

    assert [radicals for radicals, _ in reachable] == [0, 2]


def test_xyz2omol_state_uses_target_bucket_dp(monkeypatch) -> None:
    calls = {"count": 0}

    original_dp = with_metals_module._group_candidates_by_target_dp

    def tracking_dp(*args, **kwargs):
        calls["count"] += 1
        return original_dp(*args, **kwargs)

    monkeypatch.setattr(with_metals_module, "_group_candidates_by_target_dp", tracking_dp)

    with_metals_module.xyz2omol_state(
        """3
LiCO
Li 0.0 0.0 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
        0,
        0,
    )

    assert calls["count"] == 1


def test_score_candidate_with_no_metal_state_uses_metal_state_fast_path(monkeypatch) -> None:
    no_metal = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    no_metal_state = ReconstructionState(
        omol=no_metal,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        metadata={
            "organic_core_score": organic_core_score(no_metal),
            "post_reinsertion_base_key": build_post_reinsertion_base_key(no_metal),
        },
    )
    candidate = make_metal_candidate_state(
        (),
        (MetalAtomPosition(1, "Li", 3, 1, 0, 0.0, 0.0, 0.0),),
        -1,
        0,
        combination_index=0,
    )
    calls = {"count": 0}

    def fake_from_metal_states(*args, **kwargs):
        calls["count"] += 1
        return 1.25

    def fail_combine(*args, **kwargs):
        raise AssertionError("metal candidates should not be combined before winner selection")

    monkeypatch.setattr(
        scoring_module,
        "combined_candidate_score_from_metal_states",
        fake_from_metal_states,
    )
    monkeypatch.setattr(with_metals_module, "combine_metal_with_omol", fail_combine)

    scored = with_metals_module._score_candidate_with_no_metal_state(candidate, no_metal_state)

    assert calls["count"] == 1
    assert scored.score == 1.25
    assert scored.combined_omol is None


def test_score_candidate_with_no_metal_state_reuses_candidate_cache(monkeypatch) -> None:
    no_metal = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    base_symmetry_penalty, charged_atom_snapshots = build_post_reinsertion_base_components(no_metal)
    no_metal_state = ReconstructionState(
        omol=no_metal,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        metadata={
            "organic_core_score": organic_core_score(no_metal),
            "post_reinsertion_base_key": build_post_reinsertion_base_key(no_metal),
            "post_reinsertion_base_symmetry_penalty": base_symmetry_penalty,
            "post_reinsertion_charged_atom_snapshots": charged_atom_snapshots,
        },
    )
    candidate = make_metal_candidate_state(
        (),
        (MetalAtomPosition(1, "Li", 3, 1, 0, 0.0, 0.0, 0.0),),
        -1,
        0,
        combination_index=0,
    )
    calls = {"count": 0}

    def fake_from_metal_states(*args, **kwargs):
        calls["count"] += 1
        return 1.25

    monkeypatch.setattr(
        scoring_module,
        "combined_candidate_score_from_metal_states",
        fake_from_metal_states,
    )

    first = with_metals_module._score_candidate_with_no_metal_state(candidate, no_metal_state)
    second = with_metals_module._score_candidate_with_no_metal_state(first, no_metal_state)

    assert calls["count"] == 1
    assert second.score == 1.25


def test_xyz2omol_state_combines_only_best_candidate(monkeypatch) -> None:
    no_metal = pybel.readstring(
        "xyz",
        """2
CO
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
    )
    base_symmetry_penalty, charged_atom_snapshots = build_post_reinsertion_base_components(no_metal)
    no_metal_state = ReconstructionState(
        omol=no_metal,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        metadata={
            "organic_core_score": organic_core_score(no_metal),
            "post_reinsertion_base_key": build_post_reinsertion_base_key(no_metal),
            "post_reinsertion_base_symmetry_penalty": base_symmetry_penalty,
            "post_reinsertion_charged_atom_snapshots": charged_atom_snapshots,
        },
    )
    candidate_a = make_metal_candidate_state(
        (),
        (MetalAtomPosition(1, "Li", 3, 1, 0, 0.0, 0.0, 0.0),),
        -1,
        0,
        combination_index=0,
    )
    candidate_b = make_metal_candidate_state(
        (),
        (MetalAtomPosition(1, "Li", 3, 2, 0, 0.0, 0.0, 0.0),),
        -2,
        0,
        combination_index=1,
    )

    monkeypatch.setattr(
        with_metals_module,
        "_group_candidates_by_target_dp",
        lambda *args, **kwargs: {
            (candidate_a.no_metal_charge_target, candidate_a.no_metal_radical_target): [candidate_a],
            (candidate_b.no_metal_charge_target, candidate_b.no_metal_radical_target): [candidate_b],
        },
    )
    monkeypatch.setattr(
        with_metals_module,
        "xyz_to_omol_no_metal_state",
        lambda *args, **kwargs: no_metal_state,
    )

    def fake_score(*args, **kwargs):
        metal_states = args[-1]
        return 5.0 if metal_states[0].valence == 1 else 1.0

    monkeypatch.setattr(
        scoring_module,
        "combined_candidate_score_from_metal_states",
        fake_score,
    )

    calls = {"count": 0, "winner_valence": None}

    def fake_combine(omol, metal_states):
        calls["count"] += 1
        calls["winner_valence"] = metal_states[0].valence
        return omol

    monkeypatch.setattr(with_metals_module, "combine_metal_with_omol", fake_combine)

    result = with_metals_module.xyz2omol_state(
        """3
LiCO
Li 0.0 0.0 0.0
C 2.0 0.0 0.0
O 3.2 0.0 0.0
""",
        0,
        0,
    )

    assert result is not None
    assert result.score == 1.0
    assert calls["count"] == 1
    assert calls["winner_valence"] == 2
    assert result.combined_omol is no_metal
