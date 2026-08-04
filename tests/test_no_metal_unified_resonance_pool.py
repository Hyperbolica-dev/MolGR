from __future__ import annotations

import pytest
from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.pipeline import reconstruct_without_metals
from molgr.fallback.state import ReconstructionState
from molgr.fallback.utils.electrons import (
    get_lone_pair_count,
    get_unpaired_electron_count,
    has_unresolved_two_electron_center,
    set_unpaired_electron_count,
    set_unresolved_two_electron_center,
)
from molgr.fallback.utils.no_metals import neighbor_radicals
from molgr.fallback.utils.no_metals import resonance as no_metal_resonance
from molgr.fallback.utils.no_metals import selection as no_metal_selection


def _oxygen_radical_pair() -> pybel.Molecule:
    obmol = ob.OBMol()
    obmol.BeginModify()
    for _ in range(2):
        atom = obmol.NewAtom()
        atom.SetAtomicNum(8)
        atom.SetFormalCharge(0)
        set_unpaired_electron_count(atom, 1)
    obmol.AddBond(1, 2, 1)
    obmol.EndModify()
    return pybel.Molecule(obmol)


def _state(omol: pybel.Molecule, *, source: str = "test") -> ReconstructionState:
    return ReconstructionState(
        omol=omol,
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("prepared_no_metal_seed",),
        metadata={"source": source},
    )


def _charge_separated_resonance(*, negative_atom_idx: int) -> ReconstructionState:
    omol = pybel.readstring("smi", "O=CO")
    obmol = omol.OBMol
    terminal_indices = (1, 3)
    positive_atom_idx = next(idx for idx in terminal_indices if idx != negative_atom_idx)
    obmol.GetAtom(negative_atom_idx).SetFormalCharge(-1)
    obmol.GetAtom(positive_atom_idx).SetFormalCharge(1)
    for bond in ob.OBMolBondIter(obmol):
        terminal_idx = (
            bond.GetEndAtomIdx() if bond.GetBeginAtomIdx() == 2 else bond.GetBeginAtomIdx()
        )
        bond.SetBondOrder(1 if terminal_idx == negative_atom_idx else 2)
    return _state(omol, source=f"negative-{negative_atom_idx}")


def test_no_metal_exact_score_ties_use_graph_key_not_enumeration_order(monkeypatch) -> None:
    lower_key = _charge_separated_resonance(negative_atom_idx=1)
    higher_key = _charge_separated_resonance(negative_atom_idx=3)
    candidates = [higher_key, lower_key]

    monkeypatch.setattr(
        reconstruct_without_metals,
        "_run_linear_preparation",
        lambda state: state,
    )
    monkeypatch.setattr(
        reconstruct_without_metals.neighbor_radicals,
        "enumerate_neighbor_radical_seeds",
        lambda state, *, exact_discrepancy=None: [state] if exact_discrepancy == 0 else [],
    )
    monkeypatch.setattr(
        reconstruct_without_metals.no_metal_resonance,
        "search_resonance_candidates",
        lambda states, **kwargs: candidates,
    )
    monkeypatch.setattr(
        reconstruct_without_metals.selection,
        "_no_metal_candidate_selection_key",
        lambda candidate, **kwargs: (0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0.0),
    )

    expected_key = no_metal_selection._no_metal_candidate_graph_tie_break_key(lower_key)
    assert expected_key < no_metal_selection._no_metal_candidate_graph_tie_break_key(higher_key)

    first = reconstruct_without_metals._run_no_metal_pipeline_from_state(lower_key)
    candidates.reverse()
    second = reconstruct_without_metals._run_no_metal_pipeline_from_state(lower_key)

    assert first is not None
    assert second is not None
    assert first.metadata["source"] == second.metadata["source"] == "negative-1"


def test_neighbor_radical_enumeration_returns_seeds_before_recovery() -> None:
    seeds = neighbor_radicals.enumerate_neighbor_radical_seeds(_state(_oxygen_radical_pair()))

    assert [seed.metadata["neighbor_radical_resolution"] for seed in seeds] == [
        "bond_order",
        "charge_separation",
        "charge_separation",
    ]
    assert [seed.metadata.get("positive_atom_idx") for seed in seeds] == [None, 1, 2]
    assert all("break_deformed_ene" not in seed.phase_history for seed in seeds)
    assert all("break_one_bond" not in seed.phase_history for seed in seeds)
    assert all("fresh_omol_charge_radical_final" not in seed.phase_history for seed in seeds)


def test_neighbor_radical_enumeration_supports_exact_charge_separation_layers() -> None:
    state = _state(_oxygen_radical_pair())

    layer_zero = neighbor_radicals.enumerate_neighbor_radical_seeds(
        state,
        exact_discrepancy=0,
    )
    layer_one = neighbor_radicals.enumerate_neighbor_radical_seeds(
        state,
        exact_discrepancy=1,
    )
    layer_two = neighbor_radicals.enumerate_neighbor_radical_seeds(
        state,
        exact_discrepancy=2,
    )

    assert [seed.metadata["neighbor_radical_resolution"] for seed in layer_zero] == ["bond_order"]
    assert [seed.metadata.get("positive_atom_idx") for seed in layer_one] == [1, 2]
    assert layer_two == []


def test_neighbor_radical_resolution_validates_before_consuming_electrons() -> None:
    omol = _oxygen_radical_pair()
    before = tuple(get_unpaired_electron_count(atom.OBAtom) for atom in omol)

    with pytest.raises(ValueError, match="positive_atom_idx"):
        neighbor_radicals._resolve_neighbor_radical_pair(
            omol,
            1,
            2,
            "charge_separation",
            3,
        )

    assert tuple(get_unpaired_electron_count(atom.OBAtom) for atom in omol) == before


def test_resonance_seed_pool_preserves_unresolved_two_electron_center() -> None:
    omol = pybel.readstring("smi", "[C]([H])[H]")
    center = omol.OBMol.GetAtom(1)
    set_unpaired_electron_count(center, 0)
    set_unresolved_two_electron_center(center, True)

    seeds = no_metal_resonance.build_resonance_seed_pool([_state(omol)])

    assert len(seeds) == 1
    center = seeds[0].omol.OBMol.GetAtom(1)
    assert (get_unpaired_electron_count(center), get_lone_pair_count(center)) == (0, 0)
    assert has_unresolved_two_electron_center(center)
    assert "assign_unresolved_two_electron_centers_for_resonance" not in seeds[0].phase_history


def test_resonance_seed_pool_preserves_charge_separated_oxygen_state(monkeypatch) -> None:
    obmol = ob.OBMol()
    oxygen = obmol.NewAtom()
    oxygen.SetAtomicNum(8)
    oxygen.SetFormalCharge(1)
    set_unpaired_electron_count(oxygen, 0)
    carbon = obmol.NewAtom()
    carbon.SetAtomicNum(6)
    carbon.SetFormalCharge(-1)
    set_unpaired_electron_count(carbon, 0)
    obmol.AddBond(1, 2, 1)

    state = _state(pybel.Molecule(obmol), source="ADEQOS-charge-separation")
    seeds = no_metal_resonance.build_resonance_seed_pool([state])

    assert seeds
    assert all(get_unpaired_electron_count(seed.omol.OBMol.GetAtom(1)) == 0 for seed in seeds)
    assert all(
        "refresh_electronic_labels_for_resonance" not in seed.phase_history for seed in seeds
    )

    monkeypatch.setattr(
        no_metal_resonance,
        "search_resonance_candidates",
        lambda states, **kwargs: [],
    )
    layer = reconstruct_without_metals._expand_resonance_layer(
        [state],
        traversal_policy=no_metal_resonance._default_resonance_traversal_policy(),
        config=None,
    )

    assert layer.resonance_seeds
    assert all(
        get_unpaired_electron_count(seed.omol.OBMol.GetAtom(1)) == 0
        for seed in layer.resonance_seeds
    )
    assert all(
        "refresh_electronic_labels_for_resonance" not in seed.phase_history
        for seed in layer.resonance_seeds
    )


def test_resonance_seed_pool_does_not_expand_multiple_unresolved_centers() -> None:
    omol = pybel.readstring("smi", "[C]([H])[H].[C]([H])[H]")
    for atom_idx in (1, 4):
        center = omol.OBMol.GetAtom(atom_idx)
        set_unpaired_electron_count(center, 0)
        set_unresolved_two_electron_center(center, True)

    seeds = no_metal_resonance.build_resonance_seed_pool([_state(omol)])

    assert len(seeds) == 1
    assert all(
        has_unresolved_two_electron_center(seeds[0].omol.OBMol.GetAtom(atom_idx))
        for atom_idx in (1, 4)
    )


def test_resonance_candidate_validation_resolves_new_unresolved_center(
    monkeypatch,
) -> None:
    omol = pybel.readstring("smi", "[C]([H])[H]")
    center = omol.OBMol.GetAtom(1)
    set_unpaired_electron_count(center, 0)
    set_unresolved_two_electron_center(center, True)
    state = _state(omol)

    def walk_root(root, *, visit, **kwargs):
        del kwargs
        visit(
            no_metal_resonance.resonance_utils.ResonanceSearchNode(
                root,
                no_metal_resonance.resonance_utils.build_resonance_state_key(root),
                0,
                0,
            )
        )

    monkeypatch.setattr(no_metal_resonance, "walk_radical_resonances", walk_root)
    monkeypatch.setattr(
        no_metal_resonance,
        "_score_reconstruction_candidate",
        lambda candidate, **kwargs: 0.0,
    )
    monkeypatch.setattr(
        no_metal_resonance,
        "_annotate_no_metal_candidate_topology",
        lambda candidate, **kwargs: None,
    )

    candidates = no_metal_resonance.search_resonance_candidates([state])

    assert len(candidates) == 1
    center = candidates[0].omol.OBMol.GetAtom(1)
    assert (get_unpaired_electron_count(center), get_lone_pair_count(center)) == (0, 1)
    assert not has_unresolved_two_electron_center(center)
    assert "resolve_unresolved_two_electron_centers_at_validation" in candidates[0].phase_history


def test_no_metal_pipeline_stops_after_first_valid_discrepancy_layer(monkeypatch) -> None:
    seed_state = _state(pybel.readstring("smi", "CC"), source="input")
    prepared_state = _state(pybel.readstring("smi", "CC"), source="prepared")
    branch_seeds = [
        _state(pybel.readstring("smi", "CC"), source="branch-0"),
        _state(pybel.readstring("smi", "C=C"), source="branch-1"),
    ]
    winner = branch_seeds[1]
    recorded: list[tuple[ReconstructionState, ...]] = []

    monkeypatch.setattr(
        reconstruct_without_metals,
        "_run_linear_preparation",
        lambda state: prepared_state,
    )
    monkeypatch.setattr(
        reconstruct_without_metals.neighbor_radicals,
        "enumerate_neighbor_radical_seeds",
        lambda state, *, exact_discrepancy=None: branch_seeds if exact_discrepancy == 0 else [],
    )

    def search(states, **kwargs):
        del kwargs
        recorded.append(tuple(states))
        return [winner]

    monkeypatch.setattr(
        reconstruct_without_metals.no_metal_resonance,
        "search_resonance_candidates",
        search,
    )
    monkeypatch.setattr(
        reconstruct_without_metals.selection,
        "_score_reconstruction_candidate",
        lambda candidate, **kwargs: 0.0,
    )
    monkeypatch.setattr(
        reconstruct_without_metals.selection,
        "_no_metal_candidate_selection_key",
        lambda candidate, **kwargs: (0.0, 0, 0.0, 0.0, 0.0, 0.0),
    )

    result = reconstruct_without_metals._run_no_metal_pipeline_from_state(seed_state)

    assert result is not None
    assert result.omol is winner.omol
    assert result.metadata["source"] == "branch-1"
    assert result.phase_history[-1] == "select_best_no_metal_candidate"
    assert recorded == [tuple(branch_seeds)]


def test_no_metal_pipeline_widens_with_one_shared_resonance_session(monkeypatch) -> None:
    seed_state = _state(pybel.readstring("smi", "CC"), source="input")
    prepared_state = _state(pybel.readstring("smi", "CC"), source="prepared")
    layer_zero = _state(pybel.readstring("smi", "CC"), source="layer-0")
    layer_one = _state(pybel.readstring("smi", "C=C"), source="layer-1")
    calls: list[tuple[int, object]] = []

    monkeypatch.setattr(
        reconstruct_without_metals,
        "_run_linear_preparation",
        lambda state: prepared_state,
    )

    def enumerate_layer(state, *, exact_discrepancy=None):
        del state
        return [layer_zero] if exact_discrepancy == 0 else [layer_one]

    monkeypatch.setattr(
        reconstruct_without_metals.neighbor_radicals,
        "enumerate_neighbor_radical_seeds",
        enumerate_layer,
    )

    def search(states, **kwargs):
        calls.append((len(states), kwargs["session"]))
        return [] if states[0] is layer_zero else [layer_one]

    monkeypatch.setattr(
        reconstruct_without_metals.no_metal_resonance,
        "search_resonance_candidates",
        search,
    )
    monkeypatch.setattr(
        reconstruct_without_metals.selection,
        "_no_metal_candidate_selection_key",
        lambda candidate, **kwargs: (0.0, 0, 0.0, 0.0, 0.0, 0, 0.0),
    )

    result = reconstruct_without_metals._run_no_metal_pipeline_from_state(seed_state)

    assert result is not None
    assert result.metadata["source"] == "layer-1"
    assert len(calls) == 2
    assert calls[0][1] is calls[1][1]


def test_resonance_candidate_deduplication_spans_all_seeds(monkeypatch) -> None:
    first = _state(pybel.readstring("smi", "[CH3]"), source="first")
    second = _state(pybel.readstring("smi", "[CH3]"), source="second")

    def walk_resonances(omol, *, visit, **kwargs):
        del kwargs
        resonance = omol.clone
        visit(
            no_metal_resonance.resonance_utils.ResonanceSearchNode(
                resonance,
                no_metal_resonance.resonance_utils.build_resonance_state_key(resonance),
                0,
                0,
            )
        )

    monkeypatch.setattr(
        no_metal_resonance,
        "walk_radical_resonances",
        walk_resonances,
    )
    monkeypatch.setattr(
        no_metal_resonance.resonance_utils,
        "process_resonance",
        lambda omol, remaining_charge: (omol, remaining_charge, False),
    )
    monkeypatch.setattr(no_metal_resonance, "validate_omol", lambda *args: True)
    monkeypatch.setattr(
        no_metal_resonance,
        "_score_reconstruction_candidate",
        lambda candidate, **kwargs: 0.0,
    )
    monkeypatch.setattr(
        no_metal_resonance,
        "_annotate_no_metal_candidate_topology",
        lambda candidate, **kwargs: None,
    )
    monkeypatch.setattr(ReconstructionState, "score", lambda self, *args, **kwargs: 0.0)

    candidates = no_metal_resonance.search_resonance_candidates([first, second])

    assert len(candidates) == 1


def test_resonance_candidate_deduplication_spans_session_calls(monkeypatch) -> None:
    state = _state(pybel.readstring("smi", "[CH3]"))

    def walk_resonances(omol, *, visit, **kwargs):
        del kwargs
        visit(
            no_metal_resonance.resonance_utils.ResonanceSearchNode(
                omol.clone,
                no_metal_resonance.resonance_utils.build_resonance_state_key(omol),
                0,
                0,
            )
        )

    monkeypatch.setattr(no_metal_resonance, "walk_radical_resonances", walk_resonances)
    monkeypatch.setattr(no_metal_resonance, "validate_omol", lambda *args: True)
    monkeypatch.setattr(
        no_metal_resonance,
        "_score_reconstruction_candidate",
        lambda candidate, **kwargs: 0.0,
    )
    monkeypatch.setattr(
        no_metal_resonance,
        "_annotate_no_metal_candidate_topology",
        lambda candidate, **kwargs: None,
    )
    session = no_metal_resonance._ResonanceSearchSession()

    first = no_metal_resonance.search_resonance_candidates([state], session=session)
    second = no_metal_resonance.search_resonance_candidates([state], session=session)

    assert len(first) == 1
    assert second == []


def test_resonance_traversal_labels_prune_only_dominated_paths() -> None:
    state_key = (("state",), 0)
    labels = {}

    assert no_metal_resonance._register_traversal_label(labels, state_key, (2, 1))
    assert not no_metal_resonance._register_traversal_label(labels, state_key, (2, 2))
    assert no_metal_resonance._register_traversal_label(labels, state_key, (1, 2))
    assert no_metal_resonance._register_traversal_label(labels, state_key, (1, 0))
    assert labels[state_key] == [(1, 0)]
