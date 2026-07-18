from __future__ import annotations

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.pipeline import reconstruct_without_metals
from molgr.fallback.state import ReconstructionState
from molgr.fallback.utils.no_metals import neighbor_radicals
from molgr.fallback.utils.no_metals import resonance as no_metal_resonance


def _oxygen_radical_pair() -> pybel.Molecule:
    obmol = ob.OBMol()
    obmol.BeginModify()
    for _ in range(2):
        atom = obmol.NewAtom()
        atom.SetAtomicNum(8)
        atom.SetFormalCharge(0)
        atom.SetSpinMultiplicity(1)
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
        reconstruct_without_metals.preparation,
        "prepare_no_metal_seed",
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
        reconstruct_without_metals.preparation,
        "prepare_no_metal_seed",
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
