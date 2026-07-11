from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import pytest
from rdkit import Chem


pytest.importorskip("openbabel")
pytest.importorskip("rdkit")

from openbabel import openbabel as ob

from molgr import _core
from molgr.config import MolGRConfig
from molgr.fallback.pipeline import reconstruct_with_metals
from molgr.fallback.utils.force_field import organic_force_field_evaluation
from molgr.fallback.utils.metals import preparation, scoring, search
from molgr.fallback.utils.no_metals import preparation as no_metal_preparation
from molgr.fallback.utils.no_metals import resonance as no_metal_resonance
from molgr.fallback.utils.no_metals import selection as no_metal_selection
from molgr.fallback.utils.organic_topology import compute_organic_topology_metrics
from molgr.interface import xyz_to_rdmol


_TMQMG_CSV = Path("/mnt/e/download/tmQMg_properties_and_targets.csv")
_TMQMG_XYZ_DIR = Path("/mnt/e/download/tmQMg_xyz/xyz")


pytestmark = pytest.mark.skipif(
    not _TMQMG_CSV.exists() or not _TMQMG_XYZ_DIR.exists(),
    reason="tmQMg local benchmark data is not available",
)


def _load_tmqmg_row(case_idx: int) -> dict[str, str]:
    with _TMQMG_CSV.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return rows[case_idx - 1]


def _load_tmqmg_xyz(case_id: str) -> str:
    return (_TMQMG_XYZ_DIR / f"{case_id}.xyz").read_text(encoding="utf-8")


def _canonical_reconstructed_smiles(
    xyz_block: str,
    *,
    backend: str,
    total_charge: int,
    total_radical_electrons: int,
    config,
) -> str:
    rdmol = xyz_to_rdmol(
        xyz_block,
        total_charge,
        total_radical_electrons + 1,
        backend=backend,
        make_dative_bonds=True,
        make_stereochemistry=True,
        config=config,
    )
    return Chem.MolToSmiles(
        Chem.RemoveHs(rdmol),
        canonical=True,
        isomericSmiles=True,
    )


def _python_scored_candidates(xyz_block: str, total_charge: int, total_radical_electrons: int):
    config = MolGRConfig()
    base_state = preparation.prepare_metal_state(
        xyz_block,
        total_charge,
        total_radical_electrons,
        config=config,
    )
    search_groups = search._build_metal_state_search_groups(
        base_state.available_valence_radical_states,
        config=config,
    )
    layered_groups = search._build_layered_metal_state_search_groups(
        search_groups,
        total_radical_electrons,
        config=config,
    )
    for layer_index, layer in enumerate(layered_groups):
        grouped_candidates = search._group_candidates_by_target_dp(
            base_state.phase_history,
            layer,
            total_charge,
            total_radical_electrons,
            config=config,
        )
        candidates = []
        for target_entry in grouped_candidates.values():
            if not target_entry:
                continue
            no_metal_state = (
                reconstruct_with_metals.reconstruct_without_metals.xyz_to_omol_no_metal_state(
                    base_state.no_metal_xyz_block,
                    target_entry[0].no_metal_charge_target,
                    target_entry[0].no_metal_radical_target,
                    config=config,
                )
            )
            if no_metal_state is None:
                continue
            for candidate in target_entry:
                candidates.append(
                    scoring._prepare_candidate_with_no_metal_state(
                        candidate,
                        no_metal_state,
                        config=config,
                    )
                )
        if candidates:
            scored = scoring.select_best_candidate(candidates, config=config)
            return layer_index, candidates, scored
    return None, [], None


def _python_no_metal_resonance_candidate_summaries(
    no_metal_xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
    *,
    config,
) -> tuple[dict[str, object], ...]:
    seed_state = no_metal_preparation._seed_state(
        no_metal_xyz_block,
        total_charge,
        total_radical_electrons,
    )
    candidates = []
    for state in no_metal_preparation._enumerate_no_metal_candidate_states(seed_state):
        candidates.extend(
            no_metal_resonance._recover_resonance_candidates(
                state,
                resonance_traversal_policy=no_metal_resonance._default_resonance_traversal_policy(
                    config
                ),
                config=config,
            )
        )

    summaries = []
    for candidate in candidates:
        no_metal_selection._score_reconstruction_candidate(candidate, config=config)
        selection_key = no_metal_selection._no_metal_candidate_selection_key(
            candidate,
            config=config,
        )
        summaries.append(
            {
                "smiles": candidate.omol.write("smi").split()[0],
                "resonance_index": candidate.metadata.get("resonance_index"),
                "score": float(candidate.metadata.get("score", 0.0)),
                "selection_key": selection_key,
            }
        )
    return tuple(summaries)


def _cpp_no_metal_resonance_candidate_summaries(
    no_metal_xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
    *,
    config,
) -> tuple[dict[str, object], ...]:
    summaries = []
    for item in _core.dev.pipeline.reconstruct_without_metals.debug_resonance_candidate_summaries(
        no_metal_xyz_block,
        total_charge,
        total_radical_electrons,
        config=config,
    ):
        selection_key = (
            -float(item["aromatic_stability_score"]),
            -int(item["aromatic_atom_count"]),
            -float(item["adjusted_max_conjugated_component_size"]),
            -float(item["adjusted_conjugated_atom_count"]),
            -float(item["adjusted_conjugated_bond_count"]),
            float(item["score"]),
        )
        summaries.append(
            {
                "smiles": item["smiles"],
                "resonance_index": item["resonance_index"],
                "score": float(item["score"]),
                "selection_key": selection_key,
            }
        )
    return tuple(summaries)


def _metal_state_signature(metal_state) -> tuple[object, ...]:
    return (
        int(metal_state.idx),
        str(metal_state.symbol),
        int(metal_state.element_idx),
        int(metal_state.valence),
        int(metal_state.radical_num),
        round(float(metal_state.position_x) * 1_000_000),
        round(float(metal_state.position_y) * 1_000_000),
        round(float(metal_state.position_z) * 1_000_000),
    )


def _metal_choice_signature(choice) -> tuple[tuple[object, ...], ...]:
    return tuple(_metal_state_signature(metal_state) for metal_state in choice)


def _metal_search_group_signature(group) -> tuple[tuple[tuple[object, ...], ...], ...]:
    return tuple(_metal_choice_signature(choice) for choice in group)


def _metal_search_layer_signature(layer) -> tuple[tuple[tuple[tuple[object, ...], ...], ...], ...]:
    return tuple(_metal_search_group_signature(group) for group in layer)


def _cpp_metal_search_group_signature(group) -> tuple[tuple[tuple[object, ...], ...], ...]:
    return tuple(tuple(tuple(metal_state) for metal_state in choice) for choice in group)


def _cpp_metal_search_layer_signature(
    layer,
) -> tuple[tuple[tuple[tuple[object, ...], ...], ...], ...]:
    return tuple(_cpp_metal_search_group_signature(group) for group in layer)


def _candidate_search_signature(candidate) -> dict[str, object]:
    return {
        "combination_index": int(candidate.metadata.get("combination_index", -1)),
        "no_metal_charge_target": int(candidate.no_metal_charge_target),
        "no_metal_radical_target": int(candidate.no_metal_radical_target),
        "metal_assignment_rank": float(candidate.metadata.get("metal_assignment_rank", 0.0)),
        "metal_states": _metal_choice_signature(candidate.metal_states),
        "phase_history": tuple(candidate.phase_history),
    }


def _python_metal_search_summary(
    xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
    *,
    config,
) -> dict[str, object]:
    base_state = preparation.prepare_metal_state(
        xyz_block,
        total_charge,
        total_radical_electrons,
        config=config,
    )
    state_search_groups = search._build_metal_state_search_groups(
        base_state.available_valence_radical_states,
        config=config,
    )
    layered_groups = search._build_layered_metal_state_search_groups(
        state_search_groups,
        total_radical_electrons,
        config=config,
    )
    target_buckets_by_layer = []
    for layer in layered_groups:
        grouped_candidates = search._group_candidates_by_target_dp(
            base_state.phase_history,
            layer,
            total_charge,
            total_radical_electrons,
            config=config,
        )
        layer_buckets = []
        for target, candidates in sorted(grouped_candidates.items()):
            layer_buckets.append(
                {
                    "target": tuple(target),
                    "candidates": tuple(
                        _candidate_search_signature(candidate) for candidate in candidates
                    ),
                }
            )
        target_buckets_by_layer.append(tuple(layer_buckets))

    return {
        "available_valence_radical_states": tuple(
            tuple(_metal_state_signature(metal_state) for metal_state in state_options)
            for state_options in base_state.available_valence_radical_states
        ),
        "base_phase_history": tuple(base_state.phase_history),
        "state_search_groups": tuple(
            _metal_search_group_signature(group) for group in state_search_groups
        ),
        "layered_state_search_groups": tuple(
            _metal_search_layer_signature(layer) for layer in layered_groups
        ),
        "target_buckets_by_layer": tuple(target_buckets_by_layer),
    }


def _cpp_candidate_search_signature(candidate: dict) -> dict[str, object]:
    return {
        "combination_index": int(candidate["combination_index"]),
        "no_metal_charge_target": int(candidate["no_metal_charge_target"]),
        "no_metal_radical_target": int(candidate["no_metal_radical_target"]),
        "metal_assignment_rank": float(candidate["metal_assignment_rank"]),
        "metal_states": tuple(tuple(item) for item in candidate["metal_states"]),
        "phase_history": tuple(candidate["phase_history"]),
    }


def _cpp_metal_search_summary(raw: dict) -> dict[str, object]:
    target_buckets_by_layer = []
    for layer_buckets in raw["target_buckets_by_layer"]:
        normalized_buckets = []
        for bucket in layer_buckets:
            normalized_buckets.append(
                {
                    "target": tuple(bucket["target"]),
                    "candidates": tuple(
                        _cpp_candidate_search_signature(candidate)
                        for candidate in bucket["candidates"]
                    ),
                }
            )
        target_buckets_by_layer.append(tuple(normalized_buckets))

    return {
        "available_valence_radical_states": tuple(
            tuple(tuple(metal_state) for metal_state in state_options)
            for state_options in raw["available_valence_radical_states"]
        ),
        "base_phase_history": tuple(raw["base_phase_history"]),
        "state_search_groups": tuple(
            _cpp_metal_search_group_signature(group) for group in raw["state_search_groups"]
        ),
        "layered_state_search_groups": tuple(
            _cpp_metal_search_layer_signature(layer) for layer in raw["layered_state_search_groups"]
        ),
        "target_buckets_by_layer": tuple(target_buckets_by_layer),
    }


def _pybel_atom_signature(omol) -> tuple[tuple[int, int, int, int, int, int, int, bool], ...]:
    atom_signature = []
    for atom in ob.OBMolAtomIter(omol.OBMol):
        atom_signature.append(
            (
                int(atom.GetAtomicNum()),
                int(atom.GetFormalCharge()),
                int(atom.GetSpinMultiplicity()),
                int(atom.GetHyb()),
                round(float(atom.GetX()) * 1_000_000),
                round(float(atom.GetY()) * 1_000_000),
                round(float(atom.GetZ()) * 1_000_000),
                bool(atom.IsAromatic()),
            )
        )
    return tuple(atom_signature)


def _pybel_bond_signature(omol) -> tuple[tuple[int, int, int, bool], ...]:
    bond_signature = []
    for bond in ob.OBMolBondIter(omol.OBMol):
        begin_idx = int(bond.GetBeginAtom().GetIdx())
        end_idx = int(bond.GetEndAtom().GetIdx())
        if begin_idx > end_idx:
            begin_idx, end_idx = end_idx, begin_idx
        bond_signature.append(
            (begin_idx, end_idx, int(bond.GetBondOrder()), bool(bond.IsAromatic()))
        )
    return tuple(sorted(bond_signature))


def _get_ptr(obmol: ob.OBMol) -> int:
    this = getattr(obmol, "this", None)
    if this is not None:
        return int(this)  # type: ignore[arg-type]
    return int(obmol)  # type: ignore[arg-type]


def _metrics_signature(metrics) -> dict[str, object]:
    return {
        "aromatic_atom_count": metrics.aromatic_atom_count,
        "aromatic_ring_count": metrics.aromatic_ring_count,
        "aromatic_stability_score": metrics.aromatic_stability_score,
        "conjugated_atom_count": metrics.conjugated_atom_count,
        "conjugated_bond_count": metrics.conjugated_bond_count,
        "max_conjugated_component_size": metrics.max_conjugated_component_size,
        "conjugated_atom_indices": tuple(metrics.conjugated_atom_indices),
    }


def _assert_cpp_python_topology_metrics_match(omol, *, config) -> None:
    py_metrics = _metrics_signature(
        compute_organic_topology_metrics(omol, config=config.organic_topology)
    )
    cpp_metrics = dict(
        _core.dev.utils.compute_organic_topology_metrics_ptr(
            _get_ptr(omol.OBMol),
            config=config,
        )
    )
    cpp_metrics["conjugated_atom_indices"] = tuple(cpp_metrics["conjugated_atom_indices"])

    for key, py_value in py_metrics.items():
        cpp_value = cpp_metrics[key]
        if isinstance(py_value, float):
            assert cpp_value == pytest.approx(py_value), key
        else:
            assert cpp_value == py_value, key


@pytest.mark.parametrize(
    "smiles",
    [
        "c1ccccc1",
        "C=CC=C",
        "[c-]1[c-][c-][cH][cH][cH]1",
        "c1cc2ccccc2cc1",
        "[n-]1cccc1",
        "C/C(=C/C(=N/c1c(C)cccc1C)/C)/[N-]c1c(C)cccc1C.Cc1ccccc1",
    ],
)
def test_cpp_python_organic_topology_metrics_match_smiles(smiles: str) -> None:
    from openbabel import pybel

    config = MolGRConfig()
    omol = pybel.readstring("smi", smiles)

    _assert_cpp_python_topology_metrics_match(omol, config=config)


def test_cpp_python_organic_topology_metrics_match_case_522_no_metal_candidate() -> None:
    row = _load_tmqmg_row(522)
    xyz_block = _load_tmqmg_xyz(row["id"])
    total_charge = int(row["charge"])
    total_radical_electrons = 0
    config = MolGRConfig()
    _layer_index, py_candidates, _py_best = _python_scored_candidates(
        xyz_block,
        total_charge,
        total_radical_electrons,
    )

    assert py_candidates
    assert py_candidates[0].no_metal_state is not None
    _assert_cpp_python_topology_metrics_match(py_candidates[0].no_metal_state.omol, config=config)


def test_cpp_python_organic_force_field_energy_matches_case_522_no_metal_candidates() -> None:
    row = _load_tmqmg_row(522)
    xyz_block = _load_tmqmg_xyz(row["id"])
    total_charge = int(row["charge"])
    total_radical_electrons = 0
    config = MolGRConfig()
    _layer_index, py_candidates, _py_best = _python_scored_candidates(
        xyz_block,
        total_charge,
        total_radical_electrons,
    )

    assert py_candidates
    for candidate in py_candidates:
        assert candidate.no_metal_state is not None
        py_energy = organic_force_field_evaluation(candidate.no_metal_state).energy_kj_mol
        cpp_energy = _core.dev.utils.organic_force_field_energy_ptr(
            _get_ptr(candidate.no_metal_state.omol.OBMol),
            config=config,
        )
        assert cpp_energy == pytest.approx(py_energy)


@pytest.mark.parametrize("case_idx", [430, 431])
def test_cpp_python_scored_candidate_parity_for_dmso_state_regressions(case_idx: int) -> None:
    row = _load_tmqmg_row(case_idx)
    xyz_block = _load_tmqmg_xyz(row["id"])
    total_charge = int(row["charge"])
    total_radical_electrons = 0
    config = MolGRConfig()
    cpp = _core.dev.pipeline.reconstruct_with_metals.debug_scored_candidate_summaries(
        xyz_block,
        total_charge,
        total_radical_electrons,
        config=config,
    )
    py_layer_index, py_candidates, py_best = _python_scored_candidates(
        xyz_block,
        total_charge,
        total_radical_electrons,
    )

    assert cpp["layer_index"] == py_layer_index
    assert len(cpp["candidates"]) == len(py_candidates)
    for cpp_item, py_item in zip(cpp["candidates"], py_candidates):
        assert cpp_item["combination_index"] == py_item.metadata.get("combination_index")
        assert cpp_item["no_metal_charge_target"] == py_item.no_metal_charge_target
        assert cpp_item["score"] == pytest.approx(py_item.score or 0.0)
        assert cpp_item["selection_key"] == py_item.metadata.get("selection_key")
        assert cpp_item["selected"] == (py_best is py_item)


def test_cpp_python_metal_search_stage_parity_case_522() -> None:
    row = _load_tmqmg_row(522)
    xyz_block = _load_tmqmg_xyz(row["id"])
    total_charge = int(row["charge"])
    total_radical_electrons = 0
    config = MolGRConfig()

    py_summary = _python_metal_search_summary(
        xyz_block,
        total_charge,
        total_radical_electrons,
        config=config,
    )
    cpp_summary = _cpp_metal_search_summary(
        _core.dev.pipeline.reconstruct_with_metals.debug_metal_search_summaries(
            xyz_block,
            total_charge,
            total_radical_electrons,
            config=config,
        )
    )

    assert cpp_summary == py_summary


def test_cpp_python_scored_candidate_parity_case_522() -> None:
    row = _load_tmqmg_row(522)
    xyz_block = _load_tmqmg_xyz(row["id"])
    total_charge = int(row["charge"])
    total_radical_electrons = 0
    config = MolGRConfig()

    _core.dev.pipeline.clear_force_field_evaluation_cache()
    _core.dev.pipeline.clear_uff_atom_typing_cache()
    cpp = _core.dev.pipeline.reconstruct_with_metals.debug_scored_candidate_summaries(
        xyz_block,
        total_charge,
        total_radical_electrons,
        config=config,
    )
    py_layer_index, py_candidates, py_best = _python_scored_candidates(
        xyz_block,
        total_charge,
        total_radical_electrons,
    )

    assert cpp["layer_index"] == py_layer_index
    assert len(cpp["candidates"]) == len(py_candidates)

    for cpp_item, py_item in zip(cpp["candidates"], py_candidates):
        assert cpp_item["combination_index"] == py_item.metadata.get("combination_index")
        assert cpp_item["no_metal_charge_target"] == py_item.no_metal_charge_target
        assert cpp_item["no_metal_radical_target"] == py_item.no_metal_radical_target
        assert py_item.no_metal_state is not None
        assert cpp_item["no_metal_total_charge"] == py_item.no_metal_state.total_charge
        assert (
            cpp_item["no_metal_total_radical_electrons"]
            == py_item.no_metal_state.total_radical_electrons
        )
        assert tuple(cpp_item["no_metal_atom_signature"]) == _pybel_atom_signature(
            py_item.no_metal_state.omol
        )
        assert tuple(cpp_item["no_metal_bond_signature"]) == _pybel_bond_signature(
            py_item.no_metal_state.omol
        )
        assert cpp_item["score"] == pytest.approx(py_item.score or 0.0)
        assert cpp_item["metal_assignment_rank"] == pytest.approx(
            float(py_item.metadata.get("metal_assignment_rank", 0.0))
        )
        assert cpp_item["organic_aromatic_ring_count"] == py_item.metadata.get(
            "organic_aromatic_ring_count"
        )
        assert cpp_item["organic_aromatic_stability_score"] == pytest.approx(
            float(py_item.metadata.get("organic_aromatic_stability_score", 0.0))
        )
        assert cpp_item["organic_charge_localization_penalty"] == pytest.approx(
            float(py_item.metadata.get("organic_charge_localization_penalty", 0.0))
        )
        assert cpp_item["organic_radical_localization_penalty"] == pytest.approx(
            float(py_item.metadata.get("organic_radical_localization_penalty", 0.0))
        )
        assert cpp_item["metal_discordance_structural_count"] == pytest.approx(
            float(py_item.metadata.get("metal_discordance_structural_count", 0.0))
        )
        assert cpp_item["metal_discordance_count"] == pytest.approx(
            float(py_item.metadata.get("metal_discordance_count", 0.0))
        )
        assert cpp_item["metal_discordance_aromatic_ring_deficit_count"] == py_item.metadata.get(
            "metal_discordance_aromatic_ring_deficit_count"
        )
        assert cpp_item["metal_discordance_aromatic_stability_deficit"] == pytest.approx(
            float(py_item.metadata.get("metal_discordance_aromatic_stability_deficit", 0.0))
        )
        assert cpp_item["passes_metal_discordance_filter"] == py_item.metadata.get(
            "passes_metal_discordance_filter"
        )
        assert cpp_item["selection_key"] == py_item.metadata.get("selection_key")

        assert cpp_item["selected"] == (py_best is py_item)


@pytest.mark.parametrize("case_idx", [113, 627, 684, 727])
def test_cpp_python_no_metal_resonance_candidate_parity_cases(case_idx: int) -> None:
    row = _load_tmqmg_row(case_idx)
    xyz_block = _load_tmqmg_xyz(row["id"])
    total_charge = int(row["charge"])
    total_radical_electrons = 0
    config = MolGRConfig()
    base_state = preparation.prepare_metal_state(
        xyz_block,
        total_charge,
        total_radical_electrons,
        config=config,
    )
    _layer_index, py_candidates, py_best = _python_scored_candidates(
        xyz_block,
        total_charge,
        total_radical_electrons,
    )

    assert py_best is not None
    assert py_candidates
    target_charge = py_best.no_metal_charge_target
    target_radicals = py_best.no_metal_radical_target

    py_summaries = _python_no_metal_resonance_candidate_summaries(
        base_state.no_metal_xyz_block,
        target_charge,
        target_radicals,
        config=config,
    )
    cpp_summaries = _cpp_no_metal_resonance_candidate_summaries(
        base_state.no_metal_xyz_block,
        target_charge,
        target_radicals,
        config=config,
    )

    assert len(cpp_summaries) == len(py_summaries)
    for cpp_item, py_item in zip(cpp_summaries, py_summaries):
        assert cpp_item["smiles"] == py_item["smiles"]
        assert cpp_item["resonance_index"] == py_item["resonance_index"]
        assert cpp_item["score"] == pytest.approx(py_item["score"])
        assert cpp_item["selection_key"] == pytest.approx(py_item["selection_key"])


@pytest.mark.parametrize("case_idx", [113, 627, 684, 727])
def test_cpp_python_final_smiles_match_resonance_representative_cases(case_idx: int) -> None:
    row = _load_tmqmg_row(case_idx)
    xyz_block = _load_tmqmg_xyz(row["id"])
    total_charge = int(row["charge"])
    total_radical_electrons = 0
    config = MolGRConfig()

    cpp_smiles = _canonical_reconstructed_smiles(
        xyz_block,
        backend="cpp",
        total_charge=total_charge,
        total_radical_electrons=total_radical_electrons,
        config=config,
    )
    py_smiles = _canonical_reconstructed_smiles(
        xyz_block,
        backend="python",
        total_charge=total_charge,
        total_radical_electrons=total_radical_electrons,
        config=config,
    )

    assert cpp_smiles == py_smiles


@pytest.mark.parametrize("case_idx", [113, 522, 627])
def test_cpp_all_accelerations_preserve_python_final_smiles(case_idx: int) -> None:
    row = _load_tmqmg_row(case_idx)
    xyz_block = _load_tmqmg_xyz(row["id"])
    total_charge = int(row["charge"])
    total_radical_electrons = 0
    baseline_config = MolGRConfig()
    cpp_config = replace(
        baseline_config,
        cpp_backend=replace(
            baseline_config.cpp_backend,
            max_threads=4,
            enable_target_bucket_parallelism=True,
            enable_candidate_scoring_parallelism=False,
            enable_uff_atom_typing_cache=True,
            enable_target_bucket_score_bundle_preheat=True,
            target_bucket_parallel_threshold=1,
            target_bucket_parallel_max_threads=1,
            candidate_score_parallel_threshold=32,
        ),
    )

    cpp_smiles = _canonical_reconstructed_smiles(
        xyz_block,
        backend="cpp",
        total_charge=total_charge,
        total_radical_electrons=total_radical_electrons,
        config=cpp_config,
    )
    py_smiles = _canonical_reconstructed_smiles(
        xyz_block,
        backend="python",
        total_charge=total_charge,
        total_radical_electrons=total_radical_electrons,
        config=baseline_config,
    )

    assert cpp_smiles == py_smiles
