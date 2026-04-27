# pyright: reportMissingImports=false

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import pytest


pytest.importorskip("openbabel")

from openbabel import pybel
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom

from molgr.config import make_default_config
from molgr.fallback.pipeline import reconstruct_without_metals as no_metal_module
from molgr.fallback.pipeline.reconstruct_without_metals import (
    xyz_to_omol_no_metal,
    xyz_to_omol_no_metal_state,
)
from molgr.fallback.stages.eliminate import eliminate_NNN
from molgr.fallback.stages.fresh import fresh_omol_charge_radical
from molgr.fallback.state import ReconstructionState
from molgr.fallback.utils.no_metals import preparation as no_metal_preparation_module
from molgr.fallback.utils.no_metals import resonance as no_metal_resonance_module
from molgr.fallback.utils.no_metals import selection as no_metal_selection_module
from molgr.fallback.utils.tools import typed_lru_cache
from molgr.utils.converter import pybel_to_rdmol
from molgr.utils.equivalence import check_equivalence


RDLogger.DisableLog("rdApp.*")  # type: ignore


def _total_charge_and_radicals(mol: Chem.Mol) -> tuple[int, int]:
    charge = 0
    radicals = 0
    for atom in mol.GetAtoms():  # pyright: ignore[reportCallIssue]
        charge += int(atom.GetFormalCharge())
        radicals += int(atom.GetNumRadicalElectrons())
    return charge, radicals


@typed_lru_cache(maxsize=128, typed=True)
def _seed_case(smiles: str) -> tuple[str, int, int]:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    mol_h = Chem.AddHs(mol)
    embed_code = rdDistGeom.EmbedMolecule(mol_h)  # pyright: ignore[reportCallIssue]
    assert int(embed_code) == 0
    charge, radicals = _total_charge_and_radicals(mol_h)
    return Chem.MolToXYZBlock(mol_h), charge, radicals


def _load_curated_smiles() -> list[str]:
    csv_path = Path(__file__).with_name("test_cases.csv")
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    curated_rows = [1, 2, 5, 10, 17]
    return [rows[idx - 1]["smiles"] for idx in curated_rows]


@pytest.mark.parametrize("smiles", _load_curated_smiles())
def test_fallback_no_metal_reconstructs_curated_cases(smiles: str) -> None:
    xyz_block, total_charge, total_radical_electrons = _seed_case(smiles)

    result = xyz_to_omol_no_metal(xyz_block, total_charge, total_radical_electrons)

    assert result is not None

    expected = Chem.MolFromSmiles(smiles)
    assert expected is not None
    equivalent, info = check_equivalence(
        expected,
        pybel_to_rdmol(result),
        use_chirality=True,
        max_resonance=100,
    )
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


def test_eliminate_nnn_negative_produces_closed_shell_azide() -> None:
    omol = pybel.readstring("smi", "[N][N][N]")

    omol, hit = fresh_omol_charge_radical(omol)

    assert hit
    assert [(atom.idx, atom.OBAtom.GetFormalCharge(), atom.OBAtom.GetSpinMultiplicity()) for atom in omol] == [
        (1, 0, 2),
        (2, 0, 1),
        (3, 0, 2),
    ]

    omol, given_charge, hit = eliminate_NNN(omol, 0, False)

    assert hit
    assert given_charge == 1
    assert [(atom.idx, atom.OBAtom.GetFormalCharge(), atom.OBAtom.GetSpinMultiplicity()) for atom in omol] == [
        (1, -1, 0),
        (2, 1, 0),
        (3, -1, 0),
    ]


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

    monkeypatch.setattr(no_metal_preparation_module, "make_connections", lambda omol: (omol, False))
    monkeypatch.setattr(no_metal_preparation_module, "pre_clean", lambda omol: (omol, False))
    monkeypatch.setattr(
        no_metal_preparation_module,
        "fresh_omol_charge_radical",
        lambda omol: (omol, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "eliminate_NNN",
        lambda omol, given_charge, positive: (omol, 7 if not positive else given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "eliminate_high_positive_charge_atoms",
        lambda omol, given_charge: (omol, given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "eliminate_CN_in_doubt",
        lambda omol, given_charge: (omol, given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "eliminate_carboxyl",
        lambda omol, given_charge: (omol, given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "eliminate_carbene_neighbor_heteroatom",
        lambda omol, given_charge: (omol, given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "clean_carbene_neighbor_unsaturated",
        lambda omol: (omol, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "clean_neighbor_radicals",
        lambda omol: (omol, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "eliminate_charge_spliting",
        lambda omol, given_charge: (omol, given_charge, False),
    )

    def record_break_deformed_ene(omol, given_charge, total_radical_electrons):
        recorded["break_deformed_ene"] = (given_charge, total_radical_electrons)
        return omol, False

    def record_break_one_bond(omol, given_charge, total_radical_electrons):
        recorded["break_one_bond"] = (given_charge, total_radical_electrons)
        return omol, given_charge, False

    monkeypatch.setattr(no_metal_preparation_module, "break_deformed_ene", record_break_deformed_ene)
    monkeypatch.setattr(no_metal_preparation_module, "break_one_bond", record_break_one_bond)

    next_state = no_metal_preparation_module._run_linear_pipeline(state)

    assert recorded["break_deformed_ene"] == (7, 3)
    assert recorded["break_one_bond"] == (7, 3)
    assert next_state.given_charge == 7


def test_fallback_no_metal_reuses_resonance_score_metadata(monkeypatch: pytest.MonkeyPatch) -> None:
    no_metal_module._run_no_metal_pipeline_cached.cache_clear()

    base_state = ReconstructionState(
        omol=pybel.readstring("smi", "CC"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz",),
    )
    weaker_candidate = ReconstructionState(
        omol=pybel.readstring("smi", "C=CC=C"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz", "validate_resonance_candidate"),
        metadata={"score": 5.0},
    )
    stronger_candidate = ReconstructionState(
        omol=pybel.readstring("smi", "C=CC=C"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz", "validate_resonance_candidate"),
        metadata={"score": 2.0},
    )

    monkeypatch.setattr(no_metal_preparation_module, "_seed_state", lambda *args, **kwargs: base_state)
    monkeypatch.setattr(no_metal_preparation_module, "_run_linear_pipeline", lambda state: state)
    monkeypatch.setattr(
        no_metal_preparation_module,
        "validate_omol",
        lambda *args, **kwargs: False,
    )
    monkeypatch.setattr(
        no_metal_resonance_module,
        "_recover_resonance_candidates",
        lambda state, **kwargs: [weaker_candidate, stronger_candidate],
    )
    monkeypatch.setattr(
        no_metal_selection_module,
        "_score_reconstruction_candidate",
        lambda candidate, **kwargs: float(candidate.metadata["score"]),
    )

    result = xyz_to_omol_no_metal_state("cached-resonance", 0, 0)

    assert result is not None
    assert result.metadata["score"] == 2.0
    assert result.phase_history[-1] == "select_best_resonance_candidate"


def test_no_metal_resonance_selection_prefers_aromatic_topology_before_force_field(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    no_metal_module._run_no_metal_pipeline_cached.cache_clear()

    base_state = ReconstructionState(
        omol=pybel.readstring("smi", "CC"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz",),
    )
    lower_force_field_candidate = ReconstructionState(
        omol=pybel.readstring("smi", "C=CC=C"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz", "validate_resonance_candidate"),
        metadata={"score": 1.0},
    )
    aromatic_candidate = ReconstructionState(
        omol=pybel.readstring("smi", "c1ccccc1"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz", "validate_resonance_candidate"),
        metadata={"score": 10.0},
    )

    monkeypatch.setattr(no_metal_preparation_module, "_seed_state", lambda *args, **kwargs: base_state)
    monkeypatch.setattr(no_metal_preparation_module, "_run_linear_pipeline", lambda state: state)
    monkeypatch.setattr(
        no_metal_preparation_module,
        "validate_omol",
        lambda *args, **kwargs: False,
    )
    monkeypatch.setattr(
        no_metal_resonance_module,
        "_recover_resonance_candidates",
        lambda state, **kwargs: [lower_force_field_candidate, aromatic_candidate],
    )
    monkeypatch.setattr(
        no_metal_selection_module,
        "_score_reconstruction_candidate",
        lambda candidate, **kwargs: float(candidate.metadata["score"]),
    )

    result = xyz_to_omol_no_metal_state("topology-first", 0, 0)

    assert result is not None
    assert result.omol.write("can").strip() == aromatic_candidate.omol.write("can").strip()
    assert result.metadata["organic_aromatic_atom_count"] == 6
    assert result.metadata["organic_topology_selection_key"][:4] == (-6, -6, -6, -6)
    assert result.phase_history[-1] == "select_best_resonance_candidate"


def test_no_metal_cached_pipeline_uses_config_in_key_and_execution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    no_metal_module._run_no_metal_pipeline_cached.cache_clear()

    base_state = ReconstructionState(
        omol=pybel.readstring("smi", "CC"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz",),
    )
    default_config = make_default_config()
    config_a = replace(
        default_config,
        resonance=replace(default_config.resonance, max_depth=2),
    )
    config_b = replace(
        default_config,
        resonance=replace(default_config.resonance, max_depth=3),
    )
    seen_configs = []

    monkeypatch.setattr(no_metal_preparation_module, "_seed_state", lambda *args, **kwargs: base_state)

    def fake_run_from_state(seed_state: ReconstructionState, *, config=None):
        assert seed_state is base_state
        seen_configs.append(config)
        return base_state

    monkeypatch.setattr(no_metal_module, "_run_no_metal_pipeline_from_state", fake_run_from_state)

    first = xyz_to_omol_no_metal_state("cfg-key", 0, 0, config=config_a)
    second = xyz_to_omol_no_metal_state("cfg-key", 0, 0, config=config_a)
    third = xyz_to_omol_no_metal_state("cfg-key", 0, 0, config=config_b)

    assert first is base_state
    assert second is base_state
    assert third is base_state
    assert seen_configs == [config_a, config_b]
