# pyright: reportMissingImports=false

from __future__ import annotations

import csv
from pathlib import Path

import pytest


pytest.importorskip("openbabel")

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom

from molgr.fallback.pipeline.reconstruct_without_metals import xyz_to_omol_no_metal_state
from molgr.fallback.state import MetalCandidateStateMachine, make_metal_candidate_state
from molgr.fallback.utils.dataclasses import MetalAtomPosition
from molgr.fallback.utils.force_field import organic_force_field_evaluation
from molgr.fallback.utils.metals.preparation import combine_metal_with_omol


RDLogger.DisableLog("rdApp.*")  # type: ignore


def _curated_case(case_idx: int) -> dict[str, object]:
    csv_path = Path(__file__).with_name("test_cases.csv")
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    smiles = rows[case_idx - 1]["smiles"]

    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    mol_h = Chem.AddHs(mol)
    embed_code = rdDistGeom.EmbedMolecule(mol_h)  # pyright: ignore[reportCallIssue]
    assert int(embed_code) == 0

    total_charge = 0
    total_radical_electrons = 0
    for atom in mol_h.GetAtoms():  # pyright: ignore[reportCallIssue]
        total_charge += int(atom.GetFormalCharge())
        total_radical_electrons += int(atom.GetNumRadicalElectrons())

    return {
        "xyz_block": Chem.MolToXYZBlock(mol_h),
        "total_charge": total_charge,
        "total_radical_electrons": total_radical_electrons,
    }


def test_reconstruction_state_uses_organic_force_field_policy() -> None:
    case = _curated_case(1)

    state = xyz_to_omol_no_metal_state(
        case["xyz_block"],
        case["total_charge"],
        case["total_radical_electrons"],
    )

    assert state is not None
    evaluation = organic_force_field_evaluation(state)

    assert state.score("organic_core") == pytest.approx(evaluation.energy_kj_mol)
    assert state.score("full") == pytest.approx(evaluation.energy_kj_mol)
    assert state.metadata["organic_core_score"] == pytest.approx(evaluation.energy_kj_mol)
    assert state.metadata["score"] == pytest.approx(evaluation.energy_kj_mol)
    assert state.metadata["force_field_requested"] == "auto"
    assert state.metadata["force_field_resolved_force_field"] == "uff"
    assert isinstance(state.metadata["force_field_score_key"], tuple)


def test_reconstruction_state_no_longer_exposes_breakdown_api() -> None:
    case = _curated_case(1)

    state = xyz_to_omol_no_metal_state(
        case["xyz_block"],
        case["total_charge"],
        case["total_radical_electrons"],
    )

    assert state is not None
    assert not hasattr(state, "score_breakdown")


def test_reconstruction_state_no_longer_exposes_rescore_api() -> None:
    case = _curated_case(1)

    state = xyz_to_omol_no_metal_state(
        case["xyz_block"],
        case["total_charge"],
        case["total_radical_electrons"],
    )

    assert state is not None
    assert not hasattr(state, "rescore")


def test_organic_force_field_evaluation_rejects_metal_containing_input() -> None:
    organic = _curated_case(1)
    state = xyz_to_omol_no_metal_state(
        organic["xyz_block"],
        organic["total_charge"],
        organic["total_radical_electrons"],
    )
    assert state is not None

    anchor = state.omol.atoms[0].OBAtom
    combined = combine_metal_with_omol(
        state.omol,
        (
            MetalAtomPosition(
                idx=1,
                symbol="Li",
                element_idx=3,
                valence=1,
                radical_num=0,
                position_x=anchor.GetX(),
                position_y=anchor.GetY(),
                position_z=anchor.GetZ() + 1.8,
            ),
        ),
    )

    with pytest.raises(ValueError, match="metal-free molecules"):
        organic_force_field_evaluation(combined)


def test_metal_candidate_uses_organic_force_field_only() -> None:
    case = _curated_case(1)

    no_metal_state = xyz_to_omol_no_metal_state(
        case["xyz_block"],
        case["total_charge"],
        case["total_radical_electrons"],
    )

    assert no_metal_state is not None
    anchor = no_metal_state.omol.atoms[0].OBAtom
    candidate = make_metal_candidate_state(
        (),
        [
            MetalAtomPosition(
                idx=1,
                symbol="Li",
                element_idx=3,
                valence=1,
                radical_num=0,
                position_x=anchor.GetX(),
                position_y=anchor.GetY(),
                position_z=anchor.GetZ() + 1.8,
            )
        ],
        no_metal_charge_target=no_metal_state.total_charge,
        no_metal_radical_target=no_metal_state.total_radical_electrons,
        combination_index=0,
    )
    machine = MetalCandidateStateMachine.from_candidate_state(candidate)
    machine.set_no_metal_state("bind_no_metal_state", no_metal_state)
    candidate = machine.freeze()

    organic_score = no_metal_state.score("organic_core")
    score = candidate.combined_score()

    assert candidate.combined_omol is None
    assert score == pytest.approx(organic_score)
    assert candidate.metadata["force_field_energy"] == pytest.approx(organic_score)
    assert candidate.metadata["score"] == pytest.approx(organic_score)
    assert candidate.metadata["force_field_requested"] == "auto"
    assert candidate.metadata["force_field_resolved_force_field"] == "uff"
    assert "selection_score_profile" not in candidate.metadata
    assert "force_field_base_energy" not in candidate.metadata
    assert "force_field_electrostatic_energy" not in candidate.metadata
    assert "force_field_electrostatic_pair_count" not in candidate.metadata
    assert "post_reinsertion_score" not in candidate.metadata
    assert not hasattr(candidate, "score_breakdown")


def test_metal_candidate_no_longer_exposes_rescore_api() -> None:
    case = _curated_case(1)

    no_metal_state = xyz_to_omol_no_metal_state(
        case["xyz_block"],
        case["total_charge"],
        case["total_radical_electrons"],
    )

    assert no_metal_state is not None
    anchor = no_metal_state.omol.atoms[0].OBAtom
    candidate = make_metal_candidate_state(
        (),
        [
            MetalAtomPosition(
                idx=1,
                symbol="Li",
                element_idx=3,
                valence=1,
                radical_num=0,
                position_x=anchor.GetX(),
                position_y=anchor.GetY(),
                position_z=anchor.GetZ() + 1.8,
            )
        ],
        no_metal_charge_target=no_metal_state.total_charge,
        no_metal_radical_target=no_metal_state.total_radical_electrons,
        combination_index=0,
    )
    machine = MetalCandidateStateMachine.from_candidate_state(candidate)
    machine.set_no_metal_state("bind_no_metal_state", no_metal_state)
    candidate = machine.freeze()

    assert not hasattr(candidate, "rescore")
