# pyright: reportCallIssue=false
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

import pytest
from rdkit import Chem

from molgr.fallback.utils.consts import NON_METAL_DICT
from molgr.interface import xyz_to_rdmol
from molgr.utils.converter import get_atom_unpaired_electrons
from molgr.utils.equivalence import check_equivalence
from scripts.reconstruction_trace import (
    TraceInputCase,
    load_review_fixture_cases,
    trace_reconstruction_case,
)


_FIXTURE_ROOT = Path(__file__).parent / "data" / "reviewed" / "tmqmg"
_MANIFEST_PATH = _FIXTURE_ROOT / "manifest.json"


def _manifest_records() -> list[dict[str, Any]]:
    payload = json.loads(_MANIFEST_PATH.read_text(encoding="utf-8"))
    assert payload["schema_version"] == 1
    return payload["fixtures"]


def test_review_fixture_manifest_pins_source_dataset() -> None:
    payload = json.loads(_MANIFEST_PATH.read_text(encoding="utf-8"))

    assert payload["source_dataset"] == {
        "name": "tmQMg",
        "properties_file": "tmQMg_properties_and_targets.csv",
        "properties_sha256": "3920c1c8f4ec81bc8e44b8d0256a7da1e36c8805c3c0adfd47e50c46e633f473",
        "publication_doi": "10.1039/D2DD00129B",
        "repository": "https://github.com/uiocompcat/tmQMg",
        "revision": "e1dc9887b8f20a217a1db6ca972d726bcbaab45b",
        "xyz_file": "tmQMg_xyz.zip",
        "xyz_sha256": "e0d15a70bcba294717cd9f9792e7fac99ef0c5c61c3a6e08dcc8a8643f53660a",
    }


def _molecule_electronic_state(mol: Chem.Mol) -> tuple[int, int]:
    total_charge = 0
    total_radical_electrons = 0
    for atom in mol.GetAtoms():  # pyright: ignore[reportCallIssue]
        total_charge += int(atom.GetFormalCharge())
        total_radical_electrons += get_atom_unpaired_electrons(atom)
    return total_charge, total_radical_electrons


def _normalize_fixture_answer(mol: Chem.Mol) -> Chem.Mol:
    """Drop coordination bonds and metal stereochemistry from an answer graph."""

    answer = Chem.RWMol(Chem.Mol(mol))
    coordination_bonds = [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for bond in answer.GetBonds()
        if bond.GetBondType() == Chem.BondType.DATIVE
    ]
    for begin_idx, end_idx in coordination_bonds:
        answer.RemoveBond(begin_idx, end_idx)
    for atom in answer.GetAtoms():
        if int(atom.GetAtomicNum()) not in NON_METAL_DICT:
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

    normalized = answer.GetMol()
    normalized.UpdatePropertyCache(strict=False)
    return normalized


def test_review_fixture_manifest_is_complete_and_unique() -> None:
    records = _manifest_records()

    assert records
    assert len({record["case_id"] for record in records}) == len(records)
    assert {record["kind"] for record in records} <= {
        "approved_graph",
        "manual_reference",
        "reference_graph",
    }
    for record in records:
        assert (_FIXTURE_ROOT / str(record["structure_file"])).is_file()
        assert int(record["spin_multiplicity"]) == int(record["total_radical_electrons"]) + 1
        assert {"reviewer", "notes", "updated_at"}.isdisjoint(record)


def test_fixture_answers_ignore_coordination_bonds_and_metal_stereochemistry() -> None:
    path = _FIXTURE_ROOT / "approved_graph" / "ABEGOD.sdf"
    source = next(
        mol
        for mol in Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
        if mol is not None
    )

    assert any(bond.GetBondType() == Chem.BondType.DATIVE for bond in source.GetBonds())
    assert any(
        atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        for atom in source.GetAtoms()
        if int(atom.GetAtomicNum()) not in NON_METAL_DICT
    )

    answer = _normalize_fixture_answer(source)

    assert all(bond.GetBondType() != Chem.BondType.DATIVE for bond in answer.GetBonds())
    assert all(
        atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
        for atom in answer.GetAtoms()
        if int(atom.GetAtomicNum()) not in NON_METAL_DICT
    )
    equivalent, info = check_equivalence(answer, source, use_chirality=False)
    assert equivalent, info.reason


def test_reviewed_approved_sdf_graphs_survive_reconstruction() -> None:
    records = [record for record in _manifest_records() if record["kind"] == "approved_graph"]

    for record in records:
        path = _FIXTURE_ROOT / str(record["structure_file"])
        expected = next(
            mol
            for mol in Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
            if mol is not None
        )
        assert expected.GetNumConformers() == 1
        total_charge = int(expected.GetProp("TOTAL_CHARGE"))
        total_radical_electrons = int(expected.GetProp("TOTAL_RADICAL_ELECTRONS"))
        spin_multiplicity = int(expected.GetProp("SPIN_MULTIPLICITY"))
        assert _molecule_electronic_state(expected) == (
            total_charge,
            total_radical_electrons,
        )
        assert spin_multiplicity == total_radical_electrons + 1
        expected_answer = _normalize_fixture_answer(expected)

        rebuilt = xyz_to_rdmol(
            Chem.MolToXYZBlock(expected),
            total_charge=total_charge,
            spin_multiplicity=spin_multiplicity,
            backend="cpp",
            make_dative_bonds=True,
            make_stereochemistry=True,
        )
        equivalent, info = check_equivalence(expected_answer, rebuilt, use_chirality=False)
        assert equivalent, f"{record['case_id']}: {info.reason}"


def test_review_manifest_is_a_live_trace_input_source() -> None:
    cases = load_review_fixture_cases(_MANIFEST_PATH)
    records = _manifest_records()

    assert {case.id for case in cases} == {str(record["case_id"]) for record in records}
    assert all(case.xyz_source == "review_fixture" for case in cases)
    assert all(case.fixture_kind for case in cases)
    assert all(case.fixture_structure_file for case in cases)


@pytest.mark.parametrize("case_id", ["ABAGAM", "ABEGOD"])
def test_review_fixture_trace_matches_approved_answer(case_id: str) -> None:
    case = next(
        case for case in load_review_fixture_cases(_MANIFEST_PATH, [case_id]) if case.id == case_id
    )
    assert isinstance(case, TraceInputCase)

    trace = trace_reconstruction_case(case, score_all_candidates=False)
    fixture_check = trace["review_fixture"]

    assert fixture_check["equivalent"] is True, fixture_check
    assert fixture_check["equivalence_method"]
    assert fixture_check["trace_smiles"]


@pytest.mark.parametrize("backend", ["cpp", "python"])
def test_reviewed_xyz_smiles_reference_graphs_survive_reconstruction(
    backend: Literal["cpp", "python"],
) -> None:
    records = [
        record
        for record in _manifest_records()
        if record["kind"] in {"manual_reference", "reference_graph"}
    ]

    for record in records:
        smiles_field = (
            "approved_smiles" if record["kind"] == "manual_reference" else "reference_smiles"
        )
        reference = Chem.MolFromSmiles(str(record[smiles_field]), sanitize=False)
        assert reference is not None
        reference.UpdatePropertyCache(strict=False)
        expected_answer = _normalize_fixture_answer(reference)
        rebuilt = xyz_to_rdmol(
            (_FIXTURE_ROOT / str(record["structure_file"])).read_text(encoding="utf-8"),
            total_charge=int(record["total_charge"]),
            spin_multiplicity=int(record["spin_multiplicity"]),
            backend=backend,
            make_dative_bonds=True,
            make_stereochemistry=True,
        )
        equivalent, info = check_equivalence(expected_answer, rebuilt, use_chirality=False)
        assert equivalent, f"{record['case_id']}: {info.reason}"
