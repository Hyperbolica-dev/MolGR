# pyright: reportCallIssue=false
from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Geometry import Point3D

from tools.molgr_review.fixture_builder import reconstruct_case_mol
from tools.molgr_review.triage_mapping import map_candidate_reference_xyz, parse_xyz_atoms


ROOT = Path(__file__).resolve().parents[3]
BASELINE_DB = ROOT / ".local/molgr_review_backups/2026-08-13_pre_manual_review/review.sqlite"
XYZ_DIR = ROOT / ".local/tmQMg/data/xyz"


def graph(
    symbols: list[str], bonds: list[tuple[int, int, Chem.BondType]], coordinates=None
) -> Chem.Mol:
    editable = Chem.RWMol()
    for symbol in symbols:
        editable.AddAtom(Chem.Atom(symbol))
    for begin, end, kind in bonds:
        editable.AddBond(begin, end, kind)
    mol = editable.GetMol()
    mol.UpdatePropertyCache(strict=False)
    if coordinates is None:
        coordinates = [(float(index), 0.0, 0.0) for index in range(len(symbols))]
    conformer = Chem.Conformer(len(symbols))
    for index, coordinate in enumerate(coordinates):
        conformer.SetAtomPosition(index, Point3D(*coordinate))
    mol.AddConformer(conformer)
    return mol


def xyz_atoms(mol: Chem.Mol):
    conformer = mol.GetConformer()
    return [
        (
            atom.GetSymbol(),
            tuple(conformer.GetAtomPosition(atom.GetIdx())),
        )
        for atom in mol.GetAtoms()
    ]


def test_unique_metal_edge_mapping_ignores_irrelevant_atom_order() -> None:
    candidate = graph(
        ["Co", "S", "N", "C"],
        [
            (0, 2, Chem.BondType.DATIVE),
            (2, 3, Chem.BondType.SINGLE),
            (3, 1, Chem.BondType.SINGLE),
        ],
    )
    reference = graph(
        ["C", "N", "S", "Co"],
        [
            (3, 1, Chem.BondType.DATIVE),
            (1, 0, Chem.BondType.SINGLE),
            (0, 2, Chem.BondType.SINGLE),
            (2, 3, Chem.BondType.DATIVE),
        ],
    )

    result = map_candidate_reference_xyz(candidate, reference, xyz_atoms(candidate))

    assert result.confidence == "unique_graph_mapping"
    assert result.mapping_signature_count == 1
    assert result.reference_to_candidate[3] == 0
    assert result.reference_to_candidate[2] == 1


def test_multiple_sulfur_correspondences_are_ambiguous() -> None:
    candidate = graph(
        ["Co", "S", "S", "C"],
        [(0, 3, Chem.BondType.SINGLE), (3, 1, Chem.BondType.SINGLE), (3, 2, Chem.BondType.SINGLE)],
    )
    reference = graph(
        ["Co", "S", "S", "C"],
        [
            (0, 3, Chem.BondType.SINGLE),
            (3, 1, Chem.BondType.SINGLE),
            (3, 2, Chem.BondType.SINGLE),
            (0, 1, Chem.BondType.DATIVE),
        ],
    )

    result = map_candidate_reference_xyz(candidate, reference, xyz_atoms(candidate))

    assert result.confidence == "ambiguous"
    assert result.mapping_signature_count == 2


def test_symmetric_oxygen_correspondences_are_ambiguous() -> None:
    candidate = graph(
        ["Co", "N", "N", "O", "O", "C"],
        [
            (0, 1, Chem.BondType.DATIVE),
            (0, 2, Chem.BondType.DATIVE),
            (1, 5, Chem.BondType.SINGLE),
            (2, 5, Chem.BondType.SINGLE),
            (3, 5, Chem.BondType.SINGLE),
            (4, 5, Chem.BondType.SINGLE),
        ],
    )
    reference = graph(
        ["Co", "N", "N", "O", "O", "C"],
        [
            (0, 1, Chem.BondType.DATIVE),
            (0, 2, Chem.BondType.DATIVE),
            (0, 3, Chem.BondType.DATIVE),
            (1, 5, Chem.BondType.SINGLE),
            (2, 5, Chem.BondType.SINGLE),
            (3, 5, Chem.BondType.SINGLE),
            (4, 5, Chem.BondType.SINGLE),
        ],
    )

    result = map_candidate_reference_xyz(candidate, reference, xyz_atoms(candidate))

    assert result.confidence == "ambiguous"
    assert result.mapping_signature_count == 2


def test_truncated_enumeration_is_never_strong() -> None:
    candidate = graph(
        ["Co", "S", "S", "C"],
        [(0, 3, Chem.BondType.SINGLE), (3, 1, Chem.BondType.SINGLE), (3, 2, Chem.BondType.SINGLE)],
    )
    reference = graph(
        ["Co", "S", "S", "C"],
        [
            (0, 3, Chem.BondType.SINGLE),
            (3, 1, Chem.BondType.SINGLE),
            (3, 2, Chem.BondType.SINGLE),
            (0, 1, Chem.BondType.DATIVE),
        ],
    )

    result = map_candidate_reference_xyz(candidate, reference, xyz_atoms(candidate), max_matches=1)

    assert result.enumeration_truncated is True
    assert result.confidence == "ambiguous"


def _real_case(case_id: str):
    if not BASELINE_DB.exists() or not (XYZ_DIR / f"{case_id}.xyz").exists():
        pytest.skip("local frozen tmQMg review data is unavailable")
    connection = sqlite3.connect(f"file:{BASELINE_DB}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        row = dict(
            connection.execute("SELECT * FROM cases WHERE case_id = ?", (case_id,)).fetchone()
        )
    finally:
        connection.close()
    candidate = reconstruct_case_mol(row, xyz_dir=XYZ_DIR)
    reference = Chem.MolFromSmiles(row["reference_smiles"])
    assert reference is not None
    xyz = parse_xyz_atoms((XYZ_DIR / f"{case_id}.xyz").read_text(encoding="utf-8"))
    return candidate, reference, xyz


def test_abegod_real_mapping_regression() -> None:
    candidate, reference, xyz = _real_case("ABEGOD")

    result = map_candidate_reference_xyz(candidate, reference, xyz)

    assert result.confidence == "unique_graph_mapping"
    assert result.mapping_signature_count == 1
    assert result.equal_best_mapping_count == 4
    assert result.reference_to_candidate[20] == 20  # common S-
    assert result.reference_to_candidate[16] == 21  # oxygen-bearing reference-only S
    signature_text = json.dumps(result.decision_relevant_signatures)
    assert '"metal_bond", [16, 21], false, true' in signature_text.lower()
    conformer = candidate.GetConformer()
    assert (
        conformer.GetAtomPosition(16) - conformer.GetAtomPosition(20)
    ).Length() == pytest.approx(2.187174678, abs=1e-9)
    assert (
        conformer.GetAtomPosition(16) - conformer.GetAtomPosition(21)
    ).Length() == pytest.approx(2.068788633, abs=1e-9)


def test_abatec_real_mapping_regression_stays_ambiguous_when_truncated() -> None:
    candidate, reference, xyz = _real_case("ABATEC")

    result = map_candidate_reference_xyz(candidate, reference, xyz)

    assert result.enumeration_truncated is True
    assert result.confidence == "ambiguous"
    conformer = candidate.GetConformer()
    assert (
        conformer.GetAtomPosition(40) - conformer.GetAtomPosition(42)
    ).Length() == pytest.approx(1.272291181, abs=1e-9)
    assert (conformer.GetAtomPosition(0) - conformer.GetAtomPosition(42)).Length() == pytest.approx(
        1.888161076, abs=1e-9
    )
