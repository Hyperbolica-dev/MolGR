from __future__ import annotations

import io
import json
import sqlite3
import sys
import threading
from http.client import HTTPConnection
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D


APP_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_DIR))

import fixture_builder  # noqa: E402
import server as review_server  # noqa: E402
from fixture_builder import (  # noqa: E402
    ACCEPT_BOTH_STATUS,
    _sdf_text,
    case_electronic_state,
    load_fixture_records,
    sync_review_fixture,
)
from server import ReviewHandler, ReviewServer  # noqa: E402

from molgr.utils.converter import (  # noqa: E402
    METAL_UNPAIRED_ELECTRONS_PROP,
    get_atom_unpaired_electrons,
)


def test_review_page_exposes_fixture_removal_action() -> None:
    html = (APP_DIR / "static" / "index.html").read_text(encoding="utf-8")
    javascript = (APP_DIR / "static" / "app.js").read_text(encoding="utf-8")
    stylesheet = (APP_DIR / "static" / "style.css").read_text(encoding="utf-8")

    assert 'id="removeFixture"' in html
    assert 'id="languageToggle"' in html
    assert 'data-i18n="languageToggle"' in html
    assert 'id="openTrace"' in html
    assert '$("removeFixture").addEventListener("click", removeCurrentFixture)' in javascript
    assert 'localStorage.setItem("moleculeReviewLanguage", state.language)' in javascript
    assert "function applyLanguage()" in javascript
    assert "function renderKetcherStatus()" in javascript
    assert 'setKetcherStatus("ketcherReady")' in javascript
    assert "window.open(" in javascript
    assert "`/trace/${encodeURIComponent(state.current.case_id)}`" in javascript
    assert 'await saveReview("needs_followup")' in javascript
    assert 'data-status="accept_both"' in html
    assert 'accept_both: "接受两者"' in javascript
    assert 'id="imageLightbox"' in html
    assert 'class="viewer-grid"' in html
    assert 'loading="lazy"' in html
    assert html.index('id="reviewControls"') < html.index('id="mainLayout"')
    assert html.index('id="mainLayout"') < html.index('id="reviewPanel"')
    assert "function openImageLightbox(box)" in javascript
    assert "dialog.showModal()" in javascript
    assert "if (workspace) workspace.scrollTop = 0;" in javascript
    assert '<code title="${escapeHtml(value)}">' in javascript
    assert (
        "grid-template-columns: clamp(260px, 34vw, var(--sidebar-width)) "
        "minmax(0, 1fr);" in stylesheet
    )
    assert "min-width: 1200px;" in stylesheet
    assert "@media (max-width: 1700px)" in stylesheet
    assert "contain: layout paint;" in stylesheet
    assert "flex-direction: column;" in stylesheet
    assert "grid-template-columns: repeat(2, minmax(0, 1fr));" in stylesheet
    assert "text-overflow: ellipsis;" in stylesheet
    assert "body {\n    display: block;" not in stylesheet
    assert 'class="render-kind' not in html
    assert "primaryKind" not in javascript
    assert "secondaryKind" not in javascript
    assert 'querySelectorAll(".render-kind")' not in javascript


def test_review_2d_render_omits_explicit_h_without_mutating_electronic_state(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    charged_radical = Chem.AddHs(Chem.MolFromSmiles("[NH3+][CH2]"))
    assert charged_radical is not None
    original_atom_count = charged_radical.GetNumAtoms()
    original_h_count = sum(
        atom.GetAtomicNum() == 1
        for atom in charged_radical.GetAtoms()  # pyright: ignore[reportCallIssue]
    )
    original_heavy_charges = [
        atom.GetFormalCharge()
        for atom in charged_radical.GetAtoms()  # pyright: ignore[reportCallIssue]
        if atom.GetAtomicNum() != 1
    ]
    original_heavy_radicals = [
        atom.GetNumRadicalElectrons()
        for atom in charged_radical.GetAtoms()  # pyright: ignore[reportCallIssue]
        if atom.GetAtomicNum() != 1
    ]
    captured: dict[str, Chem.Mol] = {}

    class SvgResult:
        data = "<svg>prepared</svg>"

    def capture_render(mol: Chem.Mol, **kwargs: object) -> SvgResult:
        captured["mol"] = Chem.Mol(mol)
        return SvgResult()

    import rdkit_dof

    monkeypatch.setattr(rdkit_dof, "MolToDofImage", capture_render)
    assert review_server._render_mol_svg(charged_radical, legend="test") == ("<svg>prepared</svg>")

    rendered = captured["mol"]
    assert original_h_count > 0
    assert charged_radical.GetNumAtoms() == original_atom_count
    assert sum(atom.GetAtomicNum() == 1 for atom in rendered.GetAtoms()) == 0  # pyright: ignore[reportCallIssue]
    assert sum(atom.GetFormalCharge() for atom in charged_radical.GetAtoms()) == sum(  # pyright: ignore[reportCallIssue]
        atom.GetFormalCharge()
        for atom in rendered.GetAtoms()  # pyright: ignore[reportCallIssue]
    )
    assert [
        atom.GetFormalCharge()
        for atom in rendered.GetAtoms()  # pyright: ignore[reportCallIssue]
    ] == original_heavy_charges
    assert [
        atom.GetNumRadicalElectrons()
        for atom in rendered.GetAtoms()  # pyright: ignore[reportCallIssue]
    ] == original_heavy_radicals
    assert sum(original_heavy_radicals) == 1


def test_review_2d_render_localizes_only_requested_heavy_atom(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule = Chem.AddHs(Chem.MolFromSmiles("CC"))
    assert molecule is not None
    captured: dict[str, object] = {}

    class SvgResult:
        data = "<svg>localized</svg>"

    def capture_render(mol: Chem.Mol, **kwargs: object) -> SvgResult:
        captured["mol"] = Chem.Mol(mol)
        captured["highlight_atoms"] = kwargs.get("highlightAtoms")
        return SvgResult()

    import rdkit_dof

    monkeypatch.setattr(rdkit_dof, "MolToDofImage", capture_render)
    review_server._render_mol_svg(
        molecule,
        legend="test",
        atom_notes={1: "XYZ #19 · H#20 attached here", 2: "H#20"},
    )

    rendered = captured["mol"]
    assert isinstance(rendered, Chem.Mol)
    assert rendered.GetNumAtoms() == 2
    assert captured["highlight_atoms"] == [1]
    assert rendered.GetAtomWithIdx(1).GetProp("atomNote") == "XYZ #19 · H#20 attached here"
    assert not molecule.GetAtomWithIdx(1).HasProp("atomNote")


def test_review_2d_hydrogen_mode_keeps_all_source_h_without_generating_h(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    base = Chem.MolFromSmiles("CN")
    assert base is not None
    molecule = Chem.AddHs(base)
    source_atoms = molecule.GetNumAtoms()
    disputed_h = next(
        atom.GetIdx()
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() == 1 and atom.GetNeighbors()[0].GetAtomicNum() == 6  # pyright: ignore[reportIndexIssue]
    )
    captured: dict[str, Chem.Mol] = {}

    class SvgResult:
        data = "<svg>hydrogen</svg>"

    def capture_render(mol: Chem.Mol, **kwargs: object) -> SvgResult:
        captured["mol"] = Chem.Mol(mol)
        return SvgResult()

    import rdkit_dof

    monkeypatch.setattr(rdkit_dof, "MolToDofImage", capture_render)
    review_server._render_mol_svg(
        molecule,
        legend="H",
        atom_notes={disputed_h: f"H · #{disputed_h}"},
        show_hydrogens=True,
    )

    assert molecule.GetNumAtoms() == source_atoms
    assert captured["mol"].GetNumAtoms() == source_atoms
    rendered_h = [
        atom
        for atom in captured["mol"].GetAtoms()  # pyright: ignore[reportCallIssue]
        if atom.GetAtomicNum() == 1
    ]
    assert len(rendered_h) == sum(
        atom.GetAtomicNum() == 1
        for atom in molecule.GetAtoms()  # pyright: ignore[reportCallIssue]
    )
    assert any(
        atom.HasProp("atomNote") and atom.GetProp("atomNote") == f"H · #{disputed_h}"
        for atom in rendered_h
    )


def test_reference_xyz_uses_source_coordinates_and_requires_reliable_mapping(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    editable = Chem.RWMol()
    for symbol in ["C", "H", "H", "H", "H"]:
        editable.AddAtom(Chem.Atom(symbol))
    for hydrogen in range(1, 5):
        editable.AddBond(0, hydrogen, Chem.BondType.SINGLE)
    candidate = editable.GetMol()
    coordinates = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, -1.0, 0.0),
    ]
    conformer = Chem.Conformer(5)
    for index, coordinate in enumerate(coordinates):
        conformer.SetAtomPosition(index, Point3D(*coordinate))
    candidate.AddConformer(conformer)
    candidate.UpdatePropertyCache(strict=False)
    reference = Chem.MolFromSmiles("C")
    assert reference is not None
    xyz_path = tmp_path / "METHANE.xyz"
    xyz_path.write_text(
        "5\nMETHANE\n"
        + "\n".join(
            f"{symbol} {x} {y} {z}"
            for symbol, (x, y, z) in zip(["C", "H", "H", "H", "H"], coordinates)
        )
        + "\n",
        encoding="utf-8",
    )
    case = {"xyz_path": str(xyz_path)}

    result, confidence, error = review_server._reference_xyz_mol(
        case,
        {"mapping_confidence": "unique_graph_mapping"},
        candidate,
        reference,
        None,
    )
    assert result is not None
    assert confidence == "unique_graph_mapping"
    assert error == ""
    assert result.GetNumAtoms() == 5
    assert result.GetNumBonds() == 4
    assert tuple(result.GetConformer().GetAtomPosition(1)) == coordinates[1]

    class AmbiguousMapping:
        confidence = "ambiguous"
        reference_to_candidate: dict[int, int] = {}

    monkeypatch.setattr(
        review_server,
        "map_candidate_reference_xyz",
        lambda *_args, **_kwargs: AmbiguousMapping(),
    )
    unavailable, ambiguous_confidence, reason = review_server._reference_xyz_mol(
        case,
        {"mapping_confidence": "unique_graph_mapping"},
        candidate,
        reference,
        None,
    )
    assert unavailable is None
    assert ambiguous_confidence == "ambiguous"
    assert reason == "atom_correspondence_not_reliable"

    class RepresentativeAmbiguousMapping:
        confidence = "ambiguous"
        reference_to_candidate = {0: 0}
        enumeration_truncated = True

    representative, representative_confidence, representative_error = (
        review_server._reference_xyz_mol(
            case,
            None,
            candidate,
            reference,
            None,
            RepresentativeAmbiguousMapping(),
        )
    )
    assert representative is not None
    assert representative_confidence == "ambiguous"
    assert representative_error == ""
    assert representative.GetNumBonds() == 4
    assert tuple(representative.GetConformer().GetAtomPosition(1)) == coordinates[1]


def test_mapped_xyz_comparison_preserves_real_edges_and_exposes_donor_pair() -> None:
    candidate_editable = Chem.RWMol()
    for symbol in ["Mo", "C", "O", "O"]:
        candidate_editable.AddAtom(Chem.Atom(symbol))
    candidate_editable.AddBond(0, 2, Chem.BondType.DATIVE)
    candidate_editable.AddBond(1, 2, Chem.BondType.SINGLE)
    candidate_editable.AddBond(1, 3, Chem.BondType.DOUBLE)
    candidate = candidate_editable.GetMol()
    conformer = Chem.Conformer(4)
    for index, point in enumerate([(0, 0, 0), (2, 0, 0), (1, 0, 0), (3, 0, 0)]):
        conformer.SetAtomPosition(index, Point3D(*point))
    candidate.AddConformer(conformer)
    candidate.UpdatePropertyCache(strict=False)

    reference_editable = Chem.RWMol()
    for symbol in ["O", "C", "O", "Mo"]:
        reference_editable.AddAtom(Chem.Atom(symbol))
    reference_editable.AddBond(3, 0, Chem.BondType.DATIVE)
    reference_editable.AddBond(1, 0, Chem.BondType.SINGLE)
    reference_editable.AddBond(1, 2, Chem.BondType.DOUBLE)
    reference = reference_editable.GetMol()
    reference.UpdatePropertyCache(strict=False)

    comparison = review_server._mapped_coordination_comparison(
        candidate,
        reference,
        {3: 0, 1: 1, 0: 2, 2: 3},
    )
    assert candidate.GetBondBetweenAtoms(0, 2) is not None
    assert candidate.GetBondBetweenAtoms(0, 3) is None
    assert reference.GetBondBetweenAtoms(3, 0) is not None
    assert reference.GetBondBetweenAtoms(3, 2) is None
    assert comparison["coordination_edges"] == [
        {
            "presence": "common",
            "metal_element": "Mo",
            "donor_element": "O",
            "candidate_metal_xyz_index": 0,
            "candidate_donor_xyz_index": 2,
            "reference_metal_atom_index": 3,
            "reference_donor_atom_index": 0,
            "distance": 1.0,
            "mapped_ligand_group": [
                {
                    "element": "O",
                    "candidate_xyz_index": 2,
                    "reference_atom_index": 0,
                    "role": "donor",
                },
                {
                    "element": "O",
                    "candidate_xyz_index": 3,
                    "reference_atom_index": 2,
                    "role": "mapped_equivalent",
                },
            ],
            "interpretation": "mapped_donor_preserved",
        }
    ]


@pytest.mark.parametrize(
    ("attributes", "expected"),
    [
        (
            {
                "timeout": False,
                "enumeration_truncated": False,
                "mapping_signature_count": 2,
                "equal_best_mapping_count": 2,
            },
            ("multiple_valid_mappings", "multiple_equally_valid_atom_mappings"),
        ),
        (
            {
                "timeout": False,
                "enumeration_truncated": True,
                "mapping_signature_count": 1,
                "equal_best_mapping_count": 8,
            },
            (
                "mapping_enumeration_truncated",
                "mapping_enumeration_truncated_before_unique_correspondence",
            ),
        ),
        (
            {
                "timeout": True,
                "enumeration_truncated": False,
                "mapping_signature_count": 1,
                "equal_best_mapping_count": 1,
            },
            ("mapping_timeout", "mapping_timeout_before_unique_correspondence"),
        ),
    ],
)
def test_reference_xyz_ambiguity_uses_existing_mapping_status(
    attributes: dict[str, object], expected: tuple[str, str]
) -> None:
    mapping = SimpleNamespace(confidence="ambiguous", **attributes)
    assert review_server._mapping_ambiguity_details(mapping) == expected


def test_reference_xyz_ambiguity_reports_decision_relevant_xyz_location() -> None:
    mapping = SimpleNamespace(
        decision_relevant_signatures=(
            (("metal_bond", (4, 9), True, False, "dative", "none"),),
            (("metal_bond", (4, 10), True, False, "dative", "none"),),
        ),
        enumeration_truncated=False,
        timeout=False,
    )
    assert review_server._mapping_ambiguity_locations(mapping) == {
        "affected_xyz_atoms": [4, 9, 10],
        "alternatives": [
            {
                "alternative": 1,
                "differences": [
                    {
                        "kind": "metal_bond",
                        "xyz_atoms": [4, 9],
                        "candidate_present": True,
                        "reference_present": False,
                        "candidate_bond": "dative",
                        "reference_bond": "none",
                    }
                ],
            },
            {
                "alternative": 2,
                "differences": [
                    {
                        "kind": "metal_bond",
                        "xyz_atoms": [4, 10],
                        "candidate_present": True,
                        "reference_present": False,
                        "candidate_bond": "dative",
                        "reference_bond": "none",
                    }
                ],
            },
        ],
        "location_proven": True,
    }


def test_review_2d_render_preserves_hydride_annotations() -> None:
    smiles = (
        "CC(C)(C)N->[RuH+]123(<-n4c(C(F)(F)F)cc(C(F)(F)F)n4[BH2-]"
        "n4nc(C(F)(F)F)cc4C(F)(F)F)<-C4=C->1CCC->2=C->3CC4"
    )
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    with_explicit_h = Chem.AddHs(mol)

    def heavy_bond_signature(value: Chem.Mol) -> list[tuple[int, int, str]]:
        heavy_indices = [
            atom.GetIdx()
            for atom in value.GetAtoms()  # pyright: ignore[reportCallIssue]
            if atom.GetAtomicNum() != 1
        ]
        positions = {atom_index: position for position, atom_index in enumerate(heavy_indices)}
        return sorted(
            (
                positions[bond.GetBeginAtomIdx()],
                positions[bond.GetEndAtomIdx()],
                str(bond.GetBondType()),
            )
            for bond in value.GetBonds()  # pyright: ignore[reportCallIssue]
            if bond.GetBeginAtomIdx() in positions and bond.GetEndAtomIdx() in positions
        )

    rendered = review_server._prepare_review_2d_mol(with_explicit_h)

    assert sum(atom.GetAtomicNum() == 1 for atom in rendered.GetAtoms()) == 0  # pyright: ignore[reportCallIssue]
    rendered_smiles = Chem.MolToSmiles(rendered, canonical=True)
    assert "[RuH+]" in rendered_smiles
    assert "[BH2-]" in rendered_smiles
    assert sum(atom.GetFormalCharge() for atom in with_explicit_h.GetAtoms()) == sum(  # pyright: ignore[reportCallIssue]
        atom.GetFormalCharge()
        for atom in rendered.GetAtoms()  # pyright: ignore[reportCallIssue]
    )
    assert heavy_bond_signature(rendered) == heavy_bond_signature(with_explicit_h)
    assert sum(atom.GetNumRadicalElectrons() for atom in rendered.GetAtoms()) == 0  # pyright: ignore[reportCallIssue]


def test_case_payload_hides_only_empty_or_exact_duplicate_organic_graphs() -> None:
    duplicate = review_server._row_dict(
        {
            "candidate_smiles": " C ",
            "candidate_organic_smiles": "C",
            "reference_smiles": "N",
            "reference_organic_smiles": " N ",
        }
    )
    assert duplicate is not None
    assert duplicate["available_render_kinds"] == ["candidate", "reference"]

    distinct = review_server._row_dict(
        {
            "candidate_smiles": "C",
            "candidate_organic_smiles": "CC",
            "reference_smiles": "N",
            "reference_organic_smiles": "NN",
        }
    )
    assert distinct is not None
    assert distinct["available_render_kinds"] == [
        "candidate",
        "reference",
        "candidate_organic",
        "reference_organic",
    ]


def test_fixture_manifest_rejects_unconfirmed_answers(tmp_path: Path) -> None:
    (tmp_path / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "fixtures": [
                    {
                        "case_id": "PENDING",
                        "kind": "pending_algorithm",
                        "structure_file": "reference_graph/PENDING.xyz",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="non-confirmed answers"):
        load_fixture_records(tmp_path)


def test_review_trace_page_renders_current_case_in_new_document(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    xyz_path = tmp_path / "TRACE.xyz"
    xyz_path.write_text(
        """1
carbon
C 0.0 0.0 0.0
""",
        encoding="utf-8",
    )
    db_path = tmp_path / "review.sqlite"
    with sqlite3.connect(db_path) as conn:
        conn.executescript((APP_DIR / "schema.sql").read_text(encoding="utf-8"))
        conn.execute(
            """
            INSERT INTO cases(
                    case_id, row_index, category, xyz_path, total_charge,
                    total_radical_electrons, spin_multiplicity,
                    reference_smiles, candidate_status, metadata_json
                )
                VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                "TRACE",
                1,
                "graph_not_equivalent",
                str(xyz_path),
                2,
                2,
                3,
                "C",
                "ok",
                json.dumps(
                    {
                        "spin_multiplicity_used": "3",
                        "total_radical_electrons_used": "2",
                    }
                ),
            ),
        )
        conn.commit()

    captured: dict[str, object] = {}

    def fake_render(
        cases: object,
        *,
        score_all_candidates: bool,
        dof_max_images: int | None,
        defer_dof_images: bool,
    ) -> str:
        captured["cases"] = cases
        captured["score_all_candidates"] = score_all_candidates
        captured["dof_max_images"] = dof_max_images
        captured["defer_dof_images"] = defer_dof_images
        return "<!doctype html><title>TRACE</title>"

    monkeypatch.setattr(review_server, "render_trace_report", fake_render)
    server = ReviewServer(("127.0.0.1", 0), ReviewHandler)
    server.db_path = db_path
    server.xyz_dir = None
    server.fixtures_dir = tmp_path / "reviewed"
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def get(path: str) -> tuple[int, str, str]:
        connection = HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        connection.request("GET", path)
        response = connection.getresponse()
        body = response.read().decode("utf-8")
        content_type = response.getheader("Content-Type") or ""
        connection.close()
        return response.status, body, content_type

    def post(path: str, payload: object) -> tuple[int, str, str]:
        connection = HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        connection.request(
            "POST",
            path,
            body=json.dumps(payload),
            headers={"Content-Type": "application/json"},
        )
        response = connection.getresponse()
        body = response.read().decode("utf-8")
        content_type = response.getheader("Content-Type") or ""
        connection.close()
        return response.status, body, content_type

    try:
        status, body, content_type = get("/trace/TRACE")
        methane = Chem.AddHs(Chem.MolFromSmiles("C"))
        dof_status, dof_body, dof_content_type = post(
            "/api/render-dof",
            {
                "render_type": "single",
                "sdf": Chem.MolToMolBlock(methane) + "\n$$$$\n",
                "legends": ["methane"],
                "size": [360, 300],
            },
        )
        invalid_dof_status, _, _ = post("/api/render-dof", {"render_type": "single"})
        missing_status, _, _ = get("/trace/MISSING")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)

    assert status == 200
    assert body == "<!doctype html><title>TRACE</title>"
    assert content_type == "text/html; charset=utf-8"
    assert captured["score_all_candidates"] is False
    assert captured["dof_max_images"] is None
    assert captured["defer_dof_images"] is True
    assert dof_status == 200
    assert dof_body.startswith("<svg")
    assert dof_content_type == "image/svg+xml; charset=utf-8"
    assert invalid_dof_status == 400
    traced_cases = captured["cases"]
    assert isinstance(traced_cases, list)
    traced_case = traced_cases[0]
    assert traced_case.id == "TRACE"
    assert traced_case.total_charge == 2
    assert traced_case.total_radical_electrons == 2
    assert traced_case.xyz_source == "review_page"
    assert missing_status == 404


def test_case_electronic_state_does_not_sanitize_reference_smiles() -> None:
    case = {
        "total_charge": 0,
        "reference_smiles": (
            "CO->[La]123456(ON(O)O)(<-OC)(<-N(O)(O1)O2)(<-N(O)(O3)O4)"
            "N1CCCCC1C1NC(C2CCCCN2)NC(C2CCCCN25)N16"
        ),
        "metadata_json": "{}",
    }

    assert case_electronic_state(case) == (0, 0, 1)


def test_review_sdf_preserves_metal_unpaired_electron_property() -> None:
    mol = Chem.MolFromSmiles("[Fe+2]", sanitize=False)
    assert mol is not None
    iron = mol.GetAtomWithIdx(0)
    iron.SetIntProp(METAL_UNPAIRED_ELECTRONS_PROP, 2)
    iron.SetNumRadicalElectrons(0)
    Chem.CreateAtomIntPropertyList(mol, METAL_UNPAIRED_ELECTRONS_PROP)

    restored = next(
        Chem.ForwardSDMolSupplier(
            io.BytesIO(_sdf_text(mol, {}).encode()),
            sanitize=False,
            removeHs=False,
        )
    )

    assert restored is not None
    restored_iron = restored.GetAtomWithIdx(0)
    assert restored_iron.GetNumRadicalElectrons() == 0
    assert get_atom_unpaired_electrons(restored_iron) == 2


def test_review_sdf_uses_case_spin_not_metal_unpaired_electrons(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mol = Chem.MolFromSmiles("[Fe+2]", sanitize=False)
    assert mol is not None
    mol.GetAtomWithIdx(0).SetIntProp(METAL_UNPAIRED_ELECTRONS_PROP, 2)
    Chem.CreateAtomIntPropertyList(mol, METAL_UNPAIRED_ELECTRONS_PROP)
    monkeypatch.setattr(fixture_builder, "reconstruct_case_mol", lambda *args, **kwargs: mol)

    xyz_path = tmp_path / "SPIN.xyz"
    xyz_path.write_text("1\nFe\nFe 0.0 0.0 0.0\n", encoding="utf-8")
    fixtures_dir = tmp_path / "fixtures"
    (fixtures_dir / "manifest.json").parent.mkdir(parents=True)
    (fixtures_dir / "manifest.json").write_text(
        json.dumps({"schema_version": 1, "fixtures": []}), encoding="utf-8"
    )
    case = {
        "case_id": "SPIN",
        "row_index": 1,
        "xyz_path": str(xyz_path),
        "total_charge": 0,
        "total_radical_electrons": 0,
        "spin_multiplicity": 1,
        "reference_smiles": "[Fe+2]",
        "metadata_json": "{}",
    }

    record = sync_review_fixture(
        case,
        {"status": ACCEPT_BOTH_STATUS},
        fixtures_dir=fixtures_dir,
    )

    assert record is not None
    assert record["total_radical_electrons"] == 0
    assert record["spin_multiplicity"] == 1
    sdf = next(
        mol
        for mol in Chem.SDMolSupplier(
            str(fixtures_dir / str(record["structure_file"])),
            sanitize=False,
            removeHs=False,
        )
        if mol is not None
    )
    assert sdf.GetProp("TOTAL_RADICAL_ELECTRONS") == "0"
    assert sdf.GetProp("SPIN_MULTIPLICITY") == "1"
    assert get_atom_unpaired_electrons(sdf.GetAtomWithIdx(0)) == 2


def test_review_fixture_sync_stores_manual_xyz_smiles_and_removes_stale_files(
    tmp_path: Path,
) -> None:
    source = Chem.AddHs(Chem.MolFromSmiles("C"))
    assert AllChem.EmbedMolecule(source, randomSeed=7) == 0  # pyright: ignore[reportAttributeAccessIssue]
    xyz_path = tmp_path / "TEST.xyz"
    xyz_path.write_text(Chem.MolToXYZBlock(source), encoding="utf-8")

    case = {
        "case_id": "TEST",
        "row_index": 1,
        "xyz_path": str(xyz_path),
        "total_charge": 0,
        "reference_smiles": "C",
        "metadata_json": json.dumps(
            {"spin_multiplicity_used": "1", "total_radical_electrons_used": "0"}
        ),
    }
    review = {
        "status": "manual_reference",
        "corrected_smiles": "[CH4]",
        "corrected_molblock": "deliberately not a MolBlock",
        "reviewer": "test",
        "notes": "approved correction",
        "updated_at": "2026-07-17T00:00:00+00:00",
    }
    fixtures_dir = tmp_path / "fixtures"
    fixtures_dir.mkdir()
    (fixtures_dir / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "source_dataset": {"name": "example", "revision": "pinned"},
                "fixtures": [],
            }
        ),
        encoding="utf-8",
    )

    record = sync_review_fixture(case, review, fixtures_dir=fixtures_dir)

    assert record is not None
    assert record["kind"] == "manual_reference"
    assert record["structure_file"] == "manual_reference/TEST.xyz"
    assert record["approved_smiles"] == "[CH4]"
    assert record["reference_smiles"] == "C"
    assert record["source"] == "corrected_smiles"
    assert record["total_charge"] == 0
    assert record["total_radical_electrons"] == 0
    assert record["spin_multiplicity"] == 1
    assert {"reviewer", "notes", "updated_at"}.isdisjoint(record)
    assert json.loads((fixtures_dir / "manifest.json").read_text())["source_dataset"] == {
        "name": "example",
        "revision": "pinned",
    }
    frozen_xyz = fixtures_dir / str(record["structure_file"])
    assert frozen_xyz.read_bytes() == xyz_path.read_bytes()

    removed = sync_review_fixture(
        case,
        {**review, "status": "needs_followup"},
        fixtures_dir=fixtures_dir,
    )
    assert removed is None
    assert json.loads((fixtures_dir / "manifest.json").read_text())["fixtures"] == []
    assert not (fixtures_dir / "manual_reference" / "TEST.xyz").exists()

    reference_record = sync_review_fixture(
        case,
        {
            **review,
            "status": "accept_reference",
            "corrected_smiles": "",
            "corrected_molblock": "",
        },
        fixtures_dir=fixtures_dir,
    )
    assert reference_record is not None
    assert reference_record["kind"] == "reference_graph"
    assert reference_record["reference_smiles"] == "C"
    assert reference_record["total_charge"] == 0
    assert reference_record["spin_multiplicity"] == 1
    frozen_xyz = fixtures_dir / str(reference_record["structure_file"])
    assert frozen_xyz.read_text(encoding="utf-8") == xyz_path.read_text(encoding="utf-8")


def test_reference_answer_wrong_status_does_not_modify_fixture(tmp_path: Path) -> None:
    fixtures_dir = tmp_path / "fixtures"
    approved_dir = fixtures_dir / "approved_graph"
    approved_dir.mkdir(parents=True)
    approved_path = approved_dir / "WRONG.xyz"
    approved_path.write_text("fixture\n", encoding="utf-8")
    record = {
        "case_id": "WRONG",
        "kind": "approved_graph",
        "structure_file": "approved_graph/WRONG.xyz",
        "row_index": 1,
        "total_charge": 0,
        "total_radical_electrons": 0,
        "spin_multiplicity": 1,
        "reference_smiles": "C",
        "approved_smiles": "C",
        "accepted_smiles": ["C"],
        "source": "molgr_reconstruction",
        "review_status": "accept_candidate",
    }
    (fixtures_dir / "manifest.json").write_text(
        json.dumps({"schema_version": 1, "fixtures": [record]}),
        encoding="utf-8",
    )

    result = sync_review_fixture(
        {"case_id": "WRONG"},
        {"status": "reference_answer_wrong"},
        fixtures_dir=fixtures_dir,
    )

    assert result == record
    assert approved_path.read_text(encoding="utf-8") == "fixture\n"
    assert json.loads((fixtures_dir / "manifest.json").read_text(encoding="utf-8"))["fixtures"] == [
        record
    ]


def test_review_api_writes_fixture_with_review_decision(tmp_path: Path) -> None:
    xyz_path = tmp_path / "API.xyz"
    xyz_path.write_text(
        """3
water
O 0.000000 0.000000 0.000000
H 0.957200 0.000000 0.000000
H -0.239987 0.927297 0.000000
""",
        encoding="utf-8",
    )
    db_path = tmp_path / "review.sqlite"
    with sqlite3.connect(db_path) as conn:
        conn.executescript((APP_DIR / "schema.sql").read_text(encoding="utf-8"))
        conn.execute(
            """
            INSERT INTO cases(
                case_id, row_index, category, xyz_path, total_charge,
                reference_smiles, candidate_status, metadata_json
            )
            VALUES(?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                "API",
                1,
                "graph_not_equivalent",
                str(xyz_path),
                0,
                "O",
                "ok",
                json.dumps(
                    {
                        "spin_multiplicity_used": "1",
                        "total_radical_electrons_used": "0",
                    }
                ),
            ),
        )
        conn.commit()

    server = ReviewServer(("127.0.0.1", 0), ReviewHandler)
    server.db_path = db_path
    server.xyz_dir = None
    server.fixtures_dir = tmp_path / "reviewed"
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def post_review(review: dict[str, str]) -> tuple[int, dict[str, Any]]:
        connection = HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        connection.request(
            "POST",
            "/api/cases/API/review",
            body=json.dumps(review),
            headers={"Content-Type": "application/json"},
        )
        response = connection.getresponse()
        response_payload = json.loads(response.read())
        connection.close()
        return response.status, response_payload

    try:
        missing_smiles_status, missing_smiles_payload = post_review(
            {
                "status": "manual_reference",
                "corrected_smiles": "",
                "corrected_molblock": "a MolBlock is not fixture-authoritative",
                "notes": "",
                "reviewer": "tester",
            }
        )
        reference_status, reference_payload = post_review(
            {
                "status": "accept_reference",
                "corrected_smiles": "",
                "corrected_molblock": "",
                "notes": "reference graph approved",
                "reviewer": "tester",
            }
        )
        both_status, both_payload = post_review(
            {
                "status": "accept_both",
                "corrected_smiles": "",
                "corrected_molblock": "",
                "notes": "candidate and reference are both acceptable",
                "reviewer": "tester",
            }
        )
        manual_status, manual_payload = post_review(
            {
                "status": "manual_reference",
                "corrected_smiles": "O",
                "corrected_molblock": "deliberately not a MolBlock",
                "notes": "manually approved graph",
                "reviewer": "tester",
            }
        )
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)

    assert missing_smiles_status == 400
    assert missing_smiles_payload == {"error": "corrected_smiles_required_for_manual_reference"}
    assert reference_status == 200
    assert reference_payload["fixture"]["kind"] == "reference_graph"
    assert both_status == 200
    assert both_payload["fixture"]["kind"] == "accepted_both"
    assert len(both_payload["fixture"]["accepted_smiles"]) == 2
    assert both_payload["fixture"]["reference_smiles"] == "O"
    assert manual_status == 200
    assert manual_payload["fixture"]["kind"] == "manual_reference"
    assert (server.fixtures_dir / manual_payload["fixture"]["structure_file"]).is_file()
    assert not (server.fixtures_dir / "reference_graph" / "API.xyz").exists()
    with sqlite3.connect(db_path) as conn:
        assert conn.execute("SELECT status FROM reviews WHERE case_id = 'API'").fetchone() == (
            "manual_reference",
        )


def test_review_api_reports_and_removes_existing_fixture(tmp_path: Path) -> None:
    xyz_path = tmp_path / "SYNC.xyz"
    xyz_path.write_text(
        """1
carbon
C 0.0 0.0 0.0
""",
        encoding="utf-8",
    )
    fixtures_dir = tmp_path / "reviewed"
    case = {
        "case_id": "SYNC",
        "row_index": 1,
        "xyz_path": str(xyz_path),
        "total_charge": 0,
        "reference_smiles": "C",
        "metadata_json": json.dumps(
            {"spin_multiplicity_used": "1", "total_radical_electrons_used": "0"}
        ),
    }
    review = {
        "status": "accept_reference",
        "corrected_smiles": "",
        "corrected_molblock": "",
        "reviewer": "tester",
        "notes": "",
        "updated_at": "2026-07-18T00:00:00+00:00",
    }

    record = sync_review_fixture(case, review, fixtures_dir=fixtures_dir)
    assert record is not None
    assert load_fixture_records(fixtures_dir)["SYNC"]["kind"] == "reference_graph"

    db_path = tmp_path / "review.sqlite"
    with sqlite3.connect(db_path) as conn:
        conn.executescript((APP_DIR / "schema.sql").read_text(encoding="utf-8"))
        conn.execute(
            """
            INSERT INTO cases(
                case_id, row_index, category, xyz_path, total_charge,
                reference_smiles, candidate_status, metadata_json
            )
            VALUES(?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                "SYNC",
                1,
                "graph_not_equivalent",
                str(xyz_path),
                0,
                "C",
                "ok",
                case["metadata_json"],
            ),
        )
        conn.commit()

    server = ReviewServer(("127.0.0.1", 0), ReviewHandler)
    server.db_path = db_path
    server.xyz_dir = None
    server.fixtures_dir = fixtures_dir
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def request(path: str, body: dict[str, str] | None = None) -> tuple[int, dict[str, Any]]:
        connection = HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        if body is None:
            connection.request("GET", path)
        else:
            connection.request(
                "POST",
                path,
                body=json.dumps(body),
                headers={"Content-Type": "application/json"},
            )
        response = connection.getresponse()
        payload = json.loads(response.read())
        connection.close()
        return response.status, payload

    try:
        status, payload = request("/api/cases/SYNC")
        assert status == 200
        assert payload["fixture"]["kind"] == "reference_graph"
        status, payload = request("/api/cases?q=SYNC")
        assert status == 200
        assert payload["items"][0]["fixture"]["kind"] == "reference_graph"

        status, payload = request(
            "/api/cases/SYNC/review",
            {
                "status": "needs_followup",
                "corrected_smiles": "",
                "corrected_molblock": "",
                "reviewer": "tester",
                "notes": "deferred",
            },
        )
        assert status == 200
        assert payload["fixture"] is None

        status, payload = request("/api/cases/SYNC")
        assert status == 200
        assert payload["fixture"] is None
        assert load_fixture_records(fixtures_dir) == {}
        assert not (fixtures_dir / "reference_graph" / "SYNC.xyz").exists()
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)


def test_review_undo_restores_unreviewed_filter_and_fixture_state(tmp_path: Path) -> None:
    xyz_path = tmp_path / "UNDO.xyz"
    xyz_path.write_text("1\ncarbon\nC 0 0 0\n", encoding="utf-8")
    fixtures_dir = tmp_path / "reviewed"
    db_path = tmp_path / "review.sqlite"
    with sqlite3.connect(db_path) as conn:
        conn.executescript((APP_DIR / "schema.sql").read_text(encoding="utf-8"))
        conn.execute(
            """
            INSERT INTO cases(case_id, row_index, xyz_path, total_charge, reference_smiles,
                              candidate_status, metadata_json)
            VALUES('UNDO', 1, ?, 0, 'C', 'ok', ?)
            """,
            (
                str(xyz_path),
                json.dumps({"spin_multiplicity_used": "1", "total_radical_electrons_used": "0"}),
            ),
        )
        conn.commit()

    server = ReviewServer(("127.0.0.1", 0), ReviewHandler)
    server.db_path = db_path
    server.xyz_dir = None
    server.fixtures_dir = fixtures_dir
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def request(path: str, *, method: str = "GET", body: dict[str, str] | None = None):
        connection = HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        headers = {"X-MolGR-Review-Session": "undo-session"}
        encoded = None
        if body is not None:
            headers["Content-Type"] = "application/json"
            encoded = json.dumps(body)
        connection.request(method, path, body=encoded, headers=headers)
        response = connection.getresponse()
        payload = json.loads(response.read())
        connection.close()
        return response.status, payload

    try:
        status, saved = request(
            "/api/cases/UNDO/review",
            method="POST",
            body={
                "status": "accept_reference",
                "reviewer": "resonance-representation",
                "notes": "evidence",
            },
        )
        assert status == 200
        assert saved["mutation"]["case_id"] == "UNDO"
        assert load_fixture_records(fixtures_dir)["UNDO"]["kind"] == "reference_graph"
        status, filtered = request("/api/cases?status=unreviewed")
        assert status == 200
        assert all(item["case_id"] != "UNDO" for item in filtered["items"])

        status, undone = request(
            "/api/review-undo",
            method="POST",
            body={"mutation_id": saved["mutation"]["mutation_id"]},
        )
        assert status == 200
        assert undone["case_id"] == "UNDO"
        assert undone["restored_review"] is None
        assert load_fixture_records(fixtures_dir) == {}
        with sqlite3.connect(db_path) as conn:
            assert conn.execute("SELECT * FROM reviews WHERE case_id='UNDO'").fetchone() is None
        status, filtered = request("/api/cases?status=unreviewed")
        assert status == 200
        assert [item["case_id"] for item in filtered["items"]] == ["UNDO"]
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)


def test_review_undo_restores_overwritten_review_and_exact_fixture(tmp_path: Path) -> None:
    xyz_path = tmp_path / "OVERWRITE.xyz"
    xyz_path.write_text("1\ncarbon\nC 0 0 0\n", encoding="utf-8")
    fixtures_dir = tmp_path / "reviewed"
    case = {
        "case_id": "OVERWRITE",
        "row_index": 2,
        "xyz_path": str(xyz_path),
        "total_charge": 0,
        "reference_smiles": "C",
        "metadata_json": json.dumps(
            {"spin_multiplicity_used": "1", "total_radical_electrons_used": "0"}
        ),
    }
    original_review = {
        "status": "accept_reference",
        "corrected_smiles": "",
        "corrected_molblock": "",
        "notes": "original evidence",
        "reviewer": "original-reason",
        "updated_at": "2026-08-17T00:00:00+00:00",
    }
    sync_review_fixture(case, original_review, fixtures_dir=fixtures_dir)
    original_manifest = (fixtures_dir / "manifest.json").read_text(encoding="utf-8")
    original_xyz = (fixtures_dir / "reference_graph/OVERWRITE.xyz").read_text(encoding="utf-8")
    db_path = tmp_path / "review.sqlite"
    with sqlite3.connect(db_path) as conn:
        conn.executescript((APP_DIR / "schema.sql").read_text(encoding="utf-8"))
        conn.execute(
            """INSERT INTO cases(case_id,row_index,xyz_path,total_charge,reference_smiles,
                                  candidate_status,metadata_json)
               VALUES('OVERWRITE',2,?,0,'C','ok',?)""",
            (str(xyz_path), case["metadata_json"]),
        )
        conn.execute(
            """INSERT INTO reviews(case_id,status,corrected_smiles,corrected_molblock,
                                    notes,reviewer,updated_at)
               VALUES('OVERWRITE',?,?,?,?,?,?)""",
            tuple(original_review[key] for key in review_server.REVIEW_COLUMNS),
        )
        conn.commit()

    server = ReviewServer(("127.0.0.1", 0), ReviewHandler)
    server.db_path = db_path
    server.xyz_dir = None
    server.fixtures_dir = fixtures_dir
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def post(path: str, body: dict[str, str]):
        connection = HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        connection.request(
            "POST",
            path,
            body=json.dumps(body),
            headers={
                "Content-Type": "application/json",
                "X-MolGR-Review-Session": "overwrite-session",
            },
        )
        response = connection.getresponse()
        payload = json.loads(response.read())
        connection.close()
        return response.status, payload

    try:
        status, saved = post(
            "/api/cases/OVERWRITE/review",
            {"status": "needs_followup", "reviewer": "new-reason", "notes": "new note"},
        )
        assert status == 200
        assert load_fixture_records(fixtures_dir) == {}
        status, undone = post("/api/review-undo", {"mutation_id": saved["mutation"]["mutation_id"]})
        assert status == 200
        assert undone["restored_review"] == original_review
        with sqlite3.connect(db_path) as conn:
            restored = conn.execute(
                "SELECT status,corrected_smiles,corrected_molblock,notes,reviewer,updated_at "
                "FROM reviews WHERE case_id='OVERWRITE'"
            ).fetchone()
        assert restored == tuple(original_review[key] for key in review_server.REVIEW_COLUMNS)
        assert (fixtures_dir / "manifest.json").read_text(encoding="utf-8") == original_manifest
        assert (fixtures_dir / "reference_graph/OVERWRITE.xyz").read_text(
            encoding="utf-8"
        ) == original_xyz
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)


def test_live_candidate_payload_is_self_consistent_and_does_not_use_snapshot_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    db_path = tmp_path / "review.sqlite"
    with sqlite3.connect(db_path) as conn:
        conn.executescript((APP_DIR / "schema.sql").read_text(encoding="utf-8"))
        conn.execute(
            """
            INSERT INTO cases(
                case_id, row_index, category, xyz_path, total_charge,
                candidate_smiles, candidate_status, metadata_json
            )
            VALUES(?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                "LIVE",
                1,
                "candidate_failed",
                str(tmp_path / "not-needed.xyz"),
                0,
                "C",
                "failed",
                "{}",
            ),
        )
        conn.execute(
            """
            INSERT INTO render_cache(case_id, kind, svg, smiles, error, generated_at)
            VALUES(?, ?, ?, ?, ?, ?)
            """,
            ("LIVE", "candidate", "<svg>stale</svg>", "N", "", "2026-07-17T00:00:00Z"),
        )
        conn.commit()

    reconstruct_calls = 0

    def reconstruct_live(*args: object, **kwargs: object) -> Chem.Mol:
        nonlocal reconstruct_calls
        reconstruct_calls += 1
        mol = Chem.MolFromSmiles("CC")
        assert mol is not None
        return Chem.AddHs(mol)

    monkeypatch.setattr(review_server, "reconstruct_case_mol", reconstruct_live)
    monkeypatch.setattr(
        review_server,
        "_render_mol_svg",
        lambda mol, *, legend, atom_notes=None, show_hydrogens=False: "<svg>live</svg>",
    )

    server = ReviewServer(("127.0.0.1", 0), ReviewHandler)
    server.db_path = db_path
    server.xyz_dir = None
    server.fixtures_dir = tmp_path / "reviewed"
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def get_json(path: str) -> dict[str, object]:
        connection = HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        connection.request("GET", path)
        response = connection.getresponse()
        assert response.status == 200
        payload = json.loads(response.read())
        connection.close()
        return payload

    try:
        case_payload = get_json("/api/cases/LIVE")
        live_payload = get_json("/api/cases/LIVE/candidate-sdf")
        first_render = get_json("/api/cases/LIVE/render?kind=candidate")
        second_render = get_json("/api/cases/LIVE/render?kind=candidate")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)

    assert case_payload["candidate_snapshot_status"] == "failed"
    assert case_payload["candidate_snapshot_smiles"] == "C"
    assert live_payload["available"] is True
    assert live_payload["source"] == "live_reconstruction"
    assert live_payload["live_candidate_status"] == "ok"
    assert live_payload["live_candidate_smiles"] == "CC"
    assert live_payload["live_candidate_smiles_exact_match"] is False
    assert live_payload["smiles"] == "CC"
    sdf_mol = Chem.MolFromMolBlock(
        str(live_payload["sdf"]).partition("$$$$")[0],
        sanitize=False,
        removeHs=False,
    )
    assert sdf_mol is not None
    assert sdf_mol.GetNumAtoms() == 8
    assert sum(atom.GetAtomicNum() == 1 for atom in sdf_mol.GetAtoms()) == 6  # pyright: ignore[reportCallIssue]
    assert live_payload["candidate_snapshot_smiles"] == "C"
    assert live_payload["live_matches_candidate_snapshot"] is False
    assert live_payload["live_candidate_equivalence_reason"] == (
        "Not equivalent: heavy-atom element counts differ."
    )
    assert live_payload["svg"] == "<svg>live</svg>"
    assert first_render["svg"] == "<svg>live</svg>"
    assert second_render["svg"] == "<svg>live</svg>"
    assert reconstruct_calls == 3


def test_reference_render_cache_is_scoped_to_renderer_version(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    db_path = tmp_path / "review.sqlite"
    with sqlite3.connect(db_path) as conn:
        conn.executescript((APP_DIR / "schema.sql").read_text(encoding="utf-8"))
        conn.execute(
            """
            INSERT INTO cases(case_id, row_index, xyz_path, reference_smiles)
            VALUES(?, ?, ?, ?)
            """,
            ("CACHE", 1, str(tmp_path / "unused.xyz"), "CC"),
        )
        conn.execute(
            """
            INSERT INTO render_cache(case_id, kind, svg, smiles, error, generated_at)
            VALUES(?, ?, ?, ?, ?, ?)
            """,
            ("CACHE", "reference", "<svg>old</svg>", "CC", "", "2026-07-17T00:00:00Z"),
        )
        conn.commit()

    render_calls = 0

    def render_current(
        mol: Chem.Mol,
        *,
        legend: str,
        atom_notes: dict[int, str] | None = None,
        show_hydrogens: bool = False,
    ) -> str:
        nonlocal render_calls
        render_calls += 1
        assert mol.GetNumAtoms() == 2
        assert sum(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()) == 0  # pyright: ignore[reportCallIssue]
        return "<svg>review-2d-v2</svg>"

    monkeypatch.setattr(review_server, "_render_mol_svg", render_current)

    server = ReviewServer(("127.0.0.1", 0), ReviewHandler)
    server.db_path = db_path
    server.xyz_dir = None
    server.fixtures_dir = tmp_path / "reviewed"
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def get_reference(mode: str = "skeleton") -> dict[str, object]:
        connection = HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        connection.request("GET", f"/api/cases/CACHE/render?kind=reference&mode={mode}")
        response = connection.getresponse()
        assert response.status == 200
        payload = json.loads(response.read())
        connection.close()
        return payload

    try:
        first = get_reference()
        second = get_reference()
        hydrogen = get_reference("hydrogen")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)

    assert first["svg"] == "<svg>review-2d-v2</svg>"
    assert second["svg"] == "<svg>review-2d-v2</svg>"
    assert hydrogen["svg"] == "<svg>review-2d-v2</svg>"
    assert render_calls == 2
    with sqlite3.connect(db_path) as conn:
        cache_kinds = {
            row[0] for row in conn.execute("SELECT kind FROM render_cache WHERE case_id = 'CACHE'")
        }
    assert cache_kinds == {
        "reference",
        "review_2d_v2:reference:skeleton",
        "review_2d_v2:reference:hydrogen",
    }


def test_case_payload_selects_candidate_snapshot_for_server_python() -> None:
    runtime_label = f"py{sys.version_info.major}{sys.version_info.minor}"
    other_label = "py310" if runtime_label != "py310" else "py38"
    payload = review_server._row_dict(
        {
            "case_id": "RUNTIME",
            "candidate_smiles": "fallback",
            "candidate_status": "fallback_status",
            "metadata_json": json.dumps(
                {
                    f"{runtime_label}_molgr_cpp_smiles": "runtime",
                    f"{runtime_label}_molgr_cpp_status": "ok",
                    f"{other_label}_molgr_cpp_smiles": "other",
                    f"{other_label}_molgr_cpp_status": "error",
                }
            ),
        }
    )

    assert payload is not None
    assert payload["candidate_snapshot_runtime"] == runtime_label
    assert payload["candidate_snapshot_smiles"] == "runtime"
    assert payload["candidate_snapshot_status"] == "ok"

    failed_payload = review_server._row_dict(
        {
            "case_id": "FAILED",
            "candidate_smiles": "other-version-result",
            "candidate_status": "ok",
            "metadata_json": json.dumps(
                {
                    f"{runtime_label}_molgr_cpp_smiles": "",
                    f"{runtime_label}_molgr_cpp_status": "error",
                }
            ),
        }
    )
    assert failed_payload is not None
    assert failed_payload["candidate_snapshot_smiles"] == ""
    assert failed_payload["candidate_snapshot_status"] == "error"
