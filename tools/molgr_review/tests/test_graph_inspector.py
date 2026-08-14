from __future__ import annotations

import json
import sqlite3
import sys
import threading
from http.client import HTTPConnection
from pathlib import Path
from typing import Any

import pytest
from rdkit import Chem


APP_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_DIR))

import server as review_server  # noqa: E402
from server import ReviewHandler, ReviewServer  # noqa: E402


def test_review_graph_payload_preserves_atom_and_directed_bond_evidence() -> None:
    mol = Chem.MolFromSmiles("N->[Co+3](<-[O-]).S#[O+].[BH2-].[RuH2]")
    assert mol is not None

    payload = review_server._review_graph_payload(mol, kind="reference", smiles="synthetic")
    atoms = {(atom["element"], atom["formal_charge"]): atom for atom in payload["atoms"]}
    dative_bonds = [bond for bond in payload["bonds"] if bond["type"] == "dative"]
    triple_bonds = [bond for bond in payload["bonds"] if bond["type"] == "triple"]

    assert payload["summary"]["total_formal_charge"] == 2
    assert payload["summary"]["total_radical_electrons"] == sum(
        atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()
    )
    assert {
        (metal["element"], metal["formal_charge"]) for metal in payload["summary"]["metals"]
    } == {
        ("Co", 3),
        ("Ru", 0),
    }
    assert atoms[("N", 0)]["implicit_h"] == 3
    assert atoms[("B", -1)]["explicit_h"] == 2
    assert atoms[("Ru", 0)]["explicit_h"] == 2
    source_oxide = next(
        atom for atom in mol.GetAtoms() if atom.GetSymbol() == "O" and atom.GetFormalCharge() == -1
    )
    assert payload["atoms"][source_oxide.GetIdx()]["radical_electrons"] == (
        source_oxide.GetNumRadicalElectrons()
    )
    assert len(dative_bonds) == 2
    assert all(bond["end_element"] == "Co" and bond["directional"] for bond in dative_bonds)
    assert len(triple_bonds) == 1
    assert {triple_bonds[0]["begin_element"], triple_bonds[0]["end_element"]} == {"S", "O"}


def test_graph_api_exposes_candidate_and_reference_without_writing_review_data(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    candidate = Chem.MolFromSmiles("[O-]->[Co+3]<-N")
    assert candidate is not None
    monkeypatch.setattr(review_server, "reconstruct_case_mol", lambda *args, **kwargs: candidate)

    db_path = tmp_path / "review.sqlite"
    reference_smiles = "S#[O+].[O-]->[Co+]"
    with sqlite3.connect(db_path) as conn:
        conn.executescript((APP_DIR / "schema.sql").read_text(encoding="utf-8"))
        conn.execute(
            """
            INSERT INTO cases(case_id, row_index, xyz_path, reference_smiles)
            VALUES('GRAPH', 1, '', ?)
            """,
            (reference_smiles,),
        )
        conn.commit()

    server = ReviewServer(("127.0.0.1", 0), ReviewHandler)
    server.db_path = db_path
    server.xyz_dir = None
    server.fixtures_dir = tmp_path / "reviewed"
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def get(path: str) -> tuple[int, dict[str, Any]]:
        connection = HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        connection.request("GET", path)
        response = connection.getresponse()
        payload = json.loads(response.read())
        connection.close()
        return response.status, payload

    try:
        candidate_status, candidate_payload = get("/api/cases/GRAPH/graph?kind=candidate")
        reference_status, reference_payload = get("/api/cases/GRAPH/graph?kind=reference")
        invalid_status, invalid_payload = get("/api/cases/GRAPH/graph?kind=organic")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)

    assert candidate_status == 200
    assert candidate_payload["summary"]["metals"] == [
        {"index": 1, "element": "Co", "formal_charge": 3}
    ]
    assert all(
        bond["end_element"] == "Co"
        for bond in candidate_payload["bonds"]
        if bond["type"] == "dative"
    )
    assert reference_status == 200
    assert reference_payload["smiles"] == reference_smiles
    assert any(bond["type"] == "triple" for bond in reference_payload["bonds"])
    assert invalid_status == 400
    assert invalid_payload == {"error": "invalid_graph_kind"}
    with sqlite3.connect(db_path) as conn:
        assert conn.execute("SELECT COUNT(*) FROM reviews").fetchone()[0] == 0
        assert conn.execute("SELECT COUNT(*) FROM render_cache").fetchone()[0] == 0


def test_reviewer_details_are_folded_separate_and_do_not_restore_organic_view() -> None:
    html = (APP_DIR / "static" / "index.html").read_text(encoding="utf-8")
    javascript = (APP_DIR / "static" / "app.js").read_text(encoding="utf-8")
    stylesheet = (APP_DIR / "static" / "style.css").read_text(encoding="utf-8")

    assert '<details id="reviewerDetails" class="reviewer-details panel">' in html
    assert '<details id="smilesGraphDetails">' in html
    assert '<details id="atomBondDetails">' in html
    assert '<details id="provenanceDetails">' in html
    assert html.index('id="secondaryVisual"') < html.index('id="reviewerDetails"')
    assert html.index('id="reviewerDetails"') < html.index('class="secondary-visuals panel"')
    assert html.index('class="secondary-visuals panel"') < html.index(
        'class="technical-details panel"'
    )
    assert html.index('class="technical-details panel"') < html.index('id="reviewPanel"')
    assert 'id="jumpToReviewerDetails"' in html
    assert '$("reviewerDetails").open' in javascript
    assert "loadGraphEvidence(state.caseRequestToken)" in javascript
    assert "function jumpToReviewerDetails()" in javascript
    assert 'details.scrollIntoView({ behavior: "smooth", block: "start" })' in javascript
    assert 'bond.directional\n            ? "→"' in javascript
    assert 'bond.type === "triple"\n                ? "≡"' in javascript
    assert 'title="${escapeHtml(tr(tooltip))}"' in javascript
    assert "qTooltip" in javascript
    assert "multiplicityTooltip" in javascript
    assert "radicalsTooltip" in javascript
    assert "formulaTooltip" in javascript
    assert 'class="render-kind' not in html
    assert ".workspace > * { flex-shrink: 0; }" in stylesheet
    assert "candidate_organic" not in html
    assert "reference_organic" not in html
