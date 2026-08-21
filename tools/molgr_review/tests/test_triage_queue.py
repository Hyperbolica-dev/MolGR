from __future__ import annotations

import csv
import json
import sqlite3
import subprocess
import sys
import threading
from http.client import HTTPConnection
from pathlib import Path


APP_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_DIR))

from server import ReviewHandler, ReviewServer, load_triage_records  # noqa: E402


def test_read_only_triage_filter_and_payload(tmp_path: Path) -> None:
    db_path = tmp_path / "review.sqlite"
    with sqlite3.connect(db_path) as connection:
        connection.executescript((APP_DIR / "schema.sql").read_text(encoding="utf-8"))
        connection.executemany(
            "INSERT INTO cases(case_id, row_index, xyz_path) VALUES(?, ?, '')",
            [("A", 1), ("B", 2), ("C", 3)],
        )
    triage_path = tmp_path / "triage.csv"
    with triage_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["case_id", "triage_bucket", "mapping_confidence", "machine_reason"],
        )
        writer.writeheader()
        writer.writerows(
            [
                {
                    "case_id": "A",
                    "triage_bucket": "strong_xyz_candidate_evidence",
                    "mapping_confidence": "unique_graph_mapping",
                    "machine_reason": "H assignment",
                },
                {
                    "case_id": "B",
                    "triage_bucket": "unknown",
                    "mapping_confidence": "ambiguous",
                    "machine_reason": "manual",
                },
            ]
        )

    server = ReviewServer(("127.0.0.1", 0), ReviewHandler)
    server.db_path = db_path
    server.xyz_dir = None
    server.fixtures_dir = tmp_path / "fixtures"
    server.triage_records = load_triage_records(triage_path)
    server.runtime_info = {}
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def get(path: str):
        connection = HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        connection.request("GET", path)
        response = connection.getresponse()
        payload = json.loads(response.read())
        connection.close()
        return response.status, payload

    try:
        status, stats = get("/api/stats")
        assert status == 200
        assert stats["triage_buckets"] == {
            "strong_xyz_candidate_evidence": 1,
            "unknown": 1,
        }
        status, listing = get(
            "/api/cases?status=unreviewed&triage_bucket=strong_xyz_candidate_evidence"
        )
        assert status == 200
        assert listing["total"] == 1
        assert listing["triage_bucket_counts"] == {
            "strong_xyz_candidate_evidence": 1,
            "unknown": 1,
        }
        assert listing["items"][0]["case_id"] == "A"
        assert listing["items"][0]["triage"]["mapping_confidence"] == "unique_graph_mapping"
        with sqlite3.connect(db_path) as connection:
            connection.execute(
                "INSERT INTO reviews(case_id, status, reviewer) VALUES('A', 'accept_candidate', 'test')"
            )
            connection.commit()
        status, remaining = get("/api/cases?status=unreviewed")
        assert status == 200
        assert remaining["triage_bucket_counts"] == {"unknown": 1}
        status, case = get("/api/cases/A")
        assert status == 200
        assert case["triage_bucket"] == "strong_xyz_candidate_evidence"
        assert case["triage"]["machine_reason"] == "H assignment"
        with sqlite3.connect(db_path) as connection:
            assert connection.execute("SELECT COUNT(*) FROM reviews").fetchone()[0] == 1
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)


def test_frontend_has_compact_triage_queue_and_shortcuts() -> None:
    html = (APP_DIR / "static/index.html").read_text(encoding="utf-8")
    javascript = (APP_DIR / "static/app.js").read_text(encoding="utf-8")
    stylesheet = (APP_DIR / "static/style.css").read_text(encoding="utf-8")

    assert 'id="triageFilter"' in html
    assert 'id="triageEvidence"' in html
    assert 'params.set("triage_bucket", triageBucket)' in javascript
    assert "renderTriageBucketOptions(data.triage_bucket_counts || {})" in javascript
    assert "function renderTriageEvidence(item)" in javascript
    assert "function metalEdgeEvidence(triage, edge)" in javascript
    assert "function hydrogenEvidence(triage, hydrogen)" in javascript
    assert "function redoxEvidence(triage)" in javascript
    assert "function triageLocalization(item)" not in javascript
    assert "function applyTriageLocalization(viewer, item)" not in javascript
    assert "viewer.addLabel" not in javascript
    assert "viewer.addLine" not in javascript
    assert "&localize=1" in javascript
    assert 'data-mode="skeleton"' in html
    assert 'data-mode="hydrogen"' in html
    assert 'data-mode="raw"' in html
    assert 'data-mode="mapped"' in html
    assert 'id="mappedComparisonNote"' in html
    assert "function loadReferenceXyz(item, token)" in javascript
    assert "function renderReferenceXyz3d(" in javascript
    assert "function applyMappedComparison(viewer, side)" not in javascript
    assert "function renderMappedComparisonNote()" in javascript
    assert 'class="reference-xyz-ambiguity"' in javascript
    assert 'mappingAmbiguityTypes.${data.mapping_ambiguity_type}' in javascript
    assert "multiple_equally_valid_atom_mappings" in javascript
    assert "function mappingAmbiguityEvidence(data)" in javascript
    assert 'tr("ambiguityLocation")' in javascript
    assert "mapping_ambiguity_locations" in javascript
    assert 'state.xyzComparisonMode === "mapped" && state.currentCandidateSdf' in javascript
    assert 'class="triage-evidence-grid"' in javascript
    assert 'class="triage-trace"' in javascript
    triage_trace_css = stylesheet.split(".triage-trace", 1)[1].split(".triage-tag", 1)[0]
    assert "text-overflow: ellipsis" not in triage_trace_css
    assert 'strong_xyz_candidate_evidence: "XYZ→候选"' in javascript
    assert 'strong_xyz_candidate_evidence: "XYZ→Candidate"' in javascript
    assert "function localizeTriageFilterOptions()" in javascript
    assert "localizeTriageFilterOptions();" in javascript
    assert "item.triage_bucket !== filteredTriageBucket" in javascript
    assert '$("reviewer").value = item.reviewer || ""' in javascript
    assert "suggestedReviewReason" not in javascript
    assert "moleculeReviewReviewer" not in javascript
    assert "function navigateCase(delta)" in javascript
    assert "function reviewShortcut(event)" in javascript
    assert 'document.addEventListener("keydown", reviewShortcut)' in javascript
    assert "await loadCase(queuedNextCaseId)" in javascript


def test_structured_triage_evidence_preserves_reviewer_values() -> None:
    javascript = (APP_DIR / "static/app.js").read_text(encoding="utf-8")
    javascript = javascript.rsplit("\ninit().catch", 1)[0]
    harness = r"""
const assert = (condition, message) => { if (!condition) throw new Error(message); };
const triage = {
  triage_bucket: "strong_xyz_reference_evidence",
  metal_elements: "Co",
  mapping_confidence: "unique_graph_mapping",
};
const metal = metalEdgeEvidence(triage, {
  elements: ["S", "Co"], candidate_atoms: [21, 16], distance: 2.0688,
  edge_present_in: "reference", reference_coordination_number: 6,
  inside_agreed_shell_range: true,
});
assert(metal.includes("Reference-only coordination"), "missing edge side");
assert(metal.includes("Co · XYZ #16 ↔ S · XYZ #21"), "wrong atom/index order");
assert(metal.includes("2.069 Å"), "distance was lost");
assert(metal.includes("Candidate</span><strong>absent"), "candidate presence wrong");
assert(metal.includes("Reference</span><strong>present"), "reference presence wrong");
assert(metal.includes("CN=6 · disputed atom in donor shell"), "compact Trace evidence missing");

const hydrogen = hydrogenEvidence(triage, {
  h_atom: 42, candidate_center_element: "B", candidate_center: 40,
  candidate_distance: 1.272, reference_center_element: "Ru", reference_center: 0,
  reference_distance: 1.888, distance_margin: 0.616,
});
assert(hydrogen.includes("H · XYZ #42"), "H index missing");
assert(hydrogen.includes("B-H 1.272 Å · Ru-H 1.888 Å"), "H distances missing");
assert(hydrogen.includes("0.616 Å"), "distance margin missing");

const redox = redoxEvidence({
  metal_elements: "Co", candidate_metal_state: "Co(+5)", reference_metal_state: "Co(+3)",
  metal_charge_delta: "2", ligand_charge_delta: "-2", candidate_ligand_charge_sum: "-4",
  reference_ligand_charge_sum: "-2", metal_coordination_diff: "[]",
  mapping_confidence: "unique_graph_mapping",
});
assert(redox.includes("Co(+5)"), "candidate metal state missing");
assert(redox.includes("Co(+3)"), "reference metal state missing");
assert(redox.includes("Δ -2"), "ligand compensation missing");
assert(triageBucketLabel("strong_xyz_candidate_evidence") === "XYZ→候选", "Chinese bucket label missing");
assert(hasHydrogenAssignment({ triage: { hydrogen_assignment_diff: JSON.stringify([{ h_atom: 20 }]) } }), "H-assignment case must default to hydrogen mode");
assert(!hasHydrogenAssignment({ triage: { hydrogen_assignment_diff: "[]" } }), "non-H case must default to skeleton mode");
state.language = "en";
assert(triageBucketLabel("strong_xyz_candidate_evidence") === "XYZ→Candidate", "English bucket label missing");
"""
    subprocess.run(
        [
            "node",
            "-e",
            "globalThis.localStorage = { getItem() { return null; }, setItem() {} };\n"
            + javascript
            + "\n"
            + harness,
        ],
        check=True,
        cwd=APP_DIR,
        text=True,
    )
