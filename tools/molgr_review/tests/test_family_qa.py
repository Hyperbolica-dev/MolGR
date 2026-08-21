from __future__ import annotations

import csv
import json
import sqlite3
import sys
import threading
from http.client import HTTPConnection
from pathlib import Path
from typing import Any


APP_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_DIR))

from server import ReviewHandler, ReviewServer, load_family_qa  # noqa: E402


def _write_queue(path: Path) -> None:
    fields = [
        "priority", "family_id", "family_size", "case_id", "calibration_relation",
        "candidate_metal_state", "reference_metal_state", "repeat_count", "mapping_source",
        "canonical_transformation",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for case_id in ("CASE_A", "CASE_B"):
            writer.writerow(
                {
                    "priority": "1", "family_id": "RF001", "family_size": "5",
                    "case_id": case_id, "calibration_relation": "R6-like",
                    "candidate_metal_state": "Mo(+2)", "reference_metal_state": "Mo(+6)",
                    "repeat_count": "2", "mapping_source": "unique_graph_mapping",
                    "canonical_transformation": json.dumps({"metal_charge_delta": 4}),
                }
            )


def test_family_qa_updates_only_pending_manifest_and_undoes(tmp_path: Path) -> None:
    db_path = tmp_path / "review.sqlite"
    with sqlite3.connect(db_path) as conn:
        conn.executescript((APP_DIR / "schema.sql").read_text(encoding="utf-8"))
        conn.executemany(
            "INSERT INTO cases(case_id, row_index, xyz_path) VALUES(?, ?, '')",
            [("CASE_A", 1), ("CASE_B", 2)],
        )
    queue = tmp_path / "representation_qa_queue.csv"
    _write_queue(queue)
    manifest_path = tmp_path / "representation_families.pending.json"
    original_manifest = {
        "approval_required": True,
        "approved": False,
        "families": [{
            "family_id": "RF001", "size": 5,
            "representatives": ["CASE_A", "CASE_B"], "members": ["CASE_A", "CASE_B"],
            "proposed_status": None, "proposed_reason": None, "qa_passed": False,
        }],
    }
    manifest_path.write_text(json.dumps(original_manifest), encoding="utf-8")
    triage = {case_id: {"case_id": case_id, "triage_bucket": "complex_multi_difference"} for case_id in ("CASE_A", "CASE_B")}

    server = ReviewServer(("127.0.0.1", 0), ReviewHandler)
    server.db_path = db_path
    server.xyz_dir = None
    server.fixtures_dir = tmp_path / "fixtures"
    server.fixtures_dir.mkdir()
    server.triage_records = triage
    server.family_qa = load_family_qa(queue, manifest_path)
    server.family_qa_manifest_path = manifest_path
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def request(method: str, path: str, body: dict[str, str] | None = None) -> tuple[int, dict[str, Any]]:
        connection = HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        headers = {"X-MolGR-Review-Session": "family-test", "Content-Type": "application/json"}
        connection.request(method, path, body=json.dumps(body) if body else None, headers=headers)
        response = connection.getresponse()
        payload = json.loads(response.read())
        connection.close()
        return response.status, payload

    try:
        status, payload = request("GET", "/api/family-qa")
        assert status == 200
        assert payload["progress"] == {
            "reviewed_families": 0, "total_families": 1,
            "approved_cases": 0, "total_cases": 5,
        }
        assert not payload["missing_join_cases"]
        status, payload = request("POST", "/api/family-qa", {
            "family_id": "RF001", "case_id": "CASE_A",
            "action": "representative_mark", "value": "matches_family",
        })
        assert status == 200
        assert payload["mutation"]["mutation_type"] == "family_qa"
        status, payload = request("POST", "/api/family-qa", {
            "family_id": "RF001", "action": "decision", "value": "approve_redox",
        })
        assert status == 200
        saved = json.loads(manifest_path.read_text(encoding="utf-8"))
        assert saved["approval_required"] is True and saved["approved"] is False
        assert saved["families"][0]["proposed_status"] == "accept_both"
        assert saved["families"][0]["proposed_reason"] == "redox-representation"
        with sqlite3.connect(db_path) as conn:
            assert conn.execute("SELECT COUNT(*) FROM reviews").fetchone()[0] == 0

        mutation_id = payload["mutation"]["mutation_id"]
        status, _ = request("POST", "/api/review-undo", {"mutation_id": mutation_id})
        assert status == 200
        restored = json.loads(manifest_path.read_text(encoding="utf-8"))
        assert "qa_decision" not in restored["families"][0]
        assert restored["families"][0]["representative_marks"] == {"CASE_A": "matches_family"}
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)


def test_family_qa_frontend_is_pending_only_and_humanizes_transformations() -> None:
    html = (APP_DIR / "static" / "index.html").read_text(encoding="utf-8")
    javascript = (APP_DIR / "static" / "app.js").read_text(encoding="utf-8")
    assert 'id="familyQueue"' in html
    assert 'id="familyQaCard"' in html
    assert 'data-family-decision="approve_resonance"' in html
    assert 'data-rep-mark="outlier_blocker"' in html
    assert "function transformationSummary(raw, repeatCount = 1)" in javascript
    assert "charge_transitions_per_unit" in javascript
    assert "bond_transitions_per_unit" in javascript
    assert 'if (!state.current || state.savingReview || activeFamily()) return;' in javascript
