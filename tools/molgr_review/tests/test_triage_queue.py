from __future__ import annotations

import csv
import json
import sqlite3
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
        assert listing["items"][0]["case_id"] == "A"
        assert listing["items"][0]["triage"]["mapping_confidence"] == "unique_graph_mapping"
        status, case = get("/api/cases/A")
        assert status == 200
        assert case["triage_bucket"] == "strong_xyz_candidate_evidence"
        assert case["triage"]["machine_reason"] == "H assignment"
        with sqlite3.connect(db_path) as connection:
            assert connection.execute("SELECT COUNT(*) FROM reviews").fetchone()[0] == 0
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)


def test_frontend_has_compact_triage_queue_and_shortcuts() -> None:
    html = (APP_DIR / "static/index.html").read_text(encoding="utf-8")
    javascript = (APP_DIR / "static/app.js").read_text(encoding="utf-8")

    assert 'id="triageFilter"' in html
    assert 'id="triageEvidence"' in html
    assert 'params.set("triage_bucket", triageBucket)' in javascript
    assert "function renderTriageEvidence(item)" in javascript
    assert "function suggestedReviewReason(item)" not in javascript
    assert 'item.reviewer || localStorage.getItem("moleculeReviewReviewer")' in javascript
    assert "function navigateCase(delta)" in javascript
    assert "function reviewShortcut(event)" in javascript
    assert 'document.addEventListener("keydown", reviewShortcut)' in javascript
    assert "await loadCase(queuedNextCaseId)" in javascript
