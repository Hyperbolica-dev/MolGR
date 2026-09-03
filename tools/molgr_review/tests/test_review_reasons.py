from __future__ import annotations

import json
import sqlite3
import sys
import threading
from http.client import HTTPConnection
from pathlib import Path
from typing import Any


APP_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_DIR))

from server import ReviewHandler, ReviewServer  # noqa: E402


def test_review_reasons_are_distinct_counted_and_refresh_after_review(tmp_path: Path) -> None:
    db_path = tmp_path / "review.sqlite"
    with sqlite3.connect(db_path) as conn:
        conn.executescript((APP_DIR / "schema.sql").read_text(encoding="utf-8"))
        conn.executemany(
            "INSERT INTO cases(case_id, row_index, xyz_path) VALUES(?, ?, '')",
            [("OLD_A", 1), ("OLD_B", 2), ("OLD_C", 3), ("EMPTY", 4), ("SPACE", 5), ("NEW", 6)],
        )
        conn.executemany(
            """
            INSERT INTO reviews(case_id, status, reviewer, notes)
            VALUES(?, 'reference_answer_wrong', ?, '')
            """,
            [
                ("OLD_A", "automated-rule-a"),
                ("OLD_B", "automated-rule-a"),
                ("OLD_C", "automated-rule-b"),
                ("EMPTY", ""),
                ("SPACE", "   "),
            ],
        )
        conn.commit()

    server = ReviewServer(("127.0.0.1", 0), ReviewHandler)
    server.db_path = db_path
    server.xyz_dir = None
    server.fixtures_dir = tmp_path / "reviewed"
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
        status, payload = request("/api/review-reasons")
        assert status == 200
        assert payload["items"] == [
            {"reviewer": "automated-rule-a", "count": 2},
            {"reviewer": "automated-rule-b", "count": 1},
        ]

        status, _ = request(
            "/api/cases/NEW/review",
            {
                "status": "reference_answer_wrong",
                "reviewer": "automated-rule-a",
                "notes": "existing reason evidence",
            },
        )
        assert status == 200
        status, payload = request("/api/review-reasons")
        assert status == 200
        assert payload["items"][0] == {"reviewer": "automated-rule-a", "count": 3}

        with sqlite3.connect(db_path) as conn:
            assert conn.execute(
                "SELECT reviewer, notes FROM reviews WHERE case_id = 'NEW'"
            ).fetchone() == ("automated-rule-a", "existing reason evidence")

        status, _ = request(
            "/api/cases/NEW/review",
            {
                "status": "reference_answer_wrong",
                "reviewer": "xyz-hydrogen-assignment",
                "notes": "new reason evidence remains free text",
            },
        )
        assert status == 200
        status, payload = request("/api/review-reasons")
        assert status == 200
        assert {item["reviewer"]: item["count"] for item in payload["items"]} == {
            "automated-rule-a": 2,
            "automated-rule-b": 1,
            "xyz-hydrogen-assignment": 1,
        }
        with sqlite3.connect(db_path) as conn:
            assert conn.execute(
                "SELECT reviewer, notes FROM reviews WHERE case_id = 'NEW'"
            ).fetchone() == (
                "xyz-hydrogen-assignment",
                "new reason evidence remains free text",
            )
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)


def test_review_reason_control_is_editable_and_saves_only_its_value() -> None:
    html = (APP_DIR / "static" / "index.html").read_text(encoding="utf-8")
    javascript = (APP_DIR / "static" / "app.js").read_text(encoding="utf-8")

    assert 'data-i18n="reviewer">审核理由</span>' in html
    assert 'id="reviewer" type="text" list="reviewReasonOptions"' in html
    assert '<datalist id="reviewReasonOptions"></datalist>' in html
    assert 'data-i18n="notes">备注 / 证据</span>' in html
    assert 'reviewer: state.requiredReviewer || $("reviewer").value' in javascript
    assert 'notes: $("notes").value' in javascript
    assert 'value="${escapeHtml(reviewer)}"' in javascript
    assert "reviewer} (${count})" in javascript
    assert "await Promise.all([loadStats(), loadReviewReasons()])" in javascript
    assert 'class="render-kind' not in html


def test_required_reviewer_accepts_exact_match_and_rejects_other_values(tmp_path: Path) -> None:
    db_path = tmp_path / "review.sqlite"
    with sqlite3.connect(db_path) as conn:
        conn.executescript((APP_DIR / "schema.sql").read_text(encoding="utf-8"))
        conn.executemany(
            "INSERT INTO cases(case_id, row_index, xyz_path) VALUES(?, ?, '')",
            [("MATCH", 1), ("OTHER", 2)],
        )
        conn.commit()

    server = ReviewServer(("127.0.0.1", 0), ReviewHandler)
    server.db_path = db_path
    server.xyz_dir = None
    server.fixtures_dir = tmp_path / "reviewed"
    server.required_reviewer = "human/manual"
    server.runtime_info = {}
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def post(case_id: str, reviewer: str) -> tuple[int, dict[str, Any]]:
        connection = HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        connection.request(
            "POST",
            f"/api/cases/{case_id}/review",
            body=json.dumps({"status": "needs_followup", "reviewer": reviewer, "notes": "test"}),
            headers={"Content-Type": "application/json"},
        )
        response = connection.getresponse()
        payload = json.loads(response.read())
        connection.close()
        return response.status, payload

    try:
        matching_status, matching_payload = post("MATCH", "human/manual")
        other_status, other_payload = post("OTHER", "auto-oxidative-addition")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)

    assert matching_status == 200
    assert matching_payload["ok"] is True
    assert other_status == 400
    assert other_payload == {
        "error": "required_reviewer_mismatch",
        "required_reviewer": "human/manual",
    }
    with sqlite3.connect(db_path) as conn:
        assert conn.execute(
            "SELECT status, reviewer FROM reviews WHERE case_id='MATCH'"
        ).fetchone() == (
            "needs_followup",
            "human/manual",
        )
        assert conn.execute("SELECT COUNT(*) FROM reviews WHERE case_id='OTHER'").fetchone()[0] == 0
