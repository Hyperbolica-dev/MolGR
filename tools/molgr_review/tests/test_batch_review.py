from __future__ import annotations

import json
import sqlite3
import subprocess
import sys
import threading
from pathlib import Path

from scripts.batch_review import db_state, load_manifest, sqlite_backup


ROOT = Path(__file__).resolve().parents[3]
SCHEMA = ROOT / "tools/molgr_review/schema.sql"
SCRIPT = ROOT / "scripts/batch_review.py"
APP_DIR = ROOT / "tools/molgr_review"
sys.path.insert(0, str(APP_DIR))

from server import ReviewHandler, ReviewServer  # noqa: E402


def test_manifest_is_dry_run_by_default(tmp_path: Path) -> None:
    manifest = tmp_path / "reviews.json"
    manifest.write_text(
        json.dumps(
            {
                "reviews": [
                    {
                        "case_id": "CASE_A",
                        "proposed_status": "reference_answer_wrong",
                        "proposed_reason": "reference-metal-scan",
                        "evidence": "reference parse failed",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    result = subprocess.run(
        ["python", str(SCRIPT), str(manifest)],
        check=True,
        capture_output=True,
        text=True,
    )

    lines = [json.loads(line) for line in result.stdout.splitlines()]
    assert lines[0] == {"mode": "dry-run", "count": 1}
    assert lines[1] == {
        "case_id": "CASE_A",
        "status": "reference_answer_wrong",
        "reason": "reference-metal-scan",
        "notes": "reference parse failed",
    }


def test_manifest_normalizes_reason_without_changing_notes(tmp_path: Path) -> None:
    manifest = tmp_path / "reviews.csv"
    manifest.write_text(
        "case_id,status,reviewer,notes\nCASE_A,accept_both,auto-oxidative-addition,specific evidence\n",
        encoding="utf-8",
    )

    assert load_manifest(manifest)[0] == {
        "case_id": "CASE_A",
        "status": "accept_both",
        "reviewer": "auto-oxidative-addition",
        "notes": "specific evidence",
        "corrected_smiles": "",
        "corrected_molblock": "",
    }


def test_sqlite_backup_is_consistent(tmp_path: Path) -> None:
    source = tmp_path / "review.sqlite"
    backup = tmp_path / "backups/review.sqlite"
    with sqlite3.connect(source) as connection:
        connection.executescript(SCHEMA.read_text(encoding="utf-8"))
        connection.execute("INSERT INTO cases(case_id, row_index, xyz_path) VALUES('A', 1, '')")
        connection.execute(
            "INSERT INTO reviews(case_id, status, reviewer) VALUES('A', 'skip', 'test')"
        )
    sqlite_backup(source, backup)

    assert db_state(source)["integrity"] == "ok"
    assert db_state(backup)["integrity"] == "ok"
    assert db_state(backup)["review_count"] == 1


def test_apply_path_uses_review_api_not_direct_review_sql() -> None:
    source = SCRIPT.read_text(encoding="utf-8")

    assert "/api/cases/" in source
    assert "/review" in source
    assert "INSERT INTO reviews" not in source
    assert "--apply" in source


def test_explicit_apply_backs_up_and_writes_through_review_api(tmp_path: Path) -> None:
    db_path = tmp_path / "review.sqlite"
    with sqlite3.connect(db_path) as connection:
        connection.executescript(SCHEMA.read_text(encoding="utf-8"))
        connection.execute("INSERT INTO cases(case_id, row_index, xyz_path) VALUES('A', 1, '')")
    manifest = tmp_path / "reviews.json"
    manifest.write_text(
        json.dumps(
            {
                "reviews": [
                    {
                        "case_id": "A",
                        "status": "reference_answer_wrong",
                        "reviewer": "reference-metal-scan",
                        "notes": "verified evidence",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    server = ReviewServer(("127.0.0.1", 0), ReviewHandler)
    server.db_path = db_path
    server.xyz_dir = None
    server.fixtures_dir = tmp_path / "fixtures"
    server.triage_records = {}
    server.runtime_info = {}
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        subprocess.run(
            [
                "python",
                str(SCRIPT),
                str(manifest),
                "--apply",
                "--server-url",
                f"http://127.0.0.1:{server.server_port}",
                "--db",
                str(db_path),
                "--backup-dir",
                str(tmp_path / "backups"),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)

    backup = next((tmp_path / "backups").glob("*.sqlite"))
    assert db_state(backup)["review_count"] == 0
    assert db_state(db_path)["review_count"] == 1
    with sqlite3.connect(db_path) as connection:
        assert connection.execute(
            "SELECT status, reviewer, notes FROM reviews WHERE case_id = 'A'"
        ).fetchone() == (
            "reference_answer_wrong",
            "reference-metal-scan",
            "verified evidence",
        )
