#!/usr/bin/env python3
"""Export molecule review decisions from SQLite."""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path
from typing import Any


APP_DIR = Path(__file__).resolve().parent
REPO_ROOT = APP_DIR.parents[1]
STATE_DIR = REPO_ROOT / ".local" / "molgr_review"
DEFAULT_DB = STATE_DIR / "review.sqlite"
DEFAULT_JSONL = STATE_DIR / "reviews.jsonl"
DEFAULT_WHITELIST = STATE_DIR / "accepted_cases.json"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB, help="SQLite database path.")
    parser.add_argument("--jsonl", type=Path, default=DEFAULT_JSONL, help="Review JSONL output.")
    parser.add_argument(
        "--whitelist",
        type=Path,
        default=DEFAULT_WHITELIST,
        help="Optional JSON summary of accepted candidate/reference decisions.",
    )
    return parser.parse_args()


def _connect(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    return conn


def main() -> None:
    args = _parse_args()
    if not args.db.exists():
        raise SystemExit(f"Database does not exist: {args.db}")

    with _connect(args.db) as conn:
        rows = conn.execute(
            """
            SELECT
                c.case_id, c.row_index, c.source, c.category, c.reference_smiles,
                c.candidate_smiles, c.candidate_organic_smiles,
                r.status, r.corrected_smiles, r.corrected_molblock,
                r.notes, r.reviewer, r.updated_at
            FROM reviews r
            JOIN cases c ON c.case_id = r.case_id
            ORDER BY c.row_index
            """
        ).fetchall()

    reviews: list[dict[str, Any]] = [dict(row) for row in rows]
    args.jsonl.parent.mkdir(parents=True, exist_ok=True)
    with args.jsonl.open("w", encoding="utf-8") as fh:
        for row in reviews:
            fh.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")

    whitelist: dict[str, dict[str, str]] = {}
    for row in reviews:
        if row["status"] == "accept_candidate":
            whitelist[row["case_id"]] = {
                "status": "accepted_reference_issue",
                "reason": row.get("notes") or "Manual review accepted the candidate graph.",
            }
    args.whitelist.write_text(
        json.dumps(whitelist, indent=2, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )

    print(f"Exported {len(reviews)} reviews to {args.jsonl}")
    print(f"Exported {len(whitelist)} accepted-reference entries to {args.whitelist}")


if __name__ == "__main__":
    main()
