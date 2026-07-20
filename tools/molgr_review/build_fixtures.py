#!/usr/bin/env python3
"""Rebuild test fixtures from all current molecule review decisions."""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path
from typing import Any

from fixture_builder import (
    prune_unreviewed_fixtures,
    sync_review_fixture,
)


APP_DIR = Path(__file__).resolve().parent
REPO_ROOT = APP_DIR.parents[1]
DEFAULT_DB = REPO_ROOT / ".local" / "molgr_review" / "review.sqlite"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB, help="SQLite review database.")
    parser.add_argument(
        "--fixtures-dir",
        type=Path,
        required=True,
        help="Output root for reviewed fixtures.",
    )
    parser.add_argument(
        "--xyz-dir",
        type=Path,
        default=None,
        help="Optional directory used to remap stored XYZ paths by basename.",
    )
    return parser.parse_args()


def _rows(conn: sqlite3.Connection) -> list[dict[str, Any]]:
    conn.row_factory = sqlite3.Row
    return [
        dict(row)
        for row in conn.execute(
            """
            SELECT c.*, r.status, r.corrected_smiles, r.corrected_molblock,
                   r.notes, r.reviewer, r.updated_at
            FROM reviews r
            JOIN cases c ON c.case_id = r.case_id
            ORDER BY c.row_index
            """
        )
    ]


def main() -> None:
    args = _parse_args()
    if not args.db.exists():
        raise SystemExit(f"Database does not exist: {args.db}")

    with sqlite3.connect(args.db) as conn:
        rows = _rows(conn)
    reviewed_ids = {str(row["case_id"]) for row in rows}
    prune_unreviewed_fixtures(args.fixtures_dir, reviewed_ids)

    created = 0
    removed = 0
    failures: list[str] = []
    for row in rows:
        try:
            fixture = sync_review_fixture(
                row,
                row,
                fixtures_dir=args.fixtures_dir,
                xyz_dir=args.xyz_dir,
            )
        except Exception as exc:  # noqa: BLE001
            failures.append(f"{row['case_id']}: {type(exc).__name__}: {exc}")
            continue
        if fixture is None:
            removed += 1
        else:
            created += 1

    print(
        f"Built {created} reviewed fixtures in {args.fixtures_dir}; {removed} non-fixture reviews"
    )
    if failures:
        raise SystemExit("Fixture generation failed:\n" + "\n".join(failures))


if __name__ == "__main__":
    main()
