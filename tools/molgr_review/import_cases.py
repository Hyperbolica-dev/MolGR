#!/usr/bin/env python3
"""Import a molecule-review queue into SQLite."""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


APP_DIR = Path(__file__).resolve().parent
REPO_ROOT = APP_DIR.parents[1]
DEFAULT_DB = REPO_ROOT / ".local" / "molgr_review" / "review.sqlite"
SCHEMA_PATH = APP_DIR / "schema.sql"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Generic review queue CSV.")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB, help="SQLite database path.")
    parser.add_argument("--reset", action="store_true", help="Delete the database first.")
    parser.add_argument(
        "--reviews-jsonl",
        type=Path,
        default=None,
        help="Optional exported review decisions to restore after importing cases.",
    )
    return parser.parse_args()


def _first(row: dict[str, str], *names: str, default: str = "") -> str:
    for name in names:
        value = row.get(name)
        if value not in (None, ""):
            return str(value)
    return default


def _int_value(row: dict[str, str], *names: str, default: int = 0) -> int:
    raw = _first(row, *names)
    if not raw:
        return default
    return int(raw)


def _normalize_row(row: dict[str, str], index: int) -> dict[str, Any]:
    case_id = _first(row, "case_id")
    if not case_id:
        raise ValueError(f"queue row {index} is missing case_id")
    charge = _int_value(row, "total_charge")
    radicals = _int_value(row, "total_radical_electrons")
    spin = _int_value(row, "spin_multiplicity", default=radicals + 1)
    if spin < 1:
        raise ValueError(f"queue row {index} has invalid spin_multiplicity: {spin}")
    xyz_path = _first(row, "xyz_path")
    if not xyz_path:
        raise ValueError(f"queue row {index} is missing xyz_path")
    category = _first(row, "category")
    reserved = {
        "case_id",
        "row_index",
        "source",
        "category",
        "xyz_path",
        "total_charge",
        "total_radical_electrons",
        "spin_multiplicity",
        "reference_smiles",
        "candidate_smiles",
        "candidate_organic_smiles",
        "reference_organic_smiles",
        "candidate_status",
    }
    metadata = {key: value for key, value in row.items() if key not in reserved and value != ""}
    return {
        "case_id": case_id,
        "row_index": _int_value(row, "row_index", default=index),
        "source": _first(row, "source", default="csv"),
        "category": category,
        "xyz_path": xyz_path,
        "total_charge": charge,
        "total_radical_electrons": radicals,
        "spin_multiplicity": spin,
        "reference_smiles": _first(row, "reference_smiles"),
        "candidate_smiles": _first(row, "candidate_smiles"),
        "candidate_organic_smiles": _first(row, "candidate_organic_smiles"),
        "reference_organic_smiles": _first(row, "reference_organic_smiles"),
        "candidate_status": _first(row, "candidate_status"),
        "metadata_json": json.dumps(metadata, ensure_ascii=False, sort_keys=True),
    }


def _prune_stale_cases(conn: sqlite3.Connection, case_ids: list[str]) -> None:
    conn.execute("DROP TABLE IF EXISTS incoming_case_ids")
    conn.execute("CREATE TEMP TABLE incoming_case_ids(case_id TEXT PRIMARY KEY)")
    conn.executemany("INSERT INTO incoming_case_ids(case_id) VALUES(?)", ((value,) for value in case_ids))
    for table in ("render_cache", "reviews", "cases"):
        conn.execute(
            f"DELETE FROM {table} WHERE case_id NOT IN (SELECT case_id FROM incoming_case_ids)"
        )


def _ensure_generic_schema(conn: sqlite3.Connection) -> None:
    """Add the generic columns when reusing a pre-migration SQLite file."""

    columns = {row[1] for row in conn.execute("PRAGMA table_info(cases)")}
    additions = {
        "source": "TEXT NOT NULL DEFAULT ''",
        "total_charge": "INTEGER NOT NULL DEFAULT 0",
        "total_radical_electrons": "INTEGER NOT NULL DEFAULT 0",
        "spin_multiplicity": "INTEGER NOT NULL DEFAULT 1",
        "candidate_smiles": "TEXT NOT NULL DEFAULT ''",
        "candidate_organic_smiles": "TEXT NOT NULL DEFAULT ''",
        "reference_organic_smiles": "TEXT NOT NULL DEFAULT ''",
    }
    for name, declaration in additions.items():
        if name not in columns:
            conn.execute(f"ALTER TABLE cases ADD COLUMN {name} {declaration}")
    if "charge" in columns:
        conn.execute("UPDATE cases SET total_charge = charge WHERE total_charge = 0 AND charge != 0")
    if "candidate_smiles_canonical" in columns:
        conn.execute(
            "UPDATE cases SET candidate_smiles = candidate_smiles_canonical "
            "WHERE candidate_smiles = '' AND candidate_smiles_canonical != ''"
        )
    conn.execute("UPDATE cases SET source = 'legacy' WHERE source = ''")
    if conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'reviews'"
    ).fetchone():
        conn.execute("UPDATE reviews SET status = 'accept_candidate' WHERE status IN ('accept_molgr', 'reference_wrong')")
        conn.execute("UPDATE reviews SET status = 'accept_reference' WHERE status = 'accept_tmqmg'")


def _restore_reviews(conn: sqlite3.Connection, path: Path) -> int:
    restored = 0
    with path.open(encoding="utf-8") as fh:
        for line_number, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            payload = json.loads(line)
            if not isinstance(payload, dict) or not payload.get("case_id"):
                raise ValueError(f"invalid review record at {path}:{line_number}")
            status = str(payload.get("status") or "")
            status = {
                "accept_molgr": "accept_candidate",
                "accept_tmqmg": "accept_reference",
                "reference_wrong": "accept_candidate",
            }.get(status, status)
            conn.execute(
                """
                INSERT INTO reviews(
                    case_id, status, corrected_smiles, corrected_molblock,
                    notes, reviewer, updated_at
                ) VALUES(?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(case_id) DO UPDATE SET
                    status = excluded.status,
                    corrected_smiles = excluded.corrected_smiles,
                    corrected_molblock = excluded.corrected_molblock,
                    notes = excluded.notes,
                    reviewer = excluded.reviewer,
                    updated_at = excluded.updated_at
                """,
                (
                    str(payload["case_id"]),
                    status,
                    str(payload.get("corrected_smiles") or ""),
                    str(payload.get("corrected_molblock") or ""),
                    str(payload.get("notes") or ""),
                    str(payload.get("reviewer") or ""),
                    str(payload.get("updated_at") or ""),
                ),
            )
            restored += 1
    return restored


def main() -> None:
    args = _parse_args()
    if not args.input.is_file():
        raise SystemExit(f"Input queue does not exist: {args.input}")
    if args.reset and args.db.exists():
        args.db.unlink()
    with args.input.open(newline="", encoding="utf-8") as fh:
        rows = [_normalize_row(row, index) for index, row in enumerate(csv.DictReader(fh), start=1)]
    if not rows:
        raise SystemExit(f"Input queue is empty: {args.input}")

    args.db.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(args.db) as conn:
        conn.execute("PRAGMA foreign_keys = ON")
        conn.executescript(SCHEMA_PATH.read_text(encoding="utf-8"))
        _ensure_generic_schema(conn)
        conn.commit()
        conn.execute("BEGIN")
        conn.execute("DELETE FROM render_cache")
        _prune_stale_cases(conn, [row["case_id"] for row in rows])
        conn.executemany(
            """
            INSERT INTO cases(
                case_id, row_index, source, category, xyz_path,
                total_charge, total_radical_electrons, spin_multiplicity,
                reference_smiles, candidate_smiles, candidate_organic_smiles,
                reference_organic_smiles, candidate_status, metadata_json
            ) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(case_id) DO UPDATE SET
                row_index = excluded.row_index,
                source = excluded.source,
                category = excluded.category,
                xyz_path = excluded.xyz_path,
                total_charge = excluded.total_charge,
                total_radical_electrons = excluded.total_radical_electrons,
                spin_multiplicity = excluded.spin_multiplicity,
                reference_smiles = excluded.reference_smiles,
                candidate_smiles = excluded.candidate_smiles,
                candidate_organic_smiles = excluded.candidate_organic_smiles,
                reference_organic_smiles = excluded.reference_organic_smiles,
                candidate_status = excluded.candidate_status,
                metadata_json = excluded.metadata_json
            """,
            [
                (
                    row["case_id"],
                    row["row_index"],
                    row["source"],
                    row["category"],
                    row["xyz_path"],
                    row["total_charge"],
                    row["total_radical_electrons"],
                    row["spin_multiplicity"],
                    row["reference_smiles"],
                    row["candidate_smiles"],
                    row["candidate_organic_smiles"],
                    row["reference_organic_smiles"],
                    row["candidate_status"],
                    row["metadata_json"],
                )
                for row in rows
            ],
        )
        conn.execute("DELETE FROM metadata")
        conn.executemany(
            "INSERT INTO metadata(key, value) VALUES(?, ?)",
            [("source_csv", str(args.input)), ("record_count", str(len(rows))), ("imported_at", datetime.now(timezone.utc).isoformat())],
        )
        restored = _restore_reviews(conn, args.reviews_jsonl) if args.reviews_jsonl else 0
        conn.commit()
    print(f"Imported {len(rows)} cases into {args.db}; restored {restored} reviews")


if __name__ == "__main__":
    main()
