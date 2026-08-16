#!/usr/bin/env python3
"""Preview or explicitly apply review manifests through the running Review API."""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.error import HTTPError
from urllib.parse import quote
from urllib.request import Request, urlopen


ALLOWED_STATUSES = {
    "accept_both",
    "accept_candidate",
    "accept_reference",
    "manual_reference",
    "reference_answer_wrong",
    "needs_followup",
    "skip",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="JSON or CSV review proposal manifest")
    parser.add_argument("--apply", action="store_true", help="Write through the running Review API")
    parser.add_argument("--server-url", default="http://127.0.0.1:8765")
    parser.add_argument("--db", type=Path, help="Review SQLite used for backup and verification")
    parser.add_argument("--backup-dir", type=Path)
    return parser.parse_args()


def load_manifest(path: Path) -> list[dict[str, str]]:
    if path.suffix.lower() == ".csv":
        with path.open(encoding="utf-8", newline="") as handle:
            source = list(csv.DictReader(handle))
    else:
        payload = json.loads(path.read_text(encoding="utf-8"))
        source = payload.get("reviews", []) if isinstance(payload, dict) else payload
    if not isinstance(source, list):
        raise ValueError("manifest_reviews_must_be_a_list")
    reviews = []
    seen = set()
    for raw in source:
        if not isinstance(raw, dict):
            raise ValueError("manifest_review_must_be_an_object")
        case_id = str(raw.get("case_id") or "").strip()
        status = str(raw.get("status") or raw.get("proposed_status") or "").strip()
        reviewer = str(
            raw.get("reviewer") or raw.get("reason") or raw.get("proposed_reason") or ""
        ).strip()
        notes = str(
            raw.get("notes") or raw.get("proposed_notes") or raw.get("evidence") or ""
        ).strip()
        if not case_id:
            raise ValueError("case_id_is_required")
        if case_id in seen:
            raise ValueError(f"duplicate_case_id:{case_id}")
        if status not in ALLOWED_STATUSES:
            raise ValueError(f"invalid_status:{case_id}:{status}")
        corrected_smiles = str(raw.get("corrected_smiles") or "").strip()
        corrected_molblock = str(raw.get("corrected_molblock") or "")
        if status == "manual_reference" and not corrected_smiles:
            raise ValueError(f"corrected_smiles_required:{case_id}")
        reviews.append(
            {
                "case_id": case_id,
                "status": status,
                "reviewer": reviewer,
                "notes": notes,
                "corrected_smiles": corrected_smiles,
                "corrected_molblock": corrected_molblock,
            }
        )
        seen.add(case_id)
    return reviews


def db_state(path: Path) -> dict[str, Any]:
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        integrity = connection.execute("PRAGMA integrity_check").fetchone()[0]
        return {
            "integrity": integrity,
            "case_count": connection.execute("SELECT COUNT(*) FROM cases").fetchone()[0],
            "review_count": connection.execute("SELECT COUNT(*) FROM reviews").fetchone()[0],
            "reviewed_case_ids": {
                row[0] for row in connection.execute("SELECT case_id FROM reviews")
            },
        }
    finally:
        connection.close()


def sqlite_backup(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    source_connection = sqlite3.connect(f"file:{source}?mode=ro", uri=True)
    destination_connection = sqlite3.connect(destination)
    try:
        source_connection.backup(destination_connection)
    finally:
        destination_connection.close()
        source_connection.close()


def post_review(server_url: str, review: dict[str, str]) -> dict[str, Any]:
    case_id = review["case_id"]
    body = json.dumps({key: value for key, value in review.items() if key != "case_id"}).encode()
    request = Request(
        f"{server_url.rstrip('/')}/api/cases/{quote(case_id, safe='')}/review",
        data=body,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urlopen(request, timeout=120) as response:
            payload = json.loads(response.read())
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"review_api_failed:{case_id}:HTTP_{exc.code}:{detail}") from exc
    if not payload.get("ok"):
        raise RuntimeError(f"review_api_failed:{case_id}:{payload}")
    return payload


def server_storage(server_url: str) -> dict[str, str]:
    with urlopen(f"{server_url.rstrip('/')}/api/stats", timeout=30) as response:
        payload = json.loads(response.read())
    storage = payload.get("storage") or {}
    return {str(key): str(value) for key, value in storage.items()}


def dry_run_payload(review: dict[str, str]) -> dict[str, str]:
    return {
        "case_id": review["case_id"],
        "status": review["status"],
        "reason": review["reviewer"],
        "notes": review["notes"],
    }


def main() -> None:
    args = parse_args()
    reviews = load_manifest(args.manifest)
    print(json.dumps({"mode": "apply" if args.apply else "dry-run", "count": len(reviews)}))
    for review in reviews:
        print(json.dumps(dry_run_payload(review), ensure_ascii=False))
    if not args.apply:
        return
    if args.db is None:
        raise SystemExit("--db is required with --apply")
    before = db_state(args.db)
    if before["integrity"] != "ok":
        raise SystemExit(f"source database integrity check failed: {before['integrity']}")
    storage = server_storage(args.server_url)
    if Path(storage.get("review_db", "")).resolve() != args.db.resolve():
        raise SystemExit(
            f"server/database mismatch: server={storage.get('review_db') or 'unknown'} "
            f"requested={args.db.resolve()}"
        )
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup_dir = args.backup_dir or args.db.parent / "batch-review-backups"
    backup_path = backup_dir / f"{args.db.stem}_{timestamp}.sqlite"
    if backup_path.exists():
        raise SystemExit(f"backup already exists: {backup_path}")
    sqlite_backup(args.db, backup_path)
    backup_state = db_state(backup_path)
    if backup_state["integrity"] != "ok" or backup_state["review_count"] != before["review_count"]:
        raise SystemExit("backup verification failed")
    print(json.dumps({"backup": str(backup_path), "review_count": before["review_count"]}))
    for review in reviews:
        post_review(args.server_url, review)
    after = db_state(args.db)
    expected_new = sum(review["case_id"] not in before["reviewed_case_ids"] for review in reviews)
    expected_count = before["review_count"] + expected_new
    if after["integrity"] != "ok" or after["review_count"] != expected_count:
        raise SystemExit(
            f"post-apply verification failed: integrity={after['integrity']} "
            f"reviews={after['review_count']} expected={expected_count}; backup={backup_path}"
        )
    print(
        json.dumps(
            {
                "applied": len(reviews),
                "review_count_before": before["review_count"],
                "review_count_after": after["review_count"],
                "integrity": after["integrity"],
                "backup": str(backup_path),
            }
        )
    )


if __name__ == "__main__":
    main()
