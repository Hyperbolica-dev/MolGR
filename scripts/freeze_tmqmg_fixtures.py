#!/usr/bin/env python3
"""Copy reviewed tmQMg cases into deterministic offline test fixtures."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any


_MANIFEST_FIELDS = (
    "case_id",
    "row_index",
    "charge",
    "metal_center",
    "reference_smiles",
    "classification",
    "reason",
    "xyz_file",
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, required=True, help="tmQMg metadata CSV.")
    parser.add_argument("--xyz-dir", type=Path, required=True, help="tmQMg XYZ directory.")
    parser.add_argument(
        "--spec",
        type=Path,
        default=Path("tests/data/tmqmg/fixture_sources.json"),
        help="Reviewed fixture selection JSON.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("tests/data/tmqmg"),
        help="Fixture output root.",
    )
    return parser.parse_args()


def _load_rows(path: Path) -> tuple[dict[str, tuple[int, dict[str, str]]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"tmQMg CSV has no header: {path}")
        rows = {row["id"]: (row_index, row) for row_index, row in enumerate(reader, start=1)}
        return rows, reader.fieldnames


def _validate_xyz(path: Path, expected_atoms: str) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError(f"XYZ fixture is incomplete: {path}")
    declared_atoms = int(lines[0].strip())
    if expected_atoms and declared_atoms != int(expected_atoms):
        raise ValueError(
            f"XYZ atom count differs from tmQMg metadata for {path.stem}: "
            f"xyz={declared_atoms}, metadata={expected_atoms}"
        )


def freeze_fixtures(
    *,
    csv_path: Path,
    xyz_dir: Path,
    spec_path: Path,
    out_dir: Path,
) -> None:
    rows_by_id, fieldnames = _load_rows(csv_path)
    required_columns = {"id", "charge", "metal_center", "smiles", "n_atoms"}
    missing = required_columns - set(fieldnames)
    if missing:
        raise ValueError("tmQMg CSV is missing columns: " + ", ".join(sorted(missing)))

    fixture_spec: dict[str, list[dict[str, Any]]] = json.loads(
        spec_path.read_text(encoding="utf-8")
    )
    seen_ids: set[str] = set()
    for group, cases in fixture_spec.items():
        group_dir = out_dir / group
        group_dir.mkdir(parents=True, exist_ok=True)
        manifest_rows: list[dict[str, object]] = []
        expected_xyz_files: set[str] = set()
        for case in cases:
            case_id = str(case["id"])
            if case_id in seen_ids:
                raise ValueError(f"Duplicate fixture id across groups: {case_id}")
            seen_ids.add(case_id)
            try:
                row_index, row = rows_by_id[case_id]
            except KeyError as exc:
                raise ValueError(f"tmQMg id is absent from metadata CSV: {case_id}") from exc

            source_xyz = xyz_dir / f"{case_id}.xyz"
            if not source_xyz.exists():
                raise FileNotFoundError(source_xyz)
            _validate_xyz(source_xyz, row["n_atoms"])
            target_xyz = group_dir / source_xyz.name
            shutil.copyfile(source_xyz, target_xyz)
            expected_xyz_files.add(target_xyz.name)
            manifest_rows.append(
                {
                    "case_id": case_id,
                    "row_index": row_index,
                    "charge": row["charge"],
                    "metal_center": row["metal_center"],
                    "reference_smiles": row["smiles"],
                    "classification": case["classification"],
                    "reason": case["reason"],
                    "xyz_file": target_xyz.name,
                }
            )

        for stale_xyz in group_dir.glob("*.xyz"):
            if stale_xyz.name not in expected_xyz_files:
                stale_xyz.unlink()

        with (group_dir / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=_MANIFEST_FIELDS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(manifest_rows)


def main() -> None:
    args = _parse_args()
    freeze_fixtures(
        csv_path=args.csv,
        xyz_dir=args.xyz_dir,
        spec_path=args.spec,
        out_dir=args.out_dir,
    )


if __name__ == "__main__":
    main()
