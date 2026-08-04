#!/usr/bin/env python3
# pyright: reportCallIssue=false
"""Validate tmQMg reference graphs against XYZ element counts.

This script checks whether the dataset reference SMILES can reproduce the
element-count formula present in the original XYZ geometry.

For each entry it:
1. reads the XYZ file and counts atoms by element symbol,
2. parses the reference SMILES with RDKit,
3. adds hydrogens to the reference graph, and
4. compares the per-element counts exactly.

The output CSV records the raw counts and mismatch details, while the summary
JSON aggregates status counts for quick triage of obviously wrong references.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any, cast


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rdkit import Chem, RDLogger


RDLogger.DisableLog("rdApp.*")  # type: ignore[arg-type]

RESULT_FIELDNAMES = (
    "row_index",
    "id",
    "metal_center",
    "charge",
    "xyz_path",
    "reference_smiles_input",
    "reference_parse_status",
    "formula_match",
    "xyz_atom_count",
    "reference_atom_count_with_h",
    "xyz_formula",
    "reference_formula_with_h",
    "mismatch_detail",
    "elapsed_seconds",
    "error",
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, required=True, help="tmQMg metadata CSV path.")
    parser.add_argument(
        "--xyz-dir",
        type=Path,
        required=True,
        help="Directory containing tmQMg XYZ files named <id>.xyz.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=100,
        help="How many CSV rows to process from the top of the file.",
    )
    parser.add_argument("--out", type=Path, required=True, help="Output CSV path.")
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=None,
        help="Optional summary JSON path. Defaults to <out>.summary.json.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=10,
        help="Print one stderr progress line every N processed entries. Use 0 to silence.",
    )
    return parser.parse_args()


def _summary_path_from_output(output_path: Path) -> Path:
    if output_path.suffix:
        return output_path.with_suffix(f"{output_path.suffix}.summary.json")
    return output_path.with_name(f"{output_path.name}.summary.json")


def _empty_result(row_index: int, row: dict[str, str], xyz_path: Path) -> dict[str, Any]:
    return {
        "row_index": row_index,
        "id": row.get("id", ""),
        "metal_center": row.get("metal_center", ""),
        "charge": row.get("charge", ""),
        "xyz_path": str(xyz_path),
        "reference_smiles_input": row.get("smiles", ""),
        "reference_parse_status": "",
        "formula_match": "",
        "xyz_atom_count": "",
        "reference_atom_count_with_h": "",
        "xyz_formula": "",
        "reference_formula_with_h": "",
        "mismatch_detail": "",
        "elapsed_seconds": "",
        "error": "",
    }


def _read_xyz_element_counts(xyz_path: Path) -> Counter[str]:
    lines = xyz_path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError(f"XYZ file is too short: {xyz_path}")
    try:
        atom_count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"Invalid XYZ atom-count header in {xyz_path}") from exc

    counts: Counter[str] = Counter()
    coordinate_lines = [line for line in lines[2:] if line.strip()]
    if len(coordinate_lines) < atom_count:
        raise ValueError(
            f"XYZ coordinate lines shorter than declared atom count: {xyz_path} "
            f"(declared={atom_count}, found={len(coordinate_lines)})"
        )
    for line in coordinate_lines[:atom_count]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Malformed XYZ line in {xyz_path}: {line!r}")
        counts[parts[0]] += 1
    return counts


def _smiles_element_counts_with_h(reference_smiles: str) -> Counter[str]:
    mol = Chem.MolFromSmiles(reference_smiles)
    if mol is None:
        # Some RDKit releases reject the unusual hypervalent metal notation
        # during sanitization even though the graph is still usable for the
        # formula audit.  Keep the unsanitized graph and calculate implicit
        # hydrogens with a non-strict property-cache update.
        mol = Chem.MolFromSmiles(reference_smiles, sanitize=False)
        if mol is None:
            raise ValueError("reference_parse_failed")
        mol.UpdatePropertyCache(strict=False)
    mol_h = Chem.AddHs(mol)
    return Counter(atom.GetSymbol() for atom in mol_h.GetAtoms())


def _formula_string(counts: Counter[str]) -> str:
    return " ".join(f"{symbol}:{counts[symbol]}" for symbol in sorted(counts))


def _mismatch_detail(
    xyz_counts: Counter[str],
    reference_counts: Counter[str],
) -> str:
    deltas: list[str] = []
    for symbol in sorted(set(xyz_counts) | set(reference_counts)):
        xyz_count = xyz_counts.get(symbol, 0)
        reference_count = reference_counts.get(symbol, 0)
        if xyz_count == reference_count:
            continue
        deltas.append(f"{symbol}:xyz={xyz_count},ref={reference_count}")
    return "; ".join(deltas)


def _process_row(
    row_index: int,
    row: dict[str, str],
    *,
    xyz_dir: Path,
) -> dict[str, Any]:
    started_at = time.perf_counter()
    xyz_path = xyz_dir / f"{row['id']}.xyz"
    result = _empty_result(row_index, row, xyz_path)

    if not xyz_path.exists():
        result["reference_parse_status"] = "xyz_missing"
        result["error"] = f"Missing XYZ file: {xyz_path}"
        result["elapsed_seconds"] = round(time.perf_counter() - started_at, 6)
        return result

    try:
        xyz_counts = _read_xyz_element_counts(xyz_path)
    except Exception as exc:  # noqa: BLE001
        result["reference_parse_status"] = f"xyz_failed:{type(exc).__name__}"
        result["error"] = str(exc)
        result["elapsed_seconds"] = round(time.perf_counter() - started_at, 6)
        return result

    result["xyz_atom_count"] = sum(xyz_counts.values())
    result["xyz_formula"] = _formula_string(xyz_counts)

    reference_smiles = row.get("smiles", "")
    if not reference_smiles:
        result["reference_parse_status"] = "missing_reference_smiles"
        result["elapsed_seconds"] = round(time.perf_counter() - started_at, 6)
        return result

    try:
        reference_counts = _smiles_element_counts_with_h(reference_smiles)
    except ValueError as exc:
        if str(exc) == "reference_parse_failed":
            result["reference_parse_status"] = "reference_parse_failed"
            result["elapsed_seconds"] = round(time.perf_counter() - started_at, 6)
            return result
        raise
    except Exception as exc:  # noqa: BLE001
        result["reference_parse_status"] = f"reference_failed:{type(exc).__name__}"
        result["error"] = str(exc)
        result["elapsed_seconds"] = round(time.perf_counter() - started_at, 6)
        return result

    result["reference_parse_status"] = "ok"
    result["reference_atom_count_with_h"] = sum(reference_counts.values())
    result["reference_formula_with_h"] = _formula_string(reference_counts)

    formula_match = xyz_counts == reference_counts
    result["formula_match"] = formula_match
    if not formula_match:
        result["mismatch_detail"] = _mismatch_detail(xyz_counts, reference_counts)

    result["elapsed_seconds"] = round(time.perf_counter() - started_at, 6)
    return result


def _counter_to_dict(counter: Counter[str]) -> dict[str, int]:
    return dict(sorted(counter.items()))


def main() -> int:
    args = _parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    summary_out = args.summary_out or _summary_path_from_output(args.out)
    summary_out.parent.mkdir(parents=True, exist_ok=True)

    parse_status_counter: Counter[str] = Counter()
    formula_match_counter: Counter[str] = Counter()
    mismatch_counter: Counter[str] = Counter()
    processed = 0
    formula_match_count = 0
    started_at = time.perf_counter()

    with args.csv.open(newline="") as input_fh, args.out.open("w", newline="") as output_fh:
        reader = csv.DictReader(input_fh)
        writer = csv.DictWriter(output_fh, fieldnames=RESULT_FIELDNAMES)
        writer.writeheader()

        for row_index, row in enumerate(reader, start=1):
            if processed >= args.limit:
                break

            result = _process_row(row_index, row, xyz_dir=args.xyz_dir)
            writer.writerow(cast(Any, result))
            processed += 1

            parse_status_counter.update([result["reference_parse_status"] or "missing"])
            formula_match_value = result["formula_match"]
            if formula_match_value != "":
                formula_match_counter.update([str(formula_match_value)])
                if formula_match_value is True:
                    formula_match_count += 1
            if result["mismatch_detail"]:
                mismatch_counter.update([result["mismatch_detail"]])

            if args.progress_every > 0 and processed % args.progress_every == 0:
                print(
                    f"[tmQMg-formula] processed {processed} rows; latest id={result['id']} "
                    f"formula_match={result['formula_match']}",
                    file=sys.stderr,
                )

    total_elapsed_seconds = time.perf_counter() - started_at
    summary = {
        "input_csv": str(args.csv),
        "xyz_dir": str(args.xyz_dir),
        "limit": args.limit,
        "processed": processed,
        "formula_match_count": formula_match_count,
        "formula_match_fraction": (formula_match_count / processed) if processed else 0.0,
        "reference_parse_status_counts": _counter_to_dict(parse_status_counter),
        "formula_match_value_counts": _counter_to_dict(formula_match_counter),
        "mismatch_detail_counts": _counter_to_dict(mismatch_counter),
        "total_elapsed_seconds": round(total_elapsed_seconds, 6),
        "output_csv": str(args.out),
    }

    with summary_out.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(
        f"Wrote {processed} tmQMg formula-check rows to {args.out} and summary to {summary_out}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
