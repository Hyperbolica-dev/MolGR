"""Build deterministic, evidence-gated review queues from a GEOM pilot result table."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
from typing import Any

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize


RDLogger.DisableLog("rdApp.*")


def _tautomer_diagnostic(reference_smiles: str, predicted_smiles: str) -> bool:
    reference = Chem.MolFromSmiles(reference_smiles)
    predicted = Chem.MolFromSmiles(predicted_smiles)
    if reference is None or predicted is None:
        return False
    enumerator = rdMolStandardize.TautomerEnumerator()
    reference_key = Chem.MolToSmiles(enumerator.Canonicalize(reference), isomericSmiles=False)
    predicted_key = Chem.MolToSmiles(enumerator.Canonicalize(predicted), isomericSmiles=False)
    return reference_key == predicted_key


def build(results_path: Path, output_path: Path, *, runtime_count: int = 20) -> None:
    with results_path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    failure_ids = {
        row["case_id"]
        for row in rows
        if row["status"] != "ok" or row["evaluator_decision"] != "equivalent"
    }
    runtime_ids = {
        row["case_id"]
        for row in sorted(rows, key=lambda item: float(item["runtime_ms"]), reverse=True)[
            :runtime_count
        ]
    }
    review_rows: list[dict[str, Any]] = []
    for row in rows:
        case_id = row["case_id"]
        if case_id not in failure_ids | runtime_ids:
            continue
        flags: list[str] = []
        if row["status"] != "ok":
            flags.append("reconstruction_failure")
        if row["evaluator_decision"] == "not_equivalent":
            flags.append("non_equivalent")
        if row["evaluator_decision"] == "inconclusive":
            flags.append("inconclusive")
        if case_id in failure_ids and int(row["heavy_atom_count"]) >= 51:
            flags.append("large_molecule_failure")
        suspected_tautomer = row["evaluator_decision"] == "not_equivalent" and _tautomer_diagnostic(
            row["reference_smiles"], row["predicted_smiles"]
        )
        if suspected_tautomer:
            flags.append("suspected_tautomer_protomer")
        if case_id in runtime_ids:
            flags.append("runtime_outlier_top20")

        graph_review = case_id in failure_ids
        if suspected_tautomer:
            reason = "tautomer canonicalization agrees; evaluator v1 remains authoritative"
        elif row["evaluator_decision"] == "inconclusive":
            reason = "evaluator v1 could not prove equivalence; representation review is blocked"
        elif row["evaluator_decision"] == "not_equivalent":
            reason = "evaluator v1 reports a strict graph mismatch"
        elif row["status"] != "ok":
            reason = "reconstruction did not produce an evaluable candidate"
        else:
            reason = "equivalent case selected only as a runtime outlier"
        signature_payload = f"{row['reference_smiles']}>>{row['predicted_smiles']}"
        review_rows.append(
            {
                "priority": ";".join(flags),
                "case_id": case_id,
                "molecule_id": row["molecule_id"],
                "conformer_id": row["conformer_id"],
                "heavy_atom_count": row["heavy_atom_count"],
                "runtime_ms": row["runtime_ms"],
                "status": row["status"],
                "evaluator_decision": row["evaluator_decision"],
                "evaluator_relation": row["evaluator_relation"],
                "reference_smiles": row["reference_smiles"],
                "predicted_smiles": row["predicted_smiles"],
                "proposed_verdict": "manual_blocked" if graph_review else "no_graph_review_needed",
                "proposed_reason": reason,
                "confidence": "medium" if suspected_tautomer else "high",
                "matched_rule": "",
                "canonical_signature": hashlib.sha256(signature_payload.encode()).hexdigest()[:20],
                "evidence_summary": (
                    f"status={row['status']}; decision={row['evaluator_decision']}; "
                    f"relation={row['evaluator_relation']}; charge={row['charge_consistent']}; "
                    f"radicals={row['radical_consistent']}"
                ),
                "blockers": (
                    "no_verified_candidate_to_xyz_mapping;no_approved_representation_rule"
                    if graph_review
                    else ""
                ),
                "needs_human_review": graph_review,
            }
        )
    review_rows.sort(
        key=lambda row: (
            row["status"] == "ok",
            row["evaluator_decision"] != "not_equivalent",
            row["evaluator_decision"] != "inconclusive",
            -int(row["heavy_atom_count"]),
            -float(row["runtime_ms"]),
            row["case_id"],
        )
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(review_rows[0]) if review_rows else []
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(review_rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--runtime-count", type=int, default=20)
    args = parser.parse_args()
    build(args.results, args.output, runtime_count=args.runtime_count)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
