"""Build read-only audit tables from the corrected frozen GEOM-Drugs outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


DIFFERENCES = {
    "geom-drugs-63faa9b9935d7e3e-c0": "diazo [N-]=[N+]=C -> neutral NN=C with hydrogen reassignment",
    "geom-drugs-f39d9697995939f4-c0": "two diazo termini reorganize; one forms a five-membered N/O ring (connectivity change)",
    "geom-drugs-bde34b203dc94766-c0": "diazo [N-]=[N+]=C -> neutral NN=C with hydrogen reassignment",
    "geom-drugs-d72399ce8d396be0-c0": "C#S-N -> C-S-N; C-S triple bond becomes single with hydrogen reassignment",
    "geom-drugs-4ef136d31a3dc647-c12": "diazo [N-]=[N+]=C -> neutral NN=C with hydrogen reassignment",
    "geom-drugs-efac12d699d63176-c2": "diazo C=[N+]=[N-] -> C=NN with hydrogen reassignment",
    "geom-drugs-3b4b29becb77f6f6-c15": "diazo [N-]=[N+]=C -> neutral NN=C with hydrogen reassignment",
    "geom-drugs-665aecefc92a9cb0-c0": "ring imine N=C/N=C -> charge-separated N-/N+ single-bond assignment",
    "geom-drugs-fbb8235a7b9a03a5-c0": "diazo C=[N+]=[N-] -> C=NN with hydrogen reassignment",
    "geom-drugs-185aef7dd76c7124-c0": "three-membered N-N-C ring opens to C(=NN) (connectivity change)",
    "geom-drugs-665d56ed02af5779-c0": "diazo C=[N+]=[N-] -> C=NN; sulfone representation also changes",
    "geom-drugs-595b6c15642ebcc7-c0": "azo/enol system cyclizes into a fused charge-separated N ring (connectivity change)",
    "geom-drugs-e69b919f025233df-c162": "diazo [N-]=[N+]=C -> neutral NN=C with hydrogen reassignment",
    "geom-drugs-fa82ea62d315f68c-c0": "diazo C(=[N+]=[N-]) -> C(=NN) with hydrogen reassignment",
    "geom-drugs-d15e9ec58dc7f35c-c0": "diazo C=[N+]=[N-] -> C=NN with hydrogen reassignment",
    "geom-drugs-fdee7a8436151c38-c0": "two ring imines reorganize to aromatic/charge-separated N-/N+ bonds",
    "geom-drugs-4dd7495a9dc742b0-c0": "diazo C=[N+]=[N-] -> C=NN with hydrogen reassignment",
}


def build(run_dir: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with (run_dir / "failures.csv").open(newline="", encoding="utf-8") as handle:
        failures = list(csv.DictReader(handle))
    radical_rows = [row for row in failures if row["radical_consistent"] == "False"]
    if {row["case_id"] for row in radical_rows} != set(DIFFERENCES):
        raise RuntimeError("frozen radical-mismatch set differs from the audited 17 cases")
    radical_columns = (
        "case_id",
        "molecule_id",
        "reference_smiles",
        "predicted_smiles",
        "evaluator_decision",
        "evaluator_relation",
        "evaluator_reason",
        "reference_raw_formal_radical_electrons",
        "predicted_raw_formal_radical_electrons",
        "reference_normalized_radical_electrons",
        "predicted_normalized_radical_electrons",
        "reference_total_charge",
        "predicted_total_charge",
        "relevant_atom_bond_differences",
        "diagnostic_family",
        "review_classification",
    )
    with (output_dir / "radical_mismatches.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=radical_columns, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        for row in radical_rows:
            writer.writerow(
                {
                    **row,
                    "reference_raw_formal_radical_electrons": 0,
                    "predicted_raw_formal_radical_electrons": 0,
                    "reference_normalized_radical_electrons": 0,
                    "predicted_normalized_radical_electrons": 2,
                    "reference_total_charge": 0,
                    "predicted_total_charge": 0,
                    "relevant_atom_bond_differences": DIFFERENCES[row["case_id"]],
                    "diagnostic_family": row["reason_family"],
                    "review_classification": "manual_blocked",
                }
            )

    with (run_dir / "review_cases.csv").open(newline="", encoding="utf-8") as handle:
        review = list(csv.DictReader(handle))
    exception_rows = [row for row in review if row["diagnostic_error"]]
    if len(exception_rows) != 9:
        raise RuntimeError(f"expected 9 secondary exceptions, found {len(exception_rows)}")
    exception_columns = (
        "case_id",
        "molecule_id",
        "reference_smiles",
        "predicted_smiles",
        "reference_preparation",
        "molgr_reconstruction",
        "rdkit_output_finalization",
        "primary_evaluator",
        "secondary_chirality_evaluator",
        "diagnostic_canonicalization",
        "primary_decision",
        "primary_relation",
        "primary_reason",
        "exception_type_message",
        "corrected_classification",
    )
    with (output_dir / "secondary_diagnostic_exceptions.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=exception_columns, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        for row in exception_rows:
            writer.writerow(
                {
                    **row,
                    "reference_preparation": "success",
                    "molgr_reconstruction": "success",
                    "rdkit_output_finalization": "success",
                    "primary_evaluator": "success",
                    "secondary_chirality_evaluator": "exception",
                    "diagnostic_canonicalization": "not reached in original runner",
                    "primary_decision": row["evaluator_decision"],
                    "primary_relation": row["evaluator_relation"],
                    "primary_reason": row["evaluator_reason"],
                    "exception_type_message": row["diagnostic_error"],
                    "corrected_classification": "secondary_chirality_diagnostic_exception",
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    build(args.run_dir, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
