from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from benchmarks.bde_db_benchmark.adapter import load_bde_cases_by_record_index
from benchmarks.bde_db_benchmark.run import BDEResult, _write_csv, run_cases


OLD_METHOD_RELATIONS = {
    "ideal": "normalized_graph_identity",
    "resonance": "resonance_equivalence",
    "carbene_zwitterion": "carbene_zwitterion_equivalence",
    "inchi_key": "identifier_equivalence",
}


@dataclass(frozen=True)
class RegressionComparison:
    case_id: str
    source_record_index: int
    comparison_category: str
    prediction_changed: bool | None
    evaluation_changed: bool | None
    old_predicted_smiles: str | None
    new_predicted_smiles: str | None
    old_equivalent: bool | None
    new_decision: str | None
    old_method: str | None
    old_relation: str | None
    new_method: str | None
    new_relation: str | None
    new_inconclusive: bool | None
    new_failure: bool
    old_formal_radical_atom_index_match: bool | None
    new_formal_radical_atom_index_match: bool | None
    formal_radical_localization_changed: bool | None
    review_required: bool
    review_reasons: str


def _optional_bool(value: str | None) -> bool | None:
    if value == "True":
        return True
    if value == "False":
        return False
    return None


def _read_old_results(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"old results contain no cases: {path}")
    indices = [int(row["source_record_index"]) for row in rows]
    if len(indices) != len(set(indices)):
        raise ValueError("old results contain duplicate source_record_index values")
    return rows


def _compare(
    old_rows: list[dict[str, str]],
    new_results: list[BDEResult],
) -> list[RegressionComparison]:
    old_by_index = {int(row["source_record_index"]): row for row in old_rows}
    comparisons: list[RegressionComparison] = []
    for result in new_results:
        old = old_by_index.get(result.source_record_index)
        if old is None:
            comparisons.append(
                RegressionComparison(
                    case_id=result.case_id,
                    source_record_index=result.source_record_index,
                    comparison_category="new_case",
                    prediction_changed=None,
                    evaluation_changed=None,
                    old_predicted_smiles=None,
                    new_predicted_smiles=result.predicted_smiles,
                    old_equivalent=None,
                    new_decision=result.evaluator_decision,
                    old_method=None,
                    old_relation=None,
                    new_method=result.equivalence_method,
                    new_relation=result.evaluator_relation,
                    new_inconclusive=result.evaluator_inconclusive,
                    new_failure=result.status != "ok",
                    old_formal_radical_atom_index_match=None,
                    new_formal_radical_atom_index_match=result.formal_radical_atom_index_match,
                    formal_radical_localization_changed=None,
                    review_required=True,
                    review_reasons="new_explicit_edge_case",
                )
            )
            continue
        old_prediction = old.get("predicted_smiles") or None
        old_bonds = old.get("predicted_bonds") or None
        prediction_changed = (
            old_prediction != result.predicted_smiles or old_bonds != result.predicted_bonds
        )
        old_equivalent = _optional_bool(old.get("equivalent"))
        old_method = old.get("equivalence_method") or None
        old_relation = OLD_METHOD_RELATIONS.get(old_method or "")
        new_equivalent = (
            True
            if result.evaluator_decision == "equivalent"
            else False
            if result.evaluator_decision == "not_equivalent"
            else None
        )
        evaluation_changed = (
            old_equivalent != new_equivalent
            or old_method != result.equivalence_method
            or old_relation != result.evaluator_relation
        )
        category = (
            f"prediction_{'changed' if prediction_changed else 'unchanged'}"
            f"_evaluation_{'changed' if evaluation_changed else 'unchanged'}"
        )
        old_radical_match = _optional_bool(old.get("radical_site_consistent"))
        new_radical_match = result.formal_radical_atom_index_match
        radical_changed = old_radical_match != new_radical_match
        reasons = []
        if prediction_changed:
            reasons.append("prediction_changed")
        if evaluation_changed:
            reasons.append("evaluation_changed")
        if result.evaluator_inconclusive:
            reasons.append("inconclusive")
        if result.evaluator_decision == "not_equivalent":
            reasons.append("non_equivalent")
        if result.status != "ok":
            reasons.append("failure")
        if radical_changed:
            reasons.append("formal_radical_localization_changed")
        comparisons.append(
            RegressionComparison(
                case_id=result.case_id,
                source_record_index=result.source_record_index,
                comparison_category=category,
                prediction_changed=prediction_changed,
                evaluation_changed=evaluation_changed,
                old_predicted_smiles=old_prediction,
                new_predicted_smiles=result.predicted_smiles,
                old_equivalent=old_equivalent,
                new_decision=result.evaluator_decision,
                old_method=old_method,
                old_relation=old_relation,
                new_method=result.equivalence_method,
                new_relation=result.evaluator_relation,
                new_inconclusive=result.evaluator_inconclusive,
                new_failure=result.status != "ok",
                old_formal_radical_atom_index_match=old_radical_match,
                new_formal_radical_atom_index_match=new_radical_match,
                formal_radical_localization_changed=radical_changed,
                review_required=bool(reasons),
                review_reasons=";".join(reasons),
            )
        )
    return comparisons


def _regression_summary(comparisons: list[RegressionComparison]) -> str:
    categories = Counter(row.comparison_category for row in comparisons)
    review = [row for row in comparisons if row.review_required]
    review_lines = "\n".join(
        f"- `{row.case_id}` (record {row.source_record_index}): {row.review_reasons}; "
        f"old/new prediction `{row.old_predicted_smiles}` / `{row.new_predicted_smiles}`; "
        f"old/new evaluation `{row.old_equivalent}:{row.old_method}` / "
        f"`{row.new_decision}:{row.new_method}:{row.new_relation}`"
        for row in review
    )
    category_names = [
        "prediction_unchanged_evaluation_unchanged",
        "prediction_unchanged_evaluation_changed",
        "prediction_changed_evaluation_unchanged",
        "prediction_changed_evaluation_changed",
        "new_case",
    ]
    counts = "\n".join(f"- {name}: {categories[name]}" for name in category_names)
    return f"""# BDE-db 101-case old-vs-new regression

The first 100 cases are replayed by exact `source_record_index` from the old result CSV.
The isolated hydrogen radical `[H]` is appended as the 101st case. Runtime is diagnostic only.

## Categories

{counts}
- new inconclusive: {sum(row.new_inconclusive is True for row in comparisons)}
- new non-equivalent: {sum(row.new_decision == "not_equivalent" for row in comparisons)}
- new failure: {sum(row.new_failure for row in comparisons)}
- formal-radical localization changed: {sum(row.formal_radical_localization_changed is True for row in comparisons)}
- manual review required: {len(review)}

## Manual review

{review_lines or "- None."}
"""


def run_regression(
    input_path: Path,
    old_results_path: Path,
    out_dir: Path,
    *,
    timeout_seconds: float | None = 5.0,
) -> tuple[list[BDEResult], list[RegressionComparison]]:
    old_rows = _read_old_results(old_results_path)
    old_indices = [int(row["source_record_index"]) for row in old_rows]
    cases, diagnostics = load_bde_cases_by_record_index(
        input_path,
        old_indices,
        required_smiles=["[H]"],
    )
    results = run_cases(
        cases,
        diagnostics,
        input_path,
        out_dir,
        timeout_seconds=timeout_seconds,
        review_limit=len(cases),
    )
    comparisons = _compare(old_rows, results)
    fields = list(RegressionComparison.__dataclass_fields__)
    _write_csv(
        out_dir / "old_vs_new.csv",
        [{field: getattr(row, field) for field in fields} for row in comparisons],
        fields,
    )
    (out_dir / "old_vs_new_summary.md").write_text(
        _regression_summary(comparisons), encoding="utf-8"
    )
    return results, comparisons


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Replay the BDE-db old 100-case regression.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--old-results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--case-timeout-seconds", type=float, default=5.0)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    run_regression(
        args.input,
        args.old_results,
        args.out,
        timeout_seconds=args.case_timeout_seconds or None,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
