from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
import time
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from rdkit import Chem

from benchmarks._timeout import CaseTimeoutError, case_timeout
from benchmarks.bde_db_benchmark.adapter import (
    EXPECTED_FILENAME,
    BDECase,
    LoadDiagnostics,
    load_bde_cases,
)
from benchmarks.smiles_xyz_benchmark.methods.base import BenchmarkMethod
from benchmarks.smiles_xyz_benchmark.methods.molgr_cpp import MolGRCppMethod
from benchmarks.smiles_xyz_benchmark.methods.postprocess import remove_hs_without_sanitize
from molgr.utils.equivalence import EquivalenceDecision, evaluate_equivalence


@dataclass(frozen=True)
class BDEResult:
    case_id: str
    source_record_index: int
    stratum: str
    reference_smiles: str
    predicted_smiles: str | None
    xyz: str
    total_charge: int
    spin_multiplicity: int
    radical_site: int | None
    reconstruction_success: bool
    status: str
    failure_kind: str | None
    error: str | None
    equivalent: bool | None
    evaluator_decision: str | None
    evaluator_relation: str | None
    evaluator_reason: str | None
    equivalence_method: str | None
    evaluator_inconclusive: bool | None
    bounded_search_attempted: bool | None
    bounded_search_limit: int | None
    bounded_search_limit_reached: bool | None
    bounded_search_exhaustive: bool | None
    bounded_search_candidate_count: int | None
    bounded_search_reference_count: int | None
    exact_smiles_match: bool | None
    charge_consistent: bool | None
    radical_electron_consistent: bool | None
    atom_order_preserved: bool | None
    atom_identity_guard_reason: str | None
    formal_radical_atom_index_match: bool | None
    reference_bonds: str
    predicted_bonds: str | None
    runtime_ms: float
    source_metadata: str


def _total_charge(mol: Chem.Mol) -> int:
    return sum(int(atom.GetFormalCharge()) for atom in mol.GetAtoms())


def _radical_electrons(mol: Chem.Mol) -> int:
    return sum(int(atom.GetNumRadicalElectrons()) for atom in mol.GetAtoms())


def _radical_sites(mol: Chem.Mol) -> list[int]:
    return [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetNumRadicalElectrons()]


def _retained_reference_atom_indices(mol: Chem.Mol) -> list[int]:
    return [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() != 1 or atom.GetDegree() == 0
    ]


def _atom_identity_mapping(
    reference_mol: Chem.Mol,
    predicted_mol: Chem.Mol,
    *,
    coordinate_tolerance: float = 1e-6,
) -> tuple[bool, dict[int, int]]:
    preserved, mapping, _ = _atom_identity_mapping_with_reason(
        reference_mol,
        predicted_mol,
        coordinate_tolerance=coordinate_tolerance,
    )
    return preserved, mapping


def _atom_identity_mapping_with_reason(
    reference_mol: Chem.Mol,
    predicted_mol: Chem.Mol,
    *,
    coordinate_tolerance: float = 1e-6,
) -> tuple[bool, dict[int, int], str]:
    reference_indices = _retained_reference_atom_indices(reference_mol)
    if predicted_mol.GetNumAtoms() != len(reference_indices):
        return False, {}, "retained atom count differs"
    if reference_mol.GetNumConformers() != 1 or predicted_mol.GetNumConformers() != 1:
        return False, {}, "expected exactly one conformer on candidate and reference"
    reference_conf = reference_mol.GetConformer()
    predicted_conf = predicted_mol.GetConformer()
    mapping: dict[int, int] = {}
    for predicted_index, reference_index in enumerate(reference_indices):
        reference_atom = reference_mol.GetAtomWithIdx(reference_index)
        predicted_atom = predicted_mol.GetAtomWithIdx(predicted_index)
        if reference_atom.GetAtomicNum() != predicted_atom.GetAtomicNum():
            return False, {}, f"element differs at retained atom position {predicted_index}"
        reference_position = reference_conf.GetAtomPosition(reference_index)
        predicted_position = predicted_conf.GetAtomPosition(predicted_index)
        distance = math.sqrt(
            (reference_position.x - predicted_position.x) ** 2
            + (reference_position.y - predicted_position.y) ** 2
            + (reference_position.z - predicted_position.z) ** 2
        )
        if distance > coordinate_tolerance:
            return False, {}, f"coordinate differs at retained atom position {predicted_index}"
        mapping[predicted_index] = reference_index
    return True, mapping, "passed"


def _bonds(mol: Chem.Mol) -> str:
    rows = [
        [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), str(bond.GetBondType())]
        for bond in mol.GetBonds()
    ]
    return json.dumps(rows, separators=(",", ":"))


def _error_result(
    case: BDECase,
    *,
    failure_kind: str,
    error: str,
    runtime_ms: float,
    predicted_smiles: str | None = None,
) -> BDEResult:
    return BDEResult(
        case_id=case.case_id,
        source_record_index=case.source_record_index,
        stratum=case.stratum,
        reference_smiles=case.reference_smiles,
        predicted_smiles=predicted_smiles,
        xyz=case.xyz,
        total_charge=case.total_charge,
        spin_multiplicity=case.spin_multiplicity,
        radical_site=case.radical_site,
        reconstruction_success=False,
        status="error",
        failure_kind=failure_kind,
        error=error,
        equivalent=None,
        evaluator_decision=None,
        evaluator_relation=None,
        evaluator_reason=None,
        equivalence_method=None,
        evaluator_inconclusive=None,
        bounded_search_attempted=None,
        bounded_search_limit=None,
        bounded_search_limit_reached=None,
        bounded_search_exhaustive=None,
        bounded_search_candidate_count=None,
        bounded_search_reference_count=None,
        exact_smiles_match=None,
        charge_consistent=None,
        radical_electron_consistent=None,
        atom_order_preserved=None,
        atom_identity_guard_reason=None,
        formal_radical_atom_index_match=None,
        reference_bonds=_bonds(case.reference_mol),
        predicted_bonds=None,
        runtime_ms=runtime_ms,
        source_metadata=json.dumps(case.source_metadata, ensure_ascii=True, sort_keys=True),
    )


def _run_case(
    case: BDECase,
    method: BenchmarkMethod,
    timeout_seconds: float | None,
) -> BDEResult:
    started = time.perf_counter()
    try:
        with case_timeout(timeout_seconds, f"molgr_cpp case {case.case_id}"):
            output = method.run(case.to_method_case())
    except CaseTimeoutError as exc:
        return _error_result(
            case,
            failure_kind="timeout",
            error=str(exc),
            runtime_ms=(time.perf_counter() - started) * 1000.0,
        )
    except Exception as exc:
        return _error_result(
            case,
            failure_kind="exception",
            error=f"method.run failed: {type(exc).__name__}: {exc}",
            runtime_ms=(time.perf_counter() - started) * 1000.0,
        )

    runtime_ms = (time.perf_counter() - started) * 1000.0
    predicted_mol = output.rdkit_mol
    if output.status != "ok" or predicted_mol is None:
        return _error_result(
            case,
            failure_kind="method_error",
            error=output.error or "molgr_cpp returned no molecule",
            runtime_ms=runtime_ms,
            predicted_smiles=output.predicted_smiles,
        )

    try:
        evaluator_reference_mol = remove_hs_without_sanitize(case.reference_mol)
        with case_timeout(timeout_seconds, f"equivalence case {case.case_id}"):
            equivalence_info = evaluate_equivalence(
                predicted_mol,
                evaluator_reference_mol,
                use_chirality=False,
                max_resonance=100,
            )
        equivalent = (
            True
            if equivalence_info.decision == EquivalenceDecision.EQUIVALENT
            else False
            if equivalence_info.decision == EquivalenceDecision.NOT_EQUIVALENT
            else None
        )
        status = "ok"
        error = None
        evaluator_failure_kind = None
        evaluator_decision = equivalence_info.decision.value
        evaluator_relation = equivalence_info.relation.value
        evaluator_reason = equivalence_info.reason
        equivalence_method = (
            equivalence_info.method.value if equivalence_info.method is not None else None
        )
        evaluator_inconclusive = equivalence_info.decision == EquivalenceDecision.INCONCLUSIVE
        bounded_search = equivalence_info.bounded_search
    except CaseTimeoutError as exc:
        equivalent = None
        status = "error"
        error = str(exc)
        evaluator_failure_kind = "evaluator_timeout"
        evaluator_decision = None
        evaluator_relation = None
        evaluator_reason = None
        equivalence_method = None
        evaluator_inconclusive = None
        bounded_search = None
    except Exception as exc:
        equivalent = None
        status = "error"
        error = f"equivalence check failed: {type(exc).__name__}: {exc}"
        evaluator_failure_kind = "evaluator_exception"
        evaluator_decision = None
        evaluator_relation = None
        evaluator_reason = None
        equivalence_method = None
        evaluator_inconclusive = None
        bounded_search = None
    atom_order_preserved, atom_mapping, atom_identity_guard_reason = (
        _atom_identity_mapping_with_reason(case.reference_mol, predicted_mol)
    )
    reference_radicals = _radical_electrons(case.reference_mol)
    predicted_radicals = _radical_electrons(predicted_mol)
    mapped_predicted_radical_sites = sorted(
        atom_mapping[index] for index in _radical_sites(predicted_mol) if index in atom_mapping
    )
    formal_radical_atom_index_match = (
        mapped_predicted_radical_sites == _radical_sites(case.reference_mol)
        if atom_order_preserved
        else None
    )
    return BDEResult(
        case_id=case.case_id,
        source_record_index=case.source_record_index,
        stratum=case.stratum,
        reference_smiles=case.reference_smiles,
        predicted_smiles=output.predicted_smiles,
        xyz=case.xyz,
        total_charge=case.total_charge,
        spin_multiplicity=case.spin_multiplicity,
        radical_site=case.radical_site,
        reconstruction_success=True,
        status=status,
        failure_kind=evaluator_failure_kind,
        error=error,
        equivalent=equivalent,
        evaluator_decision=evaluator_decision,
        evaluator_relation=evaluator_relation,
        evaluator_reason=evaluator_reason,
        equivalence_method=equivalence_method,
        evaluator_inconclusive=evaluator_inconclusive,
        bounded_search_attempted=(bounded_search.attempted if bounded_search else None),
        bounded_search_limit=(bounded_search.limit if bounded_search else None),
        bounded_search_limit_reached=(bounded_search.limit_reached if bounded_search else None),
        bounded_search_exhaustive=(bounded_search.exhaustive if bounded_search else None),
        bounded_search_candidate_count=(bounded_search.mol1_count if bounded_search else None),
        bounded_search_reference_count=(bounded_search.mol2_count if bounded_search else None),
        exact_smiles_match=case.reference_smiles == output.predicted_smiles,
        charge_consistent=_total_charge(predicted_mol) == case.total_charge,
        radical_electron_consistent=predicted_radicals == reference_radicals,
        atom_order_preserved=atom_order_preserved,
        atom_identity_guard_reason=atom_identity_guard_reason,
        formal_radical_atom_index_match=formal_radical_atom_index_match,
        reference_bonds=_bonds(case.reference_mol),
        predicted_bonds=_bonds(predicted_mol),
        runtime_ms=runtime_ms,
        source_metadata=json.dumps(case.source_metadata, ensure_ascii=True, sort_keys=True),
    )


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _percentile(values: list[float], percentile: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    index = (len(ordered) - 1) * percentile
    low = int(index)
    high = min(low + 1, len(ordered) - 1)
    weight = index - low
    return ordered[low] * (1.0 - weight) + ordered[high] * weight


def _review_reason(result: BDEResult) -> tuple[int, str]:
    if not result.reconstruction_success:
        return 0, "reconstruction_failure"
    if result.evaluator_inconclusive:
        return 1, "inconclusive"
    if result.equivalent is False:
        return 2, "non_equivalent"
    if result.equivalence_method == "resonance":
        return 3, "resonance_equivalent"
    if result.formal_radical_atom_index_match is False:
        return 4, "formal_radical_atom_index_mismatch"
    if result.charge_consistent is False or result.radical_electron_consistent is False:
        return 5, "charge_or_radical_electron_mismatch"
    if result.exact_smiles_match is False:
        return 6, "exact_smiles_mismatch"
    return 7, "representative_success"


def _select_review_cases(
    results: list[BDEResult], *, limit: int = 100
) -> list[tuple[str, BDEResult]]:
    ranked = sorted(
        ((_review_reason(result), result) for result in results),
        key=lambda item: (item[0][0], item[1].source_record_index, item[1].case_id),
    )
    return [(reason, result) for (_, reason), result in ranked[:limit]]


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "UNKNOWN"


def _summary(
    results: list[BDEResult],
    diagnostics: LoadDiagnostics,
    input_path: Path,
    setup: dict[str, Any],
    review_cases: list[tuple[str, BDEResult]],
) -> str:
    reconstruction_successes = [result for result in results if result.reconstruction_success]
    reconstruction_failures = [result for result in results if not result.reconstruction_success]
    equivalent = sum(result.equivalent is True for result in results)
    decisions = Counter(
        result.evaluator_decision for result in results if result.evaluator_decision is not None
    )
    relations = Counter(
        result.evaluator_relation for result in results if result.evaluator_relation is not None
    )
    equivalence_methods = Counter(
        result.equivalence_method for result in results if result.equivalence_method is not None
    )
    runtimes = [result.runtime_ms for result in results]
    strata = {
        name: sum(result.stratum == name for result in results)
        for name in sorted(diagnostics.strata_seen)
    }
    review_lines = "\n".join(
        f"- `{row.case_id}` [{reason}] ({row.stratum}): reference `{row.reference_smiles}`; "
        f"prediction `{row.predicted_smiles}`; error `{row.error or ''}`"
        for reason, row in review_cases[:20]
    )
    is_official_filename = input_path.name == EXPECTED_FILENAME
    data_status = "ready for this run" if is_official_filename else "non-official test input"
    smoke_status = (
        "complete for this run" if is_official_filename else "not an official-data smoke test"
    )
    return f"""# BDE-db {len(results)}-case smoke-test summary

This report is a pilot diagnostic and is not a final accuracy claim.

- Input: `{input_path}`
- Expected official filename: `{EXPECTED_FILENAME}`
- Implementation git SHA: `{setup["git_sha"]}`
- Method initialization: {setup["method_initialization_ms"]:.3f} ms
- Warm-up case: `{setup["warmup_case_id"]}`; status `{setup["warmup_status"]}`; {setup["warmup_ms"]:.3f} ms

## Counts

- SDF records scanned: {diagnostics.scanned_records}
- Eligible records: {diagnostics.eligible_records}
- Selected records: {diagnostics.selected_records}
- Reconstruction success: {len(reconstruction_successes)}
- Reconstruction failures: {len(reconstruction_failures)}
- Resonance-aware graph equivalence: {equivalent}
- Evaluator decisions: `{json.dumps(dict(sorted(decisions.items())), sort_keys=True)}`
- Evaluator relations: `{json.dumps(dict(sorted(relations.items())), sort_keys=True)}`
- Equivalence methods: `{json.dumps(dict(sorted(equivalence_methods.items())), sort_keys=True)}`
- Formal radical-electron agreement: {sum(result.radical_electron_consistent is True for result in results)}
- Exact formal-radical atom-index agreement: {sum(result.formal_radical_atom_index_match is True for result in results)}
- Atom-order preservation guard passed: {sum(result.atom_order_preserved is True for result in results)}
- Charge agreement: {sum(result.charge_consistent is True for result in results)}
- Diagnostic exact SMILES: {sum(result.exact_smiles_match is True for result in results)}
- Runtime p50/p95/max: {_percentile(runtimes, 0.50):.3f} / {_percentile(runtimes, 0.95):.3f} / {max(runtimes, default=0.0):.3f} ms
- Timeouts: {sum(result.failure_kind == "timeout" for result in results)}
- Ordinary exceptions: {sum(result.failure_kind == "exception" for result in results)}
- Method-returned errors: {sum(result.failure_kind == "method_error" for result in results)}
- Evaluator timeouts: {sum(result.failure_kind == "evaluator_timeout" for result in results)}
- Evaluator exceptions: {sum(result.failure_kind == "evaluator_exception" for result in results)}
- Selected strata: `{json.dumps(strata, sort_keys=True)}`
- Loader failures encountered: {len(diagnostics.failures)}

## Readiness questions

1. Open-shell suitability: yes for the released singlet closed-shell and doublet monoradical scope; it does not test higher-spin or multiradical reconstruction.
2. Charge source: SDF atom formal charges. The release contains one charged record, the doublet oxygen radical `[O-]`; all other records are neutral.
3. Multiplicity source: all released records encode exactly zero or one formal radical electron, supporting singlet or doublet respectively. There is no explicit multiplicity field and no evidence for higher multiplicities in this release.
4. Reference graph independence: yes; official code constructs the graph from source SMILES with RDKit and replaces only conformer coordinates with Gaussian-optimized coordinates.
5. Potential unfair structures: 6,208 SDF/SMILES stereo disagreements make chirality-aware scoring unfair; explicit-H radical representations also require graph normalization. Main equivalence therefore ignores chirality while exact isomeric SMILES remains diagnostic.
6. Before a full run: review the prioritized cases, then redesign execution and output for chunked streaming so cases, XYZ strings, bond dumps, metadata, and results are not all retained in memory.

## Manual-review cases

{review_lines or "- No cases available."}

## Component status

- PROVENANCE_READY: complete
- ADAPTER_READY: complete
- DATA_LOCAL_READY: {data_status}
- SMOKE_TEST_READY: {smoke_status}
"""


def run(
    input_path: Path,
    out_dir: Path,
    *,
    limit: int = 100,
    start: int | None = None,
    end: int | None = None,
    seed: int = 0,
    timeout_seconds: float | None = 5.0,
    review_limit: int = 100,
) -> list[BDEResult]:
    cases, diagnostics = load_bde_cases(
        input_path,
        limit=limit,
        start=start,
        end=end,
        seed=seed,
    )
    return run_cases(
        cases,
        diagnostics,
        input_path,
        out_dir,
        timeout_seconds=timeout_seconds,
        review_limit=review_limit,
    )


def run_cases(
    cases: list[BDECase],
    diagnostics: LoadDiagnostics,
    input_path: Path,
    out_dir: Path,
    *,
    timeout_seconds: float | None = 5.0,
    review_limit: int = 100,
) -> list[BDEResult]:
    initialization_started = time.perf_counter()
    method = MolGRCppMethod()
    method_initialization_ms = (time.perf_counter() - initialization_started) * 1000.0
    warmup_case = next(
        (case for case in cases if case.stratum == "closed_shell"),
        cases[0] if cases else None,
    )
    warmup_started = time.perf_counter()
    warmup_status = "skipped"
    warmup_error = None
    if warmup_case is not None:
        try:
            with case_timeout(timeout_seconds, f"warm-up case {warmup_case.case_id}"):
                warmup_output = method.run(warmup_case.to_method_case())
            warmup_status = warmup_output.status
            warmup_error = warmup_output.error
        except Exception as exc:
            warmup_status = "error"
            warmup_error = f"{type(exc).__name__}: {exc}"
    warmup_ms = (time.perf_counter() - warmup_started) * 1000.0
    results = [_run_case(case, method, timeout_seconds) for case in cases]
    out_dir.mkdir(parents=True, exist_ok=True)
    result_rows = [asdict(result) for result in results]
    result_fields = list(BDEResult.__dataclass_fields__)
    _write_csv(out_dir / "results.csv", result_rows, result_fields)
    failure_rows = [row for row in result_rows if row["status"] != "ok"]
    failure_rows.extend(
        {
            **dict.fromkeys(result_fields, ""),
            "source_record_index": failure["record_index"],
            "status": "loader_error",
            "error": failure["error"],
        }
        for failure in diagnostics.failures
    )
    _write_csv(out_dir / "failures.csv", failure_rows, result_fields)
    review_cases = _select_review_cases(results, limit=review_limit)
    review_rows = [{"review_priority": reason, **asdict(result)} for reason, result in review_cases]
    _write_csv(
        out_dir / "review_cases.csv",
        review_rows,
        ["review_priority", *result_fields],
    )
    setup = {
        "git_sha": _git_sha(),
        "method_initialization_ms": method_initialization_ms,
        "warmup_case_id": warmup_case.case_id if warmup_case is not None else None,
        "warmup_status": warmup_status,
        "warmup_error": warmup_error,
        "warmup_ms": warmup_ms,
    }
    (out_dir / "run_metadata.json").write_text(
        json.dumps(setup, indent=2, sort_keys=True), encoding="utf-8"
    )
    (out_dir / "summary.md").write_text(
        _summary(results, diagnostics, input_path, setup, review_cases), encoding="utf-8"
    )
    return results


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the BDE-db MolGR benchmark pilot.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--case-timeout-seconds", type=float, default=5.0)
    parser.add_argument("--review-limit", type=int, default=100)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    run(
        args.input,
        args.out,
        limit=args.limit,
        start=args.start,
        end=args.end,
        seed=args.seed,
        timeout_seconds=args.case_timeout_seconds or None,
        review_limit=args.review_limit,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
