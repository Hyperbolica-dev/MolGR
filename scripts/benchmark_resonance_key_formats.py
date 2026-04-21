from __future__ import annotations

import argparse
import csv
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Sequence


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from openbabel import pybel

from molgr.fallback.pipeline.reconstruct_with_metals import prepare_metal_state
from molgr.fallback.pipeline.reconstruct_without_metals import (
    _DEFAULT_RESONANCE_TRAVERSAL_POLICY,
    _run_linear_pipeline,
    _seed_state,
)
from molgr.fallback.pipeline.resonance import (
    build_resonance_state_key,
    get_radical_resonances,
    process_resonance,
)
from scripts.molgr_cases_molfile import load_molfile_cases
from scripts.molgr_cases_smiles_csv import load_smiles_csv_cases


_DEFAULT_FORMATS = (
    "mol",
    "sdf",
    "mdl",
    "mol2",
    "can",
    "smi",
    "inchi",
    "inchikey",
    "cml",
    "cdjson",
    "mna",
)


@dataclass(frozen=True)
class CaseSpec:
    case_id: str
    xyz_block: str
    total_charge: int
    total_radical_electrons: int


@dataclass(frozen=True)
class ProcessedStateRecord:
    case_id: str
    resonance_index: int
    reference_key: object
    omol: pybel.Molecule


@dataclass(frozen=True)
class SchemeResult:
    scheme: str
    total_states: int
    success_states: int
    error_states: int
    distinct_reference_keys: int
    distinct_scheme_keys: int
    split_groups: int
    split_states: int
    collision_groups: int
    collision_states: int
    avg_us_per_call: Optional[float]
    avg_output_len: Optional[float]
    max_output_len: Optional[int]
    exact_partition_match: bool


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark candidate resonance processed-state key builders. "
            "Each scheme is checked against the current build_resonance_state_key "
            "for partition equivalence, then timed independently."
        )
    )
    parser.add_argument(
        "--smiles-input",
        type=Path,
        default=Path("tests/test_cases.csv"),
        help="SMILES CSV used to generate benchmark XYZ cases.",
    )
    parser.add_argument(
        "--smiles-limit",
        type=int,
        default=None,
        help="Optional cap for SMILES benchmark cases.",
    )
    parser.add_argument(
        "--molfile-input",
        action="append",
        type=Path,
        default=[Path("tests/data/sdf/MoNNMo.sdf")],
        help=(
            "Optional molfile/SDF inputs. Metals are stripped via prepare_metal_state "
            "before running the no-metal resonance pipeline. Can be passed multiple times."
        ),
    )
    parser.add_argument(
        "--molfile-limit",
        type=int,
        default=None,
        help="Optional cap for molfile cases loaded from each molfile input path.",
    )
    parser.add_argument(
        "--policy",
        choices=("default", "none"),
        default="default",
        help=(
            "Resonance traversal policy. 'default' matches current production pruning. "
            "'none' enumerates without the traversal policy for a broader state set."
        ),
    )
    parser.add_argument(
        "--formats",
        type=str,
        default=",".join(_DEFAULT_FORMATS),
        help="Comma-separated pybel output formats to benchmark.",
    )
    parser.add_argument(
        "--time-repeats",
        type=int,
        default=5,
        help="How many full passes over the collected state set to use for timing.",
    )
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=None,
        help="Optional CSV path for the summary table.",
    )
    parser.add_argument(
        "--details-out",
        type=Path,
        default=None,
        help="Optional CSV path for split/collision details.",
    )
    return parser.parse_args()


def _iter_smiles_case_specs(input_path: Path, limit: Optional[int]) -> list[CaseSpec]:
    specs: list[CaseSpec] = []
    for case in load_smiles_csv_cases(input_path=input_path, limit=limit):
        if case.get("provider_error"):
            continue
        xyz_block = case.get("xyz_block")
        total_charge = case.get("total_charge")
        total_radical_electrons = case.get("total_radical_electrons")
        case_idx = case.get("case_idx")
        if not isinstance(xyz_block, str) or not isinstance(total_charge, int):
            continue
        if not isinstance(total_radical_electrons, int) or not isinstance(case_idx, int):
            continue
        specs.append(
            CaseSpec(
                case_id=f"smiles:{case_idx}",
                xyz_block=xyz_block,
                total_charge=total_charge,
                total_radical_electrons=total_radical_electrons,
            )
        )
    return specs


def _iter_molfile_case_specs(input_path: Path, limit: Optional[int]) -> list[CaseSpec]:
    specs: list[CaseSpec] = []
    for case in load_molfile_cases(input_path=input_path, limit=limit):
        if case.get("provider_error"):
            continue
        xyz_block = case.get("xyz_block")
        total_charge = case.get("total_charge")
        total_radical_electrons = case.get("total_radical_electrons")
        case_idx = case.get("case_idx")
        ground_truth_path = case.get("ground_truth_path")
        if not isinstance(xyz_block, str) or not isinstance(total_charge, int):
            continue
        if not isinstance(total_radical_electrons, int) or not isinstance(case_idx, int):
            continue
        prepared = prepare_metal_state(
            xyz_block,
            total_charge,
            total_radical_electrons,
        )
        label = input_path.name
        if isinstance(ground_truth_path, str) and ground_truth_path:
            label = Path(ground_truth_path).name
        specs.append(
            CaseSpec(
                case_id=f"molfile:{case_idx}:{label}",
                xyz_block=prepared.no_metal_xyz_block,
                total_charge=prepared.total_charge,
                total_radical_electrons=prepared.total_radical_electrons,
            )
        )
    return specs


def _collect_case_specs(args: argparse.Namespace) -> list[CaseSpec]:
    specs = _iter_smiles_case_specs(args.smiles_input, args.smiles_limit)
    for molfile_input in args.molfile_input:
        specs.extend(_iter_molfile_case_specs(molfile_input, args.molfile_limit))
    return specs


def _collect_processed_states(
    case_specs: Sequence[CaseSpec],
    *,
    policy_name: str,
) -> tuple[list[ProcessedStateRecord], list[dict[str, object]]]:
    records: list[ProcessedStateRecord] = []
    case_rows: list[dict[str, object]] = []
    traversal_policy = _DEFAULT_RESONANCE_TRAVERSAL_POLICY if policy_name == "default" else None

    for case_spec in case_specs:
        state = _seed_state(
            case_spec.xyz_block,
            case_spec.total_charge,
            case_spec.total_radical_electrons,
        )
        state = _run_linear_pipeline(state)
        resonances = get_radical_resonances(
            state.omol,
            traversal_policy=traversal_policy,
        )

        reference_keys = set()
        for resonance_index, resonance in enumerate(resonances):
            processed_omol, _processed_charge, _hit = process_resonance(
                resonance,
                state.given_charge,
            )
            reference_key = build_resonance_state_key(processed_omol)
            reference_keys.add(reference_key)
            records.append(
                ProcessedStateRecord(
                    case_id=case_spec.case_id,
                    resonance_index=resonance_index,
                    reference_key=reference_key,
                    omol=processed_omol,
                )
            )

        case_rows.append(
            {
                "case_id": case_spec.case_id,
                "resonance_count": len(resonances),
                "distinct_reference_keys": len(reference_keys),
            }
        )

    return records, case_rows


def _format_builder(format_name: str) -> Callable[[pybel.Molecule], object]:
    return lambda omol: omol.write(format_name)


def _evaluate_scheme(
    scheme: str,
    builder: Callable[[pybel.Molecule], object],
    records: Sequence[ProcessedStateRecord],
    *,
    time_repeats: int,
) -> tuple[SchemeResult, list[dict[str, object]]]:
    outputs: list[tuple[ProcessedStateRecord, object]] = []
    details: list[dict[str, object]] = []
    output_lengths: list[int] = []

    for record in records:
        try:
            key = builder(record.omol)
        except Exception as exc:
            details.append(
                {
                    "scheme": scheme,
                    "issue_type": "error",
                    "case_id": record.case_id,
                    "resonance_index": record.resonance_index,
                    "detail": f"{type(exc).__name__}: {exc}",
                }
            )
            continue
        outputs.append((record, key))
        if isinstance(key, str):
            output_lengths.append(len(key))

    ref_to_scheme: dict[object, set[object]] = defaultdict(set)
    scheme_to_ref: dict[object, set[object]] = defaultdict(set)
    ref_state_counts: Counter[object] = Counter()
    scheme_state_counts: Counter[object] = Counter()
    for record, scheme_key in outputs:
        ref_to_scheme[record.reference_key].add(scheme_key)
        scheme_to_ref[scheme_key].add(record.reference_key)
        ref_state_counts[record.reference_key] += 1
        scheme_state_counts[scheme_key] += 1

    split_groups = sum(1 for scheme_keys in ref_to_scheme.values() if len(scheme_keys) > 1)
    split_states = sum(
        ref_state_counts[reference_key]
        for reference_key, scheme_keys in ref_to_scheme.items()
        if len(scheme_keys) > 1
    )
    collision_groups = sum(
        1 for reference_keys in scheme_to_ref.values() if len(reference_keys) > 1
    )
    collision_states = sum(
        scheme_state_counts[scheme_key]
        for scheme_key, reference_keys in scheme_to_ref.items()
        if len(reference_keys) > 1
    )

    for reference_key, scheme_keys in ref_to_scheme.items():
        if len(scheme_keys) <= 1:
            continue
        details.append(
            {
                "scheme": scheme,
                "issue_type": "split",
                "case_id": "",
                "resonance_index": "",
                "detail": f"reference_key_split_count={len(scheme_keys)} states={ref_state_counts[reference_key]}",
            }
        )
    for scheme_key, reference_keys in scheme_to_ref.items():
        if len(reference_keys) <= 1:
            continue
        detail_text = (
            f"collision_ref_count={len(reference_keys)} states={scheme_state_counts[scheme_key]}"
        )
        if isinstance(scheme_key, str):
            detail_text = f"{detail_text} serialized_prefix={scheme_key[:120]!r}"
        details.append(
            {
                "scheme": scheme,
                "issue_type": "collision",
                "case_id": "",
                "resonance_index": "",
                "detail": detail_text,
            }
        )

    avg_us_per_call: Optional[float] = None
    if len(outputs) == len(records) and records:
        for record in records:
            builder(record.omol)
        started = time.perf_counter()
        for _ in range(max(time_repeats, 0)):
            for record in records:
                builder(record.omol)
        elapsed = time.perf_counter() - started
        total_calls = len(records) * max(time_repeats, 0)
        avg_us_per_call = (elapsed / total_calls) * 1_000_000.0 if total_calls else 0.0

    avg_output_len = (
        (sum(output_lengths) / len(output_lengths))
        if output_lengths and len(output_lengths) == len(outputs)
        else None
    )
    max_output_len = max(output_lengths) if output_lengths else None
    result = SchemeResult(
        scheme=scheme,
        total_states=len(records),
        success_states=len(outputs),
        error_states=len(records) - len(outputs),
        distinct_reference_keys=len({record.reference_key for record in records}),
        distinct_scheme_keys=len(scheme_to_ref),
        split_groups=split_groups,
        split_states=split_states,
        collision_groups=collision_groups,
        collision_states=collision_states,
        avg_us_per_call=avg_us_per_call,
        avg_output_len=avg_output_len,
        max_output_len=max_output_len,
        exact_partition_match=(
            len(outputs) == len(records) and split_groups == 0 and collision_groups == 0
        ),
    )
    return result, details


def _write_summary_csv(path: Path, results: Sequence[SchemeResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "scheme",
                "total_states",
                "success_states",
                "error_states",
                "distinct_reference_keys",
                "distinct_scheme_keys",
                "split_groups",
                "split_states",
                "collision_groups",
                "collision_states",
                "avg_us_per_call",
                "avg_output_len",
                "max_output_len",
                "exact_partition_match",
            ],
        )
        writer.writeheader()
        for result in results:
            writer.writerow(
                {
                    "scheme": result.scheme,
                    "total_states": result.total_states,
                    "success_states": result.success_states,
                    "error_states": result.error_states,
                    "distinct_reference_keys": result.distinct_reference_keys,
                    "distinct_scheme_keys": result.distinct_scheme_keys,
                    "split_groups": result.split_groups,
                    "split_states": result.split_states,
                    "collision_groups": result.collision_groups,
                    "collision_states": result.collision_states,
                    "avg_us_per_call": result.avg_us_per_call,
                    "avg_output_len": result.avg_output_len,
                    "max_output_len": result.max_output_len,
                    "exact_partition_match": result.exact_partition_match,
                }
            )


def _write_details_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["scheme", "issue_type", "case_id", "resonance_index", "detail"],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    args = _parse_args()
    case_specs = _collect_case_specs(args)
    if not case_specs:
        raise ValueError("No usable benchmark cases were loaded")

    records, case_rows = _collect_processed_states(case_specs, policy_name=args.policy)
    if not records:
        raise ValueError("No processed resonance states were collected")

    format_names = [name.strip() for name in args.formats.split(",") if name.strip()]
    results: list[SchemeResult] = []
    detail_rows: list[dict[str, object]] = []

    current_result, current_details = _evaluate_scheme(
        "current_tuple",
        build_resonance_state_key,
        records,
        time_repeats=args.time_repeats,
    )
    results.append(current_result)
    detail_rows.extend(current_details)

    for format_name in format_names:
        result, details = _evaluate_scheme(
            f"write:{format_name}",
            _format_builder(format_name),
            records,
            time_repeats=args.time_repeats,
        )
        results.append(result)
        detail_rows.extend(details)

    results.sort(
        key=lambda result: (
            not result.exact_partition_match,
            result.avg_us_per_call if result.avg_us_per_call is not None else float("inf"),
            result.scheme,
        )
    )

    print(
        f"cases={len(case_specs)} collected_states={len(records)} "
        f"distinct_reference_keys={len({record.reference_key for record in records})} "
        f"policy={args.policy}"
    )
    print("case_id,resonance_count,distinct_reference_keys")
    for row in case_rows:
        print(f"{row['case_id']},{row['resonance_count']},{row['distinct_reference_keys']}")
    print()
    print(
        "scheme,total_states,success_states,error_states,distinct_reference_keys,"
        "distinct_scheme_keys,split_groups,split_states,collision_groups,collision_states,"
        "avg_us_per_call,avg_output_len,max_output_len,exact_partition_match"
    )
    for result in results:
        avg_us_text = "" if result.avg_us_per_call is None else f"{result.avg_us_per_call:.3f}"
        avg_len_text = "" if result.avg_output_len is None else f"{result.avg_output_len:.1f}"
        max_len_text = "" if result.max_output_len is None else str(result.max_output_len)
        print(
            f"{result.scheme},{result.total_states},{result.success_states},{result.error_states},"
            f"{result.distinct_reference_keys},{result.distinct_scheme_keys},{result.split_groups},"
            f"{result.split_states},{result.collision_groups},{result.collision_states},"
            f"{avg_us_text},{avg_len_text},{max_len_text},{result.exact_partition_match}"
        )

    if args.summary_out is not None:
        _write_summary_csv(args.summary_out, results)
    if args.details_out is not None:
        _write_details_csv(args.details_out, detail_rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
