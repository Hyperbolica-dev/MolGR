from __future__ import annotations

import argparse
import csv
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Sequence


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from openbabel import pybel

from molgr.fallback import xyz2omol
from molgr.fallback.pipeline.reconstruct_with_metals import prepare_metal_state
from molgr.fallback.pipeline.reconstruct_without_metals import (
    _DEFAULT_RESONANCE_TRAVERSAL_POLICY,
    _recover_resonance_candidates,
    _run_linear_pipeline,
    _seed_state,
    xyz_to_omol_no_metal_state,
)
from molgr.fallback.utils.scoring import (
    _compute_organic_core_score,
    _compute_post_reinsertion_score,
    build_combined_metal_state_key,
    build_post_reinsertion_base_components,
    uncached_omol_score,
)
from scripts.molgr_cases_molfile import load_molfile_cases
from scripts.molgr_cases_smiles_csv import load_smiles_csv_cases


@dataclass(frozen=True)
class CaseSpec:
    case_id: str
    xyz_block: str
    total_charge: int
    total_radical_electrons: int


@dataclass(frozen=True)
class EvaluationResult:
    dataset: str
    scheme: str
    total_samples: int
    distinct_signatures: int
    distinct_keys: int
    collision_groups: int
    collision_samples: int
    avg_us_per_call: float
    exact_partition_match: bool


def _serializer_name(serializer: Callable[[pybel.Molecule], str]) -> str:
    prefix = "_serialize_"
    name = serializer.__name__
    if name.startswith(prefix):
        return name[len(prefix) :]
    return name


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate whether alternate score-key serializers preserve the required cache "
            "partitions for fallback scoring."
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
            "Optional molfile/SDF inputs. For no-metal scoring samples, metals are stripped via "
            "prepare_metal_state before running the no-metal pipeline. Can be passed multiple times."
        ),
    )
    parser.add_argument(
        "--molfile-limit",
        type=int,
        default=None,
        help="Optional cap for molfile cases loaded from each molfile input path.",
    )
    parser.add_argument(
        "--time-repeats",
        type=int,
        default=10,
        help="How many full passes over each sample set to use for timing.",
    )
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=None,
        help="Optional CSV path for the summary table.",
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
        label = input_path.name
        if isinstance(ground_truth_path, str) and ground_truth_path:
            label = Path(ground_truth_path).name
        specs.append(
            CaseSpec(
                case_id=f"molfile:{case_idx}:{label}",
                xyz_block=xyz_block,
                total_charge=total_charge,
                total_radical_electrons=total_radical_electrons,
            )
        )
    return specs


def _collect_case_specs(args: argparse.Namespace) -> list[CaseSpec]:
    specs = _iter_smiles_case_specs(args.smiles_input, args.smiles_limit)
    for molfile_input in args.molfile_input:
        specs.extend(_iter_molfile_case_specs(molfile_input, args.molfile_limit))
    return specs


def _collect_no_metal_samples(case_specs: Sequence[CaseSpec]) -> list[pybel.Molecule]:
    samples: list[pybel.Molecule] = []
    for case_spec in case_specs:
        prepared = prepare_metal_state(
            case_spec.xyz_block,
            case_spec.total_charge,
            case_spec.total_radical_electrons,
        )
        state = xyz_to_omol_no_metal_state(
            prepared.no_metal_xyz_block,
            prepared.total_charge,
            prepared.total_radical_electrons,
        )
        if state is not None:
            samples.append(state.omol)

        seed_state = _seed_state(
            prepared.no_metal_xyz_block,
            prepared.total_charge,
            prepared.total_radical_electrons,
        )
        linear_state = _run_linear_pipeline(seed_state)
        for candidate in _recover_resonance_candidates(
            linear_state,
            resonance_traversal_policy=_DEFAULT_RESONANCE_TRAVERSAL_POLICY,
        ):
            samples.append(candidate.omol)
    return samples


def _collect_full_samples(case_specs: Sequence[CaseSpec]) -> list[pybel.Molecule]:
    samples: list[pybel.Molecule] = []
    for case_spec in case_specs:
        omol = xyz2omol(
            case_spec.xyz_block,
            total_charge=case_spec.total_charge,
            total_radical_electrons=case_spec.total_radical_electrons,
        )
        if omol is not None:
            samples.append(omol)
    return samples


def _serialize_mol(omol: pybel.Molecule) -> str:
    return omol.write("mol")


def _serialize_molreport(omol: pybel.Molecule) -> str:
    return omol.write("molreport")


def _organic_signature(omol: pybel.Molecule) -> object:
    return (
        _compute_organic_core_score(omol),
        *build_post_reinsertion_base_components(omol),
    )


def _full_score_signature(omol: pybel.Molecule) -> object:
    return uncached_omol_score(omol)


def _post_pair_signature(omol: pybel.Molecule) -> object:
    return _compute_post_reinsertion_score(omol)


def _evaluate_dataset(
    dataset: str,
    samples: Sequence[pybel.Molecule],
    serializer: Callable[[pybel.Molecule], str],
    signature_builder: Callable[[pybel.Molecule], object],
    *,
    pair_with_metal_key: bool,
    time_repeats: int,
) -> EvaluationResult:
    entries: list[tuple[object, object]] = []
    for omol in samples:
        serialized = serializer(omol)
        if pair_with_metal_key:
            key = (serialized, build_combined_metal_state_key(omol))
        else:
            key = serialized
        entries.append((key, signature_builder(omol)))

    signatures_by_key: dict[object, set[object]] = {}
    for key, signature in entries:
        signatures_by_key.setdefault(key, set()).add(signature)

    collision_groups = 0
    collision_samples = 0
    for key, signatures in signatures_by_key.items():
        if len(signatures) <= 1:
            continue
        collision_groups += 1
        collision_samples += sum(1 for entry_key, _signature in entries if entry_key == key)

    started = time.perf_counter()
    for _ in range(max(time_repeats, 0)):
        for omol in samples:
            serializer(omol)
    elapsed = time.perf_counter() - started
    total_calls = len(samples) * max(time_repeats, 0)
    avg_us_per_call = (elapsed / total_calls) * 1_000_000.0 if total_calls else 0.0

    return EvaluationResult(
        dataset=dataset,
        scheme=_serializer_name(serializer),
        total_samples=len(samples),
        distinct_signatures=len({signature for _key, signature in entries}),
        distinct_keys=len(signatures_by_key),
        collision_groups=collision_groups,
        collision_samples=collision_samples,
        avg_us_per_call=avg_us_per_call,
        exact_partition_match=(collision_groups == 0),
    )


def _score_key_quantization_check(serializer: Callable[[pybel.Molecule], str]) -> bool:
    first = pybel.readstring(
        "xyz",
        """2
CO
C 2.0000004 0.0 0.0
O 3.2000004 0.0 0.0
""",
    )
    second = pybel.readstring(
        "xyz",
        """2
CO
C 2.00000049 0.0 0.0
O 3.20000049 0.0 0.0
""",
    )
    return serializer(first) == serializer(second)


def _write_summary_csv(path: Path, results: Sequence[EvaluationResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "dataset",
                "scheme",
                "total_samples",
                "distinct_signatures",
                "distinct_keys",
                "collision_groups",
                "collision_samples",
                "avg_us_per_call",
                "exact_partition_match",
            ],
        )
        writer.writeheader()
        for result in results:
            writer.writerow(
                {
                    "dataset": result.dataset,
                    "scheme": result.scheme,
                    "total_samples": result.total_samples,
                    "distinct_signatures": result.distinct_signatures,
                    "distinct_keys": result.distinct_keys,
                    "collision_groups": result.collision_groups,
                    "collision_samples": result.collision_samples,
                    "avg_us_per_call": result.avg_us_per_call,
                    "exact_partition_match": result.exact_partition_match,
                }
            )


def main() -> int:
    args = _parse_args()
    case_specs = _collect_case_specs(args)
    if not case_specs:
        raise ValueError("No usable benchmark cases were loaded")

    no_metal_samples = _collect_no_metal_samples(case_specs)
    full_samples = _collect_full_samples(case_specs)

    serializers = (_serialize_mol, _serialize_molreport)
    results: list[EvaluationResult] = []
    for serializer in serializers:
        results.append(
            _evaluate_dataset(
                "organic_base_key",
                no_metal_samples,
                serializer,
                _organic_signature,
                pair_with_metal_key=False,
                time_repeats=args.time_repeats,
            )
        )
        results.append(
            _evaluate_dataset(
                "full_score_key",
                full_samples,
                serializer,
                _full_score_signature,
                pair_with_metal_key=False,
                time_repeats=args.time_repeats,
            )
        )
        results.append(
            _evaluate_dataset(
                "direct_post_key_pair",
                full_samples,
                serializer,
                _post_pair_signature,
                pair_with_metal_key=True,
                time_repeats=args.time_repeats,
            )
        )

    results.sort(
        key=lambda result: (
            result.dataset,
            not result.exact_partition_match,
            result.avg_us_per_call,
            result.scheme,
        )
    )

    print(
        f"case_count={len(case_specs)} no_metal_samples={len(no_metal_samples)} "
        f"full_samples={len(full_samples)}"
    )
    print("scheme,quantized_six_decimal_pair")
    for serializer in serializers:
        print(f"{_serializer_name(serializer)},{_score_key_quantization_check(serializer)}")
    print()
    print(
        "dataset,scheme,total_samples,distinct_signatures,distinct_keys,collision_groups,"
        "collision_samples,avg_us_per_call,exact_partition_match"
    )
    for result in results:
        print(
            f"{result.dataset},{result.scheme},{result.total_samples},{result.distinct_signatures},"
            f"{result.distinct_keys},{result.collision_groups},{result.collision_samples},"
            f"{result.avg_us_per_call:.3f},{result.exact_partition_match}"
        )

    if args.summary_out is not None:
        _write_summary_csv(args.summary_out, results)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
