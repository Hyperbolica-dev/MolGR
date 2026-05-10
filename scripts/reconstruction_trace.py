#!/usr/bin/env python3
"""Trace MolGR reconstruction for arbitrary XYZ/electronic-state inputs.

This script intentionally uses the Python fallback internals. The public C++
backend does not expose every intermediate metal candidate, while the fallback
path is the semantic reference used for parity.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
import sys
import time
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, cast


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import RDLogger

from molgr.config import MolGRConfig
from molgr.fallback.pipeline import reconstruct_without_metals
from molgr.fallback.stages.break_bond import break_deformed_ene, break_one_bond
from molgr.fallback.stages.clean import (
    clean_carbene_neighbor_unsaturated,
    clean_neighbor_radicals,
    clean_resonances,
)
from molgr.fallback.stages.eliminate import (
    eliminate_carbene_neighbor_heteroatom,
    eliminate_carboxyl,
    eliminate_charge_spliting,
    eliminate_CN_in_doubt,
    eliminate_high_positive_charge_atoms,
    eliminate_NNN,
)
from molgr.fallback.stages.fresh import fresh_omol_charge_radical
from molgr.fallback.stages.preprocess import make_connections, pre_clean, validate_omol
from molgr.fallback.state import MetalCandidateState, OmolStateMachine, ReconstructionState
from molgr.fallback.utils import resonance as resonance_utils
from molgr.fallback.utils.dataclasses import MetalAtomPosition
from molgr.fallback.utils.metals import preparation, scoring, search
from molgr.fallback.utils.no_metals import preparation as no_metal_preparation
from molgr.fallback.utils.no_metals import resonance as no_metal_resonance
from molgr.fallback.utils.no_metals import selection as no_metal_selection


RDLogger.DisableLog("rdApp.*")  # type: ignore[arg-type]


@dataclasses.dataclass(frozen=True)
class TraceInputCase:
    """One coordinate/electronic-state input to trace."""

    id: str
    xyz_block: str
    total_charge: int
    total_radical_electrons: int
    xyz_path: Path | None = None
    xyz_source: str = ""


_SCORE_DETAIL_PREFIXES = (
    "analysis_",
    "force_field_",
    "metal_",
    "organic_",
    "passes_",
    "selection_",
)
_SCORE_DETAIL_KEYS = {
    "combination_index",
    "score",
}
_IMPORTANT_SCORE_KEYS = (
    "metal_discordance_count",
    "metal_discordance_structural_count",
    "metal_discordance_aromatic_ring_deficit_count",
    "metal_discordance_max_aromatic_ring_count",
    "metal_discordance_inner_visible_diradical_count",
    "metal_discordance_outer_or_invisible_adjacent_double_charge_count",
    "metal_discordance_outer_or_invisible_adjacent_same_sign_double_charge_count",
    "metal_discordance_outer_or_invisible_adjacent_opposite_sign_double_charge_count",
    "metal_discordance_outer_or_invisible_adjacent_unknown_metal_sign_double_charge_count",
    "metal_discordance_inner_visible_adjacent_carbanion_pair_count",
    "metal_discordance_inner_visible_conjugated_carbanion_pair_count",
    "metal_discordance_inner_visible_same_sign_charge_count",
    "metal_discordance_negative_metal_count",
    "metal_discordance_zero_valent_metals_with_organic_cation_count",
    "metal_discordance_negative_metal_outer_sphere_cation_exception",
    "metal_discordance_negative_metal_positive_metal_counterion_exception",
    "passes_metal_discordance_filter",
    "score",
    "force_field_energy",
    "organic_aromatic_ring_loss",
    "organic_max_conjugated_component_loss",
    "organic_charge_localization_penalty",
    "organic_aromatic_atom_loss",
    "organic_conjugated_atom_loss",
    "organic_radical_localization_penalty",
    "selection_key",
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--xyz",
        dest="xyz_paths",
        action="append",
        type=Path,
        default=[],
        help="XYZ file to trace. May be repeated.",
    )
    parser.add_argument(
        "--xyz-block",
        dest="xyz_blocks",
        action="append",
        default=[],
        help="Inline XYZ block to trace. May be repeated.",
    )
    parser.add_argument(
        "--stdin",
        action="store_true",
        help="Read one XYZ block from stdin.",
    )
    parser.add_argument(
        "--id",
        dest="ids",
        action="append",
        default=[],
        help=(
            "Optional case id/label. May be repeated or comma-separated. Labels are assigned "
            "to --xyz, --xyz-block, and --stdin inputs in order."
        ),
    )
    parser.add_argument(
        "--case-json",
        type=Path,
        default=None,
        help=(
            "Optional JSON case file. The file may contain one object or a list of objects with "
            "id, xyz_path or xyz_block, total_charge/charge, and total_radical_electrons or "
            "spin_multiplicity."
        ),
    )
    parser.add_argument(
        "--total-charge",
        "--charge",
        dest="total_charge",
        type=int,
        default=0,
        help="Total molecular charge for direct XYZ inputs. Default: 0.",
    )
    radical_group = parser.add_mutually_exclusive_group()
    radical_group.add_argument(
        "--total-radical-electrons",
        type=int,
        default=None,
        help="Total radical-electron count for direct XYZ inputs. Default: 0.",
    )
    radical_group.add_argument(
        "--spin-multiplicity",
        type=int,
        default=None,
        help="Spin multiplicity for direct XYZ inputs. Converted to radical electrons as M - 1.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Optional output path. Defaults to stdout.",
    )
    parser.add_argument(
        "--format",
        choices=("auto", "markdown", "json"),
        default="auto",
        help=(
            "Output format. In auto mode, .json writes JSON; every other output path and stdout "
            "write a Markdown report."
        ),
    )
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="JSON indentation when --format json is used. Use 0 for compact JSON.",
    )
    parser.add_argument(
        "--no-score-all-candidates",
        action="store_true",
        help=(
            "Do not compute analysis metal-candidate metrics for candidates filtered out by "
            "discordance. Production metadata is still reported."
        ),
    )
    args = parser.parse_args()
    if args.total_radical_electrons is not None and args.total_radical_electrons < 0:
        parser.error("--total-radical-electrons must be >= 0")
    if args.spin_multiplicity is not None and args.spin_multiplicity < 1:
        parser.error("--spin-multiplicity must be >= 1")
    return args


def split_repeated_values(raw_values: Sequence[str]) -> list[str]:
    values: list[str] = []
    for raw_value in raw_values:
        values.extend(item.strip() for item in raw_value.split(",") if item.strip())
    return values


def _radicals_from_spin_multiplicity(spin_multiplicity: int) -> int:
    if spin_multiplicity < 1:
        raise ValueError("spin_multiplicity must be >= 1")
    return spin_multiplicity - 1


def _direct_total_radicals(args: argparse.Namespace) -> int:
    if args.total_radical_electrons is not None:
        return int(args.total_radical_electrons)
    if args.spin_multiplicity is not None:
        return _radicals_from_spin_multiplicity(int(args.spin_multiplicity))
    return 0


def _require_int_field(
    raw_case: dict[str, Any], keys: Sequence[str], *, default: int | None
) -> int:
    for key in keys:
        value = raw_case.get(key)
        if value is not None and value != "":
            return int(value)
    if default is not None:
        return default
    raise ValueError(f"JSON case is missing one of: {', '.join(keys)}")


def _json_case_to_input_case(
    raw_case: dict[str, Any],
    *,
    base_dir: Path,
    fallback_index: int,
) -> TraceInputCase:
    case_id = str(raw_case.get("id") or raw_case.get("label") or f"case_{fallback_index}")
    total_charge = _require_int_field(raw_case, ("total_charge", "charge"), default=0)
    if raw_case.get("total_radical_electrons") is not None or raw_case.get("radicals") is not None:
        total_radicals = _require_int_field(
            raw_case,
            ("total_radical_electrons", "radicals"),
            default=None,
        )
    elif raw_case.get("spin_multiplicity") is not None:
        total_radicals = _radicals_from_spin_multiplicity(int(raw_case["spin_multiplicity"]))
    else:
        total_radicals = 0

    xyz_block = raw_case.get("xyz_block") or raw_case.get("xyz")
    xyz_path: Path | None = None
    xyz_source = "json_xyz_block"
    if xyz_block is None:
        raw_path = raw_case.get("xyz_path") or raw_case.get("path")
        if raw_path is None:
            raise ValueError(f"JSON case {case_id!r} is missing xyz_block or xyz_path")
        xyz_path = Path(str(raw_path))
        if not xyz_path.is_absolute():
            xyz_path = base_dir / xyz_path
        xyz_block = xyz_path.read_text(encoding="utf-8")
        xyz_source = "json_xyz_path"

    return TraceInputCase(
        id=case_id,
        xyz_block=str(xyz_block),
        total_charge=total_charge,
        total_radical_electrons=total_radicals,
        xyz_path=xyz_path,
        xyz_source=xyz_source,
    )


def _load_json_cases(path: Path) -> list[TraceInputCase]:
    raw_payload = json.loads(path.read_text(encoding="utf-8"))
    raw_cases = raw_payload if isinstance(raw_payload, list) else [raw_payload]
    cases: list[TraceInputCase] = []
    for case_index, raw_case in enumerate(raw_cases, start=1):
        if not isinstance(raw_case, dict):
            raise ValueError(f"JSON case #{case_index} must be an object")
        cases.append(
            _json_case_to_input_case(
                raw_case,
                base_dir=path.parent,
                fallback_index=case_index,
            )
        )
    return cases


def _collect_input_cases(args: argparse.Namespace) -> list[TraceInputCase]:
    cases: list[TraceInputCase] = []
    if args.case_json is not None:
        cases.extend(_load_json_cases(args.case_json))

    direct_total_radicals = _direct_total_radicals(args)
    labels = split_repeated_values(args.ids)
    label_index = 0

    def next_label(default: str) -> str:
        nonlocal label_index
        if label_index >= len(labels):
            return default
        label = labels[label_index]
        label_index += 1
        return label

    direct_input_index = 1
    for xyz_path in args.xyz_paths:
        cases.append(
            TraceInputCase(
                id=next_label(xyz_path.stem),
                xyz_block=xyz_path.read_text(encoding="utf-8"),
                total_charge=int(args.total_charge),
                total_radical_electrons=direct_total_radicals,
                xyz_path=xyz_path,
                xyz_source="xyz_path",
            )
        )
        direct_input_index += 1

    for xyz_block in args.xyz_blocks:
        cases.append(
            TraceInputCase(
                id=next_label(f"inline_{direct_input_index}"),
                xyz_block=str(xyz_block),
                total_charge=int(args.total_charge),
                total_radical_electrons=direct_total_radicals,
                xyz_source="xyz_block",
            )
        )
        direct_input_index += 1

    if args.stdin:
        cases.append(
            TraceInputCase(
                id=next_label("stdin"),
                xyz_block=sys.stdin.read(),
                total_charge=int(args.total_charge),
                total_radical_electrons=direct_total_radicals,
                xyz_source="stdin",
            )
        )

    if label_index < len(labels):
        raise ValueError("more --id labels were provided than direct coordinate inputs")
    if not cases:
        raise ValueError("provide at least one of --xyz, --xyz-block, --stdin, or --case-json")
    return cases


def _first_smiles_token(raw_smiles: str) -> str:
    tokens = raw_smiles.strip().split()
    return tokens[0] if tokens else ""


def _write_omol_smiles(omol: pybel.Molecule, output_format: str) -> str:
    try:
        return _first_smiles_token(cast(str, omol.write(output_format)))
    except (OSError, ValueError, RuntimeError):
        return ""


def _omol_smiles_detail(omol: pybel.Molecule) -> dict[str, str]:
    return {
        "smiles": _write_omol_smiles(omol, "smi"),
        "canonical_smiles": _write_omol_smiles(omol, "can"),
    }


def _omol_charge_radical_summary(omol: pybel.Molecule) -> dict[str, Any]:
    obmol = cast(ob.OBMol, omol.OBMol)
    formal_charge = 0
    spin_multiplicity_sum = 0
    spin_multiplicity_singlet_sum = 0
    atom_charge_counts: dict[str, int] = {}
    atom_spin_counts: dict[str, int] = {}
    for atom_iter in ob.OBMolAtomIter(obmol):
        atom = cast(ob.OBAtom, atom_iter)
        charge = int(atom.GetFormalCharge())
        spin = int(atom.GetSpinMultiplicity())
        formal_charge += charge
        spin_multiplicity_sum += spin
        spin_multiplicity_singlet_sum += spin % 2
        if charge != 0:
            atom_charge_counts[f"{atom.GetAtomicNum()}:{charge:+d}"] = (
                atom_charge_counts.get(f"{atom.GetAtomicNum()}:{charge:+d}", 0) + 1
            )
        if spin != 0:
            atom_spin_counts[f"{atom.GetAtomicNum()}:r{spin}"] = (
                atom_spin_counts.get(f"{atom.GetAtomicNum()}:r{spin}", 0) + 1
            )

    return {
        "atom_count": int(obmol.NumAtoms()),
        "bond_count": int(obmol.NumBonds()),
        "formal_charge_sum": formal_charge,
        "spin_multiplicity_sum": spin_multiplicity_sum,
        "spin_multiplicity_singlet_sum": spin_multiplicity_singlet_sum,
        "charged_atom_counts": atom_charge_counts,
        "radical_atom_counts": atom_spin_counts,
    }


def _omol_state_snapshot(
    omol: pybel.Molecule,
    *,
    given_charge: int,
    target_charge: int,
    target_radical_electrons: int,
) -> dict[str, Any]:
    summary = _omol_charge_radical_summary(omol)
    summary.update(_omol_smiles_detail(omol))
    summary["given_charge"] = int(given_charge)
    summary["valid_for_target"] = validate_omol(
        omol,
        target_charge,
        target_radical_electrons,
    )
    return summary


def _metal_state_to_dict(metal_state: MetalAtomPosition) -> dict[str, Any]:
    return {
        "idx": int(metal_state.idx),
        "symbol": metal_state.symbol,
        "element_idx": int(metal_state.element_idx),
        "valence": int(metal_state.valence),
        "radical_num": int(metal_state.radical_num),
        "position": [
            float(metal_state.position_x),
            float(metal_state.position_y),
            float(metal_state.position_z),
        ],
    }


def _jsonable(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if math.isfinite(value):
            return value
        if math.isnan(value):
            return "NaN"
        return "Infinity" if value > 0.0 else "-Infinity"
    if isinstance(value, Path):
        return str(value)
    if dataclasses.is_dataclass(value):
        return _jsonable(dataclasses.asdict(value))
    if isinstance(value, dict):
        return {str(_jsonable(key)): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, set):
        return [_jsonable(item) for item in sorted(value, key=repr)]
    return repr(value)


def _format_float(value: float) -> str:
    if not math.isfinite(value):
        return _jsonable(value)
    if abs(value) >= 1000.0 or (0.0 < abs(value) < 0.001):
        return f"{value:.6e}"
    return f"{value:.6g}"


def _format_value(value: Any) -> str:
    value = _jsonable(value)
    if value is None:
        return ""
    if isinstance(value, bool):
        return "是" if value else "否"
    if isinstance(value, float):
        return _format_float(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, str):
        return value
    if isinstance(value, (list, tuple)):
        if all(not isinstance(item, (dict, list, tuple)) for item in value):
            return "[" + ", ".join(_format_value(item) for item in value) + "]"
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    if isinstance(value, dict):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return str(value)


def _markdown_cell(value: Any) -> str:
    text = _format_value(value)
    if text == "":
        return ""
    return text.replace("\\", "\\\\").replace("|", "\\|").replace("\n", "<br>")


def _markdown_table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> str:
    if not rows:
        return "_无_"
    lines = [
        "| " + " | ".join(_markdown_cell(header) for header in headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_markdown_cell(cell) for cell in row) + " |")
    return "\n".join(lines)


def _markdown_kv_table(items: Sequence[tuple[str, Any]]) -> str:
    return _markdown_table(("字段", "值"), items)


def _resolve_output_format(args: argparse.Namespace) -> str:
    if args.format != "auto":
        return cast(str, args.format)
    if args.out is not None and args.out.suffix.lower() == ".json":
        return "json"
    return "markdown"


def _score_detail_value(candidate: dict[str, Any], key: str, default: Any = "") -> Any:
    details = cast(Dict[str, Any], candidate.get("score_details", {}))
    metadata = cast(Dict[str, Any], candidate.get("metadata", {}))
    if key in details:
        return details[key]
    return metadata.get(key, default)


def _short_hash(value: Any) -> str:
    return hashlib.sha1(repr(value).encode("utf-8")).hexdigest()[:12]


def _metal_state_label(metal_state: dict[str, Any]) -> str:
    valence = int(metal_state.get("valence", 0))
    valence_label = f"{valence:+d}"
    return (
        f"#{metal_state.get('idx', '')} {metal_state.get('symbol', '')}"
        f"({valence_label}, r{metal_state.get('radical_num', 0)})"
    )


def _metal_states_label(metal_states: Sequence[dict[str, Any]]) -> str:
    return "<br>".join(_metal_state_label(metal_state) for metal_state in metal_states)


def _candidate_short_label(candidate: dict[str, Any]) -> str:
    layer_index = candidate.get("search_layer_index", "")
    combination_index = candidate.get("combination_index", "")
    selected = "选中, " if candidate.get("selected") else ""
    return f"{selected}L{layer_index}/C{combination_index}"


def _candidate_score_rows(candidate: dict[str, Any]) -> list[tuple[str, Any]]:
    details = cast(Dict[str, Any], candidate.get("score_details", {}))
    keys: list[str] = [key for key in _IMPORTANT_SCORE_KEYS if key in details]
    keys.extend(sorted(key for key in details if key not in set(keys)))
    return [(key, details[key]) for key in keys]


def _render_case_summary(case: dict[str, Any]) -> str:
    selected = case.get("selected_candidate")
    selected_label = ""
    if case.get("trace_kind") == "no_metal":
        no_metal_trace = cast(Dict[str, Any], case.get("no_metal_trace", {}))
        if no_metal_trace.get("status") == "direct_valid":
            selected_label = "direct"
        elif no_metal_trace.get("selected_candidate"):
            selected_label = "resonance"
    elif isinstance(selected, dict):
        selected_label = (
            f"L{selected.get('score_details', {}).get('search_layer_index', case.get('search', {}).get('selected_layer_index', ''))}"
            f"/C{selected.get('combination_index', '')}"
        )
    rows: list[tuple[str, Any]] = [
        ("id", case.get("id", "")),
        ("状态", case.get("status", "")),
        ("重建类型", case.get("trace_kind", "")),
    ]
    if case.get("row_index", "") != "":
        rows.append(("CSV 行号", case.get("row_index", "")))
    rows.extend(
        [
            ("总电荷", case.get("charge", "")),
            ("总自由基电子数", case.get("total_radical_electrons", "")),
            ("自旋多重度", case.get("spin_multiplicity", "")),
        ]
    )
    if case.get("spin_source", "") != "":
        rows.append(("自旋来源", case.get("spin_source", "")))
    rows.extend(
        [
            ("XYZ 路径", case.get("xyz_path", "")),
            ("XYZ 来源", case.get("xyz_source", "")),
        ]
    )
    if case.get("reference_smiles", "") != "":
        rows.append(("参考 SMILES", case.get("reference_smiles", "")))
    rows.extend(
        [
            (
                "金属原子数",
                cast(Dict[str, Any], case.get("base_state", {})).get("metal_atom_count", ""),
            ),
            (
                "生产选择层",
                cast(Dict[str, Any], case.get("search", {})).get("selected_layer_index", ""),
            ),
            ("生产候选数", case.get("production_candidate_count", "")),
            ("全部已评分候选数", case.get("candidate_count", "")),
            ("选中候选", selected_label),
            ("耗时秒", case.get("elapsed_seconds", "")),
        ]
    )
    return _markdown_kv_table(rows)


def _render_available_metal_states(case: dict[str, Any]) -> str:
    base_state = cast(Dict[str, Any], case.get("base_state", {}))
    state_groups = cast(List[Any], base_state.get("available_metal_states_by_site", []))
    rows: list[list[Any]] = []
    for site_index, state_options in enumerate(state_groups):
        if not isinstance(state_options, list):
            continue
        for option_index, metal_state in enumerate(state_options):
            if not isinstance(metal_state, dict):
                continue
            rows.append(
                [
                    site_index,
                    option_index,
                    metal_state.get("idx", ""),
                    metal_state.get("symbol", ""),
                    metal_state.get("valence", ""),
                    metal_state.get("radical_num", ""),
                    _format_value(metal_state.get("position", [])),
                ]
            )
    return _markdown_table(
        ("位点", "选项", "原子序号", "元素", "价态", "自由基数", "坐标"),
        rows,
    )


def _render_layer_table(case: dict[str, Any]) -> str:
    search_summary = cast(Dict[str, Any], case.get("search", {}))
    layer_summaries = cast(List[Any], search_summary.get("layer_summaries", []))
    rows: list[list[Any]] = []
    for layer in layer_summaries:
        if not isinstance(layer, dict):
            continue
        rows.append(
            [
                layer.get("layer_index", ""),
                layer.get("status", ""),
                layer.get("production_selected_layer", ""),
                layer.get("state_group_count", ""),
                layer.get("state_options_per_group", ""),
                layer.get("target_bucket_count", ""),
                layer.get("candidate_count", ""),
                layer.get("prepared_candidate_count", ""),
            ]
        )
    return _markdown_table(
        (
            "层",
            "状态",
            "生产选择层",
            "金属组数",
            "每组候选数",
            "目标桶数",
            "枚举候选数",
            "已评分候选数",
        ),
        rows,
    )


def _render_target_bucket_table(case: dict[str, Any]) -> str:
    search_summary = cast(Dict[str, Any], case.get("search", {}))
    layer_summaries = cast(List[Any], search_summary.get("layer_summaries", []))
    rows: list[list[Any]] = []
    for layer in layer_summaries:
        if not isinstance(layer, dict):
            continue
        for target_bucket in cast(List[Any], layer.get("target_buckets", [])):
            if not isinstance(target_bucket, dict):
                continue
            target = cast(Dict[str, Any], target_bucket.get("target", {}))
            organic_part = cast(Dict[str, Any], target_bucket.get("organic_part", {}))
            rows.append(
                [
                    layer.get("layer_index", ""),
                    target.get("no_metal_charge", ""),
                    target.get("no_metal_radical_electrons", ""),
                    target_bucket.get("status", ""),
                    target_bucket.get("candidate_count", ""),
                    target_bucket.get("prepared_candidate_count", ""),
                    organic_part.get("canonical_smiles", organic_part.get("smiles", "")),
                    target_bucket.get("no_metal_score", ""),
                ]
            )
    return _markdown_table(
        (
            "层",
            "有机目标电荷",
            "有机自由基电子",
            "状态",
            "候选数",
            "已评分数",
            "有机部分 canonical SMILES",
            "有机力场分",
        ),
        rows,
    )


def _state_summary_cells(state: dict[str, Any]) -> list[Any]:
    return [
        state.get("canonical_smiles", state.get("smiles", "")),
        state.get("formal_charge_sum", ""),
        state.get("spin_multiplicity_sum", ""),
        state.get("spin_multiplicity_singlet_sum", ""),
        state.get("given_charge", ""),
        state.get("valid_for_target", ""),
        state.get("charged_atom_counts", ""),
        state.get("radical_atom_counts", ""),
    ]


def _render_no_metal_linear_steps(no_metal_trace: dict[str, Any]) -> str:
    rows: list[list[Any]] = []
    for step in cast(List[Any], no_metal_trace.get("linear_steps", [])):
        if not isinstance(step, dict):
            continue
        state = cast(Dict[str, Any], step.get("state", {}))
        rows.append(
            [
                step.get("step_index", ""),
                step.get("phase", ""),
                step.get("kind", ""),
                step.get("hit", ""),
                step.get("omol_revision", ""),
                *_state_summary_cells(state),
            ]
        )
    return _markdown_table(
        (
            "步骤",
            "阶段",
            "类型",
            "命中",
            "修订",
            "canonical SMILES",
            "形式电荷",
            "自由基和",
            "自由基奇偶和",
            "剩余电荷预算",
            "匹配目标",
            "带电原子",
            "自由基原子",
        ),
        rows,
    )


def _render_no_metal_resonance_candidates(no_metal_trace: dict[str, Any]) -> str:
    resonance = cast(Dict[str, Any], no_metal_trace.get("resonance", {}))
    rows: list[list[Any]] = []
    for candidate in cast(List[Any], resonance.get("candidates", [])):
        if not isinstance(candidate, dict):
            continue
        state = cast(Dict[str, Any], candidate.get("state", {}))
        rows.append(
            [
                candidate.get("resonance_index", ""),
                candidate.get("raw_state_key_hash", ""),
                candidate.get("processed_state_key_hash", ""),
                candidate.get("process_resonance_hit", ""),
                candidate.get("duplicate_processed_state", ""),
                candidate.get("valid_for_target", ""),
                candidate.get("score", ""),
                candidate.get("organic_topology_selection_key", ""),
                state.get("canonical_smiles", state.get("smiles", "")),
                state.get("formal_charge_sum", ""),
                state.get("spin_multiplicity_sum", ""),
                state.get("spin_multiplicity_singlet_sum", ""),
                state.get("charged_atom_counts", ""),
                state.get("radical_atom_counts", ""),
                candidate.get("score_error", ""),
            ]
        )
    return _markdown_table(
        (
            "共振序号",
            "raw key",
            "processed key",
            "process 命中",
            "重复",
            "匹配目标",
            "分数",
            "选择键",
            "canonical SMILES",
            "形式电荷",
            "自由基和",
            "自由基奇偶和",
            "带电原子",
            "自由基原子",
            "评分错误",
        ),
        rows,
    )


def _render_no_metal_trace(no_metal_trace: dict[str, Any]) -> str:
    target = cast(Dict[str, Any], no_metal_trace.get("target", {}))
    selected = cast(Dict[str, Any], no_metal_trace.get("selected_candidate", {}))
    selected_state = cast(Dict[str, Any], selected.get("state", {}))
    lines = [
        _markdown_kv_table(
            (
                ("状态", no_metal_trace.get("status", "")),
                ("目标电荷", target.get("total_charge", "")),
                ("目标自由基电子", target.get("total_radical_electrons", "")),
                ("线性后直接有效", no_metal_trace.get("direct_validation", "")),
                ("选中 canonical SMILES", selected_state.get("canonical_smiles", "")),
                ("选中分数", selected.get("score", "")),
                ("选中选择键", selected.get("organic_topology_selection_key", "")),
            )
        ),
        "",
        "**线性阶段**",
        _render_no_metal_linear_steps(no_metal_trace),
    ]
    resonance = cast(Dict[str, Any], no_metal_trace.get("resonance", {}))
    if resonance:
        lines.extend(
            [
                "",
                "**共振枚举**",
                _markdown_kv_table(
                    (
                        ("遍历策略", resonance.get("traversal_policy", "")),
                        ("最大深度", resonance.get("max_depth", "")),
                        ("候选数", resonance.get("candidate_count", "")),
                        ("有效去重候选数", resonance.get("valid_unique_candidate_count", "")),
                    )
                ),
                _render_no_metal_resonance_candidates(no_metal_trace),
            ]
        )
    direct_candidate = cast(Dict[str, Any], no_metal_trace.get("direct_candidate", {}))
    if direct_candidate:
        lines.extend(
            [
                "",
                "**直接候选**",
                _markdown_kv_table(
                    (
                        ("clean_resonances 命中", direct_candidate.get("clean_resonances_hit", "")),
                        ("评分错误", direct_candidate.get("score_error", "")),
                        ("分数", direct_candidate.get("score", "")),
                        ("选择键", direct_candidate.get("organic_topology_selection_key", "")),
                    )
                ),
            ]
        )
    return "\n".join(lines)


def _render_no_metal_traces(case: dict[str, Any]) -> str:
    search_summary = cast(Dict[str, Any], case.get("search", {}))
    layer_summaries = cast(List[Any], search_summary.get("layer_summaries", []))
    sections: list[str] = []
    for layer in layer_summaries:
        if not isinstance(layer, dict):
            continue
        for target_bucket in cast(List[Any], layer.get("target_buckets", [])):
            if not isinstance(target_bucket, dict):
                continue
            no_metal_trace = target_bucket.get("no_metal_trace")
            if not isinstance(no_metal_trace, dict):
                continue
            target = cast(Dict[str, Any], target_bucket.get("target", {}))
            summary = (
                f"L{layer.get('layer_index', '')} | "
                f"Q={target.get('no_metal_charge', '')}, "
                f"R={target.get('no_metal_radical_electrons', '')} | "
                f"{no_metal_trace.get('status', '')}"
            )
            sections.append(
                "\n".join(
                    [
                        f"<details><summary>{_markdown_cell(summary)}</summary>",
                        "",
                        _render_no_metal_trace(no_metal_trace),
                        "",
                        "</details>",
                    ]
                )
            )
    return "\n\n".join(sections) if sections else "_无_"


def _render_candidate_overview(case: dict[str, Any]) -> str:
    candidates = cast(List[Any], case.get("candidates", []))
    rows: list[list[Any]] = []
    for raw_candidate in candidates:
        if not isinstance(raw_candidate, dict):
            continue
        target = cast(Dict[str, Any], raw_candidate.get("target", {}))
        organic_part = cast(Dict[str, Any], raw_candidate.get("organic_part", {}))
        rows.append(
            [
                "是" if raw_candidate.get("selected") else "",
                raw_candidate.get("search_layer_index", ""),
                raw_candidate.get("combination_index", ""),
                "是" if raw_candidate.get("in_production_selection_layer") else "",
                _metal_states_label(
                    cast(List[Dict[str, Any]], raw_candidate.get("metal_states", []))
                ),
                organic_part.get("canonical_smiles", organic_part.get("smiles", "")),
                target.get("no_metal_charge", ""),
                target.get("no_metal_radical_electrons", ""),
                _score_detail_value(raw_candidate, "metal_discordance_count", ""),
                _score_detail_value(
                    raw_candidate,
                    "metal_discordance_structural_count",
                    "",
                ),
                _score_detail_value(
                    raw_candidate,
                    "metal_discordance_aromatic_ring_deficit_count",
                    "",
                ),
                _score_detail_value(
                    raw_candidate,
                    "metal_discordance_inner_visible_diradical_count",
                    "",
                ),
                _score_detail_value(
                    raw_candidate,
                    "metal_discordance_outer_or_invisible_adjacent_double_charge_count",
                    "",
                ),
                _score_detail_value(
                    raw_candidate,
                    "metal_discordance_outer_or_invisible_adjacent_same_sign_double_charge_count",
                    "",
                ),
                _score_detail_value(
                    raw_candidate,
                    "metal_discordance_outer_or_invisible_adjacent_opposite_sign_double_charge_count",
                    "",
                ),
                _score_detail_value(
                    raw_candidate,
                    "metal_discordance_inner_visible_adjacent_carbanion_pair_count",
                    "",
                ),
                _score_detail_value(
                    raw_candidate,
                    "metal_discordance_inner_visible_conjugated_carbanion_pair_count",
                    "",
                ),
                _score_detail_value(
                    raw_candidate,
                    "metal_discordance_inner_visible_same_sign_charge_count",
                    "",
                ),
                _score_detail_value(
                    raw_candidate,
                    "metal_discordance_negative_metal_count",
                    "",
                ),
                _score_detail_value(
                    raw_candidate,
                    "metal_discordance_zero_valent_metals_with_organic_cation_count",
                    "",
                ),
                _score_detail_value(
                    raw_candidate,
                    "metal_discordance_negative_metal_outer_sphere_cation_exception",
                    "",
                ),
                _score_detail_value(
                    raw_candidate,
                    "metal_discordance_negative_metal_positive_metal_counterion_exception",
                    "",
                ),
                _score_detail_value(raw_candidate, "score", ""),
            ]
        )
    return _markdown_table(
        (
            "选中",
            "层",
            "组合",
            "生产层",
            "金属状态",
            "有机 canonical SMILES",
            "有机电荷",
            "有机自由基",
            "失谐",
            "结构失谐",
            "芳环损失",
            "内圈可见双自由基",
            "外圈/不可见邻位双电荷",
            "外圈/不可见相对金属同号双电荷",
            "外圈/不可见相对金属异号双电荷",
            "内圈可见相邻同号碳离子",
            "内圈可见共轭同号碳离子",
            "内圈可见同号电荷",
            "负价金属",
            "零价金属有机阳离子",
            "负价金属外圈H+例外",
            "负价金属阳离子金属例外",
            "有机分",
        ),
        rows,
    )


def _render_candidate_details(case: dict[str, Any]) -> str:
    candidates = cast(List[Any], case.get("candidates", []))
    sections: list[str] = []
    for raw_candidate in candidates:
        if not isinstance(raw_candidate, dict):
            continue
        organic_part = cast(Dict[str, Any], raw_candidate.get("organic_part", {}))
        summary = (
            f"{_candidate_short_label(raw_candidate)} | "
            f"失谐 {_format_value(_score_detail_value(raw_candidate, 'metal_discordance_count', ''))} | "
            f"{organic_part.get('canonical_smiles', organic_part.get('smiles', ''))}"
        )
        lines = [f"<details><summary>{_markdown_cell(summary)}</summary>", ""]
        lines.append("**金属状态**")
        metal_rows = []
        for metal_state in cast(List[Any], raw_candidate.get("metal_states", [])):
            if not isinstance(metal_state, dict):
                continue
            metal_rows.append(
                [
                    metal_state.get("idx", ""),
                    metal_state.get("symbol", ""),
                    metal_state.get("valence", ""),
                    metal_state.get("radical_num", ""),
                    metal_state.get("position", ""),
                ]
            )
        lines.append(_markdown_table(("原子序号", "元素", "价态", "自由基数", "坐标"), metal_rows))
        lines.extend(["", "**分数明细**"])
        lines.append(_markdown_table(("分数项", "值"), _candidate_score_rows(raw_candidate)))

        production_metadata = cast(Dict[str, Any], raw_candidate.get("production_metadata", {}))
        if production_metadata:
            production_rows = [
                (key, value)
                for key, value in production_metadata.items()
                if key in _IMPORTANT_SCORE_KEYS or key.startswith(("passes_", "selection_"))
            ]
            if production_rows:
                lines.extend(["", "**生产选择元数据**"])
                lines.append(_markdown_table(("字段", "值"), production_rows))

        lines.extend(["", "</details>"])
        sections.append("\n".join(lines))
    return "\n\n".join(sections) if sections else "_无_"


def _render_case_report(case: dict[str, Any]) -> str:
    lines = [f"## {case.get('id', 'unknown')}", ""]
    if case.get("status") == "error":
        lines.append(_markdown_kv_table((("状态", "error"), ("错误", case.get("error", "")))))
        return "\n".join(lines)
    if case.get("status") == "missing_csv_row":
        lines.append(_markdown_kv_table((("状态", "missing_csv_row"),)))
        return "\n".join(lines)

    if case.get("trace_kind") == "no_metal":
        lines.extend(
            [
                "### 基本信息",
                _render_case_summary(case),
                "",
                "### 无金属重建过程",
                _render_no_metal_trace(cast(Dict[str, Any], case.get("no_metal_trace", {}))),
            ]
        )
        return "\n".join(lines)

    lines.extend(
        [
            "### 基本信息",
            _render_case_summary(case),
            "",
            "### 金属位点候选",
            _render_available_metal_states(case),
            "",
            "### 搜索层",
            _render_layer_table(case),
            "",
            "### 有机目标桶",
            _render_target_bucket_table(case),
            "",
            "### 无金属重建过程",
            _render_no_metal_traces(case),
            "",
            "### 候选总览",
            _render_candidate_overview(case),
            "",
            "### 候选分数明细",
            _render_candidate_details(case),
        ]
    )
    return "\n".join(lines)


def _render_markdown_report(output: dict[str, Any]) -> str:
    input_summary = cast(Dict[str, Any], output.get("input", {}))
    source = input_summary.get("source", "")
    title = "# tmQMg 重建轨迹分析" if source == "tmQMg" else "# MolGR 重建轨迹分析"
    input_rows: list[tuple[str, Any]] = [
        ("来源", source),
    ]
    if input_summary.get("csv", "") != "":
        input_rows.append(("CSV", input_summary.get("csv", "")))
    if input_summary.get("xyz_dir", "") != "":
        input_rows.append(("XYZ 目录", input_summary.get("xyz_dir", "")))
    input_rows.extend(
        [
            ("id", input_summary.get("ids", [])),
            ("默认总电荷", input_summary.get("total_charge", "")),
            ("默认总自由基电子数", input_summary.get("total_radical_electrons", "")),
        ]
    )
    if input_summary.get("spin_source", "") != "":
        input_rows.append(("自旋来源", input_summary.get("spin_source", "")))
    input_rows.append(("样本数", output.get("case_count", "")))
    lines = [
        title,
        "",
        "## 输入",
        _markdown_kv_table(input_rows),
    ]
    for case in cast(List[Any], output.get("cases", [])):
        if isinstance(case, dict):
            lines.extend(["", _render_case_report(case)])
    return "\n".join(lines).rstrip() + "\n"


def _score_details(metadata: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in metadata.items()
        if key in _SCORE_DETAIL_KEYS or key.startswith(_SCORE_DETAIL_PREFIXES)
    }


def _candidate_identity(
    candidate: MetalCandidateState,
    *,
    fallback_candidate_index: int = 0,
) -> tuple[int, int]:
    return (
        int(candidate.metadata.get("search_layer_index", -1)),
        int(candidate.metadata.get("combination_index", fallback_candidate_index)),
    )


def _copy_metadata_by_identity(
    candidates: Sequence[MetalCandidateState],
) -> dict[tuple[int, int], dict[str, Any]]:
    copied: dict[tuple[int, int], dict[str, Any]] = {}
    for candidate_index, candidate in enumerate(candidates):
        copied[_candidate_identity(candidate, fallback_candidate_index=candidate_index)] = deepcopy(
            candidate.metadata
        )
    return copied


def _record_no_metal_stage(
    steps: list[dict[str, Any]],
    machine: OmolStateMachine,
    *,
    phase: str,
    hit: bool | None,
    target_charge: int,
    target_radical_electrons: int,
    kind: str = "stage",
) -> None:
    steps.append(
        {
            "step_index": len(steps),
            "phase": phase,
            "kind": kind,
            "hit": hit,
            "omol_revision": int(machine.omol_revision),
            "phase_history_length": len(machine.phase_history),
            "state": _omol_state_snapshot(
                machine.omol,
                given_charge=machine.given_charge,
                target_charge=target_charge,
                target_radical_electrons=target_radical_electrons,
            ),
        }
    )


def _run_no_metal_linear_trace(
    seed_state: ReconstructionState,
) -> tuple[ReconstructionState, list[dict[str, Any]]]:
    machine = OmolStateMachine.from_reconstruction_state(seed_state)
    steps: list[dict[str, Any]] = []
    target_charge = seed_state.total_charge
    target_radicals = seed_state.total_radical_electrons
    _record_no_metal_stage(
        steps,
        machine,
        phase="read_xyz",
        hit=None,
        target_charge=target_charge,
        target_radical_electrons=target_radicals,
        kind="seed",
    )

    for phase, stage, args in (
        ("make_connections", make_connections, ()),
        ("pre_clean", pre_clean, ()),
        ("fresh_omol_charge_radical_initial", fresh_omol_charge_radical, ()),
    ):
        hit = machine.run_omol_stage(phase, stage, *args)
        _record_no_metal_stage(
            steps,
            machine,
            phase=phase,
            hit=hit,
            target_charge=target_charge,
            target_radical_electrons=target_radicals,
        )

    initial_given_charge = target_charge - sum(
        cast(ob.OBAtom, atom.OBAtom).GetFormalCharge() for atom in machine.omol.atoms
    )
    machine.set_given_charge("initialize_charge_budget", initial_given_charge)
    _record_no_metal_stage(
        steps,
        machine,
        phase="initialize_charge_budget",
        hit=None,
        target_charge=target_charge,
        target_radical_electrons=target_radicals,
        kind="charge_budget",
    )

    charge_stages: tuple[tuple[str, Any, tuple[Any, ...]], ...] = (
        ("eliminate_NNN_negative", eliminate_NNN, (False,)),
        ("eliminate_high_positive_charge_atoms", eliminate_high_positive_charge_atoms, ()),
        ("eliminate_CN_in_doubt", eliminate_CN_in_doubt, ()),
        ("eliminate_NNN_positive", eliminate_NNN, (True,)),
        ("eliminate_carboxyl", eliminate_carboxyl, ()),
    )
    for phase, stage, args in charge_stages:
        hit = machine.run_omol_charge_stage(phase, stage, *args)
        _record_no_metal_stage(
            steps,
            machine,
            phase=phase,
            hit=hit,
            target_charge=target_charge,
            target_radical_electrons=target_radicals,
        )

    hit = machine.run_omol_stage(
        "clean_carbene_neighbor_unsaturated_first",
        clean_carbene_neighbor_unsaturated,
    )
    _record_no_metal_stage(
        steps,
        machine,
        phase="clean_carbene_neighbor_unsaturated_first",
        hit=hit,
        target_charge=target_charge,
        target_radical_electrons=target_radicals,
    )
    hit = machine.run_omol_charge_stage(
        "eliminate_carbene_neighbor_heteroatom",
        eliminate_carbene_neighbor_heteroatom,
    )
    _record_no_metal_stage(
        steps,
        machine,
        phase="eliminate_carbene_neighbor_heteroatom",
        hit=hit,
        target_charge=target_charge,
        target_radical_electrons=target_radicals,
    )
    for phase, stage in (
        ("clean_neighbor_radicals", clean_neighbor_radicals),
        ("clean_carbene_neighbor_unsaturated_second", clean_carbene_neighbor_unsaturated),
    ):
        hit = machine.run_omol_stage(phase, stage)
        _record_no_metal_stage(
            steps,
            machine,
            phase=phase,
            hit=hit,
            target_charge=target_charge,
            target_radical_electrons=target_radicals,
        )
    hit = machine.run_omol_charge_stage("eliminate_charge_spliting", eliminate_charge_spliting)
    _record_no_metal_stage(
        steps,
        machine,
        phase="eliminate_charge_spliting",
        hit=hit,
        target_charge=target_charge,
        target_radical_electrons=target_radicals,
    )
    hit = machine.run_omol_stage(
        "break_deformed_ene",
        break_deformed_ene,
        machine.given_charge,
        target_radicals,
        5.0,
    )
    _record_no_metal_stage(
        steps,
        machine,
        phase="break_deformed_ene",
        hit=hit,
        target_charge=target_charge,
        target_radical_electrons=target_radicals,
    )
    hit = machine.run_omol_charge_stage("break_one_bond", break_one_bond, target_radicals)
    _record_no_metal_stage(
        steps,
        machine,
        phase="break_one_bond",
        hit=hit,
        target_charge=target_charge,
        target_radical_electrons=target_radicals,
    )
    hit = machine.run_omol_stage("fresh_omol_charge_radical_final", fresh_omol_charge_radical)
    _record_no_metal_stage(
        steps,
        machine,
        phase="fresh_omol_charge_radical_final",
        hit=hit,
        target_charge=target_charge,
        target_radical_electrons=target_radicals,
    )

    return machine.freeze_like(seed_state), steps


def _no_metal_candidate_trace(candidate: ReconstructionState, *, selected: bool) -> dict[str, Any]:
    selection_key = no_metal_selection._no_metal_candidate_selection_key(candidate)
    return {
        "selected": selected,
        "phase_history": list(candidate.phase_history),
        "state": _omol_state_snapshot(
            candidate.omol,
            given_charge=candidate.given_charge,
            target_charge=candidate.total_charge,
            target_radical_electrons=candidate.total_radical_electrons,
        ),
        "score": candidate.metadata.get("score"),
        "force_field_energy": candidate.metadata.get("force_field_energy"),
        "organic_topology_selection_key": selection_key,
        "metadata": _jsonable(candidate.metadata),
    }


def _trace_no_metal_reconstruction(
    xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
    *,
    config: MolGRConfig | None = None,
) -> dict[str, Any]:
    seed_state = no_metal_preparation._seed_state(
        xyz_block,
        total_charge,
        total_radical_electrons,
    )
    trace: dict[str, Any] = {
        "target": {
            "total_charge": total_charge,
            "total_radical_electrons": total_radical_electrons,
        },
        "status": "pending",
    }
    if seed_state.total_radical_electrons < 0:
        trace["status"] = "invalid_negative_radicals"
        return trace

    linear_state, linear_steps = _run_no_metal_linear_trace(seed_state)
    trace["linear_steps"] = linear_steps
    trace["linear_result"] = {
        "phase_history": list(linear_state.phase_history),
        "state": _omol_state_snapshot(
            linear_state.omol,
            given_charge=linear_state.given_charge,
            target_charge=total_charge,
            target_radical_electrons=total_radical_electrons,
        ),
    }

    direct_valid = validate_omol(
        linear_state.omol,
        total_charge,
        total_radical_electrons,
    )
    trace["direct_validation"] = direct_valid
    if direct_valid:
        result_machine = OmolStateMachine.from_reconstruction_state(linear_state)
        result_machine.annotate("validate_direct_candidate")
        clean_hit = result_machine.run_omol_stage("clean_resonances", clean_resonances)
        direct_candidate = result_machine.freeze_like(linear_state)
        direct_trace: dict[str, Any] = {
            "clean_resonances_hit": clean_hit,
            "score_error": None,
        }
        try:
            no_metal_selection._score_reconstruction_candidate(direct_candidate, config=config)
            no_metal_selection._annotate_no_metal_candidate_topology(direct_candidate)
        except ValueError as exc:
            direct_trace["score_error"] = str(exc)
            trace["status"] = "direct_score_error"
            trace["direct_candidate"] = direct_trace
            return trace
        direct_trace.update(_no_metal_candidate_trace(direct_candidate, selected=True))
        trace["status"] = "direct_valid"
        trace["direct_candidate"] = direct_trace
        trace["selected_candidate"] = direct_trace
        return trace

    traversal_policy = no_metal_resonance._default_resonance_traversal_policy(config)
    resonance_max_depth = no_metal_resonance._resonance_max_depth(config)
    resonance_iterable = no_metal_resonance.get_radical_resonances(
        linear_state.omol,
        max_depth=resonance_max_depth,
        traversal_policy=traversal_policy,
    )

    resonance_candidates: list[ReconstructionState] = []
    resonance_reports: list[dict[str, Any]] = []
    seen_processed_states = set()
    base_machine = OmolStateMachine.from_reconstruction_state(linear_state)
    for resonance_index, resonance in enumerate(resonance_iterable):
        raw_state_key = resonance_utils.build_resonance_state_key(resonance)
        candidate_machine = base_machine.branch("branch_resonance_candidate", omol=resonance)
        process_hit = candidate_machine.run_omol_charge_stage(
            "process_resonance",
            resonance_utils.process_resonance,
        )
        processed_state_key = candidate_machine.get_cached_omol_value(
            "resonance_state_key",
            resonance_utils.build_processed_resonance_key,
        )
        duplicate = processed_state_key in seen_processed_states
        if not duplicate:
            seen_processed_states.add(processed_state_key)
        valid = validate_omol(
            candidate_machine.omol,
            total_charge,
            total_radical_electrons,
        )
        report: dict[str, Any] = {
            "resonance_index": resonance_index,
            "raw_state_key_hash": _short_hash(raw_state_key),
            "processed_state_key_hash": _short_hash(processed_state_key),
            "process_resonance_hit": process_hit,
            "duplicate_processed_state": duplicate,
            "valid_for_target": valid,
            "state": _omol_state_snapshot(
                candidate_machine.omol,
                given_charge=candidate_machine.given_charge,
                target_charge=total_charge,
                target_radical_electrons=total_radical_electrons,
            ),
        }
        if duplicate or not valid:
            resonance_reports.append(report)
            continue

        candidate_machine.annotate("validate_resonance_candidate", resonance_index=resonance_index)
        candidate = candidate_machine.freeze_like(linear_state)
        try:
            if config is None:
                candidate.score("organic_core")
            else:
                candidate.score("organic_core", config=config)
            no_metal_selection._score_reconstruction_candidate(candidate, config=config)
            no_metal_selection._annotate_no_metal_candidate_topology(candidate)
            report["score"] = candidate.metadata.get("score")
            report["organic_topology_selection_key"] = (
                no_metal_selection._no_metal_candidate_selection_key(candidate)
            )
            report["phase_history"] = list(candidate.phase_history)
            report["metadata"] = _jsonable(candidate.metadata)
            resonance_candidates.append(candidate)
        except ValueError as exc:
            report["score_error"] = str(exc)
        resonance_reports.append(report)

    trace["resonance"] = {
        "traversal_policy": type(traversal_policy).__name__,
        "max_depth": resonance_max_depth,
        "candidate_count": len(resonance_reports),
        "valid_unique_candidate_count": len(resonance_candidates),
        "candidates": resonance_reports,
    }
    if not resonance_candidates:
        trace["status"] = "no_valid_resonance_candidate"
        return trace

    best_candidate: Optional[ReconstructionState] = None
    best_selection_key: Optional[tuple[int, int, int, int, float]] = None
    for candidate in resonance_candidates:
        selection_key = no_metal_selection._no_metal_candidate_selection_key(candidate)
        if best_selection_key is not None and selection_key >= best_selection_key:
            continue
        best_selection_key = selection_key
        best_candidate = candidate

    if best_candidate is None:
        trace["status"] = "no_best_resonance_candidate"
        return trace

    result_machine = OmolStateMachine.from_reconstruction_state(best_candidate)
    result_machine.annotate("select_best_resonance_candidate")
    selected = result_machine.freeze_like(best_candidate)
    trace["status"] = "resonance_selected"
    trace["selected_candidate"] = _no_metal_candidate_trace(selected, selected=True)
    return trace


def _annotate_analysis_scores_for_all_candidates(
    candidates: Sequence[MetalCandidateState],
    *,
    config: MolGRConfig | None = None,
) -> dict[str, Any]:
    if not candidates:
        return {}

    scoring._annotate_candidate_set_discordance_features(candidates)
    for candidate in candidates:
        scoring._annotate_selected_candidate_metrics(candidate, config=config)

    for candidate in candidates:
        score_value = float(cast(float, candidate.score))
        candidate.metadata["analysis_selection_key_all_candidates"] = (
            float(candidate.metadata.get("metal_discordance_count", 0)),
            score_value,
            int(candidate.metadata.get("combination_index", 0)),
        )

    return {}


def _candidate_report(
    candidate: MetalCandidateState,
    *,
    candidate_index: int,
    selected_candidate_identity: Optional[tuple[int, int]],
    production_metadata_by_identity: dict[tuple[int, int], dict[str, Any]],
) -> dict[str, Any]:
    combination_index = int(candidate.metadata.get("combination_index", candidate_index))
    candidate_identity = _candidate_identity(
        candidate,
        fallback_candidate_index=candidate_index,
    )
    no_metal_state = candidate.no_metal_state
    no_metal: dict[str, Any] = {}
    if no_metal_state is not None:
        no_metal = {
            "target_charge": int(no_metal_state.total_charge),
            "target_radical_electrons": int(no_metal_state.total_radical_electrons),
            "given_charge": int(no_metal_state.given_charge),
            "phase_history": list(no_metal_state.phase_history),
            "smiles": _omol_smiles_detail(no_metal_state.omol),
            "score": no_metal_state.metadata.get("score"),
            "metadata": _jsonable(no_metal_state.metadata),
        }

    metadata = candidate.metadata
    production_metadata = production_metadata_by_identity.get(candidate_identity, {})
    return {
        "candidate_index": candidate_index,
        "candidate_identity": {
            "search_layer_index": candidate_identity[0],
            "combination_index": candidate_identity[1],
        },
        "search_layer_index": candidate_identity[0],
        "combination_index": combination_index,
        "selected": selected_candidate_identity == candidate_identity,
        "in_production_selection_layer": bool(production_metadata),
        "target": {
            "no_metal_charge": int(candidate.no_metal_charge_target),
            "no_metal_radical_electrons": int(candidate.no_metal_radical_target),
        },
        "metal_states": [
            _metal_state_to_dict(metal_state) for metal_state in candidate.metal_states
        ],
        "organic_part": {
            "smiles": no_metal.get("smiles", {}).get("smiles", ""),
            "canonical_smiles": no_metal.get("smiles", {}).get("canonical_smiles", ""),
        },
        "score": _jsonable(candidate.score),
        "score_details": _jsonable(_score_details(metadata)),
        "metadata": _jsonable(metadata),
        "production_metadata": _jsonable(production_metadata),
        "no_metal_state": no_metal,
        "phase_history": list(candidate.phase_history),
    }


def _target_summary(
    target: tuple[int, int],
    candidates: Sequence[MetalCandidateState],
) -> dict[str, Any]:
    return {
        "target": {
            "no_metal_charge": int(target[0]),
            "no_metal_radical_electrons": int(target[1]),
        },
        "candidate_count": len(candidates),
        "prepared_candidate_count": 0,
        "status": "pending",
    }


def _trace_candidates(
    xyz_block: str,
    *,
    total_charge: int,
    total_radical_electrons: int,
    score_all_candidates: bool,
    config: MolGRConfig | None = None,
) -> dict[str, Any]:
    base_state = preparation.prepare_metal_state(
        xyz_block,
        total_charge,
        total_radical_electrons,
        config=config,
    )
    state_search_groups = search._build_metal_state_search_groups(
        base_state.available_valence_radical_states,
        config=config,
    )
    layered_state_search_groups = search._build_layered_metal_state_search_groups(
        state_search_groups,
        total_radical_electrons,
        config=config,
    )

    layer_summaries: list[dict[str, Any]] = []
    all_scored_candidates: list[MetalCandidateState] = []
    production_scored_candidates: list[MetalCandidateState] = []
    scored_candidates_by_layer: dict[int, list[MetalCandidateState]] = {}
    selected_layer_index: Optional[int] = None

    for layer_index, available_valence_radical_states in enumerate(layered_state_search_groups):
        grouped_candidates = search._group_candidates_by_target_dp(
            base_state.phase_history,
            available_valence_radical_states,
            total_charge,
            total_radical_electrons,
            config=config,
        )
        layer_summary: dict[str, Any] = {
            "layer_index": layer_index,
            "state_group_count": len(available_valence_radical_states),
            "state_options_per_group": [
                len(state_options) for state_options in available_valence_radical_states
            ],
            "target_bucket_count": len(grouped_candidates),
            "candidate_count": sum(len(candidates) for candidates in grouped_candidates.values()),
            "prepared_candidate_count": 0,
            "target_buckets": [],
            "status": "no_candidate_targets" if not grouped_candidates else "pending",
        }
        if not grouped_candidates:
            layer_summaries.append(layer_summary)
            continue

        current_layer_scored_candidates: list[MetalCandidateState] = []
        for target, candidates in grouped_candidates.items():
            target_summary = _target_summary(target, candidates)
            if not candidates:
                target_summary["status"] = "empty_candidate_bucket"
                layer_summary["target_buckets"].append(target_summary)
                continue

            prototype = candidates[0]
            try:
                target_summary["no_metal_trace"] = _trace_no_metal_reconstruction(
                    base_state.no_metal_xyz_block,
                    prototype.no_metal_charge_target,
                    prototype.no_metal_radical_target,
                    config=config,
                )
            except Exception as exc:
                target_summary["no_metal_trace"] = {
                    "status": "trace_error",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            try:
                no_metal_state = reconstruct_without_metals.xyz_to_omol_no_metal_state(
                    base_state.no_metal_xyz_block,
                    prototype.no_metal_charge_target,
                    prototype.no_metal_radical_target,
                    config=config,
                )
            except (OSError, ValueError) as exc:
                target_summary["status"] = "no_metal_reconstruction_error"
                target_summary["error"] = str(exc)
                layer_summary["target_buckets"].append(target_summary)
                continue

            if no_metal_state is None:
                target_summary["status"] = "no_metal_reconstruction_none"
                layer_summary["target_buckets"].append(target_summary)
                continue

            target_summary["status"] = "prepared"
            target_summary["organic_part"] = _omol_smiles_detail(no_metal_state.omol)
            target_summary["no_metal_phase_history"] = list(no_metal_state.phase_history)

            for candidate in candidates:
                try:
                    scored_candidate = scoring._prepare_candidate_with_no_metal_state(
                        candidate,
                        no_metal_state,
                        config=config,
                    )
                except ValueError as exc:
                    target_summary.setdefault("candidate_errors", []).append(
                        {
                            "combination_index": candidate.metadata.get("combination_index"),
                            "error": str(exc),
                        }
                    )
                    continue
                if cast(Optional[float], scored_candidate.score) is None:
                    continue
                scored_candidate.metadata["search_layer_index"] = layer_index
                current_layer_scored_candidates.append(scored_candidate)

            target_summary["no_metal_score"] = no_metal_state.metadata.get("score")
            target_summary["prepared_candidate_count"] = len(
                [
                    candidate
                    for candidate in current_layer_scored_candidates
                    if candidate.no_metal_charge_target == target[0]
                    and candidate.no_metal_radical_target == target[1]
                ]
            )
            layer_summary["target_buckets"].append(target_summary)

        layer_summary["prepared_candidate_count"] = len(current_layer_scored_candidates)
        if current_layer_scored_candidates:
            scored_candidates_by_layer[layer_index] = current_layer_scored_candidates
            all_scored_candidates.extend(current_layer_scored_candidates)
            if selected_layer_index is None:
                selected_layer_index = layer_index
                production_scored_candidates = current_layer_scored_candidates
                layer_summary["status"] = "production_scoring_layer"
                layer_summary["production_selected_layer"] = True
            else:
                layer_summary["status"] = "additional_scored_layer"
                layer_summary["production_selected_layer"] = False
            layer_summaries.append(layer_summary)
            continue

        layer_summary["status"] = "no_scored_candidates"
        layer_summary["production_selected_layer"] = False
        layer_summaries.append(layer_summary)

    selected_candidate = None
    production_metadata_by_identity: dict[tuple[int, int], dict[str, Any]] = {}
    analysis_score_context: dict[str, Any] = {}
    selected_candidate_identity: Optional[tuple[int, int]] = None
    if production_scored_candidates:
        selected_candidate = scoring.select_best_candidate(
            production_scored_candidates, config=config
        )
        production_metadata_by_identity = _copy_metadata_by_identity(production_scored_candidates)
        if selected_candidate is not None:
            selected_candidate_identity = _candidate_identity(
                selected_candidate,
                fallback_candidate_index=-1,
            )

        if score_all_candidates:
            analysis_score_context = {
                str(layer_index): _annotate_analysis_scores_for_all_candidates(
                    layer_candidates,
                    config=config,
                )
                for layer_index, layer_candidates in sorted(scored_candidates_by_layer.items())
            }

    if score_all_candidates and analysis_score_context:
        for layer_summary in layer_summaries:
            layer_context = analysis_score_context.get(str(layer_summary["layer_index"]))
            if layer_context is not None:
                layer_summary["analysis_score_context"] = _jsonable(layer_context)

    candidate_reports = [
        _candidate_report(
            candidate,
            candidate_index=candidate_index,
            selected_candidate_identity=selected_candidate_identity,
            production_metadata_by_identity=production_metadata_by_identity,
        )
        for candidate_index, candidate in enumerate(all_scored_candidates)
    ]

    selected_summary = None
    if selected_candidate_identity is not None:
        selected_summary = next(
            (
                {
                    "candidate_index": candidate["candidate_index"],
                    "combination_index": candidate["combination_index"],
                    "target": candidate["target"],
                    "metal_states": candidate["metal_states"],
                    "organic_part": candidate["organic_part"],
                    "score": candidate["score"],
                    "score_details": candidate["score_details"],
                }
                for candidate in candidate_reports
                if (
                    candidate["candidate_identity"]["search_layer_index"],
                    candidate["candidate_identity"]["combination_index"],
                )
                == selected_candidate_identity
            ),
            None,
        )

    return {
        "status": "ok" if production_scored_candidates else "no_scored_metal_candidates",
        "base_state": {
            "metal_atom_count": int(base_state.metadata.get("metal_atom_count", 0)),
            "available_metal_states_by_site": [
                [_metal_state_to_dict(metal_state) for metal_state in state_options]
                for state_options in base_state.available_valence_radical_states
            ],
            "no_metal_xyz_block": base_state.no_metal_xyz_block,
            "phase_history": list(base_state.phase_history),
            "metadata": _jsonable(base_state.metadata),
        },
        "search": {
            "state_search_group_count": len(state_search_groups),
            "layer_count": len(layered_state_search_groups),
            "selected_layer_index": selected_layer_index,
            "layer_summaries": _jsonable(layer_summaries),
        },
        "analysis": {
            "score_all_candidates": score_all_candidates,
            "score_context": _jsonable(analysis_score_context),
            "note": (
                "production_metadata records the exact metadata after select_best_candidate; "
                "score_details also includes analysis_* metrics when score-all-candidates is enabled."
            ),
        },
        "candidate_count": len(candidate_reports),
        "production_candidate_count": len(production_scored_candidates),
        "selected_candidate": selected_summary,
        "candidates": candidate_reports,
    }


def _metal_atom_count_from_xyz(xyz_block: str) -> int:
    omol = pybel.readstring("xyz", xyz_block)
    return sum(1 for atom in omol.atoms if atom.OBAtom.IsMetal())


def trace_reconstruction_case(
    input_case: TraceInputCase,
    *,
    score_all_candidates: bool,
    config: MolGRConfig | None = None,
) -> dict[str, Any]:
    started = time.perf_counter()
    total_charge = int(input_case.total_charge)
    total_radicals = int(input_case.total_radical_electrons)
    metal_atom_count = _metal_atom_count_from_xyz(input_case.xyz_block)

    if metal_atom_count == 0:
        no_metal_trace = _trace_no_metal_reconstruction(
            input_case.xyz_block,
            total_charge,
            total_radicals,
            config=config,
        )
        trace = {
            "status": no_metal_trace.get("status", "unknown"),
            "trace_kind": "no_metal",
            "base_state": {
                "metal_atom_count": 0,
                "phase_history": ["read_xyz"],
            },
            "no_metal_trace": no_metal_trace,
            "selected_candidate": no_metal_trace.get("selected_candidate"),
            "candidate_count": 1 if no_metal_trace.get("selected_candidate") else 0,
            "production_candidate_count": 1 if no_metal_trace.get("selected_candidate") else 0,
        }
    else:
        trace = _trace_candidates(
            input_case.xyz_block,
            total_charge=total_charge,
            total_radical_electrons=total_radicals,
            score_all_candidates=score_all_candidates,
            config=config,
        )
        trace["trace_kind"] = "metal"

    trace.update(
        {
            "id": input_case.id,
            "charge": total_charge,
            "xyz_path": str(input_case.xyz_path) if input_case.xyz_path is not None else "",
            "xyz_source": input_case.xyz_source,
            "total_radical_electrons": total_radicals,
            "spin_multiplicity": total_radicals + 1,
            "elapsed_seconds": time.perf_counter() - started,
        }
    )
    return trace


def trace_reconstruction_cases(
    input_cases: Sequence[TraceInputCase],
    *,
    score_all_candidates: bool = True,
    config: MolGRConfig | None = None,
) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for input_case in input_cases:
        started = time.perf_counter()
        try:
            cases.append(
                trace_reconstruction_case(
                    input_case,
                    score_all_candidates=score_all_candidates,
                    config=config,
                )
            )
        except Exception as exc:
            cases.append(
                {
                    "id": input_case.id,
                    "status": "error",
                    "trace_kind": "unknown",
                    "charge": input_case.total_charge,
                    "total_radical_electrons": input_case.total_radical_electrons,
                    "spin_multiplicity": input_case.total_radical_electrons + 1,
                    "xyz_path": str(input_case.xyz_path) if input_case.xyz_path is not None else "",
                    "xyz_source": input_case.xyz_source,
                    "error": f"{type(exc).__name__}: {exc}",
                    "elapsed_seconds": time.perf_counter() - started,
                }
            )
    return cases


def main() -> int:
    args = _parse_args()
    config: MolGRConfig | None = None
    input_cases = _collect_input_cases(args)
    cases = trace_reconstruction_cases(
        input_cases,
        score_all_candidates=not args.no_score_all_candidates,
        config=config,
    )

    output = {
        "input": {
            "source": "generic",
            "ids": [case.id for case in input_cases],
            "total_charge": args.total_charge,
            "total_radical_electrons": _direct_total_radicals(args),
        },
        "case_count": len(cases),
        "cases": cases,
    }
    output_format = _resolve_output_format(args)
    if output_format == "json":
        output_text = json.dumps(
            _jsonable(output),
            ensure_ascii=False,
            indent=None if args.indent == 0 else args.indent,
            allow_nan=False,
        )
    else:
        output_text = _render_markdown_report(_jsonable(output))

    if args.out is None:
        print(output_text, end="" if output_text.endswith("\n") else "\n")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(
            output_text if output_text.endswith("\n") else output_text + "\n", encoding="utf-8"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
