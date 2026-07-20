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
import html
import json
import math
import re
import sys
import threading
import time
from contextlib import contextmanager
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, cast


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem, RDLogger

from molgr.config import MolGRConfig
from molgr.fallback.pipeline import reconstruct_without_metals
from molgr.fallback.stages.preprocess import validate_omol
from molgr.fallback.state import (
    MetalCandidateState,
    OmolStateMachine,
    ReconstructionState,
)
from molgr.fallback.utils.dataclasses import MetalAtomPosition
from molgr.fallback.utils.metals import preparation, scoring, search
from molgr.fallback.utils.no_metals import preparation as no_metal_preparation
from molgr.fallback.utils.no_metals import selection as no_metal_selection
from molgr.utils.converter import pybel_to_rdmol
from molgr.utils.equivalence import check_equivalence


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
    fixture_kind: str = ""
    fixture_structure_file: str = ""
    expected_smiles: str = ""


@dataclasses.dataclass
class DofRenderContext:
    """Runtime state for optional rdkit-dof image rendering."""

    image_dir: Path
    display_base_dir: Path | None
    image_format: str = "svg"
    max_images: int = 1000
    image_size: tuple[int, int] = (360, 300)
    grid_sub_img_size: tuple[int, int] = (320, 260)
    grid_mols_per_row: int = 3
    grid_max_mols: int = 24
    image_count: int = 0
    skipped_count: int = 0
    errors: list[dict[str, Any]] = dataclasses.field(default_factory=list)

    @property
    def use_svg(self) -> bool:
        return self.image_format == "svg"


_OMOL_PHASE_OBSERVER_LOCK = threading.RLock()


@contextmanager
def _observe_omol_phases(observer: Any) -> Iterator[None]:
    """Temporarily observe state-machine phases for this trace request only."""

    original_methods = {
        name: getattr(OmolStateMachine, name)
        for name in ("run_omol_stage", "run_omol_charge_stage", "set_given_charge", "annotate")
    }

    def wrap(name: str) -> Any:
        original = original_methods[name]

        def wrapped(machine: Any, *args: Any, **kwargs: Any) -> Any:
            result = original(machine, *args, **kwargs)
            phase = kwargs.get("phase") if "phase" in kwargs else (args[0] if args else None)
            if phase is not None:
                observer(machine, str(phase))
            return result

        return wrapped

    with _OMOL_PHASE_OBSERVER_LOCK:
        for name in original_methods:
            setattr(OmolStateMachine, name, wrap(name))
        try:
            yield
        finally:
            for name, original in original_methods.items():
                setattr(OmolStateMachine, name, original)


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
    "metal_discordance_aromatic_stability_deficit",
    "metal_discordance_max_aromatic_stability_score",
    "metal_discordance_aromatic_ring_deficit_count",
    "metal_discordance_max_aromatic_ring_count",
    "organic_aromatic_stability_score",
    "organic_aromatic_ring_count",
    "organic_aromatic_atom_count",
    "organic_conjugated_atom_count",
    "organic_conjugated_bond_count",
    "organic_max_conjugated_component_size",
    "metal_discordance_inner_visible_diradical_count",
    "metal_discordance_outer_or_invisible_adjacent_double_charge_count",
    "metal_discordance_outer_or_invisible_adjacent_same_sign_double_charge_count",
    "metal_discordance_outer_or_invisible_adjacent_opposite_sign_double_charge_count",
    "metal_discordance_outer_or_invisible_adjacent_unknown_metal_sign_double_charge_count",
    "metal_discordance_inner_visible_adjacent_carbanion_pair_count",
    "metal_discordance_inner_visible_conjugated_carbanion_pair_count",
    "metal_discordance_inner_visible_same_sign_charge_count",
    "metal_discordance_negative_metal_count",
    "metal_discordance_negative_metal_penalty",
    "metal_discordance_zero_valent_metals_with_organic_cation_count",
    "metal_discordance_nonnegative_metal_unsaturated_organic_cation_count",
    "metal_discordance_negative_metal_outer_sphere_cation_exception",
    "metal_discordance_negative_metal_positive_metal_counterion_exception",
    "passes_metal_discordance_filter",
    "score",
    "force_field_energy",
    "organic_charge_localization_penalty",
    "organic_radical_localization_penalty",
    "selection_key",
    "analysis_selection_key_all_candidates",
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
        "--review-fixtures-manifest",
        type=Path,
        default=None,
        help=(
            "Load trace inputs directly from a reviewed fixture manifest. SDF coordinates and "
            "XYZ electronic states are resolved relative to the manifest."
        ),
    )
    parser.add_argument(
        "--fixture-id",
        dest="fixture_ids",
        action="append",
        default=[],
        help=(
            "Optional reviewed fixture id. May be repeated or comma-separated. "
            "Defaults to all manifest fixtures."
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
        choices=("auto", "json", "html"),
        default="auto",
        help=(
            "Output format. In auto mode, .json writes JSON; every other output path and stdout "
            "write an HTML report."
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
    parser.add_argument(
        "--render-dof-images",
        action="store_true",
        help=(
            "Render reconstructed molecule graphs with rdkit-dof. SVG output is embedded in "
            "the report data; non-SVG output is written as image files."
        ),
    )
    parser.add_argument(
        "--dof-image-dir",
        type=Path,
        default=None,
        help=(
            "Directory for rdkit-dof images. Defaults to <report>.dof-images for --out, or "
            "molgr_trace_dof_images when writing to stdout."
        ),
    )
    parser.add_argument(
        "--dof-image-format",
        choices=("svg", "png"),
        default="svg",
        help="rdkit-dof image format. Default: svg.",
    )
    parser.add_argument(
        "--dof-max-images",
        type=int,
        default=1000,
        help="Maximum number of individual rdkit-dof images to write. Default: 1000.",
    )
    parser.add_argument(
        "--dof-image-size",
        default="360x300",
        help="Single-molecule rdkit-dof image size as WIDTHxHEIGHT. Default: 360x300.",
    )
    parser.add_argument(
        "--dof-grid-max-mols",
        type=int,
        default=24,
        help="Maximum molecules per comparison grid image. Default: 24.",
    )
    parser.add_argument(
        "--dof-grid-mols-per-row",
        type=int,
        default=3,
        help="Molecules per row for rdkit-dof comparison grids. Default: 3.",
    )
    parser.add_argument(
        "--dof-grid-sub-img-size",
        default="320x260",
        help="Grid sub-image size as WIDTHxHEIGHT. Default: 320x260.",
    )
    args = parser.parse_args()
    if args.total_radical_electrons is not None and args.total_radical_electrons < 0:
        parser.error("--total-radical-electrons must be >= 0")
    if args.spin_multiplicity is not None and args.spin_multiplicity < 1:
        parser.error("--spin-multiplicity must be >= 1")
    if args.dof_max_images < 0:
        parser.error("--dof-max-images must be >= 0")
    if args.dof_grid_max_mols < 0:
        parser.error("--dof-grid-max-mols must be >= 0")
    if args.dof_grid_mols_per_row < 1:
        parser.error("--dof-grid-mols-per-row must be >= 1")
    try:
        args.dof_image_size = _parse_size(args.dof_image_size)
        args.dof_grid_sub_img_size = _parse_size(args.dof_grid_sub_img_size)
    except ValueError as exc:
        parser.error(str(exc))
    return args


def _parse_size(raw_size: str) -> tuple[int, int]:
    match = re.fullmatch(r"\s*(\d+)\s*[xX,]\s*(\d+)\s*", raw_size)
    if match is None:
        raise ValueError(f"invalid size {raw_size!r}; expected WIDTHxHEIGHT")
    width, height = int(match.group(1)), int(match.group(2))
    if width <= 0 or height <= 0:
        raise ValueError("image dimensions must be positive")
    return width, height


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


def _default_dof_image_dir(args: argparse.Namespace) -> Path:
    if args.dof_image_dir is not None:
        return args.dof_image_dir
    if args.out is not None:
        return args.out.with_suffix(args.out.suffix + ".dof-images")
    return Path("molgr_trace_dof_images")


def _make_dof_render_context(args: argparse.Namespace) -> DofRenderContext | None:
    if not args.render_dof_images and _resolve_output_format(args) != "html":
        return None
    if _resolve_output_format(args) == "html":
        args.dof_image_format = "svg"
    return DofRenderContext(
        image_dir=_default_dof_image_dir(args),
        display_base_dir=args.out.parent if args.out is not None else None,
        image_format=args.dof_image_format,
        max_images=int(args.dof_max_images),
        image_size=args.dof_image_size,
        grid_sub_img_size=args.dof_grid_sub_img_size,
        grid_mols_per_row=int(args.dof_grid_mols_per_row),
        grid_max_mols=int(args.dof_grid_max_mols),
    )


def dof_rendering_summary(render_context: DofRenderContext) -> dict[str, Any]:
    storage = "embedded" if render_context.use_svg else "files"
    return {
        "image_dir": "" if render_context.use_svg else str(render_context.image_dir),
        "storage": storage,
        "format": render_context.image_format,
        "image_count": render_context.image_count,
        "skipped_count": render_context.skipped_count,
        "max_images": render_context.max_images,
        "errors": render_context.errors,
    }


def _safe_filename_part(value: Any) -> str:
    text = str(value).strip().replace("/", "_")
    text = re.sub(r"[^A-Za-z0-9_.+-]+", "_", text)
    text = text.strip("._")
    return text or "item"


def _display_image_path(path: Path, render_context: DofRenderContext) -> str:
    if render_context.display_base_dir is None:
        return str(path)
    try:
        return str(path.relative_to(render_context.display_base_dir))
    except ValueError:
        return str(path)


def _dof_image_record(
    path: Path,
    *,
    render_context: DofRenderContext,
    label: str,
    kind: str,
) -> dict[str, Any]:
    return {
        "kind": kind,
        "label": label,
        "path": str(path),
        "display_path": _display_image_path(path, render_context),
        "format": render_context.image_format,
    }


def _inline_dof_svg_record(
    svg_fragment: str,
    *,
    render_context: DofRenderContext,
    label: str,
    kind: str,
) -> dict[str, Any]:
    return {
        "kind": kind,
        "label": label,
        "format": render_context.image_format,
        "svg_fragment": svg_fragment,
    }


def _svg_fragment_from_dof_image(image: Any) -> str:
    if isinstance(image, str):
        svg_text = image
    elif hasattr(image, "data"):
        svg_text = str(image.data)
    elif hasattr(image, "_repr_svg_"):
        svg_text = str(image._repr_svg_())
    else:
        svg_text = str(image)
    svg_start = svg_text.find("<svg")
    if svg_start >= 0:
        svg_text = svg_text[svg_start:]
    return svg_text


def _reserve_dof_render_slot(
    *,
    render_context: DofRenderContext,
    label: str,
    kind: str,
    animation: bool = False,
) -> int | dict[str, Any]:
    if render_context.image_count < render_context.max_images:
        index = render_context.image_count
        render_context.image_count += 1
        return index
    render_context.skipped_count += 1
    skipped = {
        "kind": kind,
        "label": label,
        "status": "skipped_limit",
        "max_images": render_context.max_images,
    }
    if animation:
        skipped["animation"] = True
    return skipped


def _render_dof_molecule(
    omol: pybel.Molecule,
    *,
    render_context: DofRenderContext | None,
    case_id: str,
    label: str,
    kind: str,
) -> dict[str, Any] | None:
    if render_context is None:
        return None
    slot = _reserve_dof_render_slot(render_context=render_context, label=label, kind=kind)
    if isinstance(slot, dict):
        return slot

    file_stem = (
        f"{_safe_filename_part(case_id)}__{slot:04d}"
        f"__{_safe_filename_part(kind)}__{_safe_filename_part(label)}"
    )
    try:
        from rdkit_dof import MolToDofImage

        rdmol = pybel_to_rdmol(omol, sanitize=False, kekulize=False)
        if render_context.use_svg:
            image = MolToDofImage(
                rdmol,
                size=render_context.image_size,
                legend=label,
                use_svg=True,
                return_image=False,
            )
            return _inline_dof_svg_record(
                _svg_fragment_from_dof_image(image),
                render_context=render_context,
                label=label,
                kind=kind,
            )
        path = render_context.image_dir / f"{file_stem}.{render_context.image_format}"
        path.parent.mkdir(parents=True, exist_ok=True)
        MolToDofImage(
            rdmol,
            size=render_context.image_size,
            legend=label,
            use_svg=False,
            return_image=False,
            filename=str(path),
        )
        return _dof_image_record(
            path,
            render_context=render_context,
            label=label,
            kind=kind,
        )
    except Exception as exc:
        error = {
            "kind": kind,
            "label": label,
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
        }
        render_context.errors.append(error)
        return error


def _render_dof_grid(
    items: Iterable[tuple[pybel.Molecule, str]],
    *,
    render_context: DofRenderContext | None,
    case_id: str,
    label: str,
    kind: str,
) -> dict[str, Any] | None:
    if render_context is None or render_context.grid_max_mols == 0:
        return None
    slot = _reserve_dof_render_slot(render_context=render_context, label=label, kind=kind)
    if isinstance(slot, dict):
        return slot
    selected_items = list(items)[: render_context.grid_max_mols]
    if not selected_items:
        render_context.image_count -= 1
        return None

    file_stem = (
        f"{_safe_filename_part(case_id)}__{slot:04d}__grid__{_safe_filename_part(kind)}"
        f"__{_safe_filename_part(label)}"
    )
    try:
        from rdkit_dof import MolsToGridDofImage

        mols = [
            pybel_to_rdmol(omol, sanitize=False, kekulize=False) for omol, _legend in selected_items
        ]
        legends = [legend for _omol, legend in selected_items]
        if render_context.use_svg:
            image = MolsToGridDofImage(
                mols,
                molsPerRow=render_context.grid_mols_per_row,
                subImgSize=render_context.grid_sub_img_size,
                legends=legends,
                use_svg=True,
                return_image=False,
            )
            return _inline_dof_svg_record(
                _svg_fragment_from_dof_image(image),
                render_context=render_context,
                label=label,
                kind=kind,
            )
        path = render_context.image_dir / f"{file_stem}.{render_context.image_format}"
        path.parent.mkdir(parents=True, exist_ok=True)
        MolsToGridDofImage(
            mols,
            molsPerRow=render_context.grid_mols_per_row,
            subImgSize=render_context.grid_sub_img_size,
            legends=legends,
            use_svg=False,
            return_image=False,
            filename=str(path),
        )
        return _dof_image_record(
            path,
            render_context=render_context,
            label=label,
            kind=kind,
        )
    except Exception as exc:
        error = {
            "kind": kind,
            "label": label,
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
        }
        render_context.errors.append(error)
        return error


def _copy_omol(omol: pybel.Molecule) -> pybel.Molecule:
    return pybel.Molecule(ob.OBMol(cast(ob.OBMol, omol.OBMol)))


def _render_dof_animation(
    items: Iterable[tuple[pybel.Molecule, str]],
    *,
    render_context: DofRenderContext | None,
    case_id: str,
    label: str,
    kind: str,
    duration: int = 650,
) -> dict[str, Any] | None:
    if render_context is None:
        return None
    slot = _reserve_dof_render_slot(
        render_context=render_context,
        label=label,
        kind=kind,
        animation=True,
    )
    if isinstance(slot, dict):
        return slot
    selected_items = list(items)
    if not selected_items:
        render_context.image_count -= 1
        return None

    try:
        from rdkit_dof import MolsToDofSvgAnimation

        mols = [
            pybel_to_rdmol(omol, sanitize=False, kekulize=False) for omol, _legend in selected_items
        ]
        legends = [legend for _omol, legend in selected_items]
        image = MolsToDofSvgAnimation(
            mols,
            size=render_context.image_size,
            legends=legends,
            duration=duration,
            loop=0,
            return_image=False,
        )
        record = _inline_dof_svg_record(
            _svg_fragment_from_dof_image(image),
            render_context=render_context,
            label=label,
            kind=kind,
        )
        record.update({"animation": True, "frame_count": len(selected_items)})
        return record
    except Exception as exc:
        error = {
            "kind": kind,
            "label": label,
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
            "animation": True,
        }
        render_context.errors.append(error)
        return error


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


_APPROVED_REVIEW_FIXTURE_KINDS = {
    "approved_graph",
    "manual_reference",
    "reference_graph",
}


def _review_fixture_xyz_block(path: Path) -> str:
    if path.suffix.lower() != ".sdf":
        return path.read_text(encoding="utf-8")
    expected = next(
        (
            mol
            for mol in Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
            if mol is not None
        ),
        None,
    )
    if expected is None:
        raise ValueError(f"review fixture SDF is unreadable: {path}")
    if expected.GetNumConformers() != 1:
        raise ValueError(f"review fixture SDF must contain one conformer: {path}")
    return Chem.MolToXYZBlock(expected)


def load_review_fixture_cases(
    manifest_path: Path,
    fixture_ids: Sequence[str] = (),
) -> list[TraceInputCase]:
    """Load trace inputs from the current reviewed fixture manifest."""

    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    raw_records = payload.get("fixtures") if isinstance(payload, dict) else None
    if not isinstance(raw_records, list):
        raise ValueError(f"invalid reviewed fixture manifest: {manifest_path}")

    requested_ids = set(fixture_ids)
    found_ids: set[str] = set()
    cases: list[TraceInputCase] = []
    for raw_record in raw_records:
        if not isinstance(raw_record, dict):
            raise ValueError(f"invalid reviewed fixture record: {raw_record!r}")
        case_id = str(raw_record.get("case_id") or "").strip()
        if not case_id:
            raise ValueError("reviewed fixture record is missing case_id")
        if requested_ids and case_id not in requested_ids:
            continue
        if case_id in found_ids:
            raise ValueError(f"duplicate reviewed fixture id: {case_id}")
        found_ids.add(case_id)

        structure_file = str(raw_record.get("structure_file") or "").strip()
        if not structure_file:
            raise ValueError(f"reviewed fixture {case_id} is missing structure_file")
        structure_path = manifest_path.parent / structure_file
        if not structure_path.is_file():
            raise FileNotFoundError(structure_path)
        fixture_kind = str(raw_record.get("kind") or "").strip()
        expected_smiles = (
            str(raw_record.get("approved_smiles") or "").strip()
            if fixture_kind in _APPROVED_REVIEW_FIXTURE_KINDS
            else ""
        )
        cases.append(
            TraceInputCase(
                id=case_id,
                xyz_block=_review_fixture_xyz_block(structure_path),
                total_charge=int(raw_record.get("total_charge") or 0),
                total_radical_electrons=int(raw_record.get("total_radical_electrons") or 0),
                xyz_path=structure_path,
                xyz_source="review_fixture",
                fixture_kind=fixture_kind,
                fixture_structure_file=structure_file,
                expected_smiles=expected_smiles,
            )
        )

    missing_ids = requested_ids - found_ids
    if missing_ids:
        raise ValueError("reviewed fixture ids not found: " + ", ".join(sorted(missing_ids)))
    if not cases:
        raise ValueError(f"reviewed fixture manifest selected no cases: {manifest_path}")
    return cases


def _collect_input_cases(args: argparse.Namespace) -> list[TraceInputCase]:
    cases: list[TraceInputCase] = []
    if args.case_json is not None:
        cases.extend(_load_json_cases(args.case_json))
    fixture_ids = split_repeated_values(args.fixture_ids)
    if args.review_fixtures_manifest is not None:
        cases.extend(load_review_fixture_cases(args.review_fixtures_manifest, fixture_ids))
    elif fixture_ids:
        raise ValueError("--fixture-id requires --review-fixtures-manifest")

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
        raise ValueError(
            "provide at least one of --xyz, --xyz-block, --stdin, --case-json, "
            "or --review-fixtures-manifest"
        )
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
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _jsonable(dataclasses.asdict(cast(Any, value)))
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


def _is_rendered_dof_image(image: Any) -> bool:
    return (
        isinstance(image, dict)
        and not image.get("status")
        and bool(image.get("svg_fragment") or image.get("display_path") or image.get("path"))
    )


def _dof_image_path_text(image: Any) -> str:
    if not isinstance(image, dict):
        return ""
    if image.get("svg_fragment"):
        return "embedded svg"
    return image.get("display_path") or image.get("path") or image.get("status", "")


def _resolve_output_format(args: argparse.Namespace) -> str:
    if args.format != "auto":
        return cast(str, args.format)
    if args.out is not None and args.out.suffix.lower() == ".json":
        return "json"
    return "html"


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
    return "; ".join(_metal_state_label(metal_state) for metal_state in metal_states)


def _html_escape(value: Any) -> str:
    return html.escape(_format_value(value), quote=True)


def _html_json(value: Any) -> str:
    return html.escape(
        json.dumps(_jsonable(value), ensure_ascii=False, allow_nan=False), quote=False
    )


def _html_script_json(value: Any) -> str:
    return (
        json.dumps(_jsonable(value), ensure_ascii=False, allow_nan=False)
        .replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
    )


def _read_svg_fragment(image: Any) -> str:
    if not _is_rendered_dof_image(image):
        return ""
    if isinstance(image.get("svg_fragment"), str):
        return cast(str, image["svg_fragment"])
    path = Path(str(image.get("path") or image.get("display_path")))
    try:
        svg_text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        svg_text = path.read_text(encoding="iso-8859-1")
    except OSError:
        return ""
    svg_start = svg_text.find("<svg")
    if svg_start >= 0:
        svg_text = svg_text[svg_start:]
    return svg_text


def _with_inline_dof_svgs(value: Any) -> Any:
    if isinstance(value, dict):
        copied = {key: _with_inline_dof_svgs(item) for key, item in value.items()}
        if _is_rendered_dof_image(copied) and "svg_fragment" not in copied:
            svg_fragment = _read_svg_fragment(copied)
            if svg_fragment:
                copied["svg_fragment"] = svg_fragment
        return copied
    if isinstance(value, list):
        return [_with_inline_dof_svgs(item) for item in value]
    return value


def _html_metric_grid(items: Sequence[tuple[str, Any]]) -> str:
    cells = []
    for label, value in items:
        if value in ("", None, [], {}):
            continue
        cells.append(
            '<div class="metric">'
            f'<span class="metric-label">{_html_escape(label)}</span>'
            f'<span class="metric-value">{_html_escape(value)}</span>'
            "</div>"
        )
    return '<div class="metrics">' + "".join(cells) + "</div>" if cells else ""


def _html_table(
    headers: Sequence[str], rows: Sequence[Sequence[Any]], *, class_name: str = ""
) -> str:
    if not rows:
        return '<p class="empty-inline">无</p>'
    class_attr = f' class="{html.escape(class_name)}"' if class_name else ""
    head = "".join(f"<th>{_html_escape(header)}</th>" for header in headers)
    body_rows = []
    for row in rows:
        body_rows.append(
            "<tr>" + "".join(f"<td>{_html_escape(cell)}</td>" for cell in row) + "</tr>"
        )
    return f"<table{class_attr}><thead><tr>{head}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def _html_kv_table(items: Sequence[tuple[str, Any]], *, class_name: str = "") -> str:
    rows = [(key, value) for key, value in items if value not in ("", None, [], {})]
    return _html_table(("字段", "值"), rows, class_name=class_name)


def _html_details(title: str, body: str, *, open_: bool = False, class_name: str = "") -> str:
    class_attr = f' class="{html.escape(class_name)}"' if class_name else ""
    open_attr = " open" if open_ else ""
    return (
        f"<details{class_attr}{open_attr}><summary>{_html_escape(title)}</summary>{body}</details>"
    )


def _html_json_block(value: Any, *, class_name: str = "json-block") -> str:
    return f'<pre class="{html.escape(class_name)}">{_html_json(value)}</pre>'


def _html_score_table(mapping: Any) -> str:
    if not isinstance(mapping, dict) or not mapping:
        return '<p class="empty-inline">无</p>'
    keys = [key for key in _IMPORTANT_SCORE_KEYS if key in mapping]
    keys.extend(sorted(key for key in mapping if key not in set(keys)))
    return _html_table(
        ("分数项", "值"), [(key, mapping[key]) for key in keys], class_name="score-table"
    )


def _candidate_title(candidate: dict[str, Any]) -> str:
    if candidate.get("trace_item_title"):
        return str(candidate["trace_item_title"])
    if "search_layer_index" in candidate or "combination_index" in candidate:
        selected = "selected " if candidate.get("selected") else ""
        return (
            f"{selected}L{candidate.get('search_layer_index', '')}/"
            f"C{candidate.get('combination_index', '')}"
        ).strip()
    if "resonance_index" in candidate:
        return f"resonance {candidate.get('resonance_index', '')}"
    return str(candidate.get("label") or candidate.get("kind") or "item")


def _collect_html_image_items(output: dict[str, Any]) -> list[dict[str, Any]]:
    items: list[dict[str, Any]] = []

    def add_item(
        *,
        case: dict[str, Any],
        title: str,
        image: Any,
        group: str,
        selected: bool = False,
        metadata: dict[str, Any] | None = None,
    ) -> None:
        if not isinstance(image, dict):
            return
        svg = _read_svg_fragment(image)
        if not svg:
            return
        item_metadata = {} if metadata is None else metadata
        items.append(
            {
                "id": f"item-{len(items)}",
                "case_id": case.get("id", ""),
                "title": title,
                "group": group,
                "selected": selected,
                "svg": svg,
                "image": image,
                "metadata": _jsonable(item_metadata),
                "score": item_metadata.get("score", ""),
                "layer": item_metadata.get("search_layer_index", ""),
                "combination": item_metadata.get("combination_index", ""),
                "smiles": item_metadata.get("canonical_smiles")
                or item_metadata.get("smiles")
                or item_metadata.get("organic_smiles", ""),
            }
        )

    for raw_case in cast(List[Any], output.get("cases", [])):
        if not isinstance(raw_case, dict):
            continue
        case = raw_case
        no_metal_trace = cast(Dict[str, Any], case.get("no_metal_trace", {}))
        selected = cast(Dict[str, Any], no_metal_trace.get("selected_candidate", {}))
        selected_state = cast(Dict[str, Any], selected.get("state", {}))
        add_item(
            case=case,
            title="selected no-metal candidate",
            image=selected_state.get("dof_image"),
            group="selected",
            selected=True,
            metadata=selected,
        )
        add_item(
            case=case,
            title="metal candidate grid",
            image=case.get("dof_candidate_grid"),
            group="grid",
            metadata={
                "candidate_count": case.get("candidate_count", ""),
                "production_candidate_count": case.get("production_candidate_count", ""),
                "selected_candidate": case.get("selected_candidate"),
            },
        )
        for candidate in cast(List[Any], case.get("candidates", [])):
            if not isinstance(candidate, dict):
                continue
            add_item(
                case=case,
                title=_candidate_title(candidate),
                image=candidate.get("dof_image"),
                group="candidate",
                selected=bool(candidate.get("selected")),
                metadata=candidate,
            )
        search_summary = cast(Dict[str, Any], case.get("search", {}))
        for layer in cast(List[Any], search_summary.get("layer_summaries", [])):
            if not isinstance(layer, dict):
                continue
            for bucket in cast(List[Any], layer.get("target_buckets", [])):
                if not isinstance(bucket, dict):
                    continue
                target = cast(Dict[str, Any], bucket.get("target", {}))
                add_item(
                    case=case,
                    title=(
                        f"organic target L{layer.get('layer_index', '')} "
                        f"Q={target.get('no_metal_charge', '')} "
                        f"R={target.get('no_metal_radical_electrons', '')}"
                    ),
                    image=bucket.get("dof_image"),
                    group="organic",
                    metadata=bucket,
                )
    return items


def _html_no_metal_trace(no_metal_trace: dict[str, Any]) -> str:
    if not no_metal_trace:
        return '<p class="empty-inline">无</p>'

    target = cast(Dict[str, Any], no_metal_trace.get("target", {}))
    selected = cast(Dict[str, Any], no_metal_trace.get("selected_candidate", {}))
    selected_state = cast(Dict[str, Any], selected.get("state", {}))
    resonance = cast(Dict[str, Any], no_metal_trace.get("resonance", {}))
    actions = cast(List[Any], no_metal_trace.get("neighbor_radical_actions", []))
    sections = [
        _html_kv_table(
            (
                ("状态", no_metal_trace.get("status", "")),
                ("目标电荷", target.get("total_charge", "")),
                ("目标自由基电子", target.get("total_radical_electrons", "")),
                ("邻位自由基动作", actions),
                ("恢复层", no_metal_trace.get("recovery_tier", 0)),
                ("共振种子序号", resonance.get("seed_index", "")),
                ("共振序号", resonance.get("resonance_index", "")),
                ("raw 共振序号", resonance.get("raw_index", "")),
                ("归一化方式", resonance.get("normalization", "")),
                ("选中 canonical SMILES", selected_state.get("canonical_smiles", "")),
                ("选中分数", selected.get("score", "")),
                ("选中选择键", selected.get("organic_topology_selection_key", "")),
            )
        )
    ]

    step_rows = [
        (step.get("step_index", ""), step.get("phase", ""), step.get("kind", ""))
        for step in cast(List[Any], no_metal_trace.get("pipeline_steps", []))
        if isinstance(step, dict)
    ]
    sections.append(
        _html_details(
            f"生产管线阶段 ({len(step_rows)})",
            _html_table(("步骤", "阶段", "类型"), step_rows),
            open_=True,
        )
    )
    if selected:
        sections.append(_html_details("选中候选完整 JSON", _html_json_block(selected)))
    if resonance:
        sections.append(_html_details("共振选择元数据", _html_json_block(resonance)))
    sections.append(_html_details("完整 no-metal trace JSON", _html_json_block(no_metal_trace)))
    return "".join(sections)


def _html_candidate_details(candidate: dict[str, Any]) -> str:
    organic_part = cast(Dict[str, Any], candidate.get("organic_part", {}))
    metal_rows = []
    for metal_state in cast(List[Any], candidate.get("metal_states", [])):
        if isinstance(metal_state, dict):
            metal_rows.append(
                (
                    metal_state.get("idx", ""),
                    metal_state.get("symbol", ""),
                    metal_state.get("valence", ""),
                    metal_state.get("radical_num", ""),
                    metal_state.get("position", ""),
                )
            )
    body = _html_kv_table(
        (
            ("候选序号", candidate.get("candidate_index", "")),
            ("候选 identity", candidate.get("candidate_identity", "")),
            ("选中", candidate.get("selected", "")),
            ("层", candidate.get("search_layer_index", "")),
            ("组合", candidate.get("combination_index", "")),
            ("生产层", candidate.get("in_production_selection_layer", "")),
            ("候选总电荷", candidate.get("candidate_total_charge", "")),
            ("目标", candidate.get("target", "")),
            ("有机 canonical SMILES", organic_part.get("canonical_smiles", "")),
            ("分数", candidate.get("score", "")),
            ("DOF 图像", _dof_image_path_text(candidate.get("dof_image"))),
        )
    )
    body += "<h4>金属状态</h4>" + _html_table(
        ("原子序号", "元素", "价态", "自由基数", "坐标"), metal_rows
    )
    body += "<h4>分数构成</h4>" + _html_score_table(candidate.get("score_details", {}))
    production_metadata = cast(Dict[str, Any], candidate.get("production_metadata", {}))
    if production_metadata:
        body += _html_details(
            "生产选择 metadata",
            _html_score_table(production_metadata)
            + _html_details("完整 production metadata JSON", _html_json_block(production_metadata)),
        )
    body += _html_details(
        "score_details JSON", _html_json_block(candidate.get("score_details", {}))
    )
    body += _html_details("organic_part JSON", _html_json_block(candidate.get("organic_part", {})))
    body += _html_details("target JSON", _html_json_block(candidate.get("target", {})))
    body += _html_details(
        "phase_history JSON", _html_json_block(candidate.get("phase_history", []))
    )
    body += _html_details(
        "完整 candidate metadata", _html_json_block(candidate.get("metadata", {}))
    )
    body += _html_details(
        "no_metal_state",
        _html_json_block(candidate.get("no_metal_state", {})),
    )
    body += _html_details("完整 candidate JSON", _html_json_block(candidate))
    return body


def _html_case_trace(case: dict[str, Any]) -> str:
    case_id = str(case.get("id", "unknown"))
    if case.get("status") == "error":
        return (
            f'<section class="case-trace" id="case-{html.escape(_safe_filename_part(case_id))}">'
            f"<h2>{_html_escape(case_id)}</h2>"
            + _html_kv_table((("状态", "error"), ("错误", case.get("error", ""))))
            + "</section>"
        )
    base_state = cast(Dict[str, Any], case.get("base_state", {}))
    selected = cast(Dict[str, Any], case.get("selected_candidate", {}))
    review_fixture = cast(Dict[str, Any], case.get("review_fixture", {}))
    summary = _html_kv_table(
        (
            ("id", case.get("id", "")),
            ("CSV 行号", case.get("row_index", "")),
            ("状态", case.get("status", "")),
            ("重建类型", case.get("trace_kind", "")),
            ("总电荷", case.get("charge", "")),
            ("总自由基电子数", case.get("total_radical_electrons", "")),
            ("自旋多重度", case.get("spin_multiplicity", "")),
            ("自旋来源", case.get("spin_source", "")),
            ("XYZ 路径", case.get("xyz_path", "")),
            ("XYZ 来源", case.get("xyz_source", "")),
            ("review fixture 类型", review_fixture.get("kind", "")),
            ("review fixture 等价", review_fixture.get("equivalent", "")),
            ("review fixture 比较原因", review_fixture.get("equivalence_reason", "")),
            ("参考 SMILES", case.get("reference_smiles", "")),
            ("金属原子数", base_state.get("metal_atom_count", "")),
            (
                "生产选择层",
                cast(Dict[str, Any], case.get("search", {})).get("selected_layer_index", ""),
            ),
            ("生产候选数", case.get("production_candidate_count", "")),
            ("全部已评分候选数", case.get("candidate_count", "")),
            ("选中候选", selected.get("combination_index", "")),
            ("耗时秒", case.get("elapsed_seconds", "")),
        )
    )
    sections = [
        f'<section class="case-trace" id="case-{html.escape(_safe_filename_part(case_id))}">',
        f"<h2>{_html_escape(case_id)}</h2>",
        _html_details("基本信息", summary, open_=True),
        _html_details("基础状态完整 JSON", _html_json_block(base_state)),
    ]
    if review_fixture:
        sections.append(
            _html_details("review fixture 同步检查", _html_json_block(review_fixture), open_=True)
        )
    if case.get("trace_kind") == "no_metal":
        sections.append(
            _html_details(
                "无金属重建过程",
                _html_no_metal_trace(cast(Dict[str, Any], case.get("no_metal_trace", {}))),
                open_=True,
            )
        )
    else:
        state_rows = []
        for site_index, state_options in enumerate(
            cast(List[Any], base_state.get("available_metal_states_by_site", []))
        ):
            if not isinstance(state_options, list):
                continue
            for option_index, metal_state in enumerate(state_options):
                if isinstance(metal_state, dict):
                    state_rows.append(
                        (
                            site_index,
                            option_index,
                            metal_state.get("idx", ""),
                            metal_state.get("symbol", ""),
                            metal_state.get("valence", ""),
                            metal_state.get("radical_num", ""),
                            metal_state.get("position", ""),
                        )
                    )
        sections.append(
            _html_details(
                "金属位点候选",
                _html_table(
                    ("位点", "选项", "原子序号", "元素", "价态", "自由基数", "坐标"),
                    state_rows,
                ),
            )
        )
        search_summary = cast(Dict[str, Any], case.get("search", {}))
        layer_rows = []
        layer_detail_blocks = []
        target_rows = []
        no_metal_sections = []
        for layer in cast(List[Any], search_summary.get("layer_summaries", [])):
            if not isinstance(layer, dict):
                continue
            layer_detail_blocks.append(
                _html_details(
                    f"Layer {layer.get('layer_index', '')} 完整 JSON",
                    _html_json_block(layer),
                )
            )
            layer_rows.append(
                (
                    layer.get("layer_index", ""),
                    layer.get("status", ""),
                    layer.get("production_selected_layer", ""),
                    layer.get("state_group_count", ""),
                    layer.get("state_options_per_group", ""),
                    layer.get("target_bucket_count", ""),
                    layer.get("candidate_count", ""),
                    layer.get("prepared_candidate_count", ""),
                    bool(layer.get("analysis_score_context")),
                )
            )
            for bucket in cast(List[Any], layer.get("target_buckets", [])):
                if not isinstance(bucket, dict):
                    continue
                target = cast(Dict[str, Any], bucket.get("target", {}))
                organic_part = cast(Dict[str, Any], bucket.get("organic_part", {}))
                target_rows.append(
                    (
                        layer.get("layer_index", ""),
                        target.get("no_metal_charge", ""),
                        target.get("no_metal_radical_electrons", ""),
                        bucket.get("status", ""),
                        bucket.get("candidate_count", ""),
                        bucket.get("prepared_candidate_count", ""),
                        organic_part.get("canonical_smiles", organic_part.get("smiles", "")),
                        bucket.get("no_metal_score", ""),
                        _dof_image_path_text(bucket.get("dof_image")),
                    )
                )
                if isinstance(bucket.get("no_metal_trace"), dict):
                    no_metal_sections.append(
                        _html_details(
                            (
                                f"L{layer.get('layer_index', '')} "
                                f"Q={target.get('no_metal_charge', '')} "
                                f"R={target.get('no_metal_radical_electrons', '')}"
                            ),
                            _html_no_metal_trace(cast(Dict[str, Any], bucket["no_metal_trace"])),
                        )
                    )
        sections.append(
            _html_details(
                "搜索层",
                _html_table(
                    (
                        "层",
                        "状态",
                        "生产选择层",
                        "金属组数",
                        "每组候选数",
                        "目标桶数",
                        "枚举候选数",
                        "已评分候选数",
                        "analysis score context",
                    ),
                    layer_rows,
                ),
                open_=True,
            )
        )
        if layer_detail_blocks:
            sections.append(_html_details("搜索层完整 trace", "".join(layer_detail_blocks)))
        sections.append(
            _html_details(
                "有机目标桶",
                _html_table(
                    (
                        "层",
                        "有机目标电荷",
                        "有机自由基电子",
                        "状态",
                        "候选数",
                        "已评分数",
                        "有机部分 canonical SMILES",
                        "有机力场分",
                        "DOF 图像",
                    ),
                    target_rows,
                ),
                open_=True,
            )
        )
        sections.append(
            _html_details("分析和评分上下文", _html_json_block(case.get("analysis", {})))
        )
        if no_metal_sections:
            sections.append(_html_details("无金属重建过程", "".join(no_metal_sections)))
        candidate_rows = []
        for candidate in cast(List[Any], case.get("candidates", [])):
            if not isinstance(candidate, dict):
                continue
            organic_part = cast(Dict[str, Any], candidate.get("organic_part", {}))
            candidate_rows.append(
                (
                    "是" if candidate.get("selected") else "",
                    candidate.get("search_layer_index", ""),
                    candidate.get("combination_index", ""),
                    "是" if candidate.get("in_production_selection_layer") else "",
                    _metal_states_label(
                        cast(List[Dict[str, Any]], candidate.get("metal_states", []))
                    ),
                    organic_part.get("canonical_smiles", organic_part.get("smiles", "")),
                    candidate.get("candidate_total_charge", ""),
                    cast(Dict[str, Any], candidate.get("target", {})).get("no_metal_charge", ""),
                    cast(Dict[str, Any], candidate.get("target", {})).get(
                        "no_metal_radical_electrons", ""
                    ),
                    *[
                        _score_detail_value(candidate, score_key, "")
                        for score_key in _IMPORTANT_SCORE_KEYS
                    ],
                    _dof_image_path_text(candidate.get("dof_image")),
                )
            )
        sections.append(
            _html_details(
                "候选总览",
                _html_table(
                    (
                        "选中",
                        "层",
                        "组合",
                        "生产层",
                        "金属状态",
                        "有机 canonical SMILES",
                        "候选总电荷",
                        "有机电荷",
                        "有机自由基",
                        *_IMPORTANT_SCORE_KEYS,
                        "DOF 图像",
                    ),
                    candidate_rows,
                    class_name="score-overview",
                ),
                open_=True,
            )
        )
        detail_blocks = []
        for candidate in cast(List[Any], case.get("candidates", [])):
            if isinstance(candidate, dict):
                detail_blocks.append(
                    _html_details(_candidate_title(candidate), _html_candidate_details(candidate))
                )
        sections.append(_html_details("候选分数和 metadata", "".join(detail_blocks), open_=True))
    sections.append(_html_details("完整 case JSON", _html_json_block(case)))
    sections.append("</section>")
    return "".join(sections)


def _render_html_trace_sections(output: dict[str, Any]) -> str:
    sections = []
    for case in cast(List[Any], output.get("cases", [])):
        if isinstance(case, dict):
            sections.append(_html_case_trace(case))
    return "".join(sections) if sections else '<div class="empty">No cases.</div>'


def _render_html_report(output: dict[str, Any]) -> str:
    items = _collect_html_image_items(output)
    trace_sections = _render_html_trace_sections(output)
    input_summary = cast(Dict[str, Any], output.get("input", {}))
    dof_rendering = cast(Dict[str, Any], output.get("dof_rendering", {}))
    metrics = _html_metric_grid(
        (
            ("source", input_summary.get("source", "")),
            ("cases", output.get("case_count", "")),
            ("ids", ", ".join(str(item) for item in cast(List[Any], input_summary.get("ids", [])))),
            ("storage", dof_rendering.get("storage", "")),
            ("images", dof_rendering.get("image_count", "")),
            ("skipped", dof_rendering.get("skipped_count", "")),
            ("errors", len(cast(List[Any], dof_rendering.get("errors", [])))),
        )
    )
    case_options = sorted({str(item.get("case_id", "")) for item in items if item.get("case_id")})
    group_options = sorted({str(item.get("group", "")) for item in items if item.get("group")})
    nav_cards: list[str] = []
    detail_cards: list[str] = []
    for index, item in enumerate(items):
        item_id = str(item["id"])
        selected_class = " is-selected" if item.get("selected") else ""
        active_class = " is-active" if index == 0 else ""
        metadata_json = _html_json(item.get("metadata", {}))
        nav_cards.append(
            f'<button class="trace-card{selected_class}{active_class}" type="button" '
            f'data-target="{html.escape(item_id)}" data-case="{_html_escape(item.get("case_id", ""))}" '
            f'data-group="{_html_escape(item.get("group", ""))}" '
            f'data-search="{_html_escape(" ".join(str(item.get(k, "")) for k in ("case_id", "title", "group", "smiles", "score")))}">'
            f'<span class="card-title">{_html_escape(item.get("title", ""))}</span>'
            f'<span class="card-meta">{_html_escape(item.get("case_id", ""))} · {_html_escape(item.get("group", ""))}</span>'
            "</button>"
        )
        detail_cards.append(
            f'<article id="{html.escape(item_id)}" class="detail-card{active_class}" '
            f'data-case="{_html_escape(item.get("case_id", ""))}" data-group="{_html_escape(item.get("group", ""))}">'
            '<div class="detail-head">'
            f"<div><h2>{_html_escape(item.get('title', ''))}</h2>"
            f"<p>{_html_escape(item.get('case_id', ''))} · {_html_escape(item.get('group', ''))}</p></div>"
            f'<span class="badge{" selected" if item.get("selected") else ""}">'
            f"{'selected' if item.get('selected') else _html_escape(item.get('group', ''))}</span>"
            "</div>"
            '<div class="svg-wrap">'
            f"{item.get('svg', '')}"
            "</div>"
            f'<pre class="metadata">{metadata_json}</pre>'
            "</article>"
        )

    cases_json = _html_json(output)
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>MolGR Trace</title>
  <style>
    :root {{ color-scheme: light; --bg:#f6f7f9; --panel:#fff; --line:#d8dee8; --text:#172033; --muted:#667085; --accent:#2563eb; --selected:#b45309; }}
    * {{ box-sizing: border-box; }}
    body {{ margin:0; font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; background:var(--bg); color:var(--text); }}
    header {{ padding:16px 20px; border-bottom:1px solid var(--line); background:var(--panel); position:sticky; top:0; z-index:10; }}
    h1 {{ margin:0 0 10px; font-size:20px; }}
    .metrics {{ display:flex; flex-wrap:wrap; gap:8px; }}
    .metric {{ border:1px solid var(--line); border-radius:6px; background:#fafbfc; padding:6px 9px; min-width:96px; }}
    .metric-label {{ display:block; font-size:11px; color:var(--muted); text-transform:uppercase; }}
    .metric-value {{ display:block; font-size:13px; font-weight:700; overflow-wrap:anywhere; }}
    main {{ display:grid; grid-template-columns:340px minmax(0,1fr); gap:14px; padding:14px; min-height:calc(100vh - 86px); }}
    aside {{ min-width:0; border:1px solid var(--line); background:var(--panel); border-radius:8px; padding:10px; align-self:start; position:sticky; top:94px; max-height:calc(100vh - 108px); display:flex; flex-direction:column; gap:10px; }}
    .filters {{ display:grid; gap:8px; }}
    .filters input, .filters select {{ width:100%; border:1px solid var(--line); border-radius:6px; padding:8px; font:inherit; background:#fff; }}
    .toggle {{ display:flex; align-items:center; gap:8px; font-size:13px; color:var(--muted); }}
    .list {{ overflow:auto; display:grid; gap:6px; padding-right:2px; }}
    .trace-card {{ width:100%; border:1px solid var(--line); border-radius:7px; background:#fff; text-align:left; padding:8px; cursor:pointer; }}
    .trace-card:hover, .trace-card.is-active {{ border-color:var(--accent); box-shadow:0 0 0 2px rgba(37,99,235,.14); }}
    .trace-card.is-selected {{ border-color:var(--selected); background:#fff7ed; }}
    .card-title {{ display:block; font-weight:800; font-size:13px; overflow-wrap:anywhere; }}
    .card-meta {{ display:block; margin-top:3px; color:var(--muted); font-size:12px; }}
    .content {{ min-width:0; display:grid; gap:14px; }}
    .section-title {{ margin:4px 0 0; font-size:15px; color:var(--muted); text-transform:uppercase; letter-spacing:.04em; }}
    .detail-card {{ display:none; border:1px solid var(--line); border-radius:8px; background:var(--panel); padding:14px; min-width:0; }}
    .detail-card.is-active {{ display:block; }}
    .detail-head {{ display:flex; justify-content:space-between; gap:10px; align-items:flex-start; margin-bottom:12px; }}
    h2 {{ margin:0; font-size:18px; }}
    .detail-head p {{ margin:4px 0 0; color:var(--muted); }}
    .badge {{ border:1px solid var(--line); border-radius:999px; padding:4px 8px; font-size:12px; color:var(--muted); white-space:nowrap; }}
    .badge.selected {{ color:#92400e; border-color:#f59e0b; background:#fffbeb; }}
    .svg-wrap {{ width:100%; overflow:auto; border:1px solid var(--line); border-radius:8px; background:#fff; padding:10px; }}
    .svg-wrap svg {{ display:block; max-width:100%; height:auto; margin:auto; }}
    .metadata {{ margin:12px 0 0; padding:10px; border:1px solid var(--line); border-radius:8px; background:#0f172a; color:#e5e7eb; overflow:auto; max-height:300px; font-size:12px; }}
    .trace-section-wrap {{ display:grid; gap:14px; }}
    .case-trace {{ border:1px solid var(--line); border-radius:8px; background:var(--panel); padding:14px; min-width:0; }}
    .case-trace > h2 {{ margin-bottom:12px; }}
    details {{ border:1px solid var(--line); border-radius:8px; background:#fff; margin:10px 0; padding:0; overflow:hidden; }}
    details > summary {{ cursor:pointer; padding:9px 11px; font-weight:800; background:#f8fafc; border-bottom:1px solid var(--line); }}
    details:not([open]) > summary {{ border-bottom:0; }}
    details > *:not(summary) {{ margin:10px; max-width:calc(100% - 20px); overflow-x:auto; }}
    h4 {{ margin:14px 10px 6px; font-size:13px; }}
    table {{ width:100%; border-collapse:collapse; font-size:12px; margin:10px 0; }}
    th, td {{ border-bottom:1px solid #e5e7eb; padding:6px 7px; text-align:left; vertical-align:top; overflow-wrap:anywhere; }}
    th {{ color:#334155; background:#f8fafc; position:sticky; top:0; z-index:1; }}
    .score-overview {{ min-width:2400px; }}
    .json-block {{ margin:10px; padding:10px; border:1px solid var(--line); border-radius:8px; background:#0f172a; color:#e5e7eb; overflow:auto; max-height:480px; font-size:12px; }}
    .empty-inline {{ margin:10px; color:var(--muted); }}
    .empty {{ border:1px dashed var(--line); border-radius:8px; padding:24px; text-align:center; color:var(--muted); background:#fff; }}
    @media (max-width:900px) {{ main {{ grid-template-columns:1fr; }} aside {{ position:static; max-height:none; }} }}
  </style>
</head>
<body>
  <header>
    <h1>MolGR Trace</h1>
    {metrics}
  </header>
  <main>
    <aside>
      <div class="filters">
        <input id="search" type="search" placeholder="Filter case, title, group, score...">
        <select id="case-filter"><option value="">All cases</option>{"".join(f'<option value="{html.escape(case_id)}">{html.escape(case_id)}</option>' for case_id in case_options)}</select>
        <select id="group-filter"><option value="">All groups</option>{"".join(f'<option value="{html.escape(group)}">{html.escape(group)}</option>' for group in group_options)}</select>
        <label class="toggle"><input id="selected-only" type="checkbox"> Selected only</label>
      </div>
      <div id="list" class="list">{"".join(nav_cards) if nav_cards else '<div class="empty">No DOF images were rendered.</div>'}</div>
    </aside>
    <section class="content">
      <h2 class="section-title">DOF images</h2>
      {"".join(detail_cards) if detail_cards else '<div class="empty">No DOF images were rendered. Increase --dof-max-images if needed.</div>'}
      <h2 class="section-title">Full trace</h2>
      <div class="trace-section-wrap">{trace_sections}</div>
    </section>
  </main>
  <script id="trace-json" type="application/json">{cases_json}</script>
  <script>
    (() => {{
      const cards = Array.from(document.querySelectorAll(".trace-card"));
      const details = Array.from(document.querySelectorAll(".detail-card"));
      const search = document.getElementById("search");
      const caseFilter = document.getElementById("case-filter");
      const groupFilter = document.getElementById("group-filter");
      const selectedOnly = document.getElementById("selected-only");

      function activate(id) {{
        cards.forEach(card => card.classList.toggle("is-active", card.dataset.target === id));
        details.forEach(detail => detail.classList.toggle("is-active", detail.id === id));
      }}

      function applyFilters() {{
        const q = (search.value || "").toLowerCase();
        const caseValue = caseFilter.value || "";
        const groupValue = groupFilter.value || "";
        let firstVisible = null;
        cards.forEach(card => {{
          const matchesSearch = !q || (card.dataset.search || "").toLowerCase().includes(q);
          const matchesCase = !caseValue || card.dataset.case === caseValue;
          const matchesGroup = !groupValue || card.dataset.group === groupValue;
          const matchesSelected = !selectedOnly.checked || card.classList.contains("is-selected");
          const visible = matchesSearch && matchesCase && matchesGroup && matchesSelected;
          card.hidden = !visible;
          if (visible && !firstVisible) firstVisible = card;
        }});
        if (firstVisible && !cards.some(card => !card.hidden && card.classList.contains("is-active"))) {{
          activate(firstVisible.dataset.target);
        }}
      }}

      cards.forEach(card => card.addEventListener("click", () => activate(card.dataset.target)));
      [search, caseFilter, groupFilter, selectedOnly].forEach(el => el && el.addEventListener("input", applyFilters));
      applyFilters();
    }})();
  </script>
</body>
</html>
"""


def _render_html_browser_report(output: dict[str, Any]) -> str:
    inline_output = _with_inline_dof_svgs(output)
    cases_json = _html_script_json(inline_output)
    return (
        """<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>MolGR Trace</title>
  <style>
    :root { color-scheme: light; --bg:#f6f7f9; --panel:#fff; --line:#d8dee8; --text:#172033; --muted:#667085; --accent:#2563eb; --selected:#b45309; }
    * { box-sizing: border-box; }
    body { margin:0; font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; background:var(--bg); color:var(--text); }
    header { padding:14px 18px; border-bottom:1px solid var(--line); background:var(--panel); position:sticky; top:0; z-index:20; }
    h1 { margin:0; font-size:20px; }
    h2 { margin:0 0 10px; font-size:15px; }
    .global-info { display:grid; grid-template-columns:minmax(260px,.8fr) minmax(220px,.65fr) minmax(420px,1.55fr); gap:12px; padding:12px; }
    .global-card { border:1px solid var(--line); border-radius:8px; background:var(--panel); padding:12px; min-width:0; overflow:auto; }
    main { display:grid; grid-template-columns:380px minmax(0,1fr); gap:12px; padding:0 12px 12px; min-height:calc(100vh - 250px); }
    aside { border:1px solid var(--line); border-radius:8px; background:var(--panel); padding:10px; align-self:start; position:sticky; top:62px; max-height:calc(100vh - 76px); display:flex; flex-direction:column; gap:10px; min-width:0; }
    #tree-search { width:100%; border:1px solid var(--line); border-radius:6px; padding:8px; font:inherit; background:#fff; }
    .tree { overflow:auto; }
    .tree ol { list-style:none; margin:0; padding-left:14px; border-left:1px solid #e5e7eb; }
    .tree > ol { padding-left:0; border-left:0; }
    .tree li { margin:4px 0; }
    .tree-row { display:flex; align-items:flex-start; gap:4px; min-width:0; }
    .tree-toggle, .tree-toggle-spacer { flex:0 0 24px; width:24px; min-width:24px; height:24px; margin-top:2px; }
    .tree-toggle { border:1px solid transparent; border-radius:6px; background:transparent; padding:0; cursor:pointer; color:var(--muted); display:inline-flex; align-items:center; justify-content:center; }
    .tree-toggle::before { content:"▸"; font-size:12px; line-height:1; transition:transform .14s ease; }
    .tree-toggle:hover, .tree-toggle:focus-visible { border-color:var(--line); background:#f8fafc; outline:none; }
    .tree li.is-expanded > .tree-row > .tree-toggle::before { transform:rotate(90deg); }
    .tree-children { list-style:none; margin:4px 0 0 12px; padding-left:14px; border-left:1px solid #e5e7eb; }
    .tree li.is-collapsed > .tree-children { display:none; }
    .tree-node { flex:1 1 auto; width:auto; min-width:0; border:1px solid transparent; border-radius:6px; background:transparent; padding:7px 8px; text-align:left; cursor:pointer; color:var(--text); }
    .tree-node:hover, .tree-node.is-active { border-color:var(--accent); background:#eff6ff; }
    .tree-node.is-selected { border-color:var(--selected); background:#fff7ed; }
    .tree-label { display:block; font-weight:800; font-size:13px; overflow-wrap:anywhere; }
    .tree-meta { display:block; margin-top:2px; color:var(--muted); font-size:11px; }
    .node-panel { display:none; border:1px solid var(--line); border-radius:8px; background:var(--panel); padding:14px; min-width:0; }
    .node-panel.is-active { display:block; }
    .panel-head { display:flex; justify-content:space-between; align-items:flex-start; gap:10px; margin-bottom:12px; }
    .panel-head h2 { margin:0; font-size:18px; }
    .panel-head p { margin:4px 0 0; color:var(--muted); }
    .badge { border:1px solid var(--line); border-radius:999px; padding:4px 8px; font-size:12px; color:var(--muted); white-space:nowrap; }
    .badge.selected { color:#92400e; border-color:#f59e0b; background:#fffbeb; }
    .image-box { border:1px solid var(--line); border-radius:8px; background:#fff; padding:10px; overflow:auto; }
    .image-box img { display:block; max-width:100%; height:auto; margin:auto; }
    .image-box svg { display:block; max-width:100%; height:auto; margin:auto; }
    .image-empty { border:1px dashed var(--line); border-radius:8px; padding:28px; text-align:center; color:var(--muted); background:#fff; }
    .panel-info { margin-top:12px; display:grid; gap:10px; }
    details { border:1px solid var(--line); border-radius:8px; background:#fff; margin:0; overflow:hidden; }
    details > summary { cursor:pointer; padding:9px 11px; font-weight:800; background:#f8fafc; border-bottom:1px solid var(--line); }
    details:not([open]) > summary { border-bottom:0; }
    table { width:100%; border-collapse:collapse; font-size:12px; margin:10px; max-width:calc(100% - 20px); }
    th, td { border-bottom:1px solid #e5e7eb; padding:6px 7px; text-align:left; vertical-align:top; overflow-wrap:anywhere; }
    th { color:#334155; background:#f8fafc; }
    pre { margin:10px; padding:10px; border:1px solid var(--line); border-radius:8px; background:#0f172a; color:#e5e7eb; overflow:auto; max-height:520px; font-size:12px; }
    .empty { border:1px dashed var(--line); border-radius:8px; padding:24px; text-align:center; color:var(--muted); background:#fff; }
    @media (max-width:1100px) { .global-info, main { grid-template-columns:1fr; } aside { position:static; max-height:none; } }
  </style>
</head>
<body>
  <header><h1>MolGR Trace</h1></header>
  <section id="global-info" class="global-info"></section>
  <main>
    <aside>
      <input id="tree-search" type="search" placeholder="Filter trace nodes...">
      <nav id="tree" class="tree"></nav>
    </aside>
    <section id="detail" class="content"></section>
  </main>
  <script id="trace-json" type="application/json">"""
        + cases_json
        + """</script>
  <script>
    (() => {
      const trace = JSON.parse(document.getElementById("trace-json").textContent);
      const scorePrefixes = ["analysis_", "force_field_", "metal_", "organic_", "passes_", "selection_"];
      const scoreKeys = new Set(["combination_index", "score"]);
      const nodes = [];
      const nodeById = new Map();
      const roots = [];

      function hasValue(value) {
        return value !== undefined && value !== null && value !== "" &&
          !(Array.isArray(value) && value.length === 0) &&
          !(typeof value === "object" && !Array.isArray(value) && Object.keys(value).length === 0);
      }

      function fmt(value) {
        if (value === undefined || value === null) return "";
        if (typeof value === "boolean") return value ? "是" : "否";
        if (typeof value === "number") return Number.isFinite(value) ? String(Number.parseFloat(value.toPrecision(6))) : String(value);
        if (typeof value === "string") return value;
        return JSON.stringify(value);
      }

      function imagePath(image) {
        if (!image || typeof image !== "object" || image.status) return "";
        return image.display_path || image.path || "";
      }

      function imageSvg(image) {
        if (!image || typeof image !== "object" || image.status) return "";
        return image.svg_fragment || "";
      }

      function stateImage(value) {
        return value && value.state && value.state.dof_image ? value.state.dof_image : null;
      }

      function scoreDetails(metadata) {
        if (!metadata || typeof metadata !== "object") return {};
        if (metadata.score_details && typeof metadata.score_details === "object") return metadata.score_details;
        const source = metadata.metadata && typeof metadata.metadata === "object" ? metadata.metadata : metadata;
        const result = {};
        Object.entries(source).forEach(([key, value]) => {
          if (scoreKeys.has(key) || scorePrefixes.some(prefix => key.startsWith(prefix))) result[key] = value;
        });
        return result;
      }

      function makeNode({label, kind, caseId, summary = [], metadata = {}, image = null, selected = false, children = []}) {
        const node = {
          id: `node-${nodes.length}`,
          label,
          kind,
          caseId,
          summary,
          metadata,
          image,
          selected,
          children,
          collapsed: undefined,
          parent: null,
          treeLi: null,
          treeToggle: null,
          treeChildren: null,
          scoreDetails: scoreDetails(metadata),
        };
        node.search = [caseId, label, kind, fmt(summary), fmt(node.scoreDetails)].join(" ");
        nodes.push(node);
        nodeById.set(node.id, node);
        return node;
      }

      function stateSummary(state) {
        state = state || {};
        return [
          ["canonical SMILES", state.canonical_smiles || state.smiles],
          ["形式电荷", state.formal_charge_sum],
          ["自由基和", state.spin_multiplicity_sum],
          ["自由基奇偶和", state.spin_multiplicity_singlet_sum],
          ["剩余电荷预算", state.given_charge],
          ["匹配目标", state.valid_for_target],
          ["带电原子", state.charged_atom_counts],
          ["自由基原子", state.radical_atom_counts],
        ];
      }

      function metalStatesLabel(states) {
        return (states || []).map(state => {
          const valence = Number(state.valence || 0);
          const sign = valence >= 0 ? `+${valence}` : String(valence);
          return `#${state.idx || ""} ${state.symbol || ""}(${sign}, r${state.radical_num || 0})`;
        }).join("; ");
      }

      function buildSelectedPathAnimation(caseId, traceData, selected, kind) {
        const animation = traceData && traceData.selected_path_animation;
        if (!animation) return null;
        return makeNode({
          label: "最终选中路径动画",
          kind,
          caseId,
          summary: [
            ["状态", animation.status || "rendered"],
            ["帧数", animation.frame_count],
            ["动画", animation.animation],
            ["路径", animation.display_path || animation.path],
            ["目标电荷", traceData.target && traceData.target.total_charge],
            ["目标自由基电子", traceData.target && traceData.target.total_radical_electrons],
            ["选中分数", selected && selected.score],
            ["选中选择键", selected && selected.organic_topology_selection_key],
          ],
          metadata: animation,
          image: animation,
          selected: true,
        });
      }

      function selectedNoMetalTrace(item) {
        const selected = item && item.selected_candidate;
        const selectedTarget = selected && selected.target;
        const selectedLayer = item && item.search && item.search.selected_layer_index;
        const layers = item && item.search && item.search.layer_summaries || [];
        if (!selectedTarget) return null;
        const selectedCharge = String(selectedTarget.no_metal_charge);
        const selectedRadicals = String(selectedTarget.no_metal_radical_electrons);
        const orderedLayers = [
          ...layers.filter(layer => String(layer.layer_index) === String(selectedLayer)),
          ...layers.filter(layer => String(layer.layer_index) !== String(selectedLayer)),
        ];
        for (const layer of orderedLayers) {
          for (const bucket of layer.target_buckets || []) {
            const target = bucket.target || {};
            if (
              String(target.no_metal_charge) === selectedCharge &&
              String(target.no_metal_radical_electrons) === selectedRadicals &&
              bucket.no_metal_trace
            ) {
              return bucket.no_metal_trace;
            }
          }
        }
        return null;
      }

      function buildNoMetal(caseId, label, traceData) {
        traceData = traceData || {};
        const target = traceData.target || {};
        const selected = traceData.selected_candidate || {};
        const selectedState = selected.state || {};
        const resonance = traceData.resonance || {};
        const children = [];
        const selectedPathAnimation = buildSelectedPathAnimation(
          caseId,
          traceData,
          selected,
          "无金属最终路径动画",
        );
        if (selectedPathAnimation) children.push(selectedPathAnimation);

        (traceData.pipeline_steps || []).forEach(step => {
          children.push(makeNode({
            label: `${step.step_index}. ${step.phase}`,
            kind: "无金属生产管线阶段",
            caseId,
            summary: [
              ["阶段", step.phase],
              ["类型", step.kind],
            ],
            metadata: step,
            image: step.dof_image,
          }));
        });

        if (Object.keys(resonance).length) {
          children.push(makeNode({
            label: "共振选择",
            kind: "无金属共振元数据",
            caseId,
            summary: [
              ["种子序号", resonance.seed_index],
              ["共振序号", resonance.resonance_index],
              ["raw 共振序号", resonance.raw_index],
              ["归一化方式", resonance.normalization],
            ],
            metadata: resonance,
          }));
        }

        if (Object.keys(selected).length) {
          children.push(makeNode({
            label: "选中 no-metal 候选",
            kind: "无金属候选",
            caseId,
            summary: [
              ["分数", selected.score],
              ["选择键", selected.organic_topology_selection_key],
              ...stateSummary(selectedState),
            ],
            metadata: selected,
            image: selectedState.dof_image,
            selected: true,
          }));
        }

        return makeNode({
          label,
          kind: "无金属重建",
          caseId,
          summary: [
            ["状态", traceData.status],
            ["目标电荷", target.total_charge],
            ["目标自由基电子", target.total_radical_electrons],
            ["邻位自由基动作", traceData.neighbor_radical_actions],
            ["恢复层", traceData.recovery_tier],
            ["归一化方式", resonance.normalization],
            ["选中 canonical SMILES", selectedState.canonical_smiles],
            ["选中分数", selected.score],
            ["选中选择键", selected.organic_topology_selection_key],
          ],
          metadata: traceData,
          image: selectedState.dof_image,
          children,
        });
      }

      function buildTree() {
        (trace.cases || []).forEach(item => {
          const caseId = String(item.id || "unknown");
          const base = item.base_state || {};
          const selected = item.selected_candidate || {};
          const children = [];
          const caseNode = makeNode({
            label: caseId,
            kind: "case",
            caseId,
            summary: [
              ["id", item.id],
              ["CSV 行号", item.row_index],
              ["状态", item.status],
              ["重建类型", item.trace_kind],
              ["总电荷", item.charge],
              ["总自由基电子数", item.total_radical_electrons],
              ["自旋多重度", item.spin_multiplicity],
              ["自旋来源", item.spin_source],
              ["XYZ 路径", item.xyz_path],
              ["XYZ 来源", item.xyz_source],
              ["review fixture 类型", item.review_fixture && item.review_fixture.kind],
              ["review fixture 等价", item.review_fixture && item.review_fixture.equivalent],
              ["review fixture 比较原因", item.review_fixture && item.review_fixture.equivalence_reason],
              ["参考 SMILES", item.reference_smiles],
              ["金属原子数", base.metal_atom_count],
              ["生产选择层", item.search && item.search.selected_layer_index],
              ["生产候选数", item.production_candidate_count],
              ["全部已评分候选数", item.candidate_count],
              ["选中候选", selected.combination_index],
              ["耗时秒", item.elapsed_seconds],
            ],
            metadata: item,
            children,
          });
          roots.push(caseNode);
          if (item.review_fixture) {
            children.push(makeNode({
              label: "review fixture 同步检查",
              kind: "fixture equivalence",
              caseId,
              summary: [
                ["fixture 类型", item.review_fixture.kind],
                ["结构文件", item.review_fixture.structure_file],
                ["等价", item.review_fixture.equivalent],
                ["方法", item.review_fixture.equivalence_method],
                ["原因", item.review_fixture.equivalence_reason],
              ],
              metadata: item.review_fixture,
            }));
          }
          if (item.trace_kind === "no_metal") {
            children.push(buildNoMetal(caseId, "无金属重建", item.no_metal_trace || {}));
            return;
          }
          const selectedPathAnimation = buildSelectedPathAnimation(
            caseId,
            item,
            selected,
            "含金属最终路径动画",
          );
          if (selectedPathAnimation) children.push(selectedPathAnimation);
          const selectedOrganicTrace = selectedNoMetalTrace(item);
          if (selectedOrganicTrace) {
            children.push(buildNoMetal(caseId, "选中无金属重建", selectedOrganicTrace));
          }
          children.push(makeNode({
            label: "基础状态",
            kind: "含金属重建",
            caseId,
            summary: [
              ["金属原子数", base.metal_atom_count],
              ["金属位点数", (base.available_metal_states_by_site || []).length],
              ["phase_history", base.phase_history],
            ],
            metadata: base,
          }));
          if (item.dof_candidate_grid) {
            children.push(makeNode({
              label: "金属候选对比",
              kind: "含金属候选",
              caseId,
              summary: [
                ["候选数", item.candidate_count],
                ["生产候选数", item.production_candidate_count],
              ],
              metadata: {
                candidate_count: item.candidate_count,
                production_candidate_count: item.production_candidate_count,
                selected_candidate: item.selected_candidate,
              },
              image: item.dof_candidate_grid,
            }));
          }
          const candidates = (item.candidates || []).filter(candidate => candidate && typeof candidate === "object");
          ((item.search && item.search.layer_summaries) || []).forEach(layer => {
            const layerIndex = layer.layer_index;
            const layerChildren = [];
            (layer.target_buckets || []).forEach(bucket => {
              const target = bucket.target || {};
              const organic = bucket.organic_part || {};
              const bucketChildren = [];
              if (bucket.no_metal_trace) {
                bucketChildren.push(buildNoMetal(caseId, "对应无金属重建", bucket.no_metal_trace));
              }
              layerChildren.push(makeNode({
                label: `Target Q=${target.no_metal_charge} R=${target.no_metal_radical_electrons}`,
                kind: "有机目标桶",
                caseId,
                summary: [
                  ["层", layerIndex],
                  ["有机目标电荷", target.no_metal_charge],
                  ["有机自由基电子", target.no_metal_radical_electrons],
                  ["状态", bucket.status],
                  ["候选数", bucket.candidate_count],
                  ["已评分数", bucket.prepared_candidate_count],
                  ["有机 canonical SMILES", organic.canonical_smiles || organic.smiles],
                  ["有机力场分", bucket.no_metal_score],
                ],
                metadata: bucket,
                image: bucket.dof_image,
                children: bucketChildren,
              }));
            });
            candidates
              .filter(candidate => String(candidate.search_layer_index) === String(layerIndex))
              .forEach(candidate => {
                const organic = candidate.organic_part || {};
                const target = candidate.target || {};
                const title = `${candidate.selected ? "selected " : ""}L${candidate.search_layer_index}/C${candidate.combination_index}`;
                layerChildren.push(makeNode({
                  label: title,
                  kind: "金属候选",
                  caseId,
                  summary: [
                    ["选中", candidate.selected],
                    ["组合", candidate.combination_index],
                    ["生产层", candidate.in_production_selection_layer],
                    ["候选总电荷", candidate.candidate_total_charge],
                    ["有机目标电荷", target.no_metal_charge],
                    ["有机自由基电子", target.no_metal_radical_electrons],
                    ["有机 canonical SMILES", organic.canonical_smiles || organic.smiles],
                    ["分数", candidate.score],
                    ["金属状态", metalStatesLabel(candidate.metal_states)],
                  ],
                  metadata: candidate,
                  image: candidate.dof_image,
                  selected: Boolean(candidate.selected),
                }));
              });
            children.push(makeNode({
              label: `Layer ${layerIndex}`,
              kind: "金属搜索层",
              caseId,
              summary: [
                ["状态", layer.status],
                ["生产选择层", layer.production_selected_layer],
                ["金属组数", layer.state_group_count],
                ["每组候选数", layer.state_options_per_group],
                ["目标桶数", layer.target_bucket_count],
                ["枚举候选数", layer.candidate_count],
                ["已评分候选数", layer.prepared_candidate_count],
                ["analysis score context", Boolean(layer.analysis_score_context)],
              ],
              metadata: layer,
              children: layerChildren,
            }));
          });
          if (item.analysis) {
            children.push(makeNode({
              label: "分析和评分上下文",
              kind: "analysis",
              caseId,
              summary: [
                ["score_all_candidates", item.analysis.score_all_candidates],
                ["note", item.analysis.note],
              ],
              metadata: item.analysis,
            }));
          }
        });
      }

      function table(rows, headers = ["字段", "值"]) {
        const tableEl = document.createElement("table");
        const thead = tableEl.createTHead();
        const headRow = thead.insertRow();
        headers.forEach(header => {
          const th = document.createElement("th");
          th.textContent = header;
          headRow.appendChild(th);
        });
        const tbody = tableEl.createTBody();
        rows.filter(([_, value]) => hasValue(value)).forEach(row => {
          const tr = tbody.insertRow();
          row.forEach(cell => {
            const td = tr.insertCell();
            td.textContent = fmt(cell);
          });
        });
        if (!tbody.rows.length) {
          const tr = tbody.insertRow();
          const td = tr.insertCell();
          td.colSpan = headers.length;
          td.textContent = "无";
        }
        return tableEl;
      }

      function details(title, content, open = false) {
        const el = document.createElement("details");
        if (open) el.open = true;
        const summary = document.createElement("summary");
        summary.textContent = title;
        el.appendChild(summary);
        el.appendChild(content);
        return el;
      }

      function pre(value) {
        const el = document.createElement("pre");
        el.textContent = JSON.stringify(value || {}, null, 2);
        return el;
      }

      function renderGlobal() {
        const root = document.getElementById("global-info");
        const input = trace.input || {};
        const dof = trace.dof_rendering || {};
        const inputCard = document.createElement("section");
        inputCard.className = "global-card";
        inputCard.innerHTML = "<h2>输入</h2>";
        inputCard.appendChild(table([
          ["来源", input.source],
          ["CSV", input.csv],
          ["XYZ 目录", input.xyz_dir],
          ["id", input.ids],
          ["默认总电荷", input.total_charge],
          ["默认总自由基电子数", input.total_radical_electrons],
          ["自旋来源", input.spin_source],
          ["样本数", trace.case_count],
        ]));
        const dofCard = document.createElement("section");
        dofCard.className = "global-card";
        dofCard.innerHTML = "<h2>DOF 渲染</h2>";
        dofCard.appendChild(table([
          ["存储", dof.storage],
          ["图片目录", dof.image_dir],
          ["格式", dof.format],
          ["图片数", dof.image_count],
          ["跳过数", dof.skipped_count],
          ["最大图片数", dof.max_images],
          ["错误数", (dof.errors || []).length],
        ]));
        const caseCard = document.createElement("section");
        caseCard.className = "global-card";
        caseCard.innerHTML = "<h2>Cases</h2>";
        const caseRows = (trace.cases || []).map(item => [
          item.id,
          item.status,
          item.trace_kind,
          item.charge,
          item.total_radical_electrons,
          item.base_state && item.base_state.metal_atom_count,
          item.search && item.search.selected_layer_index,
          item.production_candidate_count,
          item.candidate_count,
          item.elapsed_seconds,
        ]);
        caseCard.appendChild(table(caseRows.map(row => ["", row]), ["", ""]));
        const caseTable = caseCard.querySelector("table");
        caseTable.innerHTML = "";
        const head = caseTable.createTHead().insertRow();
        ["id", "状态", "类型", "电荷", "自由基", "金属数", "生产层", "生产候选", "全部候选", "耗时秒"].forEach(text => {
          const th = document.createElement("th");
          th.textContent = text;
          head.appendChild(th);
        });
        const body = caseTable.createTBody();
        caseRows.forEach(row => {
          const tr = body.insertRow();
          row.forEach(cell => {
            const td = tr.insertCell();
            td.textContent = fmt(cell);
          });
        });
        root.append(inputCard, dofCard, caseCard);
      }

      function syncTreeNode(node, forcedOpen = false) {
        const hasChildren = node.children.length > 0;
        const isOpen = hasChildren && (!node.collapsed || forcedOpen);
        if (node.treeLi) {
          node.treeLi.classList.toggle("is-expanded", hasChildren && isOpen);
          node.treeLi.classList.toggle("is-collapsed", hasChildren && !isOpen);
        }
        if (node.treeToggle) {
          node.treeToggle.setAttribute("aria-expanded", hasChildren ? String(isOpen) : "false");
          node.treeToggle.disabled = !hasChildren;
        }
        if (node.treeChildren) {
          node.treeChildren.hidden = hasChildren && !isOpen;
        }
      }

      function renderTreeNode(node, depth = 0, parent = null) {
        const li = document.createElement("li");
        node.parent = parent;
        if (node.children.length && typeof node.collapsed !== "boolean") {
          node.collapsed = depth >= 2;
        }
        const row = document.createElement("div");
        row.className = "tree-row";
        if (node.children.length) {
          const toggle = document.createElement("button");
          toggle.type = "button";
          toggle.className = "tree-toggle";
          toggle.setAttribute("aria-label", "折叠或展开节点");
          toggle.addEventListener("click", event => {
            event.preventDefault();
            event.stopPropagation();
            node.collapsed = !node.collapsed;
            syncTreeNode(node);
            applyTreeFilter();
          });
          row.appendChild(toggle);
          node.treeToggle = toggle;
        } else {
          const spacer = document.createElement("span");
          spacer.className = "tree-toggle-spacer";
          spacer.setAttribute("aria-hidden", "true");
          row.appendChild(spacer);
        }
        const button = document.createElement("button");
        button.type = "button";
        button.className = `tree-node${node.selected ? " is-selected" : ""}`;
        button.dataset.target = node.id;
        button.dataset.search = node.search || "";
        const label = document.createElement("span");
        label.className = "tree-label";
        label.textContent = node.label;
        const meta = document.createElement("span");
        meta.className = "tree-meta";
        meta.textContent = node.kind;
        button.append(label, meta);
        button.addEventListener("click", () => activate(node.id));
        row.appendChild(button);
        li.appendChild(row);
        if (node.children.length) {
          const ol = document.createElement("ol");
          ol.className = "tree-children";
          node.children.forEach(child => ol.appendChild(renderTreeNode(child, depth + 1, node)));
          li.appendChild(ol);
          node.treeChildren = ol;
        }
        node.treeLi = li;
        syncTreeNode(node);
        return li;
      }

      function renderPanel(node) {
        const panel = document.createElement("article");
        panel.id = node.id;
        panel.className = "node-panel";
        const head = document.createElement("div");
        head.className = "panel-head";
        const titleBox = document.createElement("div");
        const title = document.createElement("h2");
        title.textContent = node.label;
        const meta = document.createElement("p");
        meta.textContent = `${node.caseId} · ${node.kind}`;
        titleBox.append(title, meta);
        head.appendChild(titleBox);
        if (node.selected) {
          const badge = document.createElement("span");
          badge.className = "badge selected";
          badge.textContent = "selected";
          head.appendChild(badge);
        }
        panel.appendChild(head);
        const svg = imageSvg(node.image);
        const imgPath = imagePath(node.image);
        if (svg) {
          const box = document.createElement("div");
          box.className = "image-box";
          box.innerHTML = svg;
          panel.appendChild(box);
        } else if (imgPath) {
          const box = document.createElement("div");
          box.className = "image-box";
          const img = document.createElement("img");
          img.src = imgPath;
          img.alt = node.label;
          box.appendChild(img);
          panel.appendChild(box);
        } else {
          const empty = document.createElement("div");
          empty.className = "image-empty";
          empty.textContent = "此节点没有可渲染的 DOF 图像";
          panel.appendChild(empty);
        }
        const info = document.createElement("section");
        info.className = "panel-info";
        info.appendChild(details("摘要", table(node.summary), true));
        if (Object.keys(node.scoreDetails || {}).length) {
          info.appendChild(details("分数构成", table(Object.entries(node.scoreDetails)), true));
        }
        info.appendChild(details("完整 JSON", pre(node.metadata)));
        panel.appendChild(info);
        return panel;
      }

      function activate(id) {
        const node = nodeById.get(id);
        let current = node ? node.parent : null;
        while (current) {
          if (current.collapsed) {
            current.collapsed = false;
            syncTreeNode(current);
          }
          current = current.parent;
        }
        document.querySelectorAll(".tree-node").forEach(button => {
          button.classList.toggle("is-active", button.dataset.target === id);
        });
        document.querySelectorAll(".node-panel").forEach(panel => {
          panel.classList.toggle("is-active", panel.id === id);
        });
        applyTreeFilter();
      }

      function applyTreeFilter() {
        const q = (document.getElementById("tree-search").value || "").trim().toLowerCase();

        function walk(node) {
          const hasChildren = node.children.length > 0;
          const selfMatches = !q || (node.search || "").toLowerCase().includes(q);
          let descendantMatches = false;
          node.children.forEach(child => {
            descendantMatches = walk(child) || descendantMatches;
          });
          const visible = selfMatches || descendantMatches || !q;
          if (node.treeLi) {
            node.treeLi.hidden = !visible;
          }
          const forcedOpen = Boolean(q && hasChildren && (selfMatches || descendantMatches));
          syncTreeNode(node, forcedOpen);
          return visible;
        }

        roots.forEach(root => walk(root));
      }

      buildTree();
      renderGlobal();
      const treeRoot = document.createElement("ol");
      roots.forEach(root => treeRoot.appendChild(renderTreeNode(root)));
      document.getElementById("tree").appendChild(treeRoot);
      nodes.forEach(node => document.getElementById("detail").appendChild(renderPanel(node)));
      if (nodes.length) activate(nodes[0].id);
      document.getElementById("tree-search").addEventListener("input", applyTreeFilter);
      applyTreeFilter();
    })();
  </script>
</body>
</html>
"""
    )


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


def _no_metal_candidate_trace(
    candidate: ReconstructionState,
    *,
    selected: bool,
    render_context: DofRenderContext | None = None,
    case_id: str = "",
    image_label: str = "no_metal_candidate",
    config: MolGRConfig | None = None,
) -> dict[str, Any]:
    selection_key = no_metal_selection._no_metal_candidate_selection_key(candidate, config=config)
    state = _omol_state_snapshot(
        candidate.omol,
        given_charge=candidate.given_charge,
        target_charge=candidate.total_charge,
        target_radical_electrons=candidate.total_radical_electrons,
    )
    image = _render_dof_molecule(
        candidate.omol,
        render_context=render_context,
        case_id=case_id,
        label=image_label,
        kind="no_metal_candidate",
    )
    if image is not None:
        state["dof_image"] = image
    return {
        "selected": selected,
        "phase_history": list(candidate.phase_history),
        "state": state,
        "score": candidate.metadata.get("score"),
        "force_field_energy": candidate.metadata.get("force_field_energy"),
        "organic_topology_selection_key": selection_key,
        "metadata": _jsonable(candidate.metadata),
    }


def _production_phase_steps(
    state: ReconstructionState,
    step_images: dict[int, dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Expose the phases recorded by the production no-metal state machine."""

    steps: list[dict[str, Any]] = []
    for step_index, phase in enumerate(state.phase_history):
        if phase.startswith("recover_"):
            kind = "recovery"
        elif phase.startswith("resolve_neighbor_radicals"):
            kind = "neighbor_radical_seed"
        elif "resonance" in phase:
            kind = "resonance"
        else:
            kind = "stage"
        step = {"step_index": step_index, "phase": phase, "kind": kind}
        if step_images is not None and step_index in step_images:
            step["dof_image"] = step_images[step_index]
        steps.append(step)
    return steps


def _trace_no_metal_reconstruction(
    xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
    *,
    render_context: DofRenderContext | None = None,
    case_id: str = "",
    config: MolGRConfig | None = None,
) -> dict[str, Any]:
    trace: dict[str, Any] = {
        "target": {
            "total_charge": total_charge,
            "total_radical_electrons": total_radical_electrons,
        },
        "status": "pending",
    }
    if total_radical_electrons < 0:
        trace["status"] = "invalid_negative_radicals"
        return trace

    seed_state = no_metal_preparation._seed_state(
        xyz_block,
        total_charge,
        total_radical_electrons,
    )
    snapshot_records: dict[int, dict[str, Any]] = {
        0: {"phase": "read_xyz", "parent_id": -1, "omol": _copy_omol(seed_state.omol)},
        1: {
            "phase": "normalize_seed_electronic_labels",
            "parent_id": 0,
            "omol": _copy_omol(seed_state.omol),
        },
    }
    seed_state.metadata["_trace_snapshot_id"] = 1
    next_snapshot_id = 2

    def observe_phase(machine: Any, phase: str) -> None:
        nonlocal next_snapshot_id
        parent_id = int(machine.metadata.get("_trace_snapshot_id", -1))
        snapshot_id = next_snapshot_id
        next_snapshot_id += 1
        machine.metadata["_trace_snapshot_id"] = snapshot_id
        snapshot_records[snapshot_id] = {
            "phase": phase,
            "parent_id": parent_id,
            "omol": _copy_omol(machine.omol),
        }

    with _observe_omol_phases(observe_phase):
        selected = reconstruct_without_metals._run_no_metal_pipeline_from_state(
            seed_state,
            config=config,
        )
    if selected is None:
        trace["status"] = "no_valid_no_metal_candidate"
        return trace

    selected_snapshot_id = int(selected.metadata.get("_trace_snapshot_id", -1))
    selected_path_records: list[dict[str, Any]] = []
    seen_snapshot_ids: set[int] = set()
    while selected_snapshot_id >= 0 and selected_snapshot_id not in seen_snapshot_ids:
        seen_snapshot_ids.add(selected_snapshot_id)
        record = snapshot_records.get(selected_snapshot_id)
        if record is None:
            break
        selected_path_records.append(record)
        selected_snapshot_id = int(record.get("parent_id", -1))
    selected_path_records.reverse()
    selected.metadata.pop("_trace_snapshot_id", None)
    step_images: dict[int, dict[str, Any]] = {}
    for step_index, record in enumerate(selected_path_records):
        image = _render_dof_molecule(
            record["omol"],
            render_context=render_context,
            case_id=case_id,
            label=f"step_{step_index}_{record['phase']}",
            kind="no_metal_pipeline_step",
        )
        if image is not None:
            step_images[step_index] = image

    trace["status"] = "selected"
    trace["pipeline_steps"] = _production_phase_steps(selected, step_images)
    trace["neighbor_radical_actions"] = list(selected.metadata.get("neighbor_radical_actions", ()))
    trace["recovery_tier"] = int(selected.metadata.get("recovery_tier", 0))
    trace["resonance"] = {
        "seed_index": selected.metadata.get("resonance_seed_index"),
        "resonance_index": selected.metadata.get("resonance_index"),
        "raw_index": selected.metadata.get("resonance_raw_index"),
        "normalization": selected.metadata.get("resonance_normalization"),
    }
    trace["selected_candidate"] = _no_metal_candidate_trace(
        selected,
        selected=True,
        render_context=render_context,
        case_id=case_id,
        image_label="selected_no_metal_candidate",
        config=config,
    )

    animation = _render_dof_animation(
        [
            (record["omol"], str(record["phase"])) for record in selected_path_records
        ],
        render_context=render_context,
        case_id=case_id,
        label="selected_no_metal_path",
        kind="selected_path_animation",
    )
    if animation is not None:
        trace["selected_path_animation"] = animation
    return trace


def _render_selected_metal_path_animation(
    selected_candidate: MetalCandidateState | None,
    *,
    no_metal_xyz_block: str,
    render_context: DofRenderContext | None,
    case_id: str,
) -> dict[str, Any] | None:
    if selected_candidate is None or selected_candidate.no_metal_state is None:
        return None
    seed_state = no_metal_preparation._seed_state(
        no_metal_xyz_block,
        selected_candidate.no_metal_charge_target,
        selected_candidate.no_metal_radical_target,
    )
    try:
        combined_omol = selected_candidate.materialize_combined_omol(
            preparation.combine_metal_with_omol
        )
    except Exception as exc:
        error = {
            "kind": "selected_metal_path_animation",
            "label": "selected_metal_path",
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
            "animation": True,
        }
        if render_context is not None:
            render_context.errors.append(error)
        return error
    return _render_dof_animation(
        [
            (_copy_omol(seed_state.omol), "normalized organic seed"),
            (_copy_omol(selected_candidate.no_metal_state.omol), "selected organic part"),
            (combined_omol, "selected metal candidate"),
        ],
        render_context=render_context,
        case_id=case_id,
        label="selected_metal_path",
        kind="selected_metal_path_animation",
    )


def _annotate_analysis_scores_for_all_candidates(
    candidates: Sequence[MetalCandidateState],
    *,
    config: MolGRConfig | None = None,
) -> dict[str, Any]:
    if not candidates:
        return {}

    scoring._annotate_candidate_set_discordance_features(candidates, config=config)
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
    render_context: DofRenderContext | None = None,
    case_id: str = "",
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
    metal_valence_sum = sum(int(metal_state.valence) for metal_state in candidate.metal_states)
    no_metal_total_charge = int(no_metal_state.total_charge) if no_metal_state is not None else None
    candidate_total_charge = (
        no_metal_total_charge + metal_valence_sum if no_metal_total_charge is not None else None
    )
    report = {
        "candidate_index": candidate_index,
        "candidate_identity": {
            "search_layer_index": candidate_identity[0],
            "combination_index": candidate_identity[1],
        },
        "search_layer_index": candidate_identity[0],
        "combination_index": combination_index,
        "selected": selected_candidate_identity == candidate_identity,
        "in_production_selection_layer": bool(production_metadata),
        "candidate_total_charge": candidate_total_charge,
        "metal_valence_sum": metal_valence_sum,
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
    if candidate.no_metal_state is not None:
        image_label = f"L{candidate_identity[0]}_C{combination_index}"
        try:
            combined_omol = candidate.materialize_combined_omol(preparation.combine_metal_with_omol)
            image = _render_dof_molecule(
                combined_omol,
                render_context=render_context,
                case_id=case_id,
                label=image_label,
                kind="metal_candidate",
            )
        except Exception as exc:
            image = {
                "kind": "metal_candidate",
                "label": image_label,
                "status": "error",
                "error": f"{type(exc).__name__}: {exc}",
            }
            if render_context is not None:
                render_context.errors.append(image)
        if image is not None:
            report["dof_image"] = image
    return report


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
    render_context: DofRenderContext | None = None,
    case_id: str = "",
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
                    render_context=render_context,
                    case_id=(
                        f"{case_id}_L{layer_index}_Q{prototype.no_metal_charge_target}"
                        f"_R{prototype.no_metal_radical_target}"
                    ),
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
            organic_image = _render_dof_molecule(
                no_metal_state.omol,
                render_context=render_context,
                case_id=case_id,
                label=(
                    f"L{layer_index}_organic_Q{prototype.no_metal_charge_target}"
                    f"_R{prototype.no_metal_radical_target}"
                ),
                kind="metal_target_organic_part",
            )
            if organic_image is not None:
                target_summary["dof_image"] = organic_image

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
            render_context=render_context,
            case_id=case_id,
        )
        for candidate_index, candidate in enumerate(all_scored_candidates)
    ]
    candidate_grid_items: list[tuple[pybel.Molecule, str]] = []
    for candidate_index, candidate in enumerate(all_scored_candidates):
        if candidate.no_metal_state is None:
            continue
        try:
            candidate_identity = _candidate_identity(
                candidate,
                fallback_candidate_index=candidate_index,
            )
            candidate_grid_items.append(
                (
                    candidate.materialize_combined_omol(preparation.combine_metal_with_omol),
                    (
                        f"{'selected ' if selected_candidate_identity == candidate_identity else ''}"
                        f"L{candidate_identity[0]}/C{candidate_identity[1]} "
                        f"score={_format_value(candidate.score)}"
                    ),
                )
            )
        except Exception:
            continue
    candidate_grid = _render_dof_grid(
        candidate_grid_items,
        render_context=render_context,
        case_id=case_id,
        label="metal_candidates",
        kind="metal_candidate_grid",
    )

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
        if selected_summary is not None and selected_candidate is not None:
            try:
                selected_omol = selected_candidate.materialize_combined_omol(
                    preparation.combine_metal_with_omol
                )
                selected_summary["graph"] = _omol_smiles_detail(selected_omol)
            except Exception as exc:
                selected_summary["graph_error"] = f"{type(exc).__name__}: {exc}"
    selected_path_animation = _render_selected_metal_path_animation(
        selected_candidate,
        no_metal_xyz_block=base_state.no_metal_xyz_block,
        render_context=render_context,
        case_id=case_id,
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
        "selected_path_animation": selected_path_animation,
        "candidates": candidate_reports,
        "dof_candidate_grid": candidate_grid,
    }


def _metal_atom_count_from_xyz(xyz_block: str) -> int:
    omol = pybel.readstring("xyz", xyz_block)
    return sum(1 for atom in omol.atoms if atom.OBAtom.IsMetal())


def _selected_trace_smiles(trace: dict[str, Any]) -> str:
    selected = trace.get("selected_candidate")
    if not isinstance(selected, dict):
        return ""
    if trace.get("trace_kind") == "no_metal":
        state = selected.get("state")
        if isinstance(state, dict):
            return str(state.get("canonical_smiles") or state.get("smiles") or "")
        return ""
    graph = selected.get("graph")
    if isinstance(graph, dict):
        return str(graph.get("canonical_smiles") or graph.get("smiles") or "")
    return ""


def _mol_from_unsanitized_smiles(smiles: str) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles, sanitize=False) if smiles.strip() else None
    if mol is not None:
        mol.UpdatePropertyCache(strict=False)
    return mol


def _review_fixture_trace_check(
    input_case: TraceInputCase,
    trace: dict[str, Any],
) -> dict[str, Any]:
    trace_smiles = _selected_trace_smiles(trace)
    check: dict[str, Any] = {
        "kind": input_case.fixture_kind,
        "structure_file": input_case.fixture_structure_file,
        "expected_smiles": input_case.expected_smiles,
        "trace_smiles": trace_smiles,
        "equivalent": None,
        "equivalence_method": "",
        "equivalence_reason": "",
    }
    if not input_case.expected_smiles:
        check["equivalence_reason"] = "fixture_has_no_approved_answer"
        return check
    if not trace_smiles:
        check["equivalent"] = False
        check["equivalence_reason"] = "trace_has_no_selected_graph"
        return check

    expected = _mol_from_unsanitized_smiles(input_case.expected_smiles)
    traced = _mol_from_unsanitized_smiles(trace_smiles)
    if expected is None:
        check["equivalent"] = False
        check["equivalence_reason"] = "fixture_approved_smiles_invalid"
        return check
    if traced is None:
        check["equivalent"] = False
        check["equivalence_reason"] = "trace_selected_smiles_invalid"
        return check

    equivalent, info = check_equivalence(expected, traced, use_chirality=False)
    check["equivalent"] = equivalent
    check["equivalence_method"] = info.method.value if info.method else ""
    check["equivalence_reason"] = info.reason
    return check


def trace_reconstruction_case(
    input_case: TraceInputCase,
    *,
    score_all_candidates: bool,
    render_context: DofRenderContext | None = None,
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
            render_context=render_context,
            case_id=input_case.id,
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
            render_context=render_context,
            case_id=input_case.id,
            config=config,
        )
        trace["trace_kind"] = "metal"

    if input_case.fixture_kind:
        trace["review_fixture"] = _review_fixture_trace_check(input_case, trace)

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
    render_context: DofRenderContext | None = None,
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
                    render_context=render_context,
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


def render_trace_report(
    input_cases: Sequence[TraceInputCase],
    *,
    score_all_candidates: bool = True,
    dof_max_images: int = 1000,
    config: MolGRConfig | None = None,
) -> str:
    """Build the browser trace report for already-resolved input cases.

    This is the shared programmatic entry point for browser-facing callers.
    Keeping report assembly here prevents the development UI from maintaining
    a second trace renderer or case representation.
    """

    render_context = DofRenderContext(
        image_dir=Path("molgr_trace_dof_images"),
        display_base_dir=None,
        image_format="svg",
        max_images=dof_max_images,
    )
    cases = trace_reconstruction_cases(
        input_cases,
        score_all_candidates=score_all_candidates,
        render_context=render_context,
        config=config,
    )
    output = {
        "input": {
            "source": "review_page",
            "ids": [case.id for case in input_cases],
        },
        "case_count": len(cases),
        "cases": cases,
        "dof_rendering": dof_rendering_summary(render_context),
    }
    return _render_html_browser_report(_jsonable(output))


def main() -> int:
    args = _parse_args()
    config: MolGRConfig | None = None
    input_cases = _collect_input_cases(args)
    render_context = _make_dof_render_context(args)
    cases = trace_reconstruction_cases(
        input_cases,
        score_all_candidates=not args.no_score_all_candidates,
        render_context=render_context,
        config=config,
    )

    output = {
        "input": {
            "source": (
                "review_fixtures" if args.review_fixtures_manifest is not None else "generic"
            ),
            "ids": [case.id for case in input_cases],
            "total_charge": args.total_charge,
            "total_radical_electrons": _direct_total_radicals(args),
        },
        "case_count": len(cases),
        "cases": cases,
    }
    if render_context is not None:
        output["dof_rendering"] = dof_rendering_summary(render_context)
    output_format = _resolve_output_format(args)
    if output_format == "json":
        output_text = json.dumps(
            _jsonable(output),
            ensure_ascii=False,
            indent=None if args.indent == 0 else args.indent,
            allow_nan=False,
        )
    else:
        output_text = _render_html_browser_report(_jsonable(output))

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
