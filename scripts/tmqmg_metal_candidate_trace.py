#!/usr/bin/env python3
# pyright: reportCallIssue=false
"""Trace tmQMg reconstruction by dataset id.

This is a tmQMg input adapter around scripts/reconstruction_trace.py. The
generic script owns the reconstruction trajectory tracing and report rendering;
this wrapper only resolves CSV rows, XYZ paths, charges, and radical counts.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any, Optional, Sequence


if __package__ in (None, ""):
    project_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(project_root))
    sys.path.insert(0, str(project_root / "src"))

from rdkit import Chem, RDLogger

from molgr.config import MolGRConfig
from scripts.reconstruction_trace import (
    DofRenderContext,
    TraceInputCase,
    _jsonable,
    _make_dof_render_context,
    _parse_size,
    _render_html_browser_report,
    _resolve_output_format,
    _with_inline_dof_svgs,
    dof_rendering_summary,
    split_repeated_values,
    trace_reconstruction_case,
)


RDLogger.DisableLog("rdApp.*")  # type: ignore[arg-type]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csv",
        type=Path,
        required=True,
        help="tmQMg metadata CSV path. Rows are selected by the 'id' column.",
    )
    parser.add_argument(
        "--xyz-dir",
        type=Path,
        required=True,
        help="Directory containing tmQMg XYZ files named <id>.xyz.",
    )
    parser.add_argument(
        "--id",
        dest="ids",
        action="append",
        required=True,
        help="tmQMg id to trace. May be repeated or comma-separated.",
    )
    parser.add_argument(
        "--spin-source",
        choices=("closed_shell", "reference_smiles"),
        default="closed_shell",
        help=(
            "How to choose total radical electrons when --total-radical-electrons is not set. "
            "Default matches the closed-shell tmQMg regression workflow."
        ),
    )
    parser.add_argument(
        "--total-radical-electrons",
        type=int,
        default=None,
        help="Explicit radical-electron count to use for every requested id.",
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


def _load_rows_by_id(path: Path) -> dict[str, tuple[int, dict[str, str]]]:
    rows: dict[str, tuple[int, dict[str, str]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "id" not in reader.fieldnames:
            raise ValueError(f"tmQMg CSV must contain an 'id' column: {path}")
        for row_index, row in enumerate(reader, start=1):
            case_id = row.get("id", "").strip()
            if case_id and case_id not in rows:
                rows[case_id] = (row_index, row)
    return rows


def _parse_int_field(row: dict[str, str], field_name: str) -> int:
    raw_value = row.get(field_name, "").strip()
    if raw_value == "":
        raise ValueError(
            f"Missing required integer field {field_name!r} for id={row.get('id', '')}"
        )
    return int(raw_value)


def _reference_smiles(row: dict[str, str]) -> str:
    for field_name in ("smiles", "SMILES", "canonical_smiles", "CanonicalSMILES"):
        value = row.get(field_name, "").strip()
        if value:
            return value
    raise ValueError(f"No reference SMILES column found for id={row.get('id', '')}")


def _radicals_from_reference_smiles(row: dict[str, str]) -> int:
    smiles = _reference_smiles(row)
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        raise ValueError(f"Could not parse reference SMILES for id={row.get('id', '')}: {smiles}")
    return int(sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()))


def _resolve_total_radicals(
    row: dict[str, str],
    *,
    spin_source: str,
    explicit_total_radicals: Optional[int],
) -> tuple[int, str]:
    if explicit_total_radicals is not None:
        return explicit_total_radicals, "explicit"
    if spin_source == "closed_shell":
        return 0, "closed_shell"
    if spin_source == "reference_smiles":
        return _radicals_from_reference_smiles(row), "reference_smiles"
    raise ValueError(f"Unsupported spin source: {spin_source!r}")


def _trace_tmqmg_case(
    row_index: int,
    row: dict[str, str],
    *,
    xyz_dir: Path,
    spin_source: str,
    explicit_total_radicals: Optional[int],
    score_all_candidates: bool,
    render_context: DofRenderContext | None,
    config: MolGRConfig | None = None,
) -> dict[str, Any]:
    started = time.perf_counter()
    case_id = row.get("id", "").strip()
    xyz_path = xyz_dir / f"{case_id}.xyz"
    total_charge = _parse_int_field(row, "charge")
    total_radicals, resolved_spin_source = _resolve_total_radicals(
        row,
        spin_source=spin_source,
        explicit_total_radicals=explicit_total_radicals,
    )
    input_case = TraceInputCase(
        id=case_id,
        xyz_block=xyz_path.read_text(encoding="utf-8"),
        total_charge=total_charge,
        total_radical_electrons=total_radicals,
        xyz_path=xyz_path,
        xyz_source="tmqmg_xyz_path",
    )
    trace = trace_reconstruction_case(
        input_case,
        score_all_candidates=score_all_candidates,
        render_context=render_context,
        config=config,
    )
    trace.update(
        {
            "row_index": row_index,
            "reference_smiles": row.get("smiles", ""),
            "spin_source": resolved_spin_source,
            "elapsed_seconds": time.perf_counter() - started,
        }
    )
    return trace


def _split_svg_fragment_lines(svg_fragment: str) -> list[str]:
    """Split embedded SVG at tag boundaries so HTML source avoids huge lines."""

    if not svg_fragment:
        return []
    return svg_fragment.replace("><", ">\n<").splitlines()


def _split_embedded_svg_fragments(value: Any) -> Any:
    if isinstance(value, dict):
        copied = {key: _split_embedded_svg_fragments(item) for key, item in value.items()}
        svg_fragment = copied.get("svg_fragment")
        if isinstance(svg_fragment, str):
            copied["svg_fragment"] = _split_svg_fragment_lines(svg_fragment)
        return copied
    if isinstance(value, list):
        return [_split_embedded_svg_fragments(item) for item in value]
    return value


def _html_script_json_pretty(value: Any) -> str:
    return (
        json.dumps(
            _jsonable(value),
            ensure_ascii=False,
            allow_nan=False,
            indent=2,
        )
        .replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
    )


def _replace_trace_json_script(html_text: str, script_json: str) -> str:
    open_tag = '<script id="trace-json" type="application/json">'
    close_tag = "</script>"
    open_index = html_text.find(open_tag)
    if open_index < 0:
        return html_text
    content_start = open_index + len(open_tag)
    close_index = html_text.find(close_tag, content_start)
    if close_index < 0:
        return html_text
    return html_text[:content_start] + "\n" + script_json + "\n  " + html_text[close_index:]


def _support_multiline_svg_fragment_arrays(html_text: str) -> str:
    return html_text.replace(
        '        return image.svg_fragment || "";',
        "        if (Array.isArray(image.svg_fragment)) {\n"
        '          return image.svg_fragment.join("\\n");\n'
        "        }\n"
        '        return image.svg_fragment || "";',
        1,
    )


def _render_tmqmg_html_report(output: dict[str, Any]) -> str:
    html_text = _render_html_browser_report(_jsonable(output))
    inline_output = _with_inline_dof_svgs(_jsonable(output))
    readable_output = _split_embedded_svg_fragments(inline_output)
    html_text = _replace_trace_json_script(
        html_text,
        _html_script_json_pretty(readable_output),
    )
    return _support_multiline_svg_fragment_arrays(html_text)


def _write_output(args: argparse.Namespace, output: dict[str, Any]) -> None:
    output_format = _resolve_output_format(args)
    if output_format == "json":
        output_text = json.dumps(
            _jsonable(output),
            ensure_ascii=False,
            indent=None if args.indent == 0 else args.indent,
            allow_nan=False,
        )
    else:
        output_text = _render_tmqmg_html_report(output)

    if args.out is None:
        print(output_text, end="" if output_text.endswith("\n") else "\n")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(
            output_text if output_text.endswith("\n") else output_text + "\n",
            encoding="utf-8",
        )


def _trace_requested_rows(
    case_ids: Sequence[str],
    rows_by_id: dict[str, tuple[int, dict[str, str]]],
    *,
    args: argparse.Namespace,
    render_context: DofRenderContext | None,
    config: MolGRConfig | None,
) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for case_id in case_ids:
        started = time.perf_counter()
        row_entry = rows_by_id.get(case_id)
        if row_entry is None:
            cases.append(
                {
                    "id": case_id,
                    "status": "missing_csv_row",
                    "trace_kind": "unknown",
                    "elapsed_seconds": time.perf_counter() - started,
                }
            )
            continue

        row_index, row = row_entry
        try:
            cases.append(
                _trace_tmqmg_case(
                    row_index,
                    row,
                    xyz_dir=args.xyz_dir,
                    spin_source=args.spin_source,
                    explicit_total_radicals=args.total_radical_electrons,
                    score_all_candidates=not args.no_score_all_candidates,
                    render_context=render_context,
                    config=config,
                )
            )
        except Exception as exc:
            cases.append(
                {
                    "row_index": row_index,
                    "id": case_id,
                    "status": "error",
                    "trace_kind": "unknown",
                    "error": f"{type(exc).__name__}: {exc}",
                    "elapsed_seconds": time.perf_counter() - started,
                }
            )
    return cases


def main() -> int:
    args = _parse_args()
    case_ids = split_repeated_values(args.ids)
    rows_by_id = _load_rows_by_id(args.csv)
    config: MolGRConfig | None = None
    render_context = _make_dof_render_context(args)
    cases = _trace_requested_rows(
        case_ids,
        rows_by_id,
        args=args,
        render_context=render_context,
        config=config,
    )

    output = {
        "input": {
            "source": "tmQMg",
            "csv": str(args.csv),
            "xyz_dir": str(args.xyz_dir),
            "ids": case_ids,
            "spin_source": args.spin_source,
            "total_radical_electrons": args.total_radical_electrons,
        },
        "case_count": len(cases),
        "cases": cases,
    }
    if render_context is not None:
        output["dof_rendering"] = dof_rendering_summary(render_context)
    _write_output(args, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
