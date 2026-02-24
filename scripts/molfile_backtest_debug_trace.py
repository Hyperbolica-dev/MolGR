"""
Author: TMJ
Date: 2026-02-24 11:52:30
LastEditors: TMJ
LastEditTime: 2026-02-24 20:33:27
Description: 请填写简介
"""

from __future__ import annotations

import argparse
from pathlib import Path


try:
    from scripts.molgr_cases_molfile import load_molfile_cases
    from scripts.molgr_trace_runner import run_trace_cases
except ModuleNotFoundError as exc:
    if exc.name != "scripts":
        raise
    from molgr_cases_molfile import load_molfile_cases
    from molgr_trace_runner import run_trace_cases


DEFAULT_OUT_ROOT = Path(".molgr_backtest_molfile_snapshots")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run molfile backtest with fallback runtime tracing"
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument(
        "--no-chirality",
        action="store_true",
        help="Disable chirality in equivalence checks",
    )
    parser.add_argument(
        "--max-resonance",
        type=int,
        default=100,
        help="Maximum resonance structures for equivalence checks",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    cases = load_molfile_cases(input_path=args.input, limit=args.limit)
    return run_trace_cases(
        cases=cases,
        out_root=args.out_root,
        use_chirality=not args.no_chirality,
        max_resonance=args.max_resonance,
    )


if __name__ == "__main__":
    raise SystemExit(main())
