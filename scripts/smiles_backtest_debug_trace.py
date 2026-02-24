from __future__ import annotations

import argparse
from pathlib import Path


try:
    from scripts.molgr_cases_smiles_csv import load_smiles_csv_cases
    from scripts.molgr_trace_runner import run_trace_cases
except ModuleNotFoundError as exc:
    if exc.name != "scripts":
        raise
    from molgr_cases_smiles_csv import load_smiles_csv_cases
    from molgr_trace_runner import run_trace_cases


DEFAULT_INPUT = Path("tests/test_cases.csv")
DEFAULT_OUT_ROOT = Path(".molgr_backtest_smiles_snapshots_full")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run SMILES backtest with fallback runtime tracing"
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--no-chirality", action="store_true")
    parser.add_argument("--max-resonance", type=int, default=100)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    cases = load_smiles_csv_cases(input_path=args.input, limit=args.limit)
    if not cases and (args.limit is None or args.limit > 0):
        cases = [
            {
                "case_idx": 1,
                "input_smiles": "",
                "xyz_block": None,
                "total_charge": None,
                "total_radical_electrons": None,
                "provider_error": "RDKit failed to parse SMILES",
            }
        ]
    return run_trace_cases(
        cases=cases,
        out_root=args.out_root,
        use_chirality=not args.no_chirality,
        max_resonance=args.max_resonance,
    )


if __name__ == "__main__":
    raise SystemExit(main())
