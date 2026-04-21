from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from molgr.fallback import xyz2omol
from molgr.fallback.pipeline import reconstruct_without_metals as no_metal_module
from molgr.fallback.pipeline import resonance as resonance_module
from molgr.fallback.pipeline.reconstruct_without_metals import xyz_to_omol_no_metal_state
from molgr.fallback.utils.scoring import omol_score_cache_clear
from molgr.utils.converter import pybel_to_rdmol
from scripts.molgr_cases_molfile import load_molfile_cases


_DEFAULT_BUILD_PROCESSED_RESONANCE_KEY = resonance_module.build_processed_resonance_key


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Profile a fallback case repeatedly. "
            "Default mode clears no-metal and scoring caches before each top-level run "
            "so py-spy reflects cold outer runs while preserving intra-run cache reuse."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("tests/data/sdf/MoNNMo.sdf"),
        help="Molfile/SDF path to profile.",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=20,
        help="How many measured fallback runs to execute.",
    )
    parser.add_argument(
        "--warmup",
        type=int,
        default=0,
        help="How many unmeasured runs to execute first.",
    )
    parser.add_argument(
        "--cache-mode",
        choices=("cold-outer", "warm"),
        default="cold-outer",
        help=(
            "Cache handling policy. "
            "'cold-outer' clears no-metal/scoring caches before each top-level run. "
            "'warm' leaves caches hot across repeated runs."
        ),
    )
    parser.add_argument(
        "--include-rdkit",
        action="store_true",
        help="Also convert the Pybel molecule to RDKit after reconstruction.",
    )
    parser.add_argument(
        "--resonance-key-mode",
        choices=("current", "write-mol"),
        default="current",
        help=(
            "Resonance processed-state key implementation to profile. "
            "'current' keeps build_processed_resonance_key as-is. "
            "'write-mol' replaces it with omol.write(\"mol\")."
        ),
    )
    return parser.parse_args()


def _load_case(input_path: Path) -> dict[str, object]:
    cases = load_molfile_cases(input_path=input_path, limit=1)
    if not cases:
        raise ValueError(f"No cases loaded from {input_path}")
    case = cases[0]
    if case.get("provider_error"):
        raise ValueError(f"Case provider failed: {case['provider_error']}")
    xyz_block = case.get("xyz_block")
    if not isinstance(xyz_block, str) or not xyz_block.strip():
        raise ValueError("Loaded case is missing xyz_block")
    return case


def _clear_runtime_caches() -> None:
    xyz_to_omol_no_metal_state.cache_clear()
    omol_score_cache_clear()


def _configure_resonance_key_mode(mode: str) -> None:
    def _write_mol_key(omol: object) -> str:
        return omol.write("mol")

    if mode == "current":
        build_key = _DEFAULT_BUILD_PROCESSED_RESONANCE_KEY
    elif mode == "write-mol":
        build_key = _write_mol_key
    else:
        raise ValueError(f"Unsupported resonance key mode: {mode}")

    no_metal_module.build_processed_resonance_key = build_key
    resonance_module.build_processed_resonance_key = build_key


def _run_once(case: dict[str, object], include_rdkit: bool) -> None:
    xyz_block = case["xyz_block"]
    total_charge = case["total_charge"]
    total_radical_electrons = case["total_radical_electrons"]
    if not isinstance(xyz_block, str):
        raise TypeError("xyz_block must be a string")
    if not isinstance(total_charge, int):
        raise TypeError("total_charge must be an int")
    if not isinstance(total_radical_electrons, int):
        raise TypeError("total_radical_electrons must be an int")

    omol = xyz2omol(
        xyz_block,
        total_charge=total_charge,
        total_radical_electrons=total_radical_electrons,
    )
    if omol is None:
        raise ValueError("fallback.xyz2omol returned None")
    if include_rdkit:
        rdmol = pybel_to_rdmol(omol)
        if rdmol is None:
            raise ValueError("pybel_to_rdmol returned None")


def main() -> int:
    args = _parse_args()
    _configure_resonance_key_mode(args.resonance_key_mode)
    case = _load_case(args.input)

    for _ in range(max(args.warmup, 0)):
        if args.cache_mode == "cold-outer":
            _clear_runtime_caches()
        _run_once(case, include_rdkit=args.include_rdkit)

    started = time.perf_counter()
    for _ in range(max(args.repeat, 0)):
        if args.cache_mode == "cold-outer":
            _clear_runtime_caches()
        _run_once(case, include_rdkit=args.include_rdkit)
    elapsed = time.perf_counter() - started

    print(
        "profiled_case=1 "
        f"input={args.input} repeat={args.repeat} warmup={args.warmup} "
        f"cache_mode={args.cache_mode} "
        f"resonance_key_mode={args.resonance_key_mode} "
        f"include_rdkit={args.include_rdkit} total_s={elapsed:.6f} "
        f"avg_ms={(elapsed / max(args.repeat, 1)) * 1000.0:.3f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
