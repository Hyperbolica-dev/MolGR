from __future__ import annotations

import argparse
import csv
import statistics
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from rdkit import Chem
from rdkit.Chem import rdDistGeom

from molgr import _core

try:
    from openbabel import openbabel as ob
except ImportError as exc:  # pragma: no cover - dev helper
    raise SystemExit("openbabel is required for stage profiling") from exc


_PIPELINE: Any = _core.pipeline
_DEV_STAGES: Any = _core.dev.stages
_DEV_RESONANCE: Any = _core.dev.pipeline.resonance


@dataclass(frozen=True)
class PreparedCase:
    smiles: str
    xyz_block: str
    total_charge: int
    total_radical_electrons: int


@dataclass(frozen=True)
class LinearStageInput:
    name: str
    mol: ob.OBMol
    given_charge: int
    runner: Callable[[ob.OBMol, int, int], int]


def _get_ptr(obmol: ob.OBMol) -> int:
    this = getattr(obmol, "this", None)
    if this is not None:
        return int(this)  # type: ignore[arg-type]
    return int(obmol)  # type: ignore[arg-type]


def _load_case_smiles(case_index: int) -> str:
    csv_path = Path("tests/test_cases.csv")
    with csv_path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if case_index < 1 or case_index > len(rows):
        raise SystemExit(f"case index out of range: {case_index}")
    return rows[case_index - 1]["smiles"]


def _prepare_case(smiles: str) -> PreparedCase:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SystemExit(f"failed to parse smiles: {smiles}")

    mol_h = Chem.AddHs(mol)
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    embed_code = rdDistGeom.EmbedMolecule(mol_h, params)
    if int(embed_code) != 0:
        raise SystemExit(f"rdkit embed failed for smiles: {smiles}")

    xyz_block = Chem.MolToXYZBlock(mol_h)
    total_charge = sum(int(atom.GetFormalCharge()) for atom in mol_h.GetAtoms())
    total_radical_electrons = sum(int(atom.GetNumRadicalElectrons()) for atom in mol_h.GetAtoms())
    return PreparedCase(
        smiles=smiles,
        xyz_block=xyz_block,
        total_charge=total_charge,
        total_radical_electrons=total_radical_electrons,
    )


def _read_xyz_obmol(xyz_block: str) -> ob.OBMol:
    conv = ob.OBConversion()
    if not conv.SetInFormat("xyz"):
        raise SystemExit("failed to set OpenBabel xyz input format")
    mol = ob.OBMol()
    if not conv.ReadString(mol, xyz_block):
        raise SystemExit("failed to parse xyz block with OpenBabel")
    return mol


def _formal_charge_sum(mol: ob.OBMol) -> int:
    return sum(atom.GetFormalCharge() for atom in ob.OBMolAtomIter(mol))


def _copy_mol(mol: ob.OBMol) -> ob.OBMol:
    return ob.OBMol(mol)


def _build_linear_stage_inputs(case: PreparedCase) -> tuple[list[LinearStageInput], ob.OBMol, int]:
    current = _read_xyz_obmol(case.xyz_block)
    given_charge = 0
    total_radical_electrons = case.total_radical_electrons
    stage_inputs: list[LinearStageInput] = []

    def add_stage(
        name: str,
        runner: Callable[[ob.OBMol, int, int], int],
    ) -> None:
        nonlocal current, given_charge
        stage_inputs.append(
            LinearStageInput(
                name=name,
                mol=_copy_mol(current),
                given_charge=given_charge,
                runner=runner,
            )
        )
        given_charge = runner(current, given_charge, total_radical_electrons)

    add_stage(
        "make_connections",
        lambda mol, gc, tr: (
            _DEV_STAGES.preprocess.make_connections_ptr(_get_ptr(mol), 1.4),
            gc,
        )[1],
    )
    add_stage(
        "pre_clean",
        lambda mol, gc, tr: (_DEV_STAGES.preprocess.pre_clean_ptr(_get_ptr(mol)), gc)[1],
    )
    add_stage(
        "fresh_omol_charge_radical_initial",
        lambda mol, gc, tr: (_DEV_STAGES.fresh.fresh_omol_charge_radical_ptr(_get_ptr(mol)), gc)[1],
    )

    given_charge = case.total_charge - _formal_charge_sum(current)

    add_stage(
        "eliminate_NNN_negative",
        lambda mol, gc, tr: _DEV_STAGES.eliminate.eliminate_nnn_ptr(_get_ptr(mol), gc, False)[0],
    )
    add_stage(
        "eliminate_high_positive_charge_atoms",
        lambda mol, gc, tr: _DEV_STAGES.eliminate.eliminate_high_positive_charge_atoms_ptr(
            _get_ptr(mol), gc
        )[0],
    )
    add_stage(
        "eliminate_CN_in_doubt",
        lambda mol, gc, tr: _DEV_STAGES.eliminate.eliminate_cn_in_doubt_ptr(_get_ptr(mol), gc)[0],
    )
    add_stage(
        "eliminate_NNN_positive",
        lambda mol, gc, tr: _DEV_STAGES.eliminate.eliminate_nnn_ptr(_get_ptr(mol), gc, True)[0],
    )
    add_stage(
        "eliminate_carboxyl",
        lambda mol, gc, tr: _DEV_STAGES.eliminate.eliminate_carboxyl_ptr(_get_ptr(mol), gc)[0],
    )
    add_stage(
        "clean_carbene_neighbor_unsaturated_first",
        lambda mol, gc, tr: (_DEV_STAGES.clean.clean_carbene_neighbor_unsaturated_ptr(_get_ptr(mol)), gc)[
            1
        ],
    )
    add_stage(
        "eliminate_carbene_neighbor_heteroatom",
        lambda mol, gc, tr: _DEV_STAGES.eliminate.eliminate_carbene_neighbor_heteroatom_ptr(
            _get_ptr(mol), gc
        )[0],
    )
    add_stage(
        "clean_neighbor_radicals",
        lambda mol, gc, tr: (_DEV_STAGES.clean.clean_neighbor_radicals_ptr(_get_ptr(mol)), gc)[1],
    )
    add_stage(
        "clean_carbene_neighbor_unsaturated_second",
        lambda mol, gc, tr: (_DEV_STAGES.clean.clean_carbene_neighbor_unsaturated_ptr(_get_ptr(mol)), gc)[
            1
        ],
    )
    add_stage(
        "eliminate_charge_spliting",
        lambda mol, gc, tr: _DEV_STAGES.eliminate.eliminate_charge_spliting_ptr(_get_ptr(mol), gc)[
            0
        ],
    )
    add_stage(
        "break_deformed_ene",
        lambda mol, gc, tr: (
            _DEV_STAGES.break_bond.break_deformed_ene_ptr(_get_ptr(mol), gc, tr, 5.0),
            gc,
        )[1],
    )
    add_stage(
        "break_one_bond",
        lambda mol, gc, tr: _DEV_STAGES.break_bond.break_one_bond_ptr(_get_ptr(mol), gc, tr)[0],
    )
    add_stage(
        "fresh_omol_charge_radical_final",
        lambda mol, gc, tr: (_DEV_STAGES.fresh.fresh_omol_charge_radical_ptr(_get_ptr(mol)), gc)[1],
    )

    return stage_inputs, current, given_charge


def _run_pipeline(case: PreparedCase, repeat: int, warmup: int) -> None:
    for _ in range(max(warmup, 0)):
        _PIPELINE.reconstruct_without_metals.xyz_to_omol_no_metal(
            case.xyz_block,
            case.total_charge,
            case.total_radical_electrons,
        )

    elapsed_ms: list[float] = []
    no_metal_ms: list[float] = []
    resonance_ms: list[float] = []
    for _ in range(max(repeat, 0)):
        started = time.perf_counter()
        result = _PIPELINE.reconstruct_without_metals.xyz_to_omol_no_metal(
            case.xyz_block,
            case.total_charge,
            case.total_radical_electrons,
        )
        elapsed_ms.append((time.perf_counter() - started) * 1000.0)
        if result is None:
            raise SystemExit("cpp no-metal pipeline returned None")
        breakdown = dict(_PIPELINE.get_last_run_timing_breakdown_ms())
        no_metal_ms.append(float(breakdown.get("no_metal_pipeline_ms", 0.0)))
        resonance_ms.append(float(breakdown.get("resonance_handling_enumeration_ms", 0.0)))

    print(f"smiles={case.smiles}")
    print(f"repeat={repeat} warmup={warmup}")
    print(f"elapsed_avg_ms={statistics.mean(elapsed_ms):.6f}")
    print(f"elapsed_min_ms={min(elapsed_ms):.6f}")
    print(f"no_metal_reducer_avg_ms={statistics.mean(no_metal_ms):.6f}")
    print(f"resonance_reducer_avg_ms={statistics.mean(resonance_ms):.6f}")


def _benchmark_stage(stage: LinearStageInput, total_radical_electrons: int, repeat: int) -> float:
    samples_ms: list[float] = []
    for _ in range(max(repeat, 0)):
        mol = _copy_mol(stage.mol)
        started = time.perf_counter()
        stage.runner(mol, stage.given_charge, total_radical_electrons)
        samples_ms.append((time.perf_counter() - started) * 1000.0)
    return statistics.mean(samples_ms)


def _benchmark_resonance(final_mol: ob.OBMol, given_charge: int, repeat: int) -> dict[str, float]:
    get_resonances_ms: list[float] = []
    process_resonance_ms: list[float] = []
    smiles_token_ms: list[float] = []

    for _ in range(max(repeat, 0)):
        mol = _copy_mol(final_mol)
        started = time.perf_counter()
        ptrs = _DEV_RESONANCE.get_radical_resonances_ptr(_get_ptr(mol))
        get_resonances_ms.append((time.perf_counter() - started) * 1000.0)
        for ptr in ptrs:
            _core.free_obmol_ptr(ptr)

    for _ in range(max(repeat, 0)):
        mol = _copy_mol(final_mol)
        started = time.perf_counter()
        out_ptr, _ = _DEV_RESONANCE.process_resonance_ptr(_get_ptr(mol), given_charge)
        process_resonance_ms.append((time.perf_counter() - started) * 1000.0)
        _core.free_obmol_ptr(out_ptr)

    for _ in range(max(repeat, 0)):
        mol = _copy_mol(final_mol)
        started = time.perf_counter()
        _DEV_RESONANCE.smiles_token_ptr(_get_ptr(mol))
        smiles_token_ms.append((time.perf_counter() - started) * 1000.0)

    return {
        "get_radical_resonances_avg_ms": statistics.mean(get_resonances_ms),
        "process_resonance_avg_ms": statistics.mean(process_resonance_ms),
        "smiles_token_avg_ms": statistics.mean(smiles_token_ms),
    }


def _run_stage_bench(case: PreparedCase, repeat: int) -> None:
    stage_inputs, final_mol, final_given_charge = _build_linear_stage_inputs(case)
    totals: list[tuple[str, float]] = [
        (stage.name, _benchmark_stage(stage, case.total_radical_electrons, repeat))
        for stage in stage_inputs
    ]
    stage_total_ms = sum(value for _, value in totals)

    print(f"smiles={case.smiles}")
    print(f"stage_repeat={repeat}")
    print("linear_stage_avg_ms")
    for name, value in totals:
        share = 0.0 if stage_total_ms == 0.0 else value / stage_total_ms * 100.0
        print(f"{name},{value:.6f},{share:.2f}%")

    resonance_metrics = _benchmark_resonance(final_mol, final_given_charge, repeat)
    print("resonance_helper_avg_ms")
    for name, value in resonance_metrics.items():
        print(f"{name},{value:.6f}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Profile C++ no-metal reconstruction hotspots.")
    parser.add_argument(
        "--mode",
        choices=("pipeline", "stages"),
        default="pipeline",
        help="pipeline: whole no-metal run, suitable for py-spy; stages: stage-level microbench.",
    )
    parser.add_argument(
        "--case-index",
        type=int,
        default=52,
        help="1-based row index in tests/test_cases.csv.",
    )
    parser.add_argument(
        "--smiles",
        type=str,
        default=None,
        help="Optional SMILES override.",
    )
    parser.add_argument("--repeat", type=int, default=200)
    parser.add_argument("--warmup", type=int, default=10)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    smiles = args.smiles if args.smiles is not None else _load_case_smiles(args.case_index)
    case = _prepare_case(smiles)

    if args.mode == "pipeline":
        _run_pipeline(case, repeat=args.repeat, warmup=args.warmup)
    else:
        _run_stage_bench(case, repeat=args.repeat)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
