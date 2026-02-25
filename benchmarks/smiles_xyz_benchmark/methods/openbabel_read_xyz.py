from __future__ import annotations

import time
from dataclasses import dataclass
from importlib import import_module
from typing import Any

from openbabel import pybel
from rdkit import Chem

from benchmarks.smiles_xyz_benchmark.methods.base import BenchmarkMethod, MethodRunOutput


@dataclass(frozen=True)
class OpenBabelReadXYZMethod(BenchmarkMethod):
    method_id: str = "openbabel_read_xyz"

    def run(self, case: dict[str, Any]) -> MethodRunOutput:
        timing_ms_breakdown: dict[str, float] = {
            "pybel_read_xyz_ms": 0.0,
            "pybel_to_rdmol_ms": 0.0,
            "postprocess_ms": 0.0,
        }

        xyz_block = case.get("xyz_block")
        if not isinstance(xyz_block, str) or not xyz_block.strip():
            return MethodRunOutput(
                status="error",
                error="missing or empty xyz_block",
                timing_ms_breakdown=timing_ms_breakdown,
            )

        pybel_read_started = time.perf_counter()
        try:
            omol = pybel.readstring("xyz", xyz_block)
        except Exception as exc:  # noqa: BLE001
            timing_ms_breakdown["pybel_read_xyz_ms"] = (
                time.perf_counter() - pybel_read_started
            ) * 1000.0
            return MethodRunOutput(
                status="error",
                error=f"pybel.readstring failed: {exc}",
                timing_ms_breakdown=timing_ms_breakdown,
            )
        timing_ms_breakdown["pybel_read_xyz_ms"] = (
            time.perf_counter() - pybel_read_started
        ) * 1000.0

        if omol is None or omol.OBMol.NumAtoms() == 0:
            return MethodRunOutput(
                status="error",
                error="pybel.readstring returned empty molecule",
                timing_ms_breakdown=timing_ms_breakdown,
            )

        pybel_to_rdmol_started = time.perf_counter()
        try:
            pybel_to_rdmol = getattr(import_module("molgr.interface"), "pybel_to_rdmol")
            rdmol = pybel_to_rdmol(omol)
            rdmol = Chem.RemoveHs(rdmol)
        except Exception as exc:  # noqa: BLE001
            timing_ms_breakdown["pybel_to_rdmol_ms"] = (
                time.perf_counter() - pybel_to_rdmol_started
            ) * 1000.0
            return MethodRunOutput(
                status="error",
                error=f"pybel_to_rdmol failed: {exc}",
                timing_ms_breakdown=timing_ms_breakdown,
            )
        timing_ms_breakdown["pybel_to_rdmol_ms"] = (
            time.perf_counter() - pybel_to_rdmol_started
        ) * 1000.0

        warnings: list[str] = []
        predicted_smiles: str | None = None
        postprocess_started = time.perf_counter()
        try:
            predicted_smiles = Chem.MolToSmiles(rdmol, canonical=True, isomericSmiles=True)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"MolToSmiles failed: {exc}")
        timing_ms_breakdown["postprocess_ms"] = (time.perf_counter() - postprocess_started) * 1000.0

        return MethodRunOutput(
            status="ok",
            error="; ".join(warnings) if warnings else None,
            predicted_smiles=predicted_smiles,
            rdkit_mol=rdmol,
            timing_ms_breakdown=timing_ms_breakdown,
        )
