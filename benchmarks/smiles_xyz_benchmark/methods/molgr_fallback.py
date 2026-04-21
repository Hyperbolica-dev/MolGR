from __future__ import annotations

import time
from dataclasses import dataclass
from importlib import import_module
from typing import Any

from benchmarks.smiles_xyz_benchmark.methods.base import BenchmarkMethod, MethodRunOutput
from benchmarks.smiles_xyz_benchmark.methods.postprocess import finalize_rdmol_with_dative_bonds


@dataclass(frozen=True)
class MolGRFallbackMethod(BenchmarkMethod):
    method_id: str = "molgr_fallback"

    def run(self, case: dict[str, Any]) -> MethodRunOutput:
        timing_ms_breakdown: dict[str, float] = {
            "fallback_xyz2omol_ms": 0.0,
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

        total_charge_value = case.get("total_charge", 0)
        try:
            total_charge = int(total_charge_value)
        except (TypeError, ValueError):
            return MethodRunOutput(
                status="error",
                error=f"invalid total_charge: {total_charge_value!r}",
                timing_ms_breakdown=timing_ms_breakdown,
            )

        total_radical_electrons_value = case.get("total_radical_electrons", 0)
        try:
            total_radical_electrons = int(total_radical_electrons_value)
        except (TypeError, ValueError):
            return MethodRunOutput(
                status="error",
                error=f"invalid total_radical_electrons: {total_radical_electrons_value!r}",
                timing_ms_breakdown=timing_ms_breakdown,
            )

        try:
            fallback = import_module("molgr.fallback")
        except Exception as exc:  # noqa: BLE001
            return MethodRunOutput(
                status="error",
                error=f"import molgr.fallback failed: {exc}",
                timing_ms_breakdown=timing_ms_breakdown,
            )

        fallback_started = time.perf_counter()
        try:
            omol = fallback.xyz2omol(
                xyz_block,
                total_charge=total_charge,
                total_radical_electrons=total_radical_electrons,
            )
        except Exception as exc:  # noqa: BLE001
            timing_ms_breakdown["fallback_xyz2omol_ms"] = (
                time.perf_counter() - fallback_started
            ) * 1000.0
            return MethodRunOutput(
                status="error",
                error=f"fallback.xyz2omol failed: {exc}",
                timing_ms_breakdown=timing_ms_breakdown,
            )
        timing_ms_breakdown["fallback_xyz2omol_ms"] = (
            time.perf_counter() - fallback_started
        ) * 1000.0

        if omol is None:
            return MethodRunOutput(
                status="error",
                error="fallback.xyz2omol returned None",
                timing_ms_breakdown=timing_ms_breakdown,
            )

        pybel_to_rdmol_started = time.perf_counter()
        try:
            pybel_to_rdmol = import_module("molgr.interface").pybel_to_rdmol
            rdkit_mol = pybel_to_rdmol(omol)
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

        if rdkit_mol is None:
            return MethodRunOutput(
                status="error",
                error="pybel_to_rdmol returned None",
                timing_ms_breakdown=timing_ms_breakdown,
            )
        postprocess_started = time.perf_counter()
        try:
            rdkit_mol, predicted_smiles = finalize_rdmol_with_dative_bonds(rdkit_mol)
        except Exception as exc:  # noqa: BLE001
            timing_ms_breakdown["postprocess_ms"] = (
                time.perf_counter() - postprocess_started
            ) * 1000.0
            return MethodRunOutput(
                status="error",
                error=f"MolToSmiles failed: {exc}",
                rdkit_mol=rdkit_mol,
                timing_ms_breakdown=timing_ms_breakdown,
            )
        timing_ms_breakdown["postprocess_ms"] = (time.perf_counter() - postprocess_started) * 1000.0

        return MethodRunOutput(
            status="ok",
            error=None,
            predicted_smiles=predicted_smiles,
            rdkit_mol=rdkit_mol,
            timing_ms_breakdown=timing_ms_breakdown,
        )
