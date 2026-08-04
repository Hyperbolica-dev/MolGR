from __future__ import annotations

import time
from dataclasses import dataclass
from importlib import import_module
from typing import Any

from rdkit import Chem

from benchmarks.smiles_xyz_benchmark.methods.base import BenchmarkMethod, MethodRunOutput
from benchmarks.smiles_xyz_benchmark.methods.postprocess import remove_hs_without_sanitize


@dataclass(frozen=True)
class MolGRFallbackMethod(BenchmarkMethod):
    method_id: str = "molgr_fallback"

    def run(self, case: dict[str, Any]) -> MethodRunOutput:
        timing_ms_breakdown: dict[str, float] = {
            "molgr_interface_ms": 0.0,
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
            xyz_to_rdmol = import_module("molgr.interface").xyz_to_rdmol
        except Exception as exc:  # noqa: BLE001
            return MethodRunOutput(
                status="error",
                error=f"import molgr interface failed: {exc}",
                timing_ms_breakdown=timing_ms_breakdown,
            )

        interface_started = time.perf_counter()
        try:
            rdkit_mol = xyz_to_rdmol(
                xyz_block,
                total_charge,
                total_radical_electrons + 1,
                backend="python",
                make_dative_bonds=True,
            )
        except Exception as exc:  # noqa: BLE001
            timing_ms_breakdown["molgr_interface_ms"] = (
                time.perf_counter() - interface_started
            ) * 1000.0
            return MethodRunOutput(
                status="error",
                error=f"molgr.interface.xyz_to_rdmol failed: {exc}",
                timing_ms_breakdown=timing_ms_breakdown,
            )
        timing_ms_breakdown["molgr_interface_ms"] = (
            time.perf_counter() - interface_started
        ) * 1000.0
        postprocess_started = time.perf_counter()
        try:
            rdkit_mol = remove_hs_without_sanitize(rdkit_mol)
            predicted_smiles = Chem.MolToSmiles(
                rdkit_mol,
                canonical=True,
                isomericSmiles=True,
            )
        except Exception as exc:  # noqa: BLE001
            timing_ms_breakdown["postprocess_ms"] = (
                time.perf_counter() - postprocess_started
            ) * 1000.0
            return MethodRunOutput(
                status="error",
                error=f"MolGR result postprocess failed: {exc}",
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
