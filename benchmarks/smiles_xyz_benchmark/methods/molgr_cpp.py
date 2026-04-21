from __future__ import annotations

import time
from dataclasses import dataclass
from importlib import import_module
from typing import Any

from benchmarks.smiles_xyz_benchmark.methods.base import BenchmarkMethod, MethodRunOutput
from benchmarks.smiles_xyz_benchmark.methods.postprocess import finalize_rdmol_with_dative_bonds


@dataclass(frozen=True)
class MolGRCppMethod(BenchmarkMethod):
    method_id: str = "molgr_cpp"

    def run(self, case: dict[str, Any]) -> MethodRunOutput:
        timing_ms_breakdown: dict[str, float] = {
            "cpp_xyz2omol_ms": 0.0,
            "cpp_no_metal_pipeline_ms": 0.0,
            "cpp_resonance_handling_enumeration_ms": 0.0,
            "cpp_metal_enumeration_combination_ms": 0.0,
            "mol_data_to_rdkit_ms": 0.0,
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
            cpp_pipeline = import_module("molgr._core.pipeline")
            cpp_xyz2omol_uncached = cpp_pipeline.reconstruct_with_metals.xyz2omol
            cpp_last_timing = getattr(
                cpp_pipeline,
                "get_last_run_timing_breakdown_ms",
                None,
            )
            mol_data_to_rdkit = import_module("molgr.interface").mol_data_to_rdkit
        except Exception as exc:  # noqa: BLE001
            return MethodRunOutput(
                status="error",
                error=f"import molgr core pipeline failed: {exc}",
                timing_ms_breakdown=timing_ms_breakdown,
            )

        def _merge_cpp_internal_timing() -> None:
            if cpp_last_timing is None:
                return
            try:
                raw = cpp_last_timing()
            except Exception:  # noqa: BLE001
                return

            timing_ms_breakdown["cpp_no_metal_pipeline_ms"] = float(
                raw.get("no_metal_pipeline_ms", 0.0)
            )
            timing_ms_breakdown["cpp_resonance_handling_enumeration_ms"] = float(
                raw.get("resonance_handling_enumeration_ms", 0.0)
            )
            timing_ms_breakdown["cpp_metal_enumeration_combination_ms"] = float(
                raw.get("metal_enumeration_combination_ms", 0.0)
            )

        cpp_started = time.perf_counter()
        try:
            mol_data = cpp_xyz2omol_uncached(
                xyz_block,
                total_charge,
                total_radical_electrons,
            )
        except Exception as exc:  # noqa: BLE001
            timing_ms_breakdown["cpp_xyz2omol_ms"] = (time.perf_counter() - cpp_started) * 1000.0
            _merge_cpp_internal_timing()
            return MethodRunOutput(
                status="error",
                error=f"molgr._core.pipeline.reconstruct_with_metals.xyz2omol failed: {exc}",
                timing_ms_breakdown=timing_ms_breakdown,
            )
        timing_ms_breakdown["cpp_xyz2omol_ms"] = (time.perf_counter() - cpp_started) * 1000.0
        _merge_cpp_internal_timing()

        if mol_data is None:
            return MethodRunOutput(
                status="error",
                error="molgr._core.pipeline.reconstruct_with_metals.xyz2omol returned None",
                timing_ms_breakdown=timing_ms_breakdown,
            )

        to_rdkit_started = time.perf_counter()
        try:
            rdkit_mol = mol_data_to_rdkit(mol_data)
        except Exception as exc:  # noqa: BLE001
            timing_ms_breakdown["mol_data_to_rdkit_ms"] = (
                time.perf_counter() - to_rdkit_started
            ) * 1000.0
            return MethodRunOutput(
                status="error",
                error=f"mol_data_to_rdkit failed: {exc}",
                timing_ms_breakdown=timing_ms_breakdown,
            )
        timing_ms_breakdown["mol_data_to_rdkit_ms"] = (
            time.perf_counter() - to_rdkit_started
        ) * 1000.0

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
