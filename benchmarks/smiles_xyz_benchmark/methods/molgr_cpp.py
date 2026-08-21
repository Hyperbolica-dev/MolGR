from __future__ import annotations

import time
from dataclasses import dataclass
from importlib import import_module
from typing import Any, Iterable

from rdkit import Chem

from benchmarks.smiles_xyz_benchmark.methods.base import BenchmarkMethod, MethodRunOutput
from benchmarks.smiles_xyz_benchmark.methods.postprocess import remove_hs_without_sanitize
from molgr.batch import ReconstructionBatchRequest, iter_xyz_to_rdmol_batch


@dataclass(frozen=True)
class MolGRCppMethod(BenchmarkMethod):
    method_id: str = "molgr_cpp"

    def run(self, case: dict[str, Any]) -> MethodRunOutput:
        timing_ms_breakdown: dict[str, float] = {
            "molgr_interface_ms": 0.0,
            "cpp_no_metal_pipeline_ms": 0.0,
            "cpp_resonance_handling_enumeration_ms": 0.0,
            "cpp_metal_enumeration_combination_ms": 0.0,
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
            molgr_interface = import_module("molgr.interface")
            xyz_to_rdmol = molgr_interface.xyz_to_rdmol
            cpp_pipeline = import_module("molgr._core.pipeline")
            cpp_last_timing = getattr(
                cpp_pipeline,
                "get_last_run_timing_breakdown_ms",
                None,
            )
        except Exception as exc:  # noqa: BLE001
            return MethodRunOutput(
                status="error",
                error=f"import molgr interface failed: {exc}",
                timing_ms_breakdown=timing_ms_breakdown,
            )

        def _merge_cpp_internal_timing() -> None:
            if cpp_last_timing is None:
                return
            try:
                raw = cpp_last_timing()
            except Exception:  # noqa: BLE001
                return

            for key, value in raw.items():
                timing_ms_breakdown[f"cpp_{key}"] = float(value)

        interface_started = time.perf_counter()
        try:
            rdkit_mol = xyz_to_rdmol(
                xyz_block,
                total_charge,
                total_radical_electrons + 1,
                backend="cpp",
                make_dative_bonds=True,
            )
        except Exception as exc:  # noqa: BLE001
            timing_ms_breakdown["molgr_interface_ms"] = (
                time.perf_counter() - interface_started
            ) * 1000.0
            _merge_cpp_internal_timing()
            return MethodRunOutput(
                status="error",
                error=f"molgr.interface.xyz_to_rdmol failed: {exc}",
                timing_ms_breakdown=timing_ms_breakdown,
            )
        timing_ms_breakdown["molgr_interface_ms"] = (
            time.perf_counter() - interface_started
        ) * 1000.0
        _merge_cpp_internal_timing()
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

    def run_batch(
        self,
        cases: Iterable[dict[str, Any]],
        *,
        max_workers: int | None = None,
    ) -> dict[int, MethodRunOutput]:
        """Run all valid cases through the native batch scheduler.

        The benchmark runner uses this entry point instead of splitting C++
        work across Python subprocesses. Post-processing remains per result,
        while native reconstruction owns the bounded worker budget.
        """

        normalized_cases = list(cases)
        outputs: dict[int, MethodRunOutput] = {}
        requests: list[ReconstructionBatchRequest] = []
        request_case_indices: list[int] = []
        for case in normalized_cases:
            case_idx = int(case["case_idx"])
            timing = {
                "molgr_interface_ms": 0.0,
                "postprocess_ms": 0.0,
            }
            xyz_block = case.get("xyz_block")
            if not isinstance(xyz_block, str) or not xyz_block.strip():
                outputs[case_idx] = MethodRunOutput(
                    status="error",
                    error="missing or empty xyz_block",
                    timing_ms_breakdown=timing,
                )
                continue
            try:
                total_charge = int(case.get("total_charge", 0))
                total_radical_electrons = int(case.get("total_radical_electrons", 0))
            except (TypeError, ValueError) as exc:
                outputs[case_idx] = MethodRunOutput(
                    status="error",
                    error=f"invalid electronic state: {exc}",
                    timing_ms_breakdown=timing,
                )
                continue
            requests.append(
                ReconstructionBatchRequest(
                    xyz_block=xyz_block,
                    total_charge=total_charge,
                    spin_multiplicity=total_radical_electrons + 1,
                )
            )
            request_case_indices.append(case_idx)

        started = time.perf_counter()
        if not requests:
            return outputs
        try:
            batch_results = iter_xyz_to_rdmol_batch(
                requests,
                backend="cpp",
                max_workers=max_workers,
                ordered=True,
            )
            for request_index, batch_result in enumerate(batch_results):
                case_idx = request_case_indices[request_index]
                timing = {
                    "molgr_interface_ms": 0.0,
                    "postprocess_ms": 0.0,
                }
                if batch_result.molecule is None:
                    diagnostics = batch_result.diagnostics
                    error = "native batch reconstruction failed"
                    if diagnostics is not None:
                        error = (
                            f"{diagnostics.code.value} at {diagnostics.stage}: "
                            f"{diagnostics.message}"
                        )
                    outputs[case_idx] = MethodRunOutput(
                        status="error",
                        error=error,
                        timing_ms_breakdown=timing,
                    )
                    continue
                postprocess_started = time.perf_counter()
                try:
                    rdkit_mol = remove_hs_without_sanitize(batch_result.molecule)
                    predicted_smiles = Chem.MolToSmiles(
                        rdkit_mol,
                        canonical=True,
                        isomericSmiles=True,
                    )
                except Exception as exc:  # noqa: BLE001
                    timing["postprocess_ms"] = (time.perf_counter() - postprocess_started) * 1000.0
                    outputs[case_idx] = MethodRunOutput(
                        status="error",
                        error=f"MolGR result postprocess failed: {exc}",
                        rdkit_mol=batch_result.molecule,
                        timing_ms_breakdown=timing,
                    )
                    continue
                timing["postprocess_ms"] = (time.perf_counter() - postprocess_started) * 1000.0
                outputs[case_idx] = MethodRunOutput(
                    status="ok",
                    predicted_smiles=predicted_smiles,
                    rdkit_mol=rdkit_mol,
                    timing_ms_breakdown=timing,
                )
        except Exception as exc:  # noqa: BLE001
            for case_idx in request_case_indices:
                outputs.setdefault(
                    case_idx,
                    MethodRunOutput(
                        status="error",
                        error=f"native batch execution failed: {exc}",
                        timing_ms_breakdown={
                            "molgr_interface_ms": (time.perf_counter() - started) * 1000.0,
                            "postprocess_ms": 0.0,
                        },
                    ),
                )
        batch_elapsed_ms = (time.perf_counter() - started) * 1000.0
        amortized_batch_ms = batch_elapsed_ms / len(requests)
        for case_idx in request_case_indices:
            output_timing = outputs[case_idx].timing_ms_breakdown
            if output_timing is None:
                continue
            output_timing["native_batch_elapsed_ms"] = batch_elapsed_ms
            # A per-case benchmark value is necessarily amortized when work is
            # scheduled by one native batch. Keep the full wall time in the
            # breakdown for throughput analysis.
            output_timing["method_ms"] = amortized_batch_ms
        return outputs
