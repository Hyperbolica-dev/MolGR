from __future__ import annotations

import contextlib
import tempfile
import time
from dataclasses import dataclass
from importlib import import_module
from pathlib import Path
from typing import Any

from rdkit import Chem

from benchmarks.smiles_xyz_benchmark.methods.base import BenchmarkMethod, MethodRunOutput


@dataclass(frozen=True)
class Cell2MolV2Method(BenchmarkMethod):
    method_id: str = "cell2mol_v2"

    def run(self, case: dict[str, Any]) -> MethodRunOutput:
        timing_ms_breakdown: dict[str, float] = {
            "write_xyz_temp_ms": 0.0,
            "cell2mol_run_ms": 0.0,
            "extract_rdkit_ms": 0.0,
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

        try:
            get_molecule = import_module("cell2mol.xyz_molecule").get_molecule
        except Exception:  # noqa: BLE001
            return MethodRunOutput(
                status="skipped",
                error="cell2mol not installed",
                timing_ms_breakdown=timing_ms_breakdown,
            )

        classes_module: Any | None = None
        original_balance_charge: Any | None = None
        try:
            try:
                classes_module = import_module("cell2mol.classes")
                new_charge_assignment = import_module("cell2mol.new_charge_assignment")
                balance_charge_impl = getattr(new_charge_assignment, "balance_charge", None)
                original = getattr(classes_module, "balance_charge", None)
                if callable(balance_charge_impl) and callable(original):
                    original_balance_charge = original

                    def _compat_balance_charge(*args: Any, **kwargs: Any) -> Any:
                        if "charges_sum" in kwargs and "input_charge" not in kwargs:
                            kwargs["input_charge"] = kwargs.pop("charges_sum")
                        return balance_charge_impl(*args, **kwargs)

                    classes_module.balance_charge = _compat_balance_charge
            except Exception:  # noqa: BLE001
                pass

            write_started = time.perf_counter()
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir) / "input.xyz"
                temp_path.write_text(xyz_block, encoding="utf-8")
                timing_ms_breakdown["write_xyz_temp_ms"] = (
                    time.perf_counter() - write_started
                ) * 1000.0

                cell2mol_started = time.perf_counter()
                mol = get_molecule(str(temp_path), "mol", total_charge, temp_dir, 0)
                timing_ms_breakdown["cell2mol_run_ms"] = (
                    time.perf_counter() - cell2mol_started
                ) * 1000.0

                extract_started = time.perf_counter()
                rdkit_mol = getattr(mol, "rdkit_obj", None)
                timing_ms_breakdown["extract_rdkit_ms"] = (
                    time.perf_counter() - extract_started
                ) * 1000.0

                if rdkit_mol is None:
                    return MethodRunOutput(
                        status="error",
                        error="cell2mol molecule missing rdkit_obj",
                        timing_ms_breakdown=timing_ms_breakdown,
                    )
                rdkit_mol = Chem.RemoveHs(rdkit_mol)
                postprocess_started = time.perf_counter()
                warnings: list[str] = []
                predicted_smiles: str | None = None
                try:
                    predicted_smiles = Chem.MolToSmiles(
                        rdkit_mol,
                        canonical=True,
                        isomericSmiles=True,
                    )
                except Exception as exc:  # noqa: BLE001
                    warnings.append(f"MolToSmiles failed: {exc}")
                timing_ms_breakdown["postprocess_ms"] = (
                    time.perf_counter() - postprocess_started
                ) * 1000.0

                return MethodRunOutput(
                    status="ok",
                    error="; ".join(warnings) if warnings else None,
                    predicted_smiles=predicted_smiles,
                    rdkit_mol=rdkit_mol,
                    timing_ms_breakdown=timing_ms_breakdown,
                )
        except Exception as exc:  # noqa: BLE001
            return MethodRunOutput(
                status="error",
                error=f"cell2mol_v2 failed: {exc}",
                timing_ms_breakdown=timing_ms_breakdown,
            )
        finally:
            if classes_module is not None and callable(original_balance_charge):
                with contextlib.suppress(Exception):
                    classes_module.balance_charge = original_balance_charge
