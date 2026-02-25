from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class MethodRunOutput:
    status: str
    error: str | None = None
    predicted_smiles: str | None = None
    rdkit_mol: Any | None = None
    equivalent: bool | None = None
    equivalence_method: str | None = None
    timing_ms_breakdown: dict[str, float] | None = None


@dataclass(frozen=True)
class BenchmarkMethod:
    method_id: str

    def run(self, case: dict[str, Any]) -> MethodRunOutput:
        del case
        return MethodRunOutput(
            status="skipped",
            error="not implemented",
            timing_ms_breakdown={"method_ms": 0.0},
        )
