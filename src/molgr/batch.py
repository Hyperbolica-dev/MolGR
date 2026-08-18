"""Streaming batch reconstruction APIs."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Iterable, Iterator, Literal, NamedTuple, Sequence, cast

from rdkit import Chem

from molgr.config import CONFIG, MolGRConfig
from molgr.diagnostics import (
    ReconstructionDiagnostics,
    ReconstructionError,
    ReconstructionFailureCode,
)
from molgr.interface import (
    _finalize_rdmol,
    _should_return_suspicious_on_reconstruction_failure,
    _suspicious_rdmol_from_input_xyz,
    xyz_to_rdmol,
)
from molgr.utils.converter import mol_data_to_rdkit

from . import _core as core


@dataclass(frozen=True)
class ReconstructionBatchRequest:
    """One XYZ reconstruction request for a batch."""

    xyz_block: str
    total_charge: int = 0
    spin_multiplicity: int = 1


class ReconstructionBatchResult(NamedTuple):
    """One completed ``(input, result, status)`` batch item.

    ``status`` is ``None`` for a successful reconstruction and contains a
    structured diagnostic for a non-fatal item failure.
    """

    input: ReconstructionBatchRequest
    result: Chem.Mol | None
    status: ReconstructionDiagnostics | None = None

    @property
    def request(self) -> ReconstructionBatchRequest:
        """Compatibility alias for the first tuple element."""

        return self.input

    @property
    def molecule(self) -> Chem.Mol | None:
        """Compatibility alias for the reconstructed molecule."""

        return self.result

    @property
    def diagnostics(self) -> ReconstructionDiagnostics | None:
        """Compatibility alias for the item status."""

        return self.status

    def as_pair(self) -> tuple[ReconstructionBatchRequest, Chem.Mol | None]:
        """Return the input request and reconstructed molecule as a pair."""

        return self.input, self.result

    def as_dict(self) -> dict[str, object]:
        """Return a self-describing representation for unstructured consumers."""

        return {
            "input": self.input,
            "result": self.result,
            "status": self.status,
        }


def _normalize_request(
    value: ReconstructionBatchRequest | Sequence[Any],
) -> ReconstructionBatchRequest:
    if isinstance(value, ReconstructionBatchRequest):
        return value
    if len(value) != 3:
        raise ValueError(
            "each batch request must contain exactly three values: "
            "xyz_block, total_charge, spin_multiplicity"
        )
    return ReconstructionBatchRequest(
        xyz_block=str(value[0]),
        total_charge=int(value[1]),
        spin_multiplicity=int(value[2]),
    )


def _diagnostics_from_suspicious_mol(
    molecule: Chem.Mol,
    *,
    backend: Literal["cpp", "python"],
) -> ReconstructionDiagnostics:
    try:
        raw = json.loads(molecule.GetProp("_MolGRReconstructionDiagnostics"))
    except Exception:
        return ReconstructionDiagnostics(
            code=ReconstructionFailureCode.NO_VALID_RECONSTRUCTION,
            stage="reconstruction",
            backend=backend,
            message="The reconstruction backend returned a suspicious molecule without diagnostics.",
        )
    return ReconstructionDiagnostics.from_mapping(raw, backend=backend)


def _batch_result_from_native(
    raw: dict[str, Any],
    request: ReconstructionBatchRequest,
    *,
    config: MolGRConfig,
    make_dative_bonds: bool,
    make_stereochemistry: bool,
    raise_on_error: bool,
) -> ReconstructionBatchResult:
    molecule_data = raw.get("molecule_data")
    raw_diagnostics = raw.get("diagnostics")
    diagnostics = (
        ReconstructionDiagnostics.from_mapping(raw_diagnostics, backend="cpp")
        if isinstance(raw_diagnostics, dict) and raw_diagnostics
        else None
    )
    if molecule_data is None:
        diagnostics = diagnostics or ReconstructionDiagnostics(
            code=ReconstructionFailureCode.NO_VALID_RECONSTRUCTION,
            stage="reconstruction",
            backend="cpp",
            message="The C++ reconstruction backend returned no molecule.",
        )
        if raise_on_error:
            raise ReconstructionError(diagnostics)
        if _should_return_suspicious_on_reconstruction_failure(config):
            molecule = _suspicious_rdmol_from_input_xyz(request.xyz_block, diagnostics)
            return ReconstructionBatchResult(request, molecule, diagnostics)
        return ReconstructionBatchResult(request, None, diagnostics)

    try:
        molecule = mol_data_to_rdkit(molecule_data)
        molecule = _finalize_rdmol(
            molecule,
            request.xyz_block,
            backend="cpp",
            make_dative_bonds=make_dative_bonds,
            make_stereochemistry=make_stereochemistry,
            config=config,
        )
    except ReconstructionError as exc:
        if raise_on_error:
            raise
        diagnostics = exc.diagnostics
        if _should_return_suspicious_on_reconstruction_failure(config):
            molecule = _suspicious_rdmol_from_input_xyz(request.xyz_block, diagnostics)
            return ReconstructionBatchResult(request, molecule, diagnostics)
        return ReconstructionBatchResult(request, None, diagnostics)
    except Exception as exc:
        diagnostics = ReconstructionDiagnostics(
            code=ReconstructionFailureCode.RDKIT_POSTPROCESS_FAILED,
            stage="rdkit.conversion",
            backend="cpp",
            message="The C++ molecule could not be converted to an RDKit molecule.",
            cause_type=type(exc).__name__,
            cause_message=str(exc),
        )
        if raise_on_error:
            raise ReconstructionError(diagnostics) from exc
        if _should_return_suspicious_on_reconstruction_failure(config):
            molecule = _suspicious_rdmol_from_input_xyz(request.xyz_block, diagnostics)
            return ReconstructionBatchResult(request, molecule, diagnostics)
        return ReconstructionBatchResult(request, None, diagnostics)
    if molecule.HasProp("_MolGRReconstructionDiagnostics"):
        diagnostics = _diagnostics_from_suspicious_mol(molecule, backend="cpp")
        if raise_on_error:
            raise ReconstructionError(diagnostics)
    return ReconstructionBatchResult(request, molecule, diagnostics)


def _iter_cpp_batch(
    requests: list[ReconstructionBatchRequest],
    *,
    config: MolGRConfig,
    max_workers: int | None,
    queue_size: int,
    ordered: bool,
    make_dative_bonds: bool,
    make_stereochemistry: bool,
    raise_on_error: bool,
) -> Iterator[ReconstructionBatchResult]:
    native_requests = [
        (request.xyz_block, request.total_charge, request.spin_multiplicity - 1)
        for request in requests
    ]
    native_iterator = core.pipeline.reconstruct_with_metals.batch_xyz2omol(
        native_requests,
        config=config,
        max_workers=0 if max_workers is None else max_workers,
        queue_size=queue_size,
        ordered=ordered,
    )
    try:
        for raw_value in native_iterator:
            raw = cast(Any, raw_value)
            request = requests[int(raw["index"])]
            yield _batch_result_from_native(
                raw,
                request,
                config=config,
                make_dative_bonds=make_dative_bonds,
                make_stereochemistry=make_stereochemistry,
                raise_on_error=raise_on_error,
            )
    finally:
        native_iterator.close()


def _iter_python_batch(
    requests: list[ReconstructionBatchRequest],
    *,
    config: MolGRConfig,
    make_dative_bonds: bool,
    make_stereochemistry: bool,
    raise_on_error: bool,
) -> Iterator[ReconstructionBatchResult]:
    for request in requests:
        try:
            molecule = xyz_to_rdmol(
                request.xyz_block,
                request.total_charge,
                request.spin_multiplicity,
                backend="python",
                make_dative_bonds=make_dative_bonds,
                make_stereochemistry=make_stereochemistry,
                config=config,
            )
        except ReconstructionError as exc:
            if raise_on_error:
                raise
            diagnostics = exc.diagnostics
            if _should_return_suspicious_on_reconstruction_failure(config):
                molecule = _suspicious_rdmol_from_input_xyz(request.xyz_block, diagnostics)
                yield ReconstructionBatchResult(request, molecule, diagnostics)
            else:
                yield ReconstructionBatchResult(request, None, diagnostics)
            continue
        diagnostics = None
        if molecule.HasProp("_MolGRReconstructionDiagnostics"):
            diagnostics = _diagnostics_from_suspicious_mol(molecule, backend="python")
            if raise_on_error:
                raise ReconstructionError(diagnostics)
        yield ReconstructionBatchResult(request, molecule, diagnostics)


def iter_xyz_to_rdmol_batch(
    requests: Iterable[ReconstructionBatchRequest | Sequence[object]],
    *,
    backend: Literal["cpp", "python"] = "cpp",
    max_workers: int | None = None,
    queue_size: int = 16,
    ordered: bool = False,
    make_dative_bonds: bool = True,
    make_stereochemistry: bool = True,
    config: MolGRConfig | None = None,
    raise_on_error: bool = False,
) -> Iterator[ReconstructionBatchResult]:
    """Stream ``(input, result, status)`` triples from a finite batch.

    The C++ backend owns the worker pool and bounded result queue. The Python
    backend intentionally remains sequential and is provided for parity and
    reference behavior without creating another Python thread pool. Every
    triple includes its original request, so ``ordered=False`` does not lose
    input/result correspondence. The input iterable is normalized once before
    native scheduling, so generators and other one-shot iterators are valid.
    """

    resolved_config = CONFIG if config is None else config
    if queue_size < 1:
        raise ValueError("queue_size must be >= 1")
    if max_workers is not None and max_workers < 1:
        raise ValueError("max_workers must be >= 1 when provided")
    if backend == "python":
        if max_workers not in (None, 1):
            raise ValueError("backend='python' only supports max_workers=1")
        normalized_requests = [_normalize_request(value) for value in requests]
        yield from _iter_python_batch(
            normalized_requests,
            config=resolved_config,
            make_dative_bonds=make_dative_bonds,
            make_stereochemistry=make_stereochemistry,
            raise_on_error=raise_on_error,
        )
        return
    if backend != "cpp":
        raise ValueError(f"Unknown backend: {backend}")
    normalized_requests = [_normalize_request(value) for value in requests]
    yield from _iter_cpp_batch(
        normalized_requests,
        config=resolved_config,
        max_workers=max_workers,
        queue_size=queue_size,
        ordered=ordered,
        make_dative_bonds=make_dative_bonds,
        make_stereochemistry=make_stereochemistry,
        raise_on_error=raise_on_error,
    )


__all__ = [
    "ReconstructionBatchRequest",
    "ReconstructionBatchResult",
    "iter_xyz_to_rdmol_batch",
]
