from __future__ import annotations

import subprocess
import sys
from dataclasses import replace

import pytest


pytest.importorskip("openbabel")
pytest.importorskip("rdkit")

from molgr.batch import (
    ReconstructionBatchRequest,
    iter_xyz_to_rdmol_batch,
)
from molgr.config import MolGRConfig, PythonInterfaceConfig
from molgr.diagnostics import ReconstructionError, ReconstructionFailureCode


_WATER_XYZ = """3
water
O 0.0000 0.0000 0.0000
H 0.7586 0.0000 0.5043
H -0.7586 0.0000 0.5043
"""


def _requests() -> list[ReconstructionBatchRequest]:
    return [
        ReconstructionBatchRequest(_WATER_XYZ),
        ReconstructionBatchRequest(_WATER_XYZ, spin_multiplicity=2),
        ReconstructionBatchRequest(_WATER_XYZ, total_charge=10),
    ]


def test_cpp_batch_streams_ordered_results_with_per_item_diagnostics() -> None:
    results = list(
        iter_xyz_to_rdmol_batch(
            _requests(),
            backend="cpp",
            max_workers=2,
            queue_size=1,
            ordered=True,
        )
    )

    assert all(isinstance(result, tuple) and len(result) == 3 for result in results)
    assert [result.input for result in results] == _requests()
    assert results[0].molecule is not None
    assert results[0].molecule.GetNumAtoms() == 3
    assert results[1].molecule is None
    assert results[1].diagnostics is not None
    assert results[1].diagnostics.code is ReconstructionFailureCode.INVALID_ELECTRONIC_TARGET
    assert results[2].molecule is None
    assert results[2].diagnostics is not None
    assert results[2].diagnostics.code is ReconstructionFailureCode.NO_VALID_ORGANIC_CANDIDATE


def test_cpp_batch_unordered_results_keep_input_correspondence() -> None:
    requests = _requests()
    results = list(
        iter_xyz_to_rdmol_batch(
            requests,
            backend="cpp",
            max_workers=2,
            queue_size=1,
            ordered=False,
        )
    )

    assert {result.input for result in results} == set(requests)
    assert all(result.input.xyz_block == _WATER_XYZ for result in results)


def test_python_batch_preserves_the_same_result_protocol() -> None:
    results = list(iter_xyz_to_rdmol_batch(_requests(), backend="python"))

    assert all(isinstance(result, tuple) and len(result) == 3 for result in results)
    assert [result.input for result in results] == _requests()
    assert results[0].molecule is not None
    assert results[0].molecule.GetNumAtoms() == 3
    assert results[1].diagnostics is not None
    assert results[1].diagnostics.code is ReconstructionFailureCode.INVALID_ELECTRONIC_TARGET
    assert results[2].diagnostics is not None
    assert results[2].diagnostics.code is ReconstructionFailureCode.NO_VALID_ORGANIC_CANDIDATE


def test_python_batch_rejects_nested_parallelism() -> None:
    with pytest.raises(ValueError, match="backend='python' only supports max_workers=1"):
        list(iter_xyz_to_rdmol_batch(_requests(), backend="python", max_workers=2))


@pytest.mark.parametrize("backend", ["cpp", "python"])
def test_batch_backends_reject_excessive_radical_electrons(backend: str) -> None:
    result = next(
        iter_xyz_to_rdmol_batch(
            [ReconstructionBatchRequest(_WATER_XYZ, spin_multiplicity=12)],
            backend=backend,
        )
    )

    assert result.molecule is None
    assert result.diagnostics is not None
    assert result.diagnostics.code is ReconstructionFailureCode.INVALID_ELECTRONIC_TARGET


def test_batch_can_return_suspicious_molecules_without_aborting() -> None:
    config = replace(
        MolGRConfig(),
        interface=PythonInterfaceConfig(reconstruction_failure_policy="return_suspicious"),
    )
    results = list(
        iter_xyz_to_rdmol_batch(
            _requests()[1:],
            backend="cpp",
            max_workers=2,
            config=config,
        )
    )

    assert len(results) == 2
    assert all(result.molecule is not None for result in results)
    assert all(
        result.molecule.GetProp("_MolGRReconstructionStatus") == "suspicious_fallback"
        for result in results
    )
    assert results[0].diagnostics is not None


def test_batch_result_can_be_consumed_as_request_molecule_pair() -> None:
    request = ReconstructionBatchRequest(_WATER_XYZ)
    result = next(iter_xyz_to_rdmol_batch([request], backend="cpp", max_workers=1))

    paired_request, molecule = result.as_pair()
    assert paired_request == request
    assert molecule is result.molecule
    assert result.as_dict() == {
        "input": request,
        "result": result.molecule,
        "status": None,
    }


def test_batch_accepts_one_shot_iterators() -> None:
    requests = _requests()
    results = list(iter_xyz_to_rdmol_batch((request for request in requests), backend="cpp"))

    assert {result.input for result in results} == set(requests)
    assert all(len(result) == 3 for result in results)


def test_batch_raise_on_error_preserves_structured_diagnostic() -> None:
    with pytest.raises(ReconstructionError) as raised:
        list(
            iter_xyz_to_rdmol_batch(
                _requests()[1:],
                backend="cpp",
                max_workers=2,
                raise_on_error=True,
            )
        )

    assert raised.value.diagnostics.code is ReconstructionFailureCode.INVALID_ELECTRONIC_TARGET


def test_cpp_batch_reuses_native_workers_across_repeated_close_paths() -> None:
    script = f"""\
from molgr import _core

xyz = {_WATER_XYZ!r}
requests = [(xyz, 0, 0)] * 60

for round_index in range(100):
    iterator = _core.pipeline.reconstruct_with_metals.batch_xyz2omol(
        requests,
        max_workers=0,
        queue_size=1,
        ordered=bool(round_index % 2),
    )
    if round_index % 2:
        first = next(iterator)
        assert first["molecule_data"] is not None
    else:
        results = list(iterator)
        assert len(results) == len(requests)
        assert all(result["molecule_data"] is not None for result in results)
    iterator.close()
print("batch-stress-ok")
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().endswith("batch-stress-ok")


def test_cpp_batches_share_workers_without_blocking_each_other() -> None:
    from molgr import _core

    requests = [(_WATER_XYZ, 0, 0)] * 60
    first = _core.pipeline.reconstruct_with_metals.batch_xyz2omol(
        requests, max_workers=8, queue_size=1, ordered=False
    )
    second = _core.pipeline.reconstruct_with_metals.batch_xyz2omol(
        requests, max_workers=8, queue_size=1, ordered=True
    )
    try:
        first_results = list(first)
        second_results = list(second)
    finally:
        first.close()
        second.close()

    assert len(first_results) == len(requests)
    assert len(second_results) == len(requests)
    assert all(result["molecule_data"] is not None for result in first_results)
    assert all(result["molecule_data"] is not None for result in second_results)
    assert [result["index"] for result in second_results] == list(range(len(requests)))
