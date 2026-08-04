from __future__ import annotations

from dataclasses import replace

import pytest


pytest.importorskip("openbabel")
pytest.importorskip("rdkit")

from molgr.config import MolGRConfig, PythonInterfaceConfig
from molgr.interface import xyz_to_rdmol


_WATER_XYZ = """3
water
O 0.0000 0.0000 0.0000
H 0.7586 0.0000 0.5043
H -0.7586 0.0000 0.5043
"""


@pytest.mark.parametrize("backend", ["python", "cpp"])
def test_even_electron_system_rejects_even_multiplicity_before_reconstruction(
    monkeypatch: pytest.MonkeyPatch,
    backend: str,
) -> None:
    backend_called = False

    def reconstruct(*args: object, **kwargs: object) -> None:
        nonlocal backend_called
        backend_called = True

    target = (
        "molgr.interface.xyz2omol"
        if backend == "python"
        else "molgr.interface.core.pipeline.reconstruct_with_metals.xyz2omol"
    )
    monkeypatch.setattr(target, reconstruct)

    with pytest.raises(
        ValueError,
        match=r"spin_multiplicity=2 is impossible for 10 total electrons; "
        r"the multiplicity must be odd",
    ):
        xyz_to_rdmol(_WATER_XYZ, spin_multiplicity=2, backend=backend)

    assert not backend_called


def test_odd_electron_system_rejects_odd_multiplicity_before_reconstruction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    backend_called = False

    def reconstruct(*args: object, **kwargs: object) -> None:
        nonlocal backend_called
        backend_called = True

    monkeypatch.setattr("molgr.interface.xyz2omol", reconstruct)

    with pytest.raises(
        ValueError,
        match=r"spin_multiplicity=1 is impossible for 9 total electrons; "
        r"the multiplicity must be even",
    ):
        xyz_to_rdmol(
            _WATER_XYZ,
            total_charge=1,
            spin_multiplicity=1,
            backend="python",
        )

    assert not backend_called


def test_odd_electron_system_accepts_even_multiplicity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    backend_called = False

    def reconstruct(*args: object, **kwargs: object) -> None:
        nonlocal backend_called
        backend_called = True

    monkeypatch.setattr("molgr.interface.xyz2omol", reconstruct)

    with pytest.raises(ValueError, match="xyz2omol failed"):
        xyz_to_rdmol(
            _WATER_XYZ,
            total_charge=1,
            spin_multiplicity=2,
            backend="python",
        )

    assert backend_called


def test_multiplicity_cannot_exceed_total_electron_limit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    backend_called = False

    def reconstruct(*args: object, **kwargs: object) -> None:
        nonlocal backend_called
        backend_called = True

    monkeypatch.setattr("molgr.interface.xyz2omol", reconstruct)

    with pytest.raises(
        ValueError,
        match=r"spin_multiplicity=12 is impossible for 10 total electrons; the maximum is 11",
    ):
        xyz_to_rdmol(_WATER_XYZ, spin_multiplicity=12, backend="python")

    assert not backend_called


def test_impossible_spin_is_not_converted_to_suspicious_fallback() -> None:
    base_config = MolGRConfig()
    config = replace(
        base_config,
        interface=PythonInterfaceConfig(reconstruction_failure_policy="return_suspicious"),
    )

    with pytest.raises(ValueError, match=r"multiplicity must be odd"):
        xyz_to_rdmol(
            _WATER_XYZ,
            spin_multiplicity=2,
            backend="python",
            config=config,
        )
