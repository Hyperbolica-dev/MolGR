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


def _force_python_backend_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("molgr.interface.xyz2omol", lambda *args, **kwargs: None)


def _force_cpp_backend_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "molgr.interface.core.pipeline.reconstruct_with_metals.xyz2omol",
        lambda *args, **kwargs: None,
    )


def test_reconstruction_failure_policy_defaults_to_raise(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _force_python_backend_failure(monkeypatch)

    with pytest.raises(ValueError, match="xyz2omol failed"):
        xyz_to_rdmol(_WATER_XYZ, backend="python")


def test_reconstruction_failure_policy_can_return_suspicious_initial_bond_perception(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _force_python_backend_failure(monkeypatch)
    base_config = MolGRConfig()
    config = replace(
        base_config,
        interface=PythonInterfaceConfig(reconstruction_failure_policy="return_suspicious"),
    )

    mol = xyz_to_rdmol(
        _WATER_XYZ,
        backend="python",
        make_dative_bonds=True,
        make_stereochemistry=True,
        config=config,
    )

    assert mol.GetProp("_MolGRReconstructionStatus") == "suspicious_fallback"
    assert mol.GetNumAtoms() == 3
    assert mol.GetNumBonds() == 2
    assert mol.GetConformer().GetAtomPosition(0).x == pytest.approx(0.0)


def test_reconstruction_failure_policy_applies_to_cpp_backend(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _force_cpp_backend_failure(monkeypatch)
    base_config = MolGRConfig()
    config = replace(
        base_config,
        interface=PythonInterfaceConfig(reconstruction_failure_policy="return_suspicious"),
    )

    mol = xyz_to_rdmol(_WATER_XYZ, backend="cpp", config=config)

    assert mol.GetProp("_MolGRReconstructionStatus") == "suspicious_fallback"
    assert mol.GetNumAtoms() == 3
    assert mol.GetNumBonds() == 2


def test_reconstruction_failure_policy_can_catch_backend_exceptions(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail(*args, **kwargs):
        raise RuntimeError("backend failed")

    monkeypatch.setattr("molgr.interface.xyz2omol", fail)
    base_config = MolGRConfig()
    config = replace(
        base_config,
        interface=PythonInterfaceConfig(reconstruction_failure_policy="return_suspicious"),
    )

    mol = xyz_to_rdmol(_WATER_XYZ, backend="python", config=config)

    assert mol.GetProp("_MolGRReconstructionStatus") == "suspicious_fallback"


def test_reconstruction_failure_policy_rejects_unknown_value(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _force_python_backend_failure(monkeypatch)
    base_config = MolGRConfig()
    config = replace(
        base_config,
        interface=PythonInterfaceConfig(reconstruction_failure_policy="return_suspicious"),
    )
    object.__setattr__(config.interface, "reconstruction_failure_policy", "unknown")

    with pytest.raises(ValueError, match="Unknown reconstruction failure policy"):
        xyz_to_rdmol(_WATER_XYZ, backend="python", config=config)
