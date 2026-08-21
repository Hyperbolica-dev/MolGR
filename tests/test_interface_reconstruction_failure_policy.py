from __future__ import annotations

from dataclasses import replace

import pytest


pytest.importorskip("openbabel")
pytest.importorskip("rdkit")

from molgr.config import MolGRConfig, PythonInterfaceConfig
from molgr.diagnostics import (
    RECONSTRUCTION_FAILURE_CODES,
    ReconstructionError,
    ReconstructionFailureCode,
)
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


def test_reconstruction_failure_codes_are_enum_members() -> None:
    assert RECONSTRUCTION_FAILURE_CODES is ReconstructionFailureCode
    assert RECONSTRUCTION_FAILURE_CODES.INVALID_XYZ.value == "INVALID_XYZ"
    assert str(RECONSTRUCTION_FAILURE_CODES.INVALID_XYZ) == "INVALID_XYZ"
    assert all(isinstance(code, ReconstructionFailureCode) for code in RECONSTRUCTION_FAILURE_CODES)


def test_reconstruction_failure_policy_defaults_to_raise(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _force_python_backend_failure(monkeypatch)

    with pytest.raises(ValueError, match="xyz2omol failed"):
        xyz_to_rdmol(_WATER_XYZ, backend="python")


def test_reconstruction_failure_exposes_structured_python_diagnostics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail(*args: object, **kwargs: object) -> None:
        del args
        diagnostics = kwargs["_diagnostics"]
        assert isinstance(diagnostics, dict)
        diagnostics.update(
            {
                "code": "NO_VALID_ORGANIC_CANDIDATE",
                "stage": "no_metal.reconstruction",
                "backend": "python",
                "message": "no organic candidate",
                "counts": {"target_buckets": 2},
                "details": {"target_charge": 3},
            }
        )
        return None

    monkeypatch.setattr("molgr.interface.xyz2omol", fail)

    with pytest.raises(ReconstructionError) as raised:
        xyz_to_rdmol(_WATER_XYZ, backend="python")

    assert raised.value.diagnostics.code is ReconstructionFailureCode.NO_VALID_ORGANIC_CANDIDATE
    assert raised.value.diagnostics.stage == "no_metal.reconstruction"
    assert raised.value.diagnostics.counts == {"target_buckets": 2}
    assert "[NO_VALID_ORGANIC_CANDIDATE]" in str(raised.value)


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


def test_suspicious_fallback_preserves_reconstruction_diagnostics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail(*args: object, **kwargs: object) -> None:
        del args
        diagnostics = kwargs["_diagnostics"]
        assert isinstance(diagnostics, dict)
        diagnostics.update(
            {
                "code": "ALL_METAL_CANDIDATES_REJECTED",
                "stage": "metal.scoring",
                "backend": "python",
                "message": "all candidates rejected",
                "counts": {"metal_candidates_scored": 4},
            }
        )
        return None

    monkeypatch.setattr("molgr.interface.xyz2omol", fail)
    config = replace(
        MolGRConfig(),
        interface=PythonInterfaceConfig(reconstruction_failure_policy="return_suspicious"),
    )

    mol = xyz_to_rdmol(_WATER_XYZ, backend="python", config=config)

    assert mol.GetProp("_MolGRReconstructionStatus") == "suspicious_fallback"
    assert mol.GetProp("_MolGRReconstructionFailureCode") == "ALL_METAL_CANDIDATES_REJECTED"
    assert mol.GetProp("_MolGRReconstructionFailureStage") == "metal.scoring"
    assert '"metal_candidates_scored": 4' in mol.GetProp("_MolGRReconstructionDiagnostics")


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


def test_cpp_reconstruction_failure_exposes_structured_diagnostics() -> None:
    with pytest.raises(ReconstructionError) as raised:
        xyz_to_rdmol(_WATER_XYZ, total_charge=10, backend="cpp")

    diagnostics = raised.value.diagnostics
    assert diagnostics.backend == "cpp"
    assert diagnostics.code is ReconstructionFailureCode.NO_VALID_ORGANIC_CANDIDATE
    assert diagnostics.stage == "no_metal.reconstruction"


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
