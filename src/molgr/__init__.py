from __future__ import annotations  # noqa: I001

import importlib.metadata

# RDKit must initialize before Open Babel's pybel on cp313 manylinux.
# isort: off
from rdkit import Chem as _rdkit_chem  # noqa: F401
from openbabel import pybel
# isort: on

from . import _core as core
from . import config
from ._core import pipeline
from .config import (
    CONFIG,
    CppBackendConfig,
    MetalRadicalInferenceConfig,
    MetalScoringConfig,
    MolGRConfig,
    OrganicTopologyConfig,
    PythonInterfaceConfig,
    ReconstructionFailurePolicy,
    ResonanceConfig,
    ResonanceTraversalScore,
)


try:
    __version__ = importlib.metadata.version("molgr")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"


pybel.ob.obErrorLog.StopLogging()


def set_log_level(level: core.LogLevel) -> None:
    """
    Set C++ backend logging level.

    Args:
        level: Can be LogLevel.DEBUG, LogLevel.INFO, etc., or int (0-4).
    """
    # Pybind11 的枚举可以直接接受 int，也可以接受枚举对象
    core.set_log_level(level)


__all__ = [
    "CONFIG",
    "CppBackendConfig",
    "MetalRadicalInferenceConfig",
    "MetalScoringConfig",
    "MolGRConfig",
    "OrganicTopologyConfig",
    "PythonInterfaceConfig",
    "ReconstructionFailurePolicy",
    "ResonanceConfig",
    "ResonanceTraversalScore",
    "config",
    "set_log_level",
    "pipeline",
]

set_log_level(core.LogLevel.WARN)
