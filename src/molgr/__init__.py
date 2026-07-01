from __future__ import annotations

import importlib.metadata

from openbabel import pybel

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


def set_log_level(level: core.LogLevel):
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

_log_level = getattr(core, "LogLevel", None)
if _log_level is not None:
    set_log_level(_log_level.WARN)
