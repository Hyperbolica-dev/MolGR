"""
Author: TMJ
Date: 2025-12-19 21:00:20
LastEditors: TMJ
LastEditTime: 2026-04-22 02:13:56
Description: 请填写简介
"""

import importlib.metadata

from openbabel import pybel

from . import _core as core
from . import config
from ._core import pipeline
from .config import (
    DEFAULT_MOLGR_CONFIG,
    CppBackendConfig,
    MetalRadicalInferenceConfig,
    MetalScoringConfig,
    MolGRConfig,
    OrganicTopologyConfig,
    PythonInterfaceConfig,
    ReconstructionFailurePolicy,
    ResonanceConfig,
    ResonanceTraversalScore,
    get_config,
    make_default_config,
    reset_config,
    set_config,
    sync_cpp_backend_default_config,
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
    "DEFAULT_MOLGR_CONFIG",
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
    "get_config",
    "make_default_config",
    "reset_config",
    "set_config",
    "sync_cpp_backend_default_config",
    "pipeline",
]

# 默认可以设为 WARN
set_log_level(core.LogLevel.WARN)
sync_cpp_backend_default_config(get_config())
