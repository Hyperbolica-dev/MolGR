"""
Author: TMJ
Date: 2025-12-19 21:00:20
LastEditors: TMJ
LastEditTime: 2026-02-27 20:49:01
Description: 请填写简介
"""

import importlib.metadata

from openbabel import pybel

from . import _core as core
from ._core import pipeline, stages, utils  # type: ignore


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
    "set_log_level",
    "pipeline",
    "stages",
    "utils",
]

# 默认可以设为 WARN
set_log_level(core.LogLevel.WARN)
