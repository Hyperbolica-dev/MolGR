"""
Author: TMJ
Date: 2025-12-19 21:00:20
LastEditors: TMJ
LastEditTime: 2025-12-19 21:54:05
Description: 请填写简介
"""
import importlib.metadata


try:
    __version__ = importlib.metadata.version("molgr")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"


from ._core import get_atom_count

__all__ = ["get_atom_count"]
