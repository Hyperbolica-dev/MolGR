"""
Author: TMJ
Date: 2025-12-19 21:00:20
LastEditors: TMJ
LastEditTime: 2025-12-26 15:54:28
Description: 请填写简介
"""

import importlib.metadata
from typing import Union

from openbabel import openbabel as ob
from openbabel import pybel

from . import _core  # type: ignore
from ._core import (  # type: ignore
    LogLevel,
    calculate_shape_quality,
    calculate_tetrahedron_volume,
    get_possible_metal_radicals,
)


try:
    __version__ = importlib.metadata.version("molgr")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"


__all__ = [
    "_core",
    "calculate_shape_quality",
    "calculate_tetrahedron_volume",
    "get_possible_metal_radicals",
]


def omol_score(mol: Union[ob.OBMol, pybel.Molecule]) -> float:
    """
    Calculate the resonance score of a molecule using the optimized C++ backend.
    Compatible with standard openbabel and pybel objects.
    """
    obmol = None

    # 1. 解包 pybel 对象
    if isinstance(mol, pybel.Molecule):
        obmol = mol.OBMol
    # 2. 处理原生 OBMol 对象
    elif isinstance(mol, ob.OBMol):
        obmol = mol
    else:
        raise TypeError(f"Input must be openbabel.OBMol or pybel.Molecule, got {type(mol)}")

    # 3. 获取 C++ 内存指针
    # SWIG 对象通常有一个 .this 属性，强转为 int 即可获得地址
    try:
        obmol_ptr = int(obmol.this) if hasattr(obmol, "this") else int(obmol)  # type: ignore
    except (ValueError, TypeError) as e:
        raise ValueError("Could not retrieve C++ memory address from OBMol object") from e

    # 4. 传递给 C++
    return _core.omol_score_from_ptr(obmol_ptr)


def set_log_level(level: LogLevel):
    """
    Set C++ backend logging level.

    Args:
        level: Can be LogLevel.DEBUG, LogLevel.INFO, etc., or int (0-4).
    """
    # Pybind11 的枚举可以直接接受 int，也可以接受枚举对象
    _core.set_log_level(level)


# 默认可以设为 WARN
set_log_level(LogLevel.WARN)
