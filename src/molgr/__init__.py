"""
Author: TMJ
Date: 2025-12-19 21:00:20
LastEditors: TMJ
LastEditTime: 2026-01-01 23:49:16
Description: 请填写简介
"""

import importlib.metadata
from typing import Union

from openbabel import openbabel as ob
from openbabel import pybel

from . import _core as core
from ._core import consts, metal, reconstruct, scoring, utils


try:
    __version__ = importlib.metadata.version("molgr")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"


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
    return scoring.omol_score_from_ptr(obmol_ptr)


def set_log_level(level: core.LogLevel):
    """
    Set C++ backend logging level.

    Args:
        level: Can be LogLevel.DEBUG, LogLevel.INFO, etc., or int (0-4).
    """
    # Pybind11 的枚举可以直接接受 int，也可以接受枚举对象
    core.set_log_level(level)


def reconstruct_to_pybel(
    xyz_block: str, total_charge: int, total_radical: int
) -> Union[pybel.Molecule, None]:
    """
    调用 C++ 核心重建，通过结构化数据中转，生成 Pybel 对象。
    """
    # 1. C++ 核心计算
    mol_ptr = reconstruct.reconstruct_from_xyz_no_metal(xyz_block, total_charge, total_radical)

    if mol_ptr == 0:
        return None

    try:
        # 2. 提取结构化数据 (MoleculeData 对象)
        # 这一步返回的是 _core.MoleculeData 实例
        mol_data = utils.extract_molecule_data(mol_ptr)

    finally:
        # 3. 释放 C++ 原始指针
        core.free_obmol_ptr(mol_ptr)

    # 4. 在 Python 侧重组
    obmol = ob.OBMol()
    obmol.BeginModify()

    # 4.1 添加原子 (使用点号访问属性，IDE 会有提示)
    # mol_data.atoms 是 _core.AtomData 对象的列表
    for atom_info in mol_data.atoms:
        a: ob.OBAtom = obmol.NewAtom()
        a.SetAtomicNum(atom_info.atomic_num)
        a.SetFormalCharge(atom_info.formal_charge)
        a.SetSpinMultiplicity(atom_info.radical_num)
        a.SetVector(atom_info.x, atom_info.y, atom_info.z)

    # 4.2 添加键
    # mol_data.bonds 是 _core.BondData 对象的列表
    for bond_info in mol_data.bonds:
        obmol.AddBond(bond_info.begin_atom_idx, bond_info.end_atom_idx, bond_info.order)

    obmol.EndModify()

    return pybel.Molecule(obmol)


# 默认可以设为 WARN
set_log_level(core.LogLevel.WARN)
