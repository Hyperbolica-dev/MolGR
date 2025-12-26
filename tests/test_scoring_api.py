"""
Author: TMJ
Date: 2025-12-26 14:06:10
LastEditors: TMJ
LastEditTime: 2025-12-26 15:45:13
Description: 请填写简介
"""

import importlib.util
import warnings

import pytest


warnings.filterwarnings("ignore", category=DeprecationWarning, module=".*importlib.*")
if importlib.util.find_spec("molgr"):
    from molgr import omol_score
if importlib.util.find_spec("openbabel"):
    from openbabel import openbabel as ob
    from openbabel import pybel


def test_omol_score_with_obmol():
    """测试直接传递 openbabel.OBMol 对象"""
    # 1. 创建转换器
    conv = ob.OBConversion()
    conv.SetInFormat("smi")

    # 2. 创建并读取分子 (苯)
    mol = ob.OBMol()
    success = conv.ReadString(mol, "c1ccccc1")
    assert success, "Failed to read SMILES"

    # 3. 调用打分函数
    score = omol_score(mol)

    # 4. 验证结果
    # 苯是对称性很高的分子，分数应该是负数（奖励）或很低
    assert isinstance(score, float)
    assert score < 0.0  # 期望它能正确计算出对称性奖励


def test_omol_score_with_pybel():
    """测试传递 pybel.Molecule 对象"""
    # 1. 使用 pybel 快捷读取
    mol = pybel.readstring("smi", "C(=O)[O+]")

    # 2. 调用打分函数
    score = omol_score(mol)

    # 3. 验证结果
    assert isinstance(score, float)
    # 只要不报错且返回浮点数，说明指针传递成功
    assert score != 0.0


def test_omol_score_consistency():
    """验证 Pybel 和 OBMol 对同一分子打分一致"""
    smi = "c1ccccc1"  # 苯

    # Pybel
    pybel_mol = pybel.readstring("smi", smi)
    score_pybel = omol_score(pybel_mol)

    # OBMol
    conv = ob.OBConversion()
    conv.SetInFormat("smi")
    ob_mol = ob.OBMol()
    conv.ReadString(ob_mol, smi)
    score_ob = omol_score(ob_mol)

    # 两者指向的底层化学结构一致，分数应相同
    assert score_pybel == pytest.approx(score_ob, abs=1e-6)


def test_omol_score_invalid_input():
    """测试非法输入"""
    with pytest.raises(TypeError):
        omol_score("This is a string, not a molecule")  # type: ignore

    with pytest.raises(TypeError):
        omol_score(123)  # type: ignore
