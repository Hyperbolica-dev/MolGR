# pyright: reportMissingImports=false

import pytest
from openbabel import pybel

from molgr import _core  # type: ignore


# 辅助函数：获取 C++ 指针
def get_ptr(mol):
    if isinstance(mol, pybel.Molecule):
        mol = mol.OBMol
    if hasattr(mol, "this"):
        return int(mol.this)
    return int(mol)


def test_metal_identification_and_strip():
    """测试金属识别与剥离"""
    # 构造一个含金属的分子：例如 CH3-Li (简化模型)
    # Li 是金属，C 是有机部分
    # Li 在原点，C 在 (2,0,0)
    xyz = """2
    Test
    Li 0.0 0.0 0.0
    C  2.0 0.0 0.0
    """
    mol = pybel.readstring("xyz", xyz)
    ptr = get_ptr(mol)

    # 1. 初始化 Handler (应该识别到 Li)
    handler = _core.pipeline.reconstruct_with_metals.MetalHandler(ptr)

    # 2. 剥离金属
    # 注意：StripMetals 会修改传入的 mol 对象
    organic_xyz = handler.strip_metals(ptr)

    # 验证 mol 中只剩下一个原子 (C)
    assert mol.OBMol.NumAtoms() == 1
    assert mol.atoms[0].atomicnum == 6  # Carbon

    # 验证返回的 XYZ 字符串只包含 1 个原子
    assert "1" in organic_xyz.split("\n")[0]
    assert "Li" not in organic_xyz


def test_generate_combinations():
    """测试金属价态组合生成"""
    # 构造含 Pd 的分子 (Pd 有多种价态)
    xyz = """1
    Pd
    Pd 0.0 0.0 0.0
    """
    mol = pybel.readstring("xyz", xyz)
    handler = _core.pipeline.reconstruct_with_metals.MetalHandler(get_ptr(mol))

    # 生成组合，限制总自由基
    # Pd 可能的价态 (从 consts.cpp): Prior: 0, 2, 4.
    # Pd(0) -> rad 0
    # Pd(2) -> rad 0 or 2 (取决于 consts 定义)
    combinations = handler.generate_combinations(total_radical_electrons=10)

    assert len(combinations) > 0
    # 检查第一个组合的结构
    first_combo = combinations[0]
    assert isinstance(first_combo, list)
    assert len(first_combo) == 1

    metal_pos = first_combo[0]
    assert metal_pos.symbol == "Pd"
    assert isinstance(metal_pos.valence, int)
    assert isinstance(metal_pos.radical_num, int)


def test_combine_and_renumber():
    """核心测试：金属回填与重排序"""
    _core.set_log_level(_core.LogLevel.DEBUG)  # 开启日志以便观察

    # 构造一个顺序敏感的分子：C - Li - O
    # 索引 (1-based): 1:C, 2:Li, 3:O
    # Li 是金属，剥离后应剩下 C 和 O
    # 回填后，Li 必须回到索引 2 的位置
    xyz = """3
    Renumber Test
    C  1.0 0.0 0.0
    Li 2.0 0.0 0.0
    O  3.0 0.0 0.0
    """
    mol = pybel.readstring("xyz", xyz)
    ptr = get_ptr(mol)

    # 记录原始地址/属性以供验证
    original_c_x = mol.atoms[0].coords[0]  # 1.0
    original_li_x = mol.atoms[1].coords[0]  # 2.0
    original_o_x = mol.atoms[2].coords[0]  # 3.0

    handler = _core.pipeline.reconstruct_with_metals.MetalHandler(ptr)

    # 1. 剥离金属 (Li 被移除)
    handler.strip_metals(ptr)
    assert mol.OBMol.NumAtoms() == 2
    # 此时 C 应该是索引 1, O 应该是索引 2 (OpenBabel 删除中间原子后，后面的会自动前移)
    assert mol.atoms[0].atomicnum == 6  # C
    assert mol.atoms[1].atomicnum == 8  # O

    # 2. 生成一个恢复方案 (原样恢复)
    combos = handler.generate_combinations(10)
    target_combo = None
    for combo in combos:
        if combo[0].symbol == "Li":
            target_combo = combo
            break

    assert target_combo is not None

    # 3. 回填金属
    _core.pipeline.reconstruct_with_metals.MetalHandler.combine_metal_with_mol(ptr, target_combo)

    # 4. 验证结果
    assert mol.OBMol.NumAtoms() == 3

    # 验证原子顺序是否恢复 (1:C, 2:Li, 3:O)
    # Atom indices are 1-based in OpenBabel API, 0-based in Pybel lists
    atom1 = mol.atoms[0]
    atom2 = mol.atoms[1]
    atom3 = mol.atoms[2]

    print(f"Atom 1: {atom1.type} at {atom1.coords}")
    print(f"Atom 2: {atom2.type} at {atom2.coords}")
    print(f"Atom 3: {atom3.type} at {atom3.coords}")

    assert atom1.atomicnum == 6  # C
    assert atom2.atomicnum == 3  # Li (应该在中间！)
    assert atom3.atomicnum == 8  # O

    # 验证坐标未乱
    assert atom1.coords[0] == pytest.approx(original_c_x)
    assert atom2.coords[0] == pytest.approx(original_li_x)
    assert atom3.coords[0] == pytest.approx(original_o_x)
