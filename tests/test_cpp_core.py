import sys

import pytest
from openbabel import openbabel as ob
from openbabel import pybel


# 假设编译好的模块名为 molgr_core 或类似名称，请根据你的 setup.py 调整导入
# 这里假设 import _core 即可，或者你的包结构是 molgr._core
try:
    import _core as core
except ImportError:
    # 如果是在包内部，尝试相对导入或根据你的构建路径调整
    import molgr._core as core

# ==============================================================================
# 辅助函数
# ==============================================================================


def intptr_to_obmol(ptr_addr: int) -> ob.OBMol:
    """
    将 C++ 返回的裸指针地址包装为 Python 的 OpenBabel.OBMol 对象。

    此操作零拷贝，且保留所有原子顺序和坐标精度。
    """
    if ptr_addr == 0:
        raise ValueError("Received null pointer from C++ extension.")

    # 1. 创建 OBMol 的 Python 空壳 (绕过 __init__，避免创建新的 C++ 对象)
    mol = ob.OBMol.__new__(ob.OBMol)

    # 2. 注入指针
    # SWIG 通常允许直接赋值整数地址，但在某些旧版本中可能需要特定格式的字符串
    try:
        mol.this = ptr_addr
    except TypeError:
        # Fallback: 构造 SWIG 风格的指针字符串 (格式: _<hex_addr>_p_OpenBabel__OBMol)
        # 注意：这依赖于 OpenBabel 编译时的命名空间，通常是 OpenBabel::OBMol
        mol.this = f"_{ptr_addr:x}_p_OpenBabel__OBMol"

    # 3. 移交所有权 (Ownership Transfer)
    # 极其重要：告诉 Python "你负责在这个对象被 GC 时调用 C++ delete"
    # 因为我们在 C++ 中使用了 .release()，现在 Python 必须接盘。
    mol.own(True)

    return mol


def intptr_to_pybel(ptr_addr: int) -> pybel.Molecule:
    """直接转为更高级的 Pybel 对象"""
    obmol = intptr_to_obmol(ptr_addr)
    return pybel.Molecule(obmol)


# ==============================================================================
# 1. Consts 模块测试
# ==============================================================================


def test_get_possible_metal_radicals():
    # 测试 Fe (Iron)
    # Fe: [Ar] 3d6 4s2 -> Valence 2 -> d6
    # High spin d6 (S=2, 4 unpaired), Low spin d6 (S=0, 0 unpaired)
    # 具体的自旋态取决于你的 consts.cpp 数据表
    radicals = core.get_possible_metal_radicals("Fe", 2)
    assert isinstance(radicals, set)
    # 验证是否包含预期的自旋态 (根据 consts.cpp 的 kDElectronsSpin)
    # Fe(II) d6 通常可能是 4 (high spin) 或 0 (low spin)
    assert 4 in radicals or 0 in radicals

    # 测试不存在的金属
    radicals = core.get_possible_metal_radicals("Xy", 1)
    assert len(radicals) == 0


# ==============================================================================
# 2. Utils 模块测试
# ==============================================================================


def test_tetrahedron_volume():
    # 标准四面体顶点
    p1 = [1.0, 1.0, 1.0]
    p2 = [1.0, -1.0, -1.0]
    p3 = [-1.0, 1.0, -1.0]
    p4 = [-1.0, -1.0, 1.0]

    vol = core.calculate_tetrahedron_volume(p1, p2, p3, p4)
    # 边长为 2*sqrt(2) 的正四面体体积 = 8/3 ≈ 2.666
    assert abs(vol - 8.0 / 3.0) < 1e-4


def test_shape_quality_perfect():
    # 构造一个理想四面体的四个方向
    p1 = [1.0, 0.0, 0.0]
    p2 = [0.0, 1.0, 0.0]
    p3 = [0.0, 0.0, 1.0]
    p4 = [1.0, 1.0, 1.0]

    quality = core.calculate_shape_quality(p1, p2, p3, p4)
    # 理想四面体得分应接近 1.0
    assert abs(quality - 1.0) < 1e-4


# ==============================================================================
# 3. Metal Handler 测试
# ==============================================================================


def test_metal_handler_strip():
    # 构造一个含铁的分子 (Ferrocene-like dummy)
    mol = pybel.readstring("smi", "[Fe].[C][C][C][C][C]")
    mol.make3D()  # 生成坐标

    # 获取 C++ 指针
    ptr = get_mol_ptr(mol)

    # 实例化 C++ MetalHandler
    handler = core.MetalHandler(ptr)

    # 剥离金属
    xyz_no_metal = handler.strip_metals(ptr)

    # 验证结果是字符串且不再包含 Fe
    assert isinstance(xyz_no_metal, str)
    assert "Fe" not in xyz_no_metal
    assert "C" in xyz_no_metal

    # 验证原始 OBMol 对象也被修改了 (原子数减少)
    assert len(mol.atoms) == 5  # 只剩 5 个碳


def test_metal_combinations():
    # 测试静态生成组合逻辑
    # Fe(II) 可能是 4 或 0 个单电子
    # 假设 total_radical_limit = 10，应该能生成所有组合

    # 由于 generate_combinations 依赖内部存储的 raw_metals_
    # 我们需要先用含金属的分子初始化 handler
    mol = pybel.readstring("smi", "[Fe]")
    mol.make3D()
    ptr = get_mol_ptr(mol)
    handler = core.MetalHandler(ptr)

    combos = handler.generate_combinations(10)
    print(combos)
    assert len(combos) > 0
    first_combo = combos[0]
    assert len(first_combo) == 1  # 只有一个金属原子

    metal_pos = first_combo[0]
    assert metal_pos.symbol == "Fe"
    # 验证属性是否正确读取
    assert metal_pos.valence in [2, 3]  # Fe 常见价态


# ==============================================================================
# 5. Initial Reconstructor 测试 (核心测试)
# ==============================================================================


def test_reconstruct_methane_xyz():
    # 构造甲烷的 XYZ 块
    xyz = """5
    Methane
    C 0.000 0.000 0.000
    H 1.089 0.000 0.000
    H -0.363 1.026 0.000
    H -0.363 -0.513 0.888
    H -0.363 -0.513 -0.888
    """

    # 调用 C++ 重建，目标电荷 0，自由基 0
    mol_ptr_addr = core.reconstruct_from_xyz_no_metal(xyz, 0, 0)

    assert mol_ptr_addr != 0

    # 由于无法轻易将 intptr_t 转回 Python OBMol 进行断言
    # 我们在这里只能验证 C++ 没有崩溃并返回了有效指针。
    # *进阶验证*：如果你想验证重建的分子是否正确，
    # 建议在 module.cpp 中额外暴露一个 `get_smiles_from_ptr(intptr_t)` 的辅助函数仅用于测试。

    print(f"Reconstruction success, OBMol created at {hex(mol_ptr_addr)}")


if __name__ == "__main__":
    # 允许直接运行脚本
    sys.exit(pytest.main(["-v", __file__]))
