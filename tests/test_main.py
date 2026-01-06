"""
Author: TMJ
Date: 2025-12-25 18:23:16
LastEditors: TMJ
LastEditTime: 2026-01-06 20:19:24
Description: 请填写简介
"""

import importlib.util
import warnings

import pytest


warnings.filterwarnings("ignore", category=DeprecationWarning, module=".*importlib.*")
if importlib.util.find_spec("molgr"):
    from molgr import _core  # type: ignore

# ==============================================================================
# 1. 测试四面体体积 (calculate_tetrahedron_volume)
# ==============================================================================


def test_calculate_tetrahedron_volume_unit():
    """测试单位四面体的体积"""
    # 构造一个简单的四面体，顶点在坐标轴上
    # 对应的行列式计算: V = 1/6 * |det(...) |
    p1 = [1.0, 0.0, 0.0]
    p2 = [0.0, 1.0, 0.0]
    p3 = [0.0, 0.0, 1.0]
    p4 = [0.0, 0.0, 0.0]  # 原点

    vol = _core.utils.calculate_tetrahedron_volume(p1, p2, p3, p4)

    # 预期体积应该是 1/6
    assert vol == pytest.approx(1.0 / 6.0, abs=1e-9)


def test_calculate_tetrahedron_volume_coplanar():
    """测试共面点（体积应为 0）"""
    # 所有点都在 Z=0 平面上
    p1 = [0.0, 0.0, 0.0]
    p2 = [1.0, 0.0, 0.0]
    p3 = [0.0, 1.0, 0.0]
    p4 = [1.0, 1.0, 0.0]

    vol = _core.utils.calculate_tetrahedron_volume(p1, p2, p3, p4)
    assert vol == pytest.approx(0.0, abs=1e-9)


def test_calculate_tetrahedron_volume_translation_invariance():
    """测试平移不变性"""
    # 将上面的单位四面体整体平移 [10, 10, 10]
    offset = 10.0
    p1 = [1.0 + offset, 0.0 + offset, 0.0 + offset]
    p2 = [0.0 + offset, 1.0 + offset, 0.0 + offset]
    p3 = [0.0 + offset, 0.0 + offset, 1.0 + offset]
    p4 = [0.0 + offset, 0.0 + offset, 0.0 + offset]

    vol = _core.utils.calculate_tetrahedron_volume(p1, p2, p3, p4)
    assert vol == pytest.approx(1.0 / 6.0, abs=1e-9)


# ==============================================================================
# 2. 测试形状质量评分 (calculate_shape_quality)
# ==============================================================================


def test_calculate_shape_quality_ideal():
    """测试理想正四面体的评分（应为 1.0）"""
    # 构造一个正四面体
    # 顶点取自立方体的交错顶点：(1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1)
    p1 = [1.0, 1.0, 1.0]
    p2 = [1.0, -1.0, -1.0]
    p3 = [-1.0, 1.0, -1.0]
    p4 = [-1.0, -1.0, 1.0]

    score = _core.utils.calculate_shape_quality(p1, p2, p3, p4)
    assert score == pytest.approx(1.0, abs=1e-6)


def test_calculate_shape_quality_flat():
    """测试扁平四面体的评分（应接近 0.0）"""
    # 共面情况
    p1 = [0.0, 0.0, 0.0]
    p2 = [1.0, 0.0, 0.0]
    p3 = [0.0, 1.0, 0.0]
    p4 = [1.0, 1.0, 0.0]  # 共面

    score = _core.utils.calculate_shape_quality(p1, p2, p3, p4)
    assert score == pytest.approx(0.0, abs=1e-9)


def test_calculate_shape_quality_distorted():
    """测试扭曲的四面体"""
    # 稍微拉长一个顶点，分数应该小于 1 但大于 0
    p1 = [1.0, 1.0, 5.0]  # 拉长 Z 轴
    p2 = [1.0, -1.0, -1.0]
    p3 = [-1.0, 1.0, -1.0]
    p4 = [-1.0, -1.0, 1.0]

    score = _core.utils.calculate_shape_quality(p1, p2, p3, p4)
    assert 0.0 < score < 1.0


# ==============================================================================
# 3. 测试金属自由基推导 (get_possible_metal_radicals)
# ==============================================================================


@pytest.mark.parametrize(
    "metal, valence, expected",
    [
        # 测试用例 1: 铁 (Fe)
        # Fe: [Ar] 3d6 4s2 (f=0, d=6, s=2, p=0)
        # Valence 2 (<= s+p): Lost 4s2 -> d6 configuration -> spins for d6 are {4, 2, 0}
        ("Fe", 2, {4, 2, 0}),
        # Valence 3 (> s+p, <= s+p+d): Lost 4s2 3d1 -> d5 configuration -> spins for d5 are {5, 3, 1}
        ("Fe", 3, {5, 3, 1}),
        # 测试用例 2: 铜 (Cu)
        # Cu: [Ar] 3d10 4s1 (f=0, d=10, s=1, p=0)
        # Valence 1 (<= s+p): Lost 4s1 -> d10 -> spins {0}
        ("Cu", 1, {0}),
        # Valence 2 (> s+p): Lost 4s1 3d1 -> d9 -> spins {1}
        ("Cu", 2, {1}),
        # 测试用例 3: 钠 (Na)
        # Na: [Ne] 3s1 (f=0, d=0, s=1, p=0)
        # Valence 1: Lost 3s1 -> Noble gas -> spins {0}
        ("Na", 1, {0}),
    ],
)
def test_get_possible_metal_radicals(metal, valence, expected):
    """参数化测试不同金属和价态的自由基情况"""
    result = _core.consts.get_possible_metal_radicals(metal, valence)
    print(f"Testing {metal} with valence {valence}: {result}")
    assert result == expected


def test_get_possible_metal_radicals_invalid():
    """测试不存在的金属"""
    result = _core.consts.get_possible_metal_radicals("UUnobtainium", 1)
    assert result == set()


def test_get_possible_metal_radicals_high_valence():
    """测试极高价态 (超过 s+p+d+f)"""
    # 例如 Li (s=1), valence=10
    # 根据 C++ 逻辑应该返回 {0} (f%2)，或者抛异常，具体看实现
    # 这里我们之前的实现是 return {f % 2}
    result = _core.consts.get_possible_metal_radicals("Li", 10)
    assert result == set()
