import importlib.util

import pytest


if importlib.util.find_spec("molgr"):
    from molgr import _core  # type: ignore


# ==============================================================================
# 1. 对称性评分 (Symmetry Penalty)
# ==============================================================================


def test_symmetry_low():
    """测试低对称性分子 (甲醇)"""
    score_1 = _core.scoring.test_symmetry_penalty("c1ccccc1[O-]")
    score_2 = _core.scoring.test_symmetry_penalty("C1=CC=C[C-]C1=O")
    assert score_2 > score_1


# ==============================================================================
# 2. 理化性质评分 (PhysChem Penalty)
# ==============================================================================
def test_physchem_charge_penalty():
    """测试不合理的电荷分布"""
    # [Na+] 是合理的，罚分应较小或为0 (取决于具体的 EN 规则)
    score_na = _core.scoring.test_physchem_penalty("[Na+]")

    # [C+] 碳正离子，极其不稳定，应该有高罚分
    score_c_plus = _core.scoring.test_physchem_penalty("[C+]")

    assert score_c_plus > score_na


def test_physchem_coulombic_repulsion():
    """测试库仑斥力"""
    # 两个正电荷相连 [NH3+][NH3+]
    # 应该触发 +15.0 的 Repulsion Penalty
    score = _core.scoring.test_physchem_penalty("[NH3+][NH3+]")
    assert score >= 15.0


# ==============================================================================
# 3. 几何偏差评分 (Geometry Deviation)
# ==============================================================================
def test_deviation_tetrahedral_ideal():
    """测试理想正四面体"""
    # 标准甲烷坐标
    xyz = """5
    Methane
    C 0.000 0.000 0.000
    H 0.629 0.629 0.629
    H -0.629 -0.629 0.629
    H 0.629 -0.629 -0.629
    H -0.629 0.629 -0.629
    """
    # C 是第1个原子
    score = _core.scoring.test_deviation_score(xyz, 1)
    # 应该是完美的，偏差接近 0
    assert score == pytest.approx(0.0, abs=1e-3)


def test_deviation_distorted():
    """测试严重扭曲的几何"""
    # 将一个 H 移到非常不合理的位置
    xyz = """5
    Distorted
    C 0.000 0.000 0.000
    H 0.629 0.629 0.629
    H -0.629 -0.629 0.629
    H 0.629 -0.629 -0.629
    H 0.100 0.100 0.100  # 这里的 H 离 C 太近且角度错误
    """
    score = _core.scoring.test_deviation_score(xyz, 1)
    assert score > 0.1  # 肯定会有显著偏差


# ==============================================================================
# 4. 总评分 (Total Score)
# ==============================================================================
def test_total_score_comparison():
    """比较好分子和坏分子的总分"""
    xyz_good = """5
    Methane
    C 0.000 0.000 0.000
    H 0.629 0.629 0.629
    H -0.629 -0.629 0.629
    H 0.629 -0.629 -0.629
    H -0.629 0.629 -0.629
    """

    # 坏分子：H2 分子但是两个 H 距离极近 (0.2A) 导致斥力
    xyz_bad = """2
    Bad H2
    H 0.000 0.000 0.000
    H 0.200 0.000 0.000
    """

    score_good = _core.scoring.test_total_score(xyz_good)
    score_bad = _core.scoring.test_total_score(xyz_bad)

    # 越低越好，所以好分子的分数应该小于坏分子
    # 注意：如果坏分子触发了极大的斥力罚分，score_bad 应该是正的大数
    # 而好分子可能有对称性奖励，甚至是负数
    assert score_good < score_bad
