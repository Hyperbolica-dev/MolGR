"""
Author: TMJ
Date: 2026-02-22 12:26:22
LastEditors: TMJ
LastEditTime: 2026-02-22 14:44:41
Description: 请填写简介
"""

import dataclasses


@dataclasses.dataclass(frozen=True)
class ElementInfo:
    symbol: str
    num_outer_electrons: int
    default_valence: int


@dataclasses.dataclass(frozen=True)
class FDSP:
    f: int
    d: int
    s: int
    p: int


@dataclasses.dataclass
class MetalAtomPosition:
    idx: int
    symbol: str
    element_idx: int
    valence: int
    radical_num: int
    position_x: float
    position_y: float
    position_z: float
