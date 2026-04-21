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


__all__ = ["ElementInfo", "FDSP", "MetalAtomPosition"]
