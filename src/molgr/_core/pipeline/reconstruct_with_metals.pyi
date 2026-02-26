"""
Fallback-aligned reconstruction helpers with metals
"""

from __future__ import annotations

import collections.abc
import typing

__all__: list[str] = [
    "MetalAtomPosition",
    "MetalHandler",
    "build_metal_states_ptr",
    "combine_metal_with_omol_ptr",
    "get_possible_metal_radicals",
    "xyz2omol",
]

class MetalAtomPosition:
    symbol: str
    def __repr__(self) -> str: ...
    @property
    def element_idx(self) -> int: ...
    @element_idx.setter
    def element_idx(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def idx(self) -> int: ...
    @idx.setter
    def idx(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def position_x(self) -> float: ...
    @position_x.setter
    def position_x(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def position_y(self) -> float: ...
    @position_y.setter
    def position_y(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def position_z(self) -> float: ...
    @position_z.setter
    def position_z(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def radical_num(self) -> int: ...
    @radical_num.setter
    def radical_num(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def valence(self) -> int: ...
    @valence.setter
    def valence(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...

class MetalHandler:
    @staticmethod
    def combine_metal_with_mol(
        mol_ptr: typing.SupportsInt | typing.SupportsIndex,
        metals: collections.abc.Sequence[MetalAtomPosition],
    ) -> None: ...
    def __init__(self, mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def generate_combinations(
        self, total_radical_electrons: typing.SupportsInt | typing.SupportsIndex
    ) -> list[list[MetalAtomPosition]]: ...
    def strip_metals(self, mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> str: ...

def build_metal_states_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    atom_idx: typing.SupportsInt | typing.SupportsIndex,
) -> list[MetalAtomPosition]: ...
def combine_metal_with_omol_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    metals: collections.abc.Sequence[MetalAtomPosition],
) -> None: ...
def get_possible_metal_radicals(
    metal: str, valence: typing.SupportsInt | typing.SupportsIndex
) -> set[int]: ...
def xyz2omol(
    xyz_block: str,
    total_charge: typing.SupportsInt | typing.SupportsIndex,
    total_radical_electrons: typing.SupportsInt | typing.SupportsIndex,
) -> typing.Any: ...
