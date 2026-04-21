"""
Development-only helpers for metal reconstruction internals
"""

from __future__ import annotations

import collections.abc
import typing

__all__: list[str] = ["MetalAtomPosition", "build_metal_states_ptr", "combine_metal_with_omol_ptr"]

class MetalAtomPosition:
    symbol: str
    def __init__(self) -> None: ...
    def __repr__(self) -> str:
        """
        Return a concise debug representation of MetalAtomPosition.
        """
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

def build_metal_states_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    atom_idx: typing.SupportsInt | typing.SupportsIndex,
) -> list[MetalAtomPosition]:
    """
    Build candidate metal states for a specific atom index on an OBMol pointer.
    """

def combine_metal_with_omol_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    metals: collections.abc.Sequence[MetalAtomPosition],
) -> None:
    """
    Combine metal states with an existing OBMol pointer in place.
    """
