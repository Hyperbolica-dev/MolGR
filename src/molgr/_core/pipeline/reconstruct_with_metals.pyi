"""
Fallback-aligned reconstruction helpers with metals
"""

from __future__ import annotations

import collections.abc
import typing

import molgr._core.utils

class MetalAtomPosition:
    symbol: str
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

class MetalHandler:
    @staticmethod
    def combine_metal_with_mol(
        mol_ptr: typing.SupportsInt | typing.SupportsIndex,
        metals: collections.abc.Sequence[MetalAtomPosition],
    ) -> None:
        """
        Apply metal positions to an OBMol pointer in place.
        """
    def __init__(self, mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Create a MetalHandler for an existing OBMol pointer.
        """
    def generate_combinations(
        self, total_radical_electrons: typing.SupportsInt | typing.SupportsIndex
    ) -> list[list[MetalAtomPosition]]:
        """
        Enumerate candidate metal attachment combinations for a radical budget.
        """
    def strip_metals(self, mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> str:
        """
        Remove metal atoms from the target molecule and return encoded metadata.
        """

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

def get_possible_metal_radicals(
    metal: str, valence: typing.SupportsInt | typing.SupportsIndex
) -> set[int]:
    """
    Get possible radical electron counts for a metal and valence.
    """

def xyz2omol(
    xyz_block: str,
    total_charge: typing.SupportsInt | typing.SupportsIndex = 0,
    total_radical_electrons: typing.SupportsInt | typing.SupportsIndex = 0,
) -> molgr._core.utils.MoleculeData:
    """
    Reconstruct molecule data from XYZ with metal-aware pipeline.
    """
