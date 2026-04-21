"""
Utilities and Data Structures
"""

from __future__ import annotations

import collections.abc
import typing

__all__: list[str] = [
    "AtomData",
    "BondData",
    "MoleculeData",
    "calculate_shape_quality",
    "calculate_tetrahedron_volume",
    "extract_molecule_data",
    "get_possible_metal_radicals",
    "molecule_data_to_obmol_ptr",
]

class AtomData:
    def __repr__(self) -> str:
        """
        Return a concise debug representation of AtomData.
        """
    @property
    def atomic_num(self) -> int: ...
    @atomic_num.setter
    def atomic_num(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def formal_charge(self) -> int: ...
    @formal_charge.setter
    def formal_charge(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def radical_num(self) -> int: ...
    @radical_num.setter
    def radical_num(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def x(self) -> float: ...
    @x.setter
    def x(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def y(self) -> float: ...
    @y.setter
    def y(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...
    @property
    def z(self) -> float: ...
    @z.setter
    def z(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None: ...

class BondData:
    def __repr__(self) -> str:
        """
        Return a concise debug representation of BondData.
        """
    @property
    def begin_atom_idx(self) -> int: ...
    @begin_atom_idx.setter
    def begin_atom_idx(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def end_atom_idx(self) -> int: ...
    @end_atom_idx.setter
    def end_atom_idx(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def order(self) -> int: ...
    @order.setter
    def order(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...

class MoleculeData:
    @property
    def atoms(self) -> list[AtomData]: ...
    @atoms.setter
    def atoms(self, arg0: collections.abc.Sequence[AtomData]) -> None: ...
    @property
    def bonds(self) -> list[BondData]: ...
    @bonds.setter
    def bonds(self, arg0: collections.abc.Sequence[BondData]) -> None: ...
    @property
    def total_charge(self) -> int: ...
    @total_charge.setter
    def total_charge(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    @property
    def total_radical_num(self) -> int: ...
    @total_radical_num.setter
    def total_radical_num(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None: ...

def calculate_shape_quality(
    p1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    p2: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    p3: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    p4: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
) -> float:
    """
    Calculate the shape quality score of a tetrahedron.

    Parameters:
        p1 (list[float]): Coordinates of the first atom.
        p2 (list[float]): Coordinates of the second atom.
        p3 (list[float]): Coordinates of the third atom.
        p4 (list[float]): Coordinates of the fourth atom.

    Returns:
        float: Quality score between 0.0 (coplanar/bad) and 1.0 (ideal).
    """

def calculate_tetrahedron_volume(
    p1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    p2: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    p3: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
    p4: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex],
) -> float:
    """
    Calculate the volume of a tetrahedron defined by 4 points.

    Parameters:
        p1 (list[float]): Coordinates of the first atom.
        p2 (list[float]): Coordinates of the second atom.
        p3 (list[float]): Coordinates of the third atom.
        p4 (list[float]): Coordinates of the fourth atom.

    Returns:
        float: The volume.
    """

def extract_molecule_data(mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> MoleculeData:
    """
    Extracts OBMol content into a structured object.
    """

def get_possible_metal_radicals(
    metal: str, valence: typing.SupportsInt | typing.SupportsIndex
) -> set[int]:
    """
    Get possible radical electron counts for a metal given its valence.

    Args:
        metal (str): The chemical symbol (e.g., "Fe").
        valence (int): The oxidation state.

    Returns:
        set[int]: A set of possible unpaired electron counts.
    """

def molecule_data_to_obmol_ptr(molecule_data: MoleculeData) -> int:
    """
    Converts MoleculeData to a newly allocated OBMol pointer. Free it with _core.free_obmol_ptr.
    """
