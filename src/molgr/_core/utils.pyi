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
]

class AtomData:
    def __repr__(self) -> str: ...
    @property
    def atomic_num(self) -> int: ...
    @atomic_num.setter
    def atomic_num(self, arg0: typing.SupportsInt) -> None: ...
    @property
    def formal_charge(self) -> int: ...
    @formal_charge.setter
    def formal_charge(self, arg0: typing.SupportsInt) -> None: ...
    @property
    def radical_num(self) -> int: ...
    @radical_num.setter
    def radical_num(self, arg0: typing.SupportsInt) -> None: ...
    @property
    def x(self) -> float: ...
    @x.setter
    def x(self, arg0: typing.SupportsFloat) -> None: ...
    @property
    def y(self) -> float: ...
    @y.setter
    def y(self, arg0: typing.SupportsFloat) -> None: ...
    @property
    def z(self) -> float: ...
    @z.setter
    def z(self, arg0: typing.SupportsFloat) -> None: ...

class BondData:
    def __repr__(self) -> str: ...
    @property
    def begin_atom_idx(self) -> int: ...
    @begin_atom_idx.setter
    def begin_atom_idx(self, arg0: typing.SupportsInt) -> None: ...
    @property
    def end_atom_idx(self) -> int: ...
    @end_atom_idx.setter
    def end_atom_idx(self, arg0: typing.SupportsInt) -> None: ...
    @property
    def order(self) -> int: ...
    @order.setter
    def order(self, arg0: typing.SupportsInt) -> None: ...

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
    def total_charge(self, arg0: typing.SupportsInt) -> None: ...
    @property
    def total_radical_num(self) -> int: ...
    @total_radical_num.setter
    def total_radical_num(self, arg0: typing.SupportsInt) -> None: ...

def calculate_shape_quality(
    p1: collections.abc.Sequence[typing.SupportsFloat],
    p2: collections.abc.Sequence[typing.SupportsFloat],
    p3: collections.abc.Sequence[typing.SupportsFloat],
    p4: collections.abc.Sequence[typing.SupportsFloat],
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
    p1: collections.abc.Sequence[typing.SupportsFloat],
    p2: collections.abc.Sequence[typing.SupportsFloat],
    p3: collections.abc.Sequence[typing.SupportsFloat],
    p4: collections.abc.Sequence[typing.SupportsFloat],
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

def extract_molecule_data(mol_ptr: typing.SupportsInt) -> MoleculeData:
    """
    Extracts OBMol content into a structured object.
    """
