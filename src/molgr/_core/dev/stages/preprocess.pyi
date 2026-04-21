"""
Fallback-aligned preprocess stage helpers
"""

from __future__ import annotations

import typing

__all__: list[str] = ["make_connections_ptr", "pre_clean_ptr", "validate_omol_ptr"]

def make_connections_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    factor: typing.SupportsFloat | typing.SupportsIndex = 1.4,
) -> bool:
    """
    Apply preprocess.make_connections to an existing OBMol.

    Args:
        mol_ptr: int address of OpenBabel::OBMol
        factor: distance factor (default matches python fallback)
    """

def pre_clean_ptr(mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> bool:
    """
    Apply preprocess.pre_clean to an existing OBMol.
    """

def validate_omol_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    total_charge: typing.SupportsInt | typing.SupportsIndex,
    total_radical: typing.SupportsInt | typing.SupportsIndex,
) -> bool:
    """
    Validate conservation of total charge and radical electrons.

    Args:
        mol_ptr: int address of OpenBabel::OBMol
        total_charge: expected total formal charge
        total_radical: expected total radical electrons
    """
