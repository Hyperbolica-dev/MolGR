"""
Fallback-aligned clean stage helpers
"""

from __future__ import annotations

import typing

__all__: list[str] = [
    "clean_carbene_neighbor_unsaturated_ptr",
    "clean_neighbor_radicals_ptr",
    "clean_resonances_ptr",
]

def clean_carbene_neighbor_unsaturated_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
) -> None:
    """
    Apply clean.clean_carbene_neighbor_unsaturated to an existing OBMol.
    """

def clean_neighbor_radicals_ptr(mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> None:
    """
    Apply clean.clean_neighbor_radicals to an existing OBMol.
    """

def clean_resonances_ptr(mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> None:
    """
    Apply clean.clean_resonances to an existing OBMol.
    """
