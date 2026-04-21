"""
Fallback-aligned reconstruction helpers with metals
"""

from __future__ import annotations

import typing

import molgr._core.utils

__all__: list[str] = ["get_possible_metal_radicals", "xyz2omol"]

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
