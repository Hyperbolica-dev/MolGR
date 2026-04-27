"""
Fallback-aligned reconstruction helpers with metals
"""

from __future__ import annotations

import typing

import molgr._core.utils

__all__: list[str] = ["xyz2omol"]

def xyz2omol(
    xyz_block: str,
    total_charge: typing.SupportsInt | typing.SupportsIndex = 0,
    total_radical_electrons: typing.SupportsInt | typing.SupportsIndex = 0,
    *,
    config: typing.Any = None,
) -> molgr._core.utils.MoleculeData:
    """
    Reconstruct molecule data from XYZ with metal-aware pipeline.
    """
