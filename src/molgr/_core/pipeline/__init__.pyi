"""
Pipeline-level helpers
"""

from __future__ import annotations

import typing

import molgr._core.utils

from . import reconstruct_with_metals, reconstruct_without_metals

__all__: list[str] = [
    "get_last_run_timing_breakdown_ms",
    "reconstruct_with_metals",
    "reconstruct_without_metals",
    "xyz2omol",
]

def get_last_run_timing_breakdown_ms() -> dict:
    """
    Return timing breakdown (milliseconds) for the most recent reconstruction run.
    """

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
