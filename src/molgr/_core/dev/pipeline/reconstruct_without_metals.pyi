"""
Development-only helpers for no-metal reconstruction internals
"""

from __future__ import annotations

import typing

import molgr.config

__all__: list[str] = [
    "debug_neighbor_radical_seeds",
    "debug_prepared_no_metal_seed",
    "debug_resonance_candidate_summaries",
]

def debug_neighbor_radical_seeds(
    xyz_block: str,
    total_charge: typing.SupportsInt | typing.SupportsIndex = 0,
    total_radical_electrons: typing.SupportsInt | typing.SupportsIndex = 0,
    exact_discrepancy: typing.SupportsInt | typing.SupportsIndex | None = None,
) -> typing.Any:
    """
    Return production C++ neighboring-radical seeds.
    """

def debug_prepared_no_metal_seed(
    xyz_block: str,
    total_charge: typing.SupportsInt | typing.SupportsIndex = 0,
    total_radical_electrons: typing.SupportsInt | typing.SupportsIndex = 0,
) -> typing.Any:
    """
    Return the production C++ prepared no-metal seed.
    """

def debug_resonance_candidate_summaries(
    xyz_block: str,
    total_charge: typing.SupportsInt | typing.SupportsIndex = 0,
    total_radical_electrons: typing.SupportsInt | typing.SupportsIndex = 0,
    *,
    config: molgr.config.MolGRConfig | None = None,
) -> list:
    """
    Return production C++ no-metal resonance candidate summaries.
    """
