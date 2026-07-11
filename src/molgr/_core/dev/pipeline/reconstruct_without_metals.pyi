"""
Development-only helpers for no-metal reconstruction internals
"""

from __future__ import annotations

import typing

__all__: list[str] = [
    "debug_linear_pipeline_state",
    "debug_linear_pipeline_trace",
    "debug_no_metal_candidate_states",
    "debug_processed_root_resonance",
    "debug_resonance_candidate_summaries",
]

def debug_no_metal_candidate_states(
    xyz_block: str,
    total_charge: typing.SupportsInt | typing.SupportsIndex = 0,
    total_radical_electrons: typing.SupportsInt | typing.SupportsIndex = 0,
) -> typing.Any:
    """
    Return C++ no-metal candidate states for parity debugging.
    """

def debug_linear_pipeline_state(
    xyz_block: str,
    total_charge: typing.SupportsInt | typing.SupportsIndex = 0,
    total_radical_electrons: typing.SupportsInt | typing.SupportsIndex = 0,
) -> typing.Any:
    """
    Return the C++ no-metal linear-pipeline state for parity debugging.
    """

def debug_linear_pipeline_trace(
    xyz_block: str,
    total_charge: typing.SupportsInt | typing.SupportsIndex = 0,
    total_radical_electrons: typing.SupportsInt | typing.SupportsIndex = 0,
) -> typing.Any:
    """
    Return per-stage C++ no-metal linear-pipeline snapshots for parity debugging.
    """

def debug_processed_root_resonance(
    xyz_block: str,
    total_charge: typing.SupportsInt | typing.SupportsIndex = 0,
    total_radical_electrons: typing.SupportsInt | typing.SupportsIndex = 0,
) -> typing.Any:
    """
    Process the linear no-metal state as a root resonance candidate for debugging.
    """

def debug_resonance_candidate_summaries(
    xyz_block: str,
    total_charge: typing.SupportsInt | typing.SupportsIndex = 0,
    total_radical_electrons: typing.SupportsInt | typing.SupportsIndex = 0,
    *,
    config: typing.Any = None,
) -> list:
    """
    Return C++ no-metal resonance candidates for parity debugging.
    """
