"""
Pipeline-level helpers
"""

from __future__ import annotations

import typing

from . import reconstruct_with_metals, reconstruct_without_metals, resonance

__all__: list[str] = [
    "get_last_run_timing_breakdown_ms",
    "get_radical_resonances_ptr",
    "process_resonance_ptr",
    "reconstruct_with_metals",
    "reconstruct_without_metals",
    "resonance",
    "smiles_token_ptr",
]

def get_last_run_timing_breakdown_ms() -> dict:
    """
    Return timing breakdown (milliseconds) for the most recent reconstruction run.
    """

def get_radical_resonances_ptr(mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> list[int]:
    """
    Get radical resonance OBMol pointers for an input OBMol pointer.
    """

def process_resonance_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    charge: typing.SupportsInt | typing.SupportsIndex,
) -> tuple:
    """
    Process resonance on an OBMol pointer and return (new_ptr, updated_charge).
    """

def smiles_token_ptr(mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> str:
    """
    Return canonical first-token SMILES string for an OBMol pointer.
    """
