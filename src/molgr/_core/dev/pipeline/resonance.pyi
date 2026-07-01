"""
Development-only helpers for resonance internals
"""

from __future__ import annotations

import typing

__all__: list[str] = [
    "get_radical_resonances_ptr",
    "get_radical_resonances_smi",
    "process_resonance_ptr",
    "process_resonance_smi",
    "smiles_token_ptr",
]

def get_radical_resonances_ptr(mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> list[int]:
    """
    Get radical resonance OBMol pointers for an input OBMol pointer.
    """

def get_radical_resonances_smi(smiles: str) -> list[str]:
    """
    Enumerate radical resonance structures from a SMILES string.
    """

def process_resonance_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    charge: typing.SupportsInt | typing.SupportsIndex,
) -> tuple:
    """
    Process resonance on an OBMol pointer and return (new_ptr, updated_charge).
    """

def process_resonance_smi(
    smiles: str, charge: typing.SupportsInt | typing.SupportsIndex
) -> tuple[str, int]:
    """
    Process one resonance step on SMILES and return updated token plus charge.
    """

def smiles_token_ptr(mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> str:
    """
    Return canonical first-token SMILES string for an OBMol pointer.
    """
