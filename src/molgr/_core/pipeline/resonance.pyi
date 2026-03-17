"""
Fallback-aligned resonance helpers
"""

from __future__ import annotations

import typing

def get_radical_resonances_smi(smiles: str) -> list[str]:
    """
    Enumerate radical resonance structures from a SMILES string.
    """

def process_resonance_smi(
    smiles: str, charge: typing.SupportsInt | typing.SupportsIndex
) -> tuple[str, int]:
    """
    Process one resonance step on SMILES and return updated token plus charge.
    """
