"""
Fallback-aligned resonance helpers
"""

from __future__ import annotations

import typing

__all__: list[str] = ["get_radical_resonances_smi", "process_resonance_smi"]

def get_radical_resonances_smi(smiles: str) -> list[str]: ...
def process_resonance_smi(
    smiles: str, charge: typing.SupportsInt | typing.SupportsIndex
) -> tuple: ...
