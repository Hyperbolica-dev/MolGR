"""
Fallback-aligned clean stage helpers
"""

from __future__ import annotations

import typing

__all__: list[str] = ["clean_resonances_ptr"]

def clean_resonances_ptr(mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> None:
    """
    Apply clean.clean_resonances to an existing OBMol.
    """
