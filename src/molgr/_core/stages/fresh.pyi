"""
Fallback-aligned fresh stage helpers
"""

from __future__ import annotations

import typing

__all__: list[str] = [
    "assign_charge_radical_for_atom_ptr",
    "assign_radical_dots_ptr",
    "fresh_omol_charge_radical_ptr",
]

def assign_charge_radical_for_atom_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    atom_idx: typing.SupportsInt | typing.SupportsIndex,
) -> None: ...
def assign_radical_dots_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    atom_idx: typing.SupportsInt | typing.SupportsIndex,
) -> int: ...
def fresh_omol_charge_radical_ptr(mol_ptr: typing.SupportsInt | typing.SupportsIndex) -> None:
    """
    Apply fresh.fresh_omol_charge_radical to an existing OBMol.
    """
