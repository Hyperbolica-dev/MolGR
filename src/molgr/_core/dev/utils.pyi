"""
Development-only utility helpers
"""

from __future__ import annotations

import typing

__all__: list[str] = [
    "test_deviation_score",
    "test_physchem_penalty",
    "test_symmetry_penalty",
    "test_total_score",
]

def test_deviation_score(
    xyz_block: str, atom_idx: typing.SupportsInt | typing.SupportsIndex
) -> float:
    """
    Calculate geometry deviation for atom (1-based index) from XYZ (For Testing)
    """

def test_physchem_penalty(smiles: str) -> float:
    """
    Calculate PhysChem penalty from SMILES (For Testing)
    """

def test_symmetry_penalty(smiles: str) -> float:
    """
    Calculate symmetry penalty from SMILES (For Testing)
    """

def test_total_score(xyz_block: str) -> float:
    """
    Calculate total OMolScore from XYZ block (For Testing)
    """
