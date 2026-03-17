"""
Fallback-aligned eliminate stage helpers
"""

from __future__ import annotations

import typing

__all__: list[str] = [
    "eliminate_1_3_dipole_ptr",
    "eliminate_carbene_neighbor_heteroatom_ptr",
    "eliminate_carboxyl_ptr",
    "eliminate_charge_spliting_ptr",
    "eliminate_cn_in_doubt_ptr",
    "eliminate_high_positive_charge_atoms_ptr",
    "eliminate_negative_charges_ptr",
    "eliminate_nnn_ptr",
    "eliminate_positive_charges_ptr",
]

def eliminate_1_3_dipole_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    given_charge: typing.SupportsInt | typing.SupportsIndex,
) -> int:
    """
    Apply eliminate.eliminate_1_3_dipole to an existing OBMol.

    Args:
        mol_ptr: int address of OpenBabel::OBMol
        given_charge: charge deficit to be updated in place and returned
    """

def eliminate_carbene_neighbor_heteroatom_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    given_charge: typing.SupportsInt | typing.SupportsIndex,
) -> int:
    """
    Apply eliminate.eliminate_carbene_neighbor_heteroatom to an existing OBMol.

    Args:
        mol_ptr: int address of OpenBabel::OBMol
        given_charge: charge deficit to be updated in place and returned
    """

def eliminate_carboxyl_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    given_charge: typing.SupportsInt | typing.SupportsIndex,
) -> int:
    """
    Apply eliminate.eliminate_carboxyl to an existing OBMol.

    Args:
        mol_ptr: int address of OpenBabel::OBMol
        given_charge: charge deficit to be updated in place and returned
    """

def eliminate_charge_spliting_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    given_charge: typing.SupportsInt | typing.SupportsIndex,
) -> int:
    """
    Apply eliminate.eliminate_charge_spliting to an existing OBMol.

    Args:
        mol_ptr: int address of OpenBabel::OBMol
        given_charge: charge deficit to be updated in place and returned
    """

def eliminate_cn_in_doubt_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    given_charge: typing.SupportsInt | typing.SupportsIndex,
) -> int:
    """
    Apply eliminate.eliminate_cn_in_doubt to an existing OBMol.

    Args:
        mol_ptr: int address of OpenBabel::OBMol
        given_charge: charge deficit to be updated in place and returned
    """

def eliminate_high_positive_charge_atoms_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    given_charge: typing.SupportsInt | typing.SupportsIndex,
) -> int:
    """
    Apply eliminate.eliminate_high_positive_charge_atoms to an existing OBMol.

    Args:
        mol_ptr: int address of OpenBabel::OBMol
        given_charge: charge deficit to be updated in place and returned
    """

def eliminate_negative_charges_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    given_charge: typing.SupportsInt | typing.SupportsIndex,
) -> int:
    """
    Apply eliminate.eliminate_negative_charges to an existing OBMol.

    Args:
        mol_ptr: int address of OpenBabel::OBMol
        given_charge: charge deficit to be updated in place and returned
    """

def eliminate_nnn_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    given_charge: typing.SupportsInt | typing.SupportsIndex,
    positive: bool = False,
) -> int:
    """
    Apply eliminate.eliminate_nnn to an existing OBMol.

    Args:
        mol_ptr: int address of OpenBabel::OBMol
        given_charge: charge deficit to be updated in place and returned
        positive: whether to run positive-direction elimination
    """

def eliminate_positive_charges_ptr(
    mol_ptr: typing.SupportsInt | typing.SupportsIndex,
    given_charge: typing.SupportsInt | typing.SupportsIndex,
) -> int:
    """
    Apply eliminate.eliminate_positive_charges to an existing OBMol.

    Args:
        mol_ptr: int address of OpenBabel::OBMol
        given_charge: charge deficit to be updated in place and returned
    """
