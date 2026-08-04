"""
Fallback-aligned no-metal reconstruction helpers
"""

from __future__ import annotations

import typing

import molgr._core.utils
import molgr.config

__all__: list[str] = ["xyz_to_omol_no_metal"]

def xyz_to_omol_no_metal(
    xyz_block: str,
    total_charge: typing.SupportsInt | typing.SupportsIndex = 0,
    total_radical_electrons: typing.SupportsInt | typing.SupportsIndex = 0,
    *,
    config: molgr.config.MolGRConfig | None = None,
) -> molgr._core.utils.MoleculeData | None:
    """
    Reconstruct molecule data from XYZ without metal handling.
    """
