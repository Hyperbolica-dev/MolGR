"""
Author: TMJ
Date: 2026-02-21 23:08:39
LastEditors: TMJ
LastEditTime: 2026-02-28 00:01:56
Description: 请填写简介
"""

from __future__ import annotations

from typing import Literal

from rdkit import Chem

from molgr.config import MolGRConfig, resolve_config
from molgr.fallback import xyz2omol
from molgr.utils.converter import mol_data_to_rdkit, pybel_to_rdmol
from molgr.utils.post_process import make_dative_bond

from . import _core as core


def xyz_to_rdmol(
    xyz_block: str,
    total_charge: int = 0,
    spin_multiplicity: int = 1,
    *,
    backend: Literal["cpp", "python"] = "cpp",
    make_dative_bonds: bool = True,
    config: MolGRConfig | None = None,
) -> Chem.Mol:
    """
    Convert XYZ block to RDKit Mol.
    """
    total_radical_electrons = spin_multiplicity - 1
    if backend == "cpp":
        resolved_config = resolve_config(config)
        moldata = core.pipeline.reconstruct_with_metals.xyz2omol(
            xyz_block,
            total_charge,
            total_radical_electrons,
            config=resolved_config,
        )
        if moldata is None:
            raise ValueError("xyz2omol failed")
        rdmol = mol_data_to_rdkit(moldata)
    elif backend == "python":
        omol = xyz2omol(
            xyz_block,
            total_charge,
            total_radical_electrons,
            config=config,
        )
        if omol is None:
            raise ValueError("xyz2omol failed")
        rdmol = pybel_to_rdmol(omol)
    else:
        raise ValueError(f"Unknown backend: {backend}")
    if make_dative_bonds:
        return make_dative_bond(rdmol)
    return rdmol
