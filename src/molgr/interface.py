"""
Author: TMJ
Date: 2026-02-21 23:08:39
LastEditors: TMJ
LastEditTime: 2026-04-28 23:25:31
Description: 请填写简介
"""

from __future__ import annotations

from typing import Literal

from openbabel import pybel
from rdkit import Chem

from molgr.config import CONFIG, MolGRConfig
from molgr.fallback import xyz2omol
from molgr.utils.converter import mol_data_to_rdkit, pybel_to_rdmol
from molgr.utils.post_process import make_dative_bond
from molgr.utils.post_process import make_stereochemistry as restore_stereochemistry

from . import _core as core


def _suspicious_rdmol_from_input_xyz(xyz_block: str) -> Chem.Mol:
    """Build an untrusted RDKit molecule from OpenBabel's initial bond perception."""

    omol = pybel.readstring("xyz", xyz_block)
    omol.OBMol.ConnectTheDots()
    omol.OBMol.PerceiveBondOrders()
    rdmol = pybel_to_rdmol(omol, sanitize=False, kekulize=False)
    rdmol.SetProp("_MolGRReconstructionStatus", "suspicious_fallback")
    return rdmol


def _should_return_suspicious_on_reconstruction_failure(config: MolGRConfig) -> bool:
    policy = config.interface.reconstruction_failure_policy
    if policy == "return_suspicious":
        return True
    if policy == "raise":
        return False
    raise ValueError(f"Unknown reconstruction failure policy: {policy!r}")


def _handle_reconstruction_failure(xyz_block: str, *, config: MolGRConfig) -> Chem.Mol:
    if _should_return_suspicious_on_reconstruction_failure(config):
        return _suspicious_rdmol_from_input_xyz(xyz_block)
    raise ValueError("xyz2omol failed")


def xyz_to_rdmol(
    xyz_block: str,
    total_charge: int = 0,
    spin_multiplicity: int = 1,
    *,
    backend: Literal["cpp", "python"] = "cpp",
    make_dative_bonds: bool = True,
    make_stereochemistry: bool = True,
    config: MolGRConfig | None = None,
) -> Chem.Mol:
    """
    Convert XYZ block to RDKit Mol.
    """
    resolved_config = CONFIG if config is None else config
    total_radical_electrons = spin_multiplicity - 1
    if backend == "cpp":
        try:
            moldata = core.pipeline.reconstruct_with_metals.xyz2omol(
                xyz_block,
                total_charge,
                total_radical_electrons,
                config=resolved_config,
            )
        except Exception:
            if _should_return_suspicious_on_reconstruction_failure(resolved_config):
                return _suspicious_rdmol_from_input_xyz(xyz_block)
            raise
        if moldata is None:
            return _handle_reconstruction_failure(xyz_block, config=resolved_config)
        rdmol = mol_data_to_rdkit(moldata)
    elif backend == "python":
        try:
            omol = xyz2omol(
                xyz_block,
                total_charge,
                total_radical_electrons,
                config=resolved_config,
            )
        except Exception:
            if _should_return_suspicious_on_reconstruction_failure(resolved_config):
                return _suspicious_rdmol_from_input_xyz(xyz_block)
            raise
        if omol is None:
            return _handle_reconstruction_failure(xyz_block, config=resolved_config)
        rdmol = pybel_to_rdmol(omol)
    else:
        raise ValueError(f"Unknown backend: {backend}")

    if make_dative_bonds:
        rdmol = make_dative_bond(rdmol, config=resolved_config)

    for atom_idx in range(rdmol.GetNumAtoms()):
        rd_atom = rdmol.GetAtomWithIdx(atom_idx)
        rd_atom.SetNoImplicit(True)
    if make_stereochemistry:
        rdmol = restore_stereochemistry(rdmol)
    return rdmol
