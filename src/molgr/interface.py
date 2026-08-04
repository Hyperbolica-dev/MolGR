"""
Author: TMJ
Date: 2026-02-21 23:08:39
LastEditors: TMJ
LastEditTime: 2026-04-28 23:25:31
Description: 请填写简介
"""

from __future__ import annotations  # noqa: I001

from typing import Literal

# Keep RDKit ahead of Open Babel: importing pybel first segfaults on cp313 manylinux.
# isort: off
from rdkit import Chem
from openbabel import pybel
# isort: on

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


def _validate_spin_multiplicity(
    xyz_block: str,
    total_charge: int,
    spin_multiplicity: int,
) -> None:
    """Reject spin states that are impossible for the input electron count."""

    if spin_multiplicity < 1:
        raise ValueError("spin_multiplicity must be >= 1")

    omol = pybel.readstring("xyz", xyz_block)
    total_electrons = sum(int(atom.atomicnum) for atom in omol.atoms) - total_charge
    if total_electrons < 0:
        raise ValueError(
            f"total_charge={total_charge} leaves a negative total electron count "
            f"({total_electrons})"
        )
    if spin_multiplicity > total_electrons + 1:
        raise ValueError(
            f"spin_multiplicity={spin_multiplicity} is impossible for "
            f"{total_electrons} total electrons; the maximum is {total_electrons + 1}"
        )
    if total_electrons % 2 != (spin_multiplicity - 1) % 2:
        required_parity = "odd" if total_electrons % 2 == 0 else "even"
        raise ValueError(
            f"spin_multiplicity={spin_multiplicity} is impossible for "
            f"{total_electrons} total electrons; the multiplicity must be {required_parity}"
        )


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
    _validate_spin_multiplicity(xyz_block, total_charge, spin_multiplicity)
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
