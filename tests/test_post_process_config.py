# pyright: reportMissingImports=false

from __future__ import annotations

from dataclasses import replace

import pytest


pytest.importorskip("rdkit")

from rdkit import Chem, RDLogger

from molgr.config import MolGRConfig
from molgr.utils.post_process import make_dative_bond


RDLogger.DisableLog("rdApp.*")  # type: ignore


def _zinc_nitrogen_pair(distance: float) -> Chem.Mol:
    mol = Chem.RWMol()
    zinc_idx = mol.AddAtom(Chem.Atom(30))
    nitrogen_idx = mol.AddAtom(Chem.Atom(7))
    conf = Chem.Conformer(2)
    conf.SetAtomPosition(zinc_idx, (0.0, 0.0, 0.0))
    conf.SetAtomPosition(nitrogen_idx, (distance, 0.0, 0.0))
    result = mol.GetMol()
    result.AddConformer(conf)
    return result


def _blocked_zinc_nitrogen_pair() -> Chem.Mol:
    mol = Chem.RWMol()
    zinc_idx = mol.AddAtom(Chem.Atom(30))
    blocker_idx = mol.AddAtom(Chem.Atom(6))
    nitrogen_idx = mol.AddAtom(Chem.Atom(7))
    conf = Chem.Conformer(3)
    conf.SetAtomPosition(zinc_idx, (0.0, 0.0, 0.0))
    conf.SetAtomPosition(blocker_idx, (1.05, 0.0, 0.0))
    conf.SetAtomPosition(nitrogen_idx, (2.10, 0.0, 0.0))
    result = mol.GetMol()
    result.AddConformer(conf)
    return result


def _nearly_blocked_zinc_nitrogen_pair() -> Chem.Mol:
    mol = Chem.RWMol()
    zinc_idx = mol.AddAtom(Chem.Atom(30))
    blocker_idx = mol.AddAtom(Chem.Atom(6))
    nitrogen_idx = mol.AddAtom(Chem.Atom(7))
    conf = Chem.Conformer(3)
    conf.SetAtomPosition(zinc_idx, (0.0, 0.0, 0.0))
    conf.SetAtomPosition(blocker_idx, (1.00, 1.00, 0.0))
    conf.SetAtomPosition(nitrogen_idx, (2.00, 0.0, 0.0))
    result = mol.GetMol()
    result.AddConformer(conf)
    return result


def _iron_alkene_pair(distance_1: float, distance_2: float) -> Chem.Mol:
    mol = Chem.RWMol()
    iron_idx = mol.AddAtom(Chem.Atom(26))
    carbon_1_idx = mol.AddAtom(Chem.Atom(6))
    carbon_2_idx = mol.AddAtom(Chem.Atom(6))
    mol.AddBond(carbon_1_idx, carbon_2_idx, Chem.BondType.DOUBLE)
    conf = Chem.Conformer(3)
    conf.SetAtomPosition(iron_idx, (0.0, 0.0, 0.0))
    conf.SetAtomPosition(carbon_1_idx, (distance_1, 0.0, 0.0))
    conf.SetAtomPosition(carbon_2_idx, (0.0, distance_2, 0.0))
    result = mol.GetMol()
    result.AddConformer(conf)
    return result


def _blocked_iron_alkene_pair() -> Chem.Mol:
    mol = Chem.RWMol()
    iron_idx = mol.AddAtom(Chem.Atom(26))
    carbon_1_idx = mol.AddAtom(Chem.Atom(6))
    carbon_2_idx = mol.AddAtom(Chem.Atom(6))
    blocker_idx = mol.AddAtom(Chem.Atom(6))
    mol.AddBond(carbon_1_idx, carbon_2_idx, Chem.BondType.DOUBLE)
    conf = Chem.Conformer(4)
    conf.SetAtomPosition(iron_idx, (0.0, 0.0, 0.0))
    conf.SetAtomPosition(carbon_1_idx, (1.80, 0.0, 0.0))
    conf.SetAtomPosition(carbon_2_idx, (0.0, 1.80, 0.0))
    conf.SetAtomPosition(blocker_idx, (0.90, 0.0, 0.0))
    result = mol.GetMol()
    result.AddConformer(conf)
    return result


def test_make_dative_bond_uses_metal_coordination_tolerance_config() -> None:
    base_config = MolGRConfig()
    tight_config = replace(
        base_config,
        metal_scoring=replace(
            base_config.metal_scoring,
            metal_coordination_extra_tolerance_angstrom=0.10,
        ),
    )
    loose_config = replace(
        base_config,
        metal_scoring=replace(
            base_config.metal_scoring,
            metal_coordination_extra_tolerance_angstrom=0.35,
        ),
    )

    tight_mol = make_dative_bond(_zinc_nitrogen_pair(2.10), config=tight_config)
    loose_mol = make_dative_bond(_zinc_nitrogen_pair(2.10), config=loose_config)

    assert tight_mol.GetBondBetweenAtoms(0, 1) is None
    loose_bond = loose_mol.GetBondBetweenAtoms(0, 1)
    assert loose_bond is not None
    assert loose_bond.GetBondType() == Chem.BondType.DATIVE


def test_make_dative_bond_requires_visible_coordination_atom() -> None:
    mol = make_dative_bond(_blocked_zinc_nitrogen_pair())

    assert mol.GetBondBetweenAtoms(0, 2) is None


def test_metal_access_radius_scale_expands_post_process_distance_cutoff() -> None:
    base_config = MolGRConfig()
    tight_config = replace(
        base_config,
        metal_scoring=replace(
            base_config.metal_scoring,
            metal_access_radius_scale=1.0,
            metal_coordination_extra_tolerance_angstrom=0.10,
        ),
    )
    loose_config = replace(
        base_config,
        metal_scoring=replace(
            base_config.metal_scoring,
            metal_access_radius_scale=1.50,
            metal_coordination_extra_tolerance_angstrom=0.10,
        ),
    )

    tight_mol = make_dative_bond(_zinc_nitrogen_pair(2.22), config=tight_config)
    loose_mol = make_dative_bond(_zinc_nitrogen_pair(2.22), config=loose_config)

    assert tight_mol.GetBondBetweenAtoms(0, 1) is None
    loose_bond = loose_mol.GetBondBetweenAtoms(0, 1)
    assert loose_bond is not None
    assert loose_bond.GetBondType() == Chem.BondType.DATIVE


def test_metal_access_radius_scale_expands_blocker_radius() -> None:
    base_config = MolGRConfig()
    unscaled_config = replace(
        base_config,
        metal_scoring=replace(
            base_config.metal_scoring,
            metal_access_radius_scale=1.0,
            metal_coordination_extra_tolerance_angstrom=0.35,
            metal_access_clearance_angstrom=0.0,
        ),
    )
    scaled_config = replace(
        base_config,
        metal_scoring=replace(
            base_config.metal_scoring,
            metal_access_radius_scale=1.50,
            metal_coordination_extra_tolerance_angstrom=0.35,
            metal_access_clearance_angstrom=0.0,
        ),
    )

    unscaled_mol = make_dative_bond(_nearly_blocked_zinc_nitrogen_pair(), config=unscaled_config)
    scaled_mol = make_dative_bond(_nearly_blocked_zinc_nitrogen_pair(), config=scaled_config)

    assert unscaled_mol.GetBondBetweenAtoms(0, 2).GetBondType() == Chem.BondType.DATIVE
    assert scaled_mol.GetBondBetweenAtoms(0, 2) is None


def test_make_dative_bond_uses_pi_dative_distance_difference_tolerance_config() -> None:
    base_config = MolGRConfig()
    tight_config = replace(
        base_config,
        metal_scoring=replace(
            base_config.metal_scoring,
            pi_dative_distance_difference_tolerance_angstrom=0.05,
        ),
    )
    loose_config = replace(
        base_config,
        metal_scoring=replace(
            base_config.metal_scoring,
            pi_dative_distance_difference_tolerance_angstrom=0.15,
        ),
    )

    tight_mol = make_dative_bond(_iron_alkene_pair(1.80, 1.90), config=tight_config)
    loose_mol = make_dative_bond(_iron_alkene_pair(1.80, 1.90), config=loose_config)

    assert tight_mol.GetBondBetweenAtoms(0, 1) is None
    assert tight_mol.GetBondBetweenAtoms(0, 2) is None
    assert loose_mol.GetBondBetweenAtoms(0, 1).GetBondType() == Chem.BondType.DATIVE
    assert loose_mol.GetBondBetweenAtoms(0, 2).GetBondType() == Chem.BondType.DATIVE


def test_make_pi_dative_bond_requires_both_atoms_visible_to_metal() -> None:
    mol = make_dative_bond(_blocked_iron_alkene_pair())

    assert mol.GetBondBetweenAtoms(0, 1) is None
    assert mol.GetBondBetweenAtoms(0, 2) is None
