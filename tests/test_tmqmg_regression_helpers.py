# pyright: reportMissingImports=false

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from rdkit import Chem


pytest.importorskip("openbabel")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.tmqmg_regression import (
    _copy_without_hydrogens_with_source_indices,
    _heavy_atom_count,
    _safe_canonical_smiles,
)


def _unsanitized_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    assert mol is not None
    mol.UpdatePropertyCache(strict=False)
    return mol


def _stereo_regression_mol() -> Chem.Mol:
    rw = Chem.RWMol()
    for symbol in ("C", "C", "C", "N"):
        rw.AddAtom(Chem.Atom(symbol))
    # Add the double bond before one substituent bond so SetStereoAtoms() would
    # fail if stereo information is copied before all neighboring bonds exist.
    rw.AddBond(0, 1, Chem.BondType.DOUBLE)
    rw.AddBond(0, 2, Chem.BondType.SINGLE)
    rw.AddBond(1, 3, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    bond = mol.GetBondBetweenAtoms(0, 1)
    assert bond is not None
    bond.SetStereoAtoms(2, 3)
    bond.SetStereo(Chem.BondStereo.STEREOZ)
    mol.UpdatePropertyCache(strict=False)
    return mol


@pytest.mark.parametrize(
    ("smiles", "expected_heavy_atom_count", "expected_smiles"),
    (
        ("[H-].[C-]#[O+]", 2, "[C-]#[O+]"),
        ("[2H]C", 1, "C"),
    ),
)
def test_heavy_atom_helpers_ignore_nonremovable_hydrogens(
    smiles: str,
    expected_heavy_atom_count: int,
    expected_smiles: str,
) -> None:
    mol = _unsanitized_mol(smiles)

    # RDKit RemoveHs() keeps these hydrogens, so it is not a valid heavy-atom count.
    old_remove_hs_count = Chem.RemoveHs(Chem.Mol(mol), sanitize=False).GetNumAtoms()

    heavy_only, source_indices = _copy_without_hydrogens_with_source_indices(mol)

    assert old_remove_hs_count > expected_heavy_atom_count
    assert _heavy_atom_count(mol) == expected_heavy_atom_count
    assert heavy_only.GetNumAtoms() == expected_heavy_atom_count
    assert len(source_indices) == expected_heavy_atom_count
    assert all(atom.GetAtomicNum() != 1 for atom in heavy_only.GetAtoms())
    assert _safe_canonical_smiles(mol) == expected_smiles


def test_copy_without_hydrogens_defers_bond_stereo_until_all_bonds_exist() -> None:
    mol = _stereo_regression_mol()

    heavy_only, source_indices = _copy_without_hydrogens_with_source_indices(mol)

    assert heavy_only.GetNumAtoms() == 4
    assert source_indices == [0, 1, 2, 3]
    bond = heavy_only.GetBondBetweenAtoms(0, 1)
    assert bond is not None
    assert bond.GetStereo() == Chem.BondStereo.STEREOZ
    assert list(bond.GetStereoAtoms()) == [2, 3]
