# pyright: reportMissingImports=false

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import pytest


pytest.importorskip("openbabel")
pytest.importorskip("rdkit")

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom

from molgr.interface import xyz_to_rdmol


RDLogger.DisableLog("rdApp.*")  # type: ignore

_EMBED_SEED = 0xC0FFEE


def _total_charge_and_radicals(mol: Chem.Mol) -> tuple[int, int]:
    charge = 0
    radicals = 0
    for atom in mol.GetAtoms():  # pyright: ignore[reportCallIssue]
        charge += int(atom.GetFormalCharge())
        radicals += int(atom.GetNumRadicalElectrons())
    return charge, radicals


def _canonical_smiles(mol: Chem.Mol) -> str:
    heavy = Chem.RemoveHs(mol)
    return Chem.MolToSmiles(
        heavy,
        canonical=True,
        isomericSmiles=True,
        allBondsExplicit=True,
        allHsExplicit=True,
    )


def _ordered_atom_signature(mol: Chem.Mol) -> tuple[tuple[Any, ...], ...]:
    return tuple(
        (
            atom.GetIdx(),
            atom.GetAtomicNum(),
            atom.GetFormalCharge(),
            atom.GetNumRadicalElectrons(),
            str(atom.GetChiralTag()),
            atom.GetIsotope(),
            atom.GetNoImplicit(),
            atom.GetIsAromatic(),
        )
        for atom in mol.GetAtoms()
    )


def _ordered_bond_signature(mol: Chem.Mol) -> tuple[tuple[Any, ...], ...]:
    return tuple(
        sorted(
            (
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                str(bond.GetBondType()),
                bond.GetIsAromatic(),
                str(bond.GetStereo()),
            )
            for bond in mol.GetBonds()
        )
    )


def _reconstruction_signature(mol: Chem.Mol) -> tuple[Any, ...]:
    return (
        _canonical_smiles(mol),
        _ordered_atom_signature(mol),
        _ordered_bond_signature(mol),
    )


def _assert_coordinates_match(cpp_mol: Chem.Mol, python_mol: Chem.Mol) -> None:
    assert cpp_mol.GetNumAtoms() == python_mol.GetNumAtoms()
    cpp_conf = cpp_mol.GetConformer()
    python_conf = python_mol.GetConformer()
    for atom_idx in range(cpp_mol.GetNumAtoms()):
        cpp_pos = cpp_conf.GetAtomPosition(atom_idx)
        python_pos = python_conf.GetAtomPosition(atom_idx)
        assert cpp_pos.x == pytest.approx(python_pos.x, abs=1e-4)
        assert cpp_pos.y == pytest.approx(python_pos.y, abs=1e-4)
        assert cpp_pos.z == pytest.approx(python_pos.z, abs=1e-4)


def _assert_backend_results_match(
    xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
    *,
    make_dative_bonds: bool,
) -> None:
    spin_multiplicity = total_radical_electrons + 1
    cpp_mol = xyz_to_rdmol(
        xyz_block,
        total_charge,
        spin_multiplicity,
        backend="cpp",
        make_dative_bonds=make_dative_bonds,
    )
    python_mol = xyz_to_rdmol(
        xyz_block,
        total_charge,
        spin_multiplicity,
        backend="python",
        make_dative_bonds=make_dative_bonds,
    )

    assert _reconstruction_signature(cpp_mol) == _reconstruction_signature(python_mol)
    _assert_coordinates_match(cpp_mol, python_mol)


def _load_smiles_backend_cases() -> list[object]:
    csv_path = Path(__file__).with_name("test_cases.csv")
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    cases: list[object] = []
    for case_idx, row in enumerate(rows, start=1):
        smiles = row["smiles"].strip()
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        mol_h = Chem.AddHs(mol)
        embed_code = rdDistGeom.EmbedMolecule(  # pyright: ignore[reportCallIssue]
            mol_h,
            randomSeed=_EMBED_SEED,
        )
        assert int(embed_code) == 0
        total_charge, total_radical_electrons = _total_charge_and_radicals(mol_h)
        cases.append(
            pytest.param(
                Chem.MolToXYZBlock(mol_h),
                total_charge,
                total_radical_electrons,
                id=f"smiles-case-{case_idx:02d}",
            )
        )
    return cases


@pytest.mark.parametrize(
    ("xyz_block", "total_charge", "total_radical_electrons"),
    _load_smiles_backend_cases(),
)
def test_cpp_and_python_backends_match_smiles_regression_cases(
    xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
) -> None:
    _assert_backend_results_match(
        xyz_block,
        total_charge,
        total_radical_electrons,
        make_dative_bonds=False,
    )
    _assert_backend_results_match(
        xyz_block,
        total_charge,
        total_radical_electrons,
        make_dative_bonds=True,
    )


def test_cpp_and_python_backends_match_monnmo_regression_case() -> None:
    mol = Chem.MolFromMolFile(
        str(Path("tests/data/sdf/MoNNMo.sdf")),
        sanitize=False,
        removeHs=False,
        strictParsing=False,
    )
    assert mol is not None
    total_charge, total_radical_electrons = _total_charge_and_radicals(mol)
    xyz_block = Chem.MolToXYZBlock(mol)

    _assert_backend_results_match(
        xyz_block,
        total_charge,
        total_radical_electrons,
        make_dative_bonds=False,
    )
    _assert_backend_results_match(
        xyz_block,
        total_charge,
        total_radical_electrons,
        make_dative_bonds=True,
    )
