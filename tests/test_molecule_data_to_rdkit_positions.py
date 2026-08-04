"""
Author: TMJ
Date: 2026-02-27 14:45:40
LastEditors: TMJ
LastEditTime: 2026-02-27 22:34:34
Description: 请填写简介
"""
# pyright: reportMissingImports=false

import io
import sys
from pathlib import Path

import pytest
from rdkit import Chem

from molgr.fallback.utils.electrons import (
    LONE_PAIR_COUNT_PROP,
    UNRESOLVED_TWO_ELECTRON_CENTER_PROP,
    set_lone_pair_count,
    set_unpaired_electron_count,
    set_unresolved_two_electron_center,
)
from molgr.utils.converter import (
    METAL_UNPAIRED_ELECTRONS_PROP,
    get_atom_lone_pair_count,
    get_atom_unpaired_electrons,
    has_atom_unresolved_two_electron_center,
    mol_data_to_rdkit,
    pybel_to_rdmol,
)


pytest.importorskip("rdkit")
pytest.importorskip("openbabel")


def _openbabel_metal_radical_molecule():
    from openbabel import openbabel as ob
    from openbabel import pybel

    obmol = ob.OBMol()
    iron = obmol.NewAtom()
    iron.SetAtomicNum(26)
    iron.SetFormalCharge(2)
    set_unpaired_electron_count(iron, 2)
    set_lone_pair_count(iron, 1)
    iron.SetVector(0.0, 0.0, 0.0)
    carbon = obmol.NewAtom()
    carbon.SetAtomicNum(6)
    set_unpaired_electron_count(carbon, 1)
    set_unresolved_two_electron_center(carbon, True)
    carbon.SetVector(3.0, 0.0, 0.0)
    return pybel.Molecule(obmol)


def _openbabel_two_coordinate_boron_anion():
    from openbabel import openbabel as ob
    from openbabel import pybel

    obmol = ob.OBMol()
    boron = obmol.NewAtom()
    boron.SetAtomicNum(5)
    boron.SetFormalCharge(-1)
    set_unpaired_electron_count(boron, 0)
    boron.SetVector(0.0, 0.0, 0.0)
    for x in (-1.3, 1.3):
        oxygen = obmol.NewAtom()
        oxygen.SetAtomicNum(8)
        oxygen.SetVector(x, 0.0, 0.0)
        obmol.AddBond(boron.GetIdx(), oxygen.GetIdx(), 1)
    return pybel.Molecule(obmol)


def _assert_metal_unpaired_electron_representation(mol) -> None:
    iron = mol.GetAtomWithIdx(0)
    carbon = mol.GetAtomWithIdx(1)

    assert iron.GetNumRadicalElectrons() == 0
    assert iron.GetIntProp(METAL_UNPAIRED_ELECTRONS_PROP) == 2
    assert get_atom_unpaired_electrons(iron) == 2
    assert get_atom_lone_pair_count(iron) == 1
    assert carbon.GetNumRadicalElectrons() == 1
    assert not carbon.HasProp(METAL_UNPAIRED_ELECTRONS_PROP)
    assert get_atom_unpaired_electrons(carbon) == 1
    assert has_atom_unresolved_two_electron_center(carbon)
    assert mol.HasProp(f"atom.iprop.{LONE_PAIR_COUNT_PROP}")
    assert mol.HasProp(f"atom.bprop.{UNRESOLVED_TWO_ELECTRON_CENTER_PROP}")


def test_rdkit_converters_store_metal_spin_outside_formal_radicals() -> None:
    from molgr import _core

    omol = _openbabel_metal_radical_molecule()
    mol_ptr = int(getattr(omol.OBMol, "this", omol.OBMol))
    mol_data = _core.utils.extract_molecule_data(mol_ptr)

    _assert_metal_unpaired_electron_representation(mol_data_to_rdkit(mol_data, kekulize=False))
    _assert_metal_unpaired_electron_representation(pybel_to_rdmol(omol, kekulize=False))


def test_metal_unpaired_electron_property_survives_sdf_round_trip() -> None:
    mol = pybel_to_rdmol(_openbabel_metal_radical_molecule(), kekulize=False)
    buffer = io.StringIO()
    writer = Chem.SDWriter(buffer)
    writer.write(mol)
    writer.flush()

    restored = next(
        Chem.ForwardSDMolSupplier(
            io.BytesIO(buffer.getvalue().encode()),
            sanitize=False,
            removeHs=False,
        )
    )
    assert restored is not None
    iron = restored.GetAtomWithIdx(0)
    assert iron.GetNumRadicalElectrons() == 0
    assert get_atom_unpaired_electrons(iron) == 2
    assert get_atom_lone_pair_count(iron) == 1
    assert has_atom_unresolved_two_electron_center(restored.GetAtomWithIdx(1))


def test_rdkit_converters_preserve_explicit_closed_shell_boron_anion() -> None:
    from molgr import _core

    omol = _openbabel_two_coordinate_boron_anion()
    mol_ptr = int(getattr(omol.OBMol, "this", omol.OBMol))
    mol_data = _core.utils.extract_molecule_data(mol_ptr)

    for converted in (mol_data_to_rdkit(mol_data), pybel_to_rdmol(omol)):
        boron = converted.GetAtomWithIdx(0)
        assert boron.GetFormalCharge() == -1
        assert boron.GetDegree() == 2
        assert boron.GetNumRadicalElectrons() == 0
        assert get_atom_unpaired_electrons(boron) == 0


def test_mol_data_to_rdkit_sets_positions_by_atom_index() -> None:
    from openbabel import pybel

    from molgr import _core  # type: ignore
    from molgr.utils.converter import mol_data_to_rdkit

    xyz_block = """3
pos
C 0.0 0.0 0.0
O 1.5 -2.0 3.0
N -4.0 5.0 -6.0
"""

    omol = pybel.readstring("xyz", xyz_block)

    def _get_ptr(obmol) -> int:
        this = getattr(obmol, "this", None)
        if this is not None:
            return int(this)  # type: ignore[arg-type]
        return int(obmol)  # type: ignore[arg-type]

    mol_ptr = _get_ptr(omol.OBMol)
    md = _core.utils.extract_molecule_data(mol_ptr)

    mol = mol_data_to_rdkit(md, sanitize=True)
    assert mol.GetNumAtoms() == 3
    conf = mol.GetConformer()

    p0 = conf.GetAtomPosition(0)
    p1 = conf.GetAtomPosition(1)
    p2 = conf.GetAtomPosition(2)

    assert p0.x == pytest.approx(0.0)
    assert p0.y == pytest.approx(0.0)
    assert p0.z == pytest.approx(0.0)

    assert p1.x == pytest.approx(1.5)
    assert p1.y == pytest.approx(-2.0)
    assert p1.z == pytest.approx(3.0)

    assert p2.x == pytest.approx(-4.0)
    assert p2.y == pytest.approx(5.0)
    assert p2.z == pytest.approx(-6.0)


@pytest.mark.parametrize("case_idx", [1])
def test_mol_data_to_rdkit_matches_fallback_pybel_conversion_for_hard_cases(case_idx: int) -> None:
    from molgr import _core  # type: ignore
    from molgr.fallback import xyz2omol
    from molgr.utils.converter import pybel_to_rdmol
    from molgr.utils.equivalence import check_equivalence

    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    from scripts.molgr_cases_smiles_csv import load_smiles_csv_cases

    cases = load_smiles_csv_cases(Path("tests/test_cases.csv"))
    case = next(row for row in cases if int(row["case_idx"]) == case_idx)

    xyz_block = str(case["xyz_block"])
    total_charge = int(case["total_charge"])
    total_radical_electrons = int(case["total_radical_electrons"])

    fallback_omol = xyz2omol(
        xyz_block,
        total_charge=total_charge,
        total_radical_electrons=total_radical_electrons,
    )
    assert fallback_omol is not None

    cpp_mol_data = _core.pipeline.reconstruct_with_metals.xyz2omol(
        xyz_block,
        total_charge,
        total_radical_electrons,
    )
    assert cpp_mol_data is not None

    fallback_rdmol = pybel_to_rdmol(fallback_omol)
    cpp_rdmol = mol_data_to_rdkit(cpp_mol_data)

    equivalent, info = check_equivalence(
        cpp_rdmol,
        fallback_rdmol,
        use_chirality=True,
        max_resonance=100,
    )
    assert equivalent, info.reason
