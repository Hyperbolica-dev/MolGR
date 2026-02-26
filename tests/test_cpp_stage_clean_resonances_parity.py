# pyright: reportMissingImports=false

from typing import Any, Iterable, List, Tuple

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr import _core  # type: ignore
from molgr.fallback.stages import clean


_stages: Any = _core.stages


def _get_ptr(obmol) -> int:
    this = getattr(obmol, "this", None)
    if this is not None:
        return int(this)  # type: ignore[arg-type]
    return int(obmol)  # type: ignore[arg-type]


def _clone_mol(mol: pybel.Molecule) -> pybel.Molecule:
    return pybel.Molecule(ob.OBMol(mol.OBMol))


def _smiles_token(mol: pybel.Molecule) -> str:
    smi = mol.write("smi")
    assert smi is not None
    return smi.split()[0]


def _state_signature(
    mol: pybel.Molecule,
) -> Tuple[Tuple[Tuple[int, int, int], ...], Tuple[Tuple[int, int, int], ...]]:
    atoms: List[Tuple[int, int, int]] = []
    for idx in range(1, mol.OBMol.NumAtoms() + 1):
        atom = mol.OBMol.GetAtom(idx)
        atoms.append((atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetSpinMultiplicity()))

    bonds: List[Tuple[int, int, int]] = []
    for bond in ob.OBMolBondIter(mol.OBMol):
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx > end_idx:
            begin_idx, end_idx = end_idx, begin_idx
        bonds.append((begin_idx, end_idx, bond.GetBondOrder()))
    bonds.sort()
    return tuple(atoms), tuple(bonds)


def _assert_atom_state_parity(py_mol: pybel.Molecule, cpp_mol: pybel.Molecule) -> None:
    assert py_mol.OBMol.NumAtoms() == cpp_mol.OBMol.NumAtoms()
    for idx in range(1, py_mol.OBMol.NumAtoms() + 1):
        py_atom = py_mol.OBMol.GetAtom(idx)
        cpp_atom = cpp_mol.OBMol.GetAtom(idx)
        assert py_atom.GetFormalCharge() == cpp_atom.GetFormalCharge()
        assert py_atom.GetSpinMultiplicity() == cpp_atom.GetSpinMultiplicity()


def _assert_key_bond_orders(
    py_mol: pybel.Molecule, cpp_mol: pybel.Molecule, bonds: Iterable[Tuple[int, int]]
) -> None:
    for begin_idx, end_idx in bonds:
        py_bond = py_mol.OBMol.GetBond(begin_idx, end_idx)
        cpp_bond = cpp_mol.OBMol.GetBond(begin_idx, end_idx)
        assert py_bond is not None and cpp_bond is not None
        assert py_bond.GetBondOrder() == cpp_bond.GetBondOrder()


def _build_multi_rule_seed_a() -> ob.OBMol:
    obmol = ob.OBMol()
    obmol.BeginModify()

    atom1 = obmol.NewAtom()
    atom1.SetAtomicNum(6)
    atom1.SetFormalCharge(-1)

    atom2 = obmol.NewAtom()
    atom2.SetAtomicNum(7)
    atom2.SetFormalCharge(1)

    atom3 = obmol.NewAtom()
    atom3.SetAtomicNum(6)
    atom3.SetFormalCharge(0)

    atom4 = obmol.NewAtom()
    atom4.SetAtomicNum(6)
    atom4.SetFormalCharge(0)

    atom5 = obmol.NewAtom()
    atom5.SetAtomicNum(6)
    atom5.SetFormalCharge(0)

    atom6 = obmol.NewAtom()
    atom6.SetAtomicNum(6)
    atom6.SetFormalCharge(-1)

    obmol.AddBond(1, 2, 2)
    obmol.AddBond(2, 3, 2)
    obmol.AddBond(4, 5, 2)
    obmol.AddBond(5, 6, 2)

    obmol.EndModify()
    return obmol


def _build_multi_rule_seed_b() -> ob.OBMol:
    obmol = ob.OBMol()
    obmol.BeginModify()

    atom1 = obmol.NewAtom()
    atom1.SetAtomicNum(6)
    atom1.SetFormalCharge(0)
    atom1.SetSpinMultiplicity(1)

    atom2 = obmol.NewAtom()
    atom2.SetAtomicNum(6)
    atom2.SetFormalCharge(0)

    atom3 = obmol.NewAtom()
    atom3.SetAtomicNum(6)
    atom3.SetFormalCharge(0)

    atom4 = obmol.NewAtom()
    atom4.SetAtomicNum(6)
    atom4.SetFormalCharge(0)
    atom4.SetSpinMultiplicity(1)

    atom5 = obmol.NewAtom()
    atom5.SetAtomicNum(7)
    atom5.SetFormalCharge(1)

    atom6 = obmol.NewAtom()
    atom6.SetAtomicNum(8)
    atom6.SetFormalCharge(-1)

    obmol.AddBond(1, 2, 1)
    obmol.AddBond(2, 3, 2)
    obmol.AddBond(3, 4, 1)
    obmol.AddBond(5, 6, 1)

    obmol.EndModify()
    return obmol


def _applied_rule_ids_in_python_order(mol: pybel.Molecule) -> List[int]:
    ordered_rules = [
        (0, clean.clean_resonances_0),
        (1, clean.clean_resonances_1),
        (2, clean.clean_resonances_2),
        (3, clean.clean_resonances_3),
        (4, clean.clean_resonances_4),
        (9, clean.clean_resonances_9),
        (5, clean.clean_resonances_5),
        (6, clean.clean_resonances_6),
        (7, clean.clean_resonances_7),
        (8, clean.clean_resonances_8),
        (9, clean.clean_resonances_9),
        (10, clean.clean_resonances_10),
        (11, clean.clean_resonances_11),
        (12, clean.clean_resonances_12),
        (13, clean.clean_resonances_13),
    ]

    probe = _clone_mol(mol)
    applied_ids: List[int] = []
    for rule_id, fn in ordered_rules:
        before = _state_signature(probe)
        probe = fn(probe)
        after = _state_signature(probe)
        if before != after:
            applied_ids.append(rule_id)
    return applied_ids


@pytest.mark.parametrize(
    ("seed_builder", "key_bonds"),
    [
        (_build_multi_rule_seed_a, [(1, 2), (2, 3), (4, 5), (5, 6)]),
        (_build_multi_rule_seed_b, [(1, 2), (2, 3), (3, 4), (5, 6)]),
    ],
)
def test_clean_resonances_cpp_matches_python(seed_builder, key_bonds) -> None:
    py_mol = pybel.Molecule(seed_builder())
    cpp_mol = pybel.Molecule(seed_builder())

    applied_rule_ids = _applied_rule_ids_in_python_order(py_mol)
    assert len(set(applied_rule_ids)) >= 2

    py_mol = clean.clean_resonances(py_mol)
    _stages.clean.clean_resonances_ptr(_get_ptr(cpp_mol.OBMol))

    _assert_atom_state_parity(py_mol, cpp_mol)
    _assert_key_bond_orders(py_mol, cpp_mol, key_bonds)
    assert _smiles_token(py_mol) == _smiles_token(cpp_mol)


def test_clean_resonances_rule12_room_guard_keeps_problematic_token() -> None:
    smiles = "C1=C2[C-]([CH]C(=O)O2)C[NH2+]1"
    expected_token = _smiles_token(pybel.readstring("smi", smiles))

    py_mol = pybel.readstring("smi", smiles)
    cpp_mol = pybel.readstring("smi", smiles)

    py_mol = clean.clean_resonances(py_mol)
    _stages.clean.clean_resonances_ptr(_get_ptr(cpp_mol.OBMol))

    assert _smiles_token(py_mol) == expected_token
    assert _smiles_token(cpp_mol) == expected_token
