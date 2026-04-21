# pyright: reportMissingImports=false

from __future__ import annotations

from typing import Any

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr import _core
from molgr.fallback.stages.break_bond import break_one_bond


def _get_ptr(obmol: ob.OBMol) -> int:
    this = getattr(obmol, "this", None)
    if this is not None:
        return int(this)  # type: ignore[arg-type]
    return int(obmol)  # type: ignore[arg-type]


def _signature(mol: pybel.Molecule) -> dict[str, Any]:
    return {
        "smiles": mol.write("smi").split()[0],
        "total_radicals": sum(atom.OBAtom.GetSpinMultiplicity() for atom in mol.atoms),
        "bonds": sorted(
            (
                min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                bond.GetBondOrder(),
            )
            for bond in ob.OBMolBondIter(mol.OBMol)
        ),
    }


def _run_cpp(smiles: str, given_charge: int, given_radical: int) -> tuple[pybel.Molecule, int]:
    mol = pybel.readstring("smi", smiles)
    charge, hit = _core.dev.stages.break_bond.break_one_bond_ptr(
        _get_ptr(mol.OBMol),
        given_charge,
        given_radical,
    )
    assert hit is True
    return mol, int(charge)


@pytest.mark.parametrize(
    ("smiles", "given_radical", "expected_signature"),
    [
        (
            "c1ccccc1",
            6,
            {
                "smiles": "[cH]1cccc[cH]1",
                "total_radicals": 2,
                "bonds": [(1, 2, 1), (1, 6, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 6, 1)],
            },
        ),
        (
            "CCC",
            5,
            {
                "smiles": "[CH3].[CH2]C",
                "total_radicals": 2,
                "bonds": [(2, 3, 1)],
            },
        ),
    ],
)
def test_break_one_bond_limits_aromatic_and_single_cases_to_one_application(
    smiles: str,
    given_radical: int,
    expected_signature: dict[str, Any],
) -> None:
    fallback_mol = pybel.readstring("smi", smiles)
    fallback_after, fallback_charge, hit = break_one_bond(
        fallback_mol,
        given_charge=0,
        given_radical=given_radical,
    )

    cpp_after, cpp_charge = _run_cpp(smiles, 0, given_radical)

    assert fallback_charge == 0
    assert cpp_charge == 0
    assert hit is True

    assert _signature(fallback_after) == expected_signature
    assert _signature(cpp_after) == expected_signature
