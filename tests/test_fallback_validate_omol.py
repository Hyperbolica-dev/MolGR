"""
Author: TMJ
Date: 2026-02-25 00:41:17
LastEditors: TMJ
LastEditTime: 2026-02-25 13:45:25
Description: 请填写简介
"""
# pyright: reportMissingImports=false

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.stages.preprocess import validate_omol


def test_validate_omol_singlet_carbene_spin_mod_2_compatibility() -> None:
    obmol = ob.OBMol()
    a = obmol.NewAtom()
    a.SetAtomicNum(6)
    a.SetFormalCharge(0)

    a.SetSpinMultiplicity(2)
    mol = pybel.Molecule(obmol)
    assert validate_omol(mol, 0, 0) is True

    a.SetSpinMultiplicity(1)
    mol = pybel.Molecule(obmol)
    assert validate_omol(mol, 0, 1) is True
