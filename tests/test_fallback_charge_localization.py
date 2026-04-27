# pyright: reportMissingImports=false
from __future__ import annotations

import pytest


pytest.importorskip("openbabel")

from openbabel import pybel

from molgr.fallback.utils.metals.scoring import _charge_localization_penalty_for_atom


def _charged_atom(smiles: str):
    mol = pybel.readstring("smi", smiles)
    charged = [atom.OBAtom for atom in mol.atoms if atom.OBAtom.GetFormalCharge() != 0]
    assert len(charged) == 1
    return mol, charged[0]


def test_two_coordinate_group15_anion_is_penalized_less_than_cation() -> None:
    cation_mol, cation = _charged_atom("[N+](C)C")
    anion_mol, anion = _charged_atom("[N-](C)C")

    cation_penalty = _charge_localization_penalty_for_atom(cation, is_conjugated=False)
    anion_penalty = _charge_localization_penalty_for_atom(anion, is_conjugated=False)

    assert cation_mol is not None and anion_mol is not None
    assert anion_penalty < cation_penalty


def test_octet_satisfied_ammonium_is_penalized_less_than_electron_deficient_nitrenium() -> None:
    ammonium_mol, ammonium = _charged_atom("[NH4+]")
    nitrenium_mol, nitrenium = _charged_atom("[N+](C)C")

    ammonium_penalty = _charge_localization_penalty_for_atom(ammonium, is_conjugated=False)
    nitrenium_penalty = _charge_localization_penalty_for_atom(nitrenium, is_conjugated=False)

    assert ammonium_mol is not None and nitrenium_mol is not None
    assert ammonium_penalty < nitrenium_penalty


def test_electron_rich_ammonium_and_borate_share_a_low_penalty_bucket() -> None:
    ammonium_mol, ammonium = _charged_atom("[NH4+]")
    borate_mol, borate = _charged_atom("[B-](F)(F)(F)F")
    nitrenium_mol, nitrenium = _charged_atom("[N+](C)C")

    ammonium_penalty = _charge_localization_penalty_for_atom(ammonium, is_conjugated=False)
    borate_penalty = _charge_localization_penalty_for_atom(borate, is_conjugated=False)
    nitrenium_penalty = _charge_localization_penalty_for_atom(nitrenium, is_conjugated=False)

    assert ammonium_mol is not None and borate_mol is not None and nitrenium_mol is not None
    assert abs(ammonium_penalty - borate_penalty) < 0.2
    assert ammonium_penalty < nitrenium_penalty
    assert borate_penalty < nitrenium_penalty


def test_boryl_and_borate_anions_are_not_treated_like_generic_carbanions() -> None:
    boryl_mol, boryl = _charged_atom("[B-](C)C")
    borate_mol, borate = _charged_atom("[B-](C)(C)(C)C")
    carbanion_mol, carbanion = _charged_atom("[C-](C)C")

    boryl_penalty = _charge_localization_penalty_for_atom(boryl, is_conjugated=False)
    borate_penalty = _charge_localization_penalty_for_atom(borate, is_conjugated=False)
    carbanion_penalty = _charge_localization_penalty_for_atom(carbanion, is_conjugated=False)

    assert boryl_mol is not None and borate_mol is not None and carbanion_mol is not None
    assert boryl_penalty < carbanion_penalty
    assert borate_penalty < carbanion_penalty
