# pyright: reportMissingImports=false
from __future__ import annotations

import pytest


pytest.importorskip("openbabel")

from openbabel import pybel

from molgr.fallback.stages.fresh import assign_radical_dots
from molgr.fallback.utils.electrons import (
    set_lone_pair_count,
    set_unpaired_electron_count,
)
from molgr.fallback.utils.metals.scoring import (
    _charge_localization_penalty_for_atom,
    _compute_organic_electronic_state_metrics,
    _radical_localization_penalty_for_atom,
)


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


def test_explicit_singlet_carbene_like_center_equals_the_triplet_penalty() -> None:
    mol = pybel.readstring("smi", "[CH2]")
    atom = mol.atoms[0].OBAtom
    assert assign_radical_dots(atom) == 2

    set_unpaired_electron_count(atom, 0)
    set_lone_pair_count(atom, 1)

    singlet_penalty = _radical_localization_penalty_for_atom(atom, is_conjugated=False)
    set_lone_pair_count(atom, 0)
    set_unpaired_electron_count(atom, 2)
    triplet_penalty = _radical_localization_penalty_for_atom(atom, is_conjugated=False)

    assert singlet_penalty == pytest.approx(5.0)
    assert triplet_penalty == pytest.approx(5.0)
    assert singlet_penalty == triplet_penalty


def test_non_carbene_active_lone_pair_is_not_a_radical_localization_penalty() -> None:
    mol = pybel.readstring("smi", "O")
    atom = mol.atoms[0].OBAtom
    set_unpaired_electron_count(atom, 0)
    set_lone_pair_count(atom, 1)

    assert _radical_localization_penalty_for_atom(atom, is_conjugated=False) == 0.0


def test_opposite_charges_cancel_within_one_connected_component() -> None:
    mol = pybel.readstring("smi", "C[n+]1nn(c([c-]1)c1ccccc1)C")
    charged_atoms = [atom.OBAtom for atom in mol.atoms if atom.formalcharge]
    raw_penalty = sum(
        _charge_localization_penalty_for_atom(atom, is_conjugated=True) for atom in charged_atoms
    )

    metrics = _compute_organic_electronic_state_metrics(mol)

    assert metrics.charge_localization_component_cancellation == pytest.approx(0.6)
    assert metrics.charge_localization_polarity_inversion_penalty == 0.0
    assert metrics.charge_localization_penalty == pytest.approx(raw_penalty - 0.6)


def test_nonconjugated_single_bond_polarity_inversion_is_penalized() -> None:
    mol = pybel.readstring("smi", "[B-](F)(F)(F)[O+]=N")

    metrics = _compute_organic_electronic_state_metrics(mol)

    assert metrics.charge_localization_component_cancellation == 0.0
    assert metrics.charge_localization_polarity_inversion_penalty == pytest.approx(0.7)
    assert metrics.charge_localization_penalty == pytest.approx(1.4)


def test_small_electronegativity_difference_does_not_define_polarity_inversion() -> None:
    mol = pybel.readstring("smi", "[B-](F)(F)(F)[P+](C)(C)C")

    metrics = _compute_organic_electronic_state_metrics(mol)

    assert metrics.charge_localization_component_cancellation == pytest.approx(0.7)
    assert metrics.charge_localization_polarity_inversion_penalty == 0.0


@pytest.mark.parametrize("smiles", ["[C-]#[O+]", "[N-]=[N+]=[N-]"])
def test_multiply_bonded_charge_pairs_are_not_polarity_inversions(smiles: str) -> None:
    mol = pybel.readstring("smi", smiles)

    metrics = _compute_organic_electronic_state_metrics(mol)

    assert metrics.charge_localization_component_cancellation > 0.0
    assert metrics.charge_localization_polarity_inversion_penalty == 0.0


def test_ring_charge_pair_is_not_treated_as_an_isolated_polarity_inversion() -> None:
    mol = pybel.readstring("smi", "[BH-]1[N+](=C)CCC1")

    metrics = _compute_organic_electronic_state_metrics(mol)

    assert metrics.charge_localization_component_cancellation > 0.0
    assert metrics.charge_localization_polarity_inversion_penalty == 0.0


def test_opposite_charges_do_not_cancel_across_disconnected_components() -> None:
    mol = pybel.readstring("smi", "[NH4+].[Cl-]")
    charged_atoms = [atom.OBAtom for atom in mol.atoms if atom.formalcharge]
    raw_penalty = sum(
        _charge_localization_penalty_for_atom(atom, is_conjugated=False) for atom in charged_atoms
    )

    metrics = _compute_organic_electronic_state_metrics(mol)

    assert metrics.charge_localization_component_cancellation == 0.0
    assert metrics.charge_localization_penalty == pytest.approx(raw_penalty)


def test_opposite_charges_do_not_cancel_across_neutral_atoms() -> None:
    mol = pybel.readstring("smi", "[O-]CC[N+](C)(C)C")

    metrics = _compute_organic_electronic_state_metrics(mol)

    assert metrics.charge_localization_component_cancellation == 0.0
