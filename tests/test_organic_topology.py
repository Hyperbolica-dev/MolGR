from __future__ import annotations

from dataclasses import replace

import pytest

from molgr.config import MolGRConfig
from molgr.fallback.utils.organic_topology import (
    _additional_atom_pi_electrons,
    compute_organic_topology_metrics,
)


def _require_pybel():
    pytest.importorskip("openbabel")
    return pytest.importorskip("openbabel.pybel")


def test_continuous_conjugation_metrics_count_only_three_bond_alternation() -> None:
    pybel = _require_pybel()
    benzene = pybel.readstring("smi", "c1ccccc1")
    diene = pybel.readstring("smi", "C=CC=C")

    benzene_metrics = compute_organic_topology_metrics(benzene)
    diene_metrics = compute_organic_topology_metrics(diene)

    assert benzene_metrics.aromatic_atom_count == 6
    assert benzene_metrics.aromatic_ring_count == 1
    assert benzene_metrics.aromatic_stability_score == pytest.approx(1.0)
    assert benzene_metrics.conjugated_atom_count == 6
    assert benzene_metrics.conjugated_bond_count == 6
    assert benzene_metrics.max_conjugated_component_size == 6

    assert diene_metrics.aromatic_atom_count == 0
    assert diene_metrics.aromatic_ring_count == 0
    assert diene_metrics.aromatic_stability_score == pytest.approx(0.0)
    assert diene_metrics.conjugated_atom_count == 4
    assert diene_metrics.conjugated_bond_count == 3
    assert diene_metrics.max_conjugated_component_size == 4


def test_charge_separated_fragment_without_bond_alternation_is_not_counted_as_conjugated() -> None:
    pybel = _require_pybel()
    alternating = pybel.readstring("smi", "[O-]/C(=C\\C(=O)C(F)(F)F)/C(F)(F)F")
    charge_separated = pybel.readstring("smi", "O=C([CH+]C(=O)C(F)(F)F)C(F)(F)F")

    alternating_metrics = compute_organic_topology_metrics(alternating)
    charge_separated_metrics = compute_organic_topology_metrics(charge_separated)

    assert charge_separated_metrics.conjugated_atom_count == 0
    assert charge_separated_metrics.conjugated_bond_count == 0
    assert charge_separated_metrics.max_conjugated_component_size == 0
    assert (
        alternating_metrics.conjugated_atom_count > charge_separated_metrics.conjugated_atom_count
    )
    assert (
        alternating_metrics.conjugated_bond_count > charge_separated_metrics.conjugated_bond_count
    )
    assert (
        alternating_metrics.max_conjugated_component_size
        > charge_separated_metrics.max_conjugated_component_size
    )


def test_high_absolute_formal_charge_sum_rejects_ring_aromaticity() -> None:
    pybel = _require_pybel()
    tolerated = pybel.readstring("smi", "[c-]1[c-][c-][cH][cH][cH]1")
    rejected = pybel.readstring("smi", "[c-]1[c-][c-][c-][cH][cH]1")

    tolerated_metrics = compute_organic_topology_metrics(tolerated)
    rejected_metrics = compute_organic_topology_metrics(rejected)

    assert tolerated_metrics.aromatic_atom_count == 6
    assert tolerated_metrics.aromatic_ring_count == 1
    assert rejected_metrics.aromatic_atom_count == 0
    assert rejected_metrics.aromatic_ring_count == 0


def test_ring_aromaticity_requires_huckel_pi_electron_count() -> None:
    pybel = _require_pybel()
    accepted = pybel.readstring("smi", "[n-]1cccc1")
    rejected = pybel.readstring("smi", "c1cc2ccccc2cc1")

    accepted_metrics = compute_organic_topology_metrics(accepted)
    rejected_metrics = compute_organic_topology_metrics(rejected)

    assert accepted_metrics.aromatic_ring_count == 1
    assert accepted_metrics.aromatic_atom_count == 5
    assert rejected_metrics.aromatic_ring_count == 1
    assert rejected_metrics.aromatic_atom_count == 6


def test_negative_ring_atom_contributes_pi_electrons_when_not_in_ring_multiple_bond() -> None:
    pybel = _require_pybel()
    omol = pybel.readstring("smi", "[N-]1C=CC=C1")
    negative_nitrogen = omol.OBMol.GetAtom(1)

    assert int(negative_nitrogen.GetAtomicNum()) == 7
    assert int(negative_nitrogen.GetFormalCharge()) == -1
    assert (
        _additional_atom_pi_electrons(
            negative_nitrogen,
            incident_to_ring_multiple_bond=False,
        )
        == 2
    )
    assert (
        _additional_atom_pi_electrons(
            negative_nitrogen,
            incident_to_ring_multiple_bond=True,
        )
        == 0
    )


def test_acasoo_like_resonance_prefers_huckel_valid_aromatic_ring_count() -> None:
    pybel = _require_pybel()
    resonance_zero = pybel.readstring(
        "smi",
        r"C\1/2=C(/C3=C(c4c(=N3)/c(=C\3/C=CC(=N3)/C(=c\3/cc/c(=C(/C(=N2)C=C1)\c1ccccc1)/[n-]3)/c1ccccc1)/c1ccccc1c4[O-])N)\c1ccccc1",
    )
    resonance_one = pybel.readstring(
        "smi",
        r"C\1/2=C(/c3c(c4c(/C(=C\5/C=CC(=N5)/C(=c\5/cc/c(=C(/C(=N2)C=C1)\c1ccccc1)/[n-]5)/c1ccccc1)/c1ccccc1C4=O)[n-]3)N)\c1ccccc1",
    )

    zero_metrics = compute_organic_topology_metrics(resonance_zero)
    one_metrics = compute_organic_topology_metrics(resonance_one)

    assert zero_metrics.aromatic_ring_count == 4
    assert zero_metrics.aromatic_atom_count == 24
    assert zero_metrics.aromatic_stability_score == pytest.approx(4.0)
    assert one_metrics.aromatic_ring_count == 5
    assert one_metrics.aromatic_atom_count == 29
    assert one_metrics.aromatic_stability_score == pytest.approx(4.63504)
    assert one_metrics.aromatic_stability_score > zero_metrics.aromatic_stability_score


def test_aromatic_stability_scores_benzene_above_heteroaromatics() -> None:
    pybel = _require_pybel()
    benzene = pybel.readstring("smi", "c1ccccc1")
    pyridine = pybel.readstring("smi", "n1ccccc1")
    pyrrole = pybel.readstring("smi", "c1cc[nH]c1")
    thiophene = pybel.readstring("smi", "c1ccsc1")

    benzene_metrics = compute_organic_topology_metrics(benzene)
    hetero_scores = [
        compute_organic_topology_metrics(mol).aromatic_stability_score
        for mol in (pyridine, pyrrole, thiophene)
    ]

    assert benzene_metrics.aromatic_stability_score == pytest.approx(1.0)
    assert all(0.0 < score < 1.0 for score in hetero_scores)


def test_aromatic_stability_uses_configured_scale_factors() -> None:
    pybel = _require_pybel()
    pyridine = pybel.readstring("smi", "n1ccccc1")
    config = replace(
        MolGRConfig().organic_topology,
        aromatic_stability_ring_size_6_factor=0.50,
        aromatic_stability_hetero_atom_penalty=0.25,
        aromatic_stability_other_ring_max_score=0.99,
    )

    metrics = compute_organic_topology_metrics(pyridine, config=config)

    assert metrics.aromatic_stability_score == pytest.approx(0.50 * 0.75)


def test_abated_like_quad_anionic_ring_is_not_counted_as_aromatic() -> None:
    pybel = _require_pybel()
    abated_like = pybel.readstring(
        "smi",
        "C/C(=C/C(=N/c1c(C)cccc1C)/C)/[N-]c1c(C)cccc1C.C[c-]1[cH-][cH-][cH-]cc1",
    )
    neutral_reference = pybel.readstring(
        "smi",
        "C/C(=C/C(=N/c1c(C)cccc1C)/C)/[N-]c1c(C)cccc1C.Cc1ccccc1",
    )

    abated_like_metrics = compute_organic_topology_metrics(abated_like)
    neutral_reference_metrics = compute_organic_topology_metrics(neutral_reference)

    assert neutral_reference_metrics.aromatic_ring_count == 3
    assert neutral_reference_metrics.aromatic_atom_count == 18
    assert abated_like_metrics.aromatic_ring_count == 2
    assert abated_like_metrics.aromatic_atom_count == 12


def test_re_dithiocarbamate_variants_do_not_get_artificial_conjugation_gap() -> None:
    pybel = _require_pybel()
    zwitterionic = pybel.readstring(
        "smi",
        "CP(c1ccccc1)c1ccccc1.CC[N+](=C([S+])[S+])CC.[C-]#[O+].[C-]#[O+].[C-]#[O+]",
    )
    anionic = pybel.readstring(
        "smi",
        "CP(c1ccccc1)c1ccccc1.CCN(C(=S)[S-])CC.[C-]#[O+].[C-]#[O+].[C-]#[O+]",
    )

    zwitterionic_metrics = compute_organic_topology_metrics(zwitterionic)
    anionic_metrics = compute_organic_topology_metrics(anionic)

    assert (
        abs(zwitterionic_metrics.conjugated_atom_count - anionic_metrics.conjugated_atom_count) <= 1
    )
    assert (
        abs(zwitterionic_metrics.conjugated_bond_count - anionic_metrics.conjugated_bond_count) <= 1
    )
    assert (
        zwitterionic_metrics.max_conjugated_component_size
        == anionic_metrics.max_conjugated_component_size
        == 6
    )


def test_au_nhc_zwitterion_variant_keeps_conjugation_component() -> None:
    pybel = _require_pybel()
    charge_separated = pybel.readstring(
        "smi",
        "Cc1cc(C)c(c(c1)C)[N+]1=C([S+]=C=C1c1ccccc1)C1=CC=CC=C1.[Cl-]",
    )
    aromatic_zwitterion = pybel.readstring(
        "smi",
        "Cc1cc(C)c(c(c1)C)n1c([s+][c-]c1c1ccccc1)c1ccccc1.[Cl-]",
    )

    charge_separated_metrics = compute_organic_topology_metrics(charge_separated)
    aromatic_zwitterion_metrics = compute_organic_topology_metrics(aromatic_zwitterion)

    assert (
        aromatic_zwitterion_metrics.max_conjugated_component_size
        >= charge_separated_metrics.max_conjugated_component_size
    )
    assert (
        aromatic_zwitterion_metrics.conjugated_atom_count
        >= charge_separated_metrics.conjugated_atom_count
    )
