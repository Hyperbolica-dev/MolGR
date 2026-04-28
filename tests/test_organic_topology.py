from __future__ import annotations

import pytest

from molgr.fallback.utils.organic_topology import compute_organic_topology_metrics


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
    assert benzene_metrics.conjugated_atom_count == 6
    assert benzene_metrics.conjugated_bond_count == 6
    assert benzene_metrics.max_conjugated_component_size == 6

    assert diene_metrics.aromatic_atom_count == 0
    assert diene_metrics.aromatic_ring_count == 0
    assert diene_metrics.conjugated_atom_count == 4
    assert diene_metrics.conjugated_bond_count == 3
    assert diene_metrics.max_conjugated_component_size == 4


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
