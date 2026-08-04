from __future__ import annotations

from pathlib import Path

import pytest


pytest.importorskip("openbabel")

from molgr import _core as core
from molgr.fallback.pipeline import reconstruct_without_metals
from molgr.fallback.utils.metals import preparation
from molgr.fallback.utils.organic_topology import compute_organic_topology_metrics


def _adugeo_xyz_block() -> str:
    return (Path(__file__).parent / "data" / "xyz" / "ADUGEO.xyz").read_text(encoding="utf-8")


def test_adugeo_copper_i_bucket_keeps_azomethine_bridge_conjugated() -> None:
    base_state = preparation.prepare_metal_state(_adugeo_xyz_block(), 1, 0)
    no_metal_state = reconstruct_without_metals.xyz_to_omol_no_metal_state(
        base_state.no_metal_xyz_block,
        total_charge=0,
        total_radical_electrons=0,
    )

    assert no_metal_state is not None
    python_metrics = compute_organic_topology_metrics(no_metal_state.omol)
    cpp_metrics = core.dev.utils.compute_organic_topology_metrics_ptr(
        int(no_metal_state.omol.OBMol.this)
    )

    assert python_metrics.conjugated_atom_count == 42
    assert python_metrics.conjugated_bond_count == 48
    assert python_metrics.max_conjugated_component_size == 42
    assert cpp_metrics["conjugated_atom_count"] == 42
    assert cpp_metrics["conjugated_bond_count"] == 48
    assert cpp_metrics["max_conjugated_component_size"] == 42
