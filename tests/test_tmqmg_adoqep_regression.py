# pyright: reportMissingImports=false

from __future__ import annotations

from pathlib import Path

import pytest
from rdkit import Chem


pytest.importorskip("openbabel")

from molgr import _core as core
from molgr.interface import xyz_to_rdmol
from molgr.utils.equivalence import EquivalenceMethod, check_equivalence
from scripts.reconstruction_trace import (
    TraceInputCase,
    _html_no_metal_trace,
    _render_html_browser_report,
    trace_reconstruction_case,
)


ADOQEP_REFERENCE_SMILES = "CC1c2ccccn2->[Cu+]23<-N=1N=C(C(=NN->2=C(C)c1ccccn->31)c1ccccc1)c1ccccc1"


def _adoqep_xyz_block() -> str:
    return (Path(__file__).parent / "data" / "xyz" / "ADOQEP.xyz").read_text(encoding="utf-8")


def test_adoqep_cpp_molecule_data_is_closed_shell() -> None:
    mol_data = core.pipeline.reconstruct_with_metals.xyz2omol(
        _adoqep_xyz_block(),
        1,
        0,
    )

    assert mol_data is not None
    assert mol_data.total_radical_num == 0
    assert [
        (atom_idx, atom.atomic_num, atom.formal_charge, atom.radical_num)
        for atom_idx, atom in enumerate(mol_data.atoms, start=1)
        if atom.radical_num
    ] == []


def test_adoqep_cpp_matches_reference_under_nonchiral_tmqmg_equivalence() -> None:
    reference = Chem.MolFromSmiles(ADOQEP_REFERENCE_SMILES)
    assert reference is not None

    reconstructed = xyz_to_rdmol(
        _adoqep_xyz_block(),
        1,
        1,
        backend="cpp",
        make_dative_bonds=True,
    )

    assert sum(atom.GetNumRadicalElectrons() for atom in reconstructed.GetAtoms()) == 0
    equivalent, info = check_equivalence(
        reference,
        reconstructed,
        use_chirality=False,
        max_resonance=100,
    )
    assert equivalent is True, info.reason
    assert info.method == EquivalenceMethod.IDEAL


def test_adoqep_trace_uses_production_no_metal_phase_history() -> None:
    trace = trace_reconstruction_case(
        TraceInputCase(
            id="ADOQEP",
            xyz_block=_adoqep_xyz_block(),
            total_charge=1,
            total_radical_electrons=0,
        ),
        score_all_candidates=False,
    )

    buckets = [
        bucket
        for layer in trace.get("search", {}).get("layer_summaries", [])
        for bucket in layer.get("target_buckets", [])
    ]
    assert buckets
    for bucket in buckets:
        no_metal_trace = bucket.get("no_metal_trace") or {}
        assert no_metal_trace.get("status") != "trace_error"
        assert no_metal_trace.get("status") == "selected"
        phases = [step["phase"] for step in no_metal_trace.get("pipeline_steps", [])]
        assert "prepare_no_metal_seed" in phases
        assert phases[-1] == "select_best_no_metal_candidate"
        assert not no_metal_trace.get("linear_branches")

        resonance = no_metal_trace.get("resonance") or {}
        assert resonance.get("normalization") in {
            "resonance_rule_normalization",
            "full_resonance_normalization",
        }

    static_html = _html_no_metal_trace(buckets[0]["no_metal_trace"])
    assert "生产管线阶段" in static_html
    assert "线性分支" not in static_html

    html = _render_html_browser_report(
        {
            "input": {
                "source": "generic",
                "ids": ["ADOQEP"],
                "total_charge": 1,
                "total_radical_electrons": 0,
            },
            "case_count": 1,
            "cases": [trace],
        }
    )
    assert "tree-toggle" in html
    assert "tree-children" in html
    assert "is-collapsed" in html
    assert "无金属生产管线阶段" in html
    assert "无金属线性分支" not in html
