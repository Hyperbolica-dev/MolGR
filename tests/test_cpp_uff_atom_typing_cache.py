# pyright: reportMissingImports=false

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest


pytest.importorskip("openbabel")
pytest.importorskip("rdkit")

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom

from molgr import _core
from molgr.config import make_default_config
from molgr.interface import xyz_to_rdmol


RDLogger.DisableLog("rdApp.*")  # type: ignore[arg-type]

_EMBED_SEED = 0xC0FFEE
_RESONANCE_CASE_INDEX = 1


def _load_csv_case(case_index: int) -> tuple[str, int, int]:
    csv_path = Path(__file__).with_name("test_cases.csv")
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    row = rows[case_index - 1]
    mol = Chem.MolFromSmiles(row["smiles"].strip())
    assert mol is not None
    mol_h = Chem.AddHs(mol)
    embed_code = rdDistGeom.EmbedMolecule(mol_h, randomSeed=_EMBED_SEED)  # pyright: ignore[reportCallIssue]
    assert int(embed_code) == 0
    total_charge = sum(int(atom.GetFormalCharge()) for atom in mol_h.GetAtoms())
    total_radical_electrons = sum(int(atom.GetNumRadicalElectrons()) for atom in mol_h.GetAtoms())
    return Chem.MolToXYZBlock(mol_h), total_charge, total_radical_electrons


def _cpp_config(*, enable_uff_atom_typing_cache: bool):
    config = make_default_config()
    return replace(
        config,
        cpp_backend=replace(
            config.cpp_backend,
            enable_uff_atom_typing_cache=enable_uff_atom_typing_cache,
        ),
    )


def _clear_cpp_scoring_caches(*, clear_uff_atom_typing_cache: bool) -> None:
    _core.dev.pipeline.clear_force_field_evaluation_cache()
    _core.dev.pipeline.clear_resonance_move_score_cache()
    if clear_uff_atom_typing_cache:
        _core.dev.pipeline.clear_uff_atom_typing_cache()


def _candidate_summaries(
    xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
    *,
    enable_uff_atom_typing_cache: bool,
) -> list[dict[str, Any]]:
    return list(
        _core.dev.pipeline.reconstruct_without_metals.debug_resonance_candidate_summaries(
            xyz_block,
            total_charge,
            total_radical_electrons,
            config=_cpp_config(enable_uff_atom_typing_cache=enable_uff_atom_typing_cache),
        )
    )


def _assert_candidate_summaries_match(
    left: list[dict[str, Any]],
    right: list[dict[str, Any]],
) -> None:
    assert len(left) == len(right)
    for left_item, right_item in zip(left, right):
        assert left_item["smiles"] == right_item["smiles"]
        assert left_item["resonance_index"] == right_item["resonance_index"]
        assert left_item["aromatic_atom_count"] == right_item["aromatic_atom_count"]
        assert (
            left_item["max_conjugated_component_size"]
            == right_item["max_conjugated_component_size"]
        )
        assert left_item["conjugated_atom_count"] == right_item["conjugated_atom_count"]
        assert left_item["conjugated_bond_count"] == right_item["conjugated_bond_count"]
        assert left_item["score"] == pytest.approx(right_item["score"], abs=1e-12)


def _reconstruction_signature(mol: Chem.Mol) -> tuple[Any, ...]:
    heavy = Chem.RemoveHs(mol)
    canonical_smiles = Chem.MolToSmiles(
        heavy,
        canonical=True,
        isomericSmiles=True,
        allBondsExplicit=True,
        allHsExplicit=True,
    )
    atom_signature = tuple(
        (
            atom.GetIdx(),
            atom.GetAtomicNum(),
            atom.GetFormalCharge(),
            atom.GetNumRadicalElectrons(),
            atom.GetIsAromatic(),
        )
        for atom in mol.GetAtoms()
    )
    bond_signature = tuple(
        sorted(
            (
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                str(bond.GetBondType()),
                bond.GetIsAromatic(),
            )
            for bond in mol.GetBonds()
        )
    )
    return canonical_smiles, atom_signature, bond_signature


def test_cpp_uff_atom_typing_cache_preserves_resonance_candidate_scores() -> None:
    xyz_block, total_charge, total_radical_electrons = _load_csv_case(_RESONANCE_CASE_INDEX)

    _clear_cpp_scoring_caches(clear_uff_atom_typing_cache=True)
    uncached = _candidate_summaries(
        xyz_block,
        total_charge,
        total_radical_electrons,
        enable_uff_atom_typing_cache=False,
    )
    assert uncached
    assert _core.dev.pipeline.get_uff_atom_typing_cache_info()["size"] == 0

    _clear_cpp_scoring_caches(clear_uff_atom_typing_cache=True)
    cached_cold = _candidate_summaries(
        xyz_block,
        total_charge,
        total_radical_electrons,
        enable_uff_atom_typing_cache=True,
    )
    cold_cache_info = _core.dev.pipeline.get_uff_atom_typing_cache_info()
    assert cold_cache_info["size"] > 0
    assert cold_cache_info["misses"] > 0

    _clear_cpp_scoring_caches(clear_uff_atom_typing_cache=False)
    cached_warm = _candidate_summaries(
        xyz_block,
        total_charge,
        total_radical_electrons,
        enable_uff_atom_typing_cache=True,
    )
    warm_cache_info = _core.dev.pipeline.get_uff_atom_typing_cache_info()
    assert warm_cache_info["hits"] > cold_cache_info["hits"]

    _assert_candidate_summaries_match(uncached, cached_cold)
    _assert_candidate_summaries_match(uncached, cached_warm)


def test_cpp_uff_atom_typing_cache_preserves_public_cpp_reconstruction_result() -> None:
    xyz_block, total_charge, total_radical_electrons = _load_csv_case(_RESONANCE_CASE_INDEX)
    spin_multiplicity = total_radical_electrons + 1

    _clear_cpp_scoring_caches(clear_uff_atom_typing_cache=True)
    uncached = xyz_to_rdmol(
        xyz_block,
        total_charge,
        spin_multiplicity,
        backend="cpp",
        config=_cpp_config(enable_uff_atom_typing_cache=False),
    )

    _clear_cpp_scoring_caches(clear_uff_atom_typing_cache=True)
    cached = xyz_to_rdmol(
        xyz_block,
        total_charge,
        spin_multiplicity,
        backend="cpp",
        config=_cpp_config(enable_uff_atom_typing_cache=True),
    )

    assert _reconstruction_signature(uncached) == _reconstruction_signature(cached)
