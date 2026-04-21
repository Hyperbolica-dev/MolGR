# pyright: reportMissingImports=false

from typing import Any, List, Sequence, Tuple

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel

from molgr import _core  # type: ignore
from molgr.fallback.pipeline.resonance import (  # type: ignore
    get_radical_resonances as py_get_radical_resonances,
)
from molgr.fallback.pipeline.resonance import process_resonance as py_process_resonance


_pipeline: Any = _core.dev.pipeline.resonance


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


def _make_seed(smiles: str, radical_atom_indices: Sequence[int]) -> pybel.Molecule:
    mol = pybel.readstring("smi", smiles)
    for idx in radical_atom_indices:
        mol.OBMol.GetAtom(idx).SetSpinMultiplicity(1)
    return mol


def _cpp_get_radical_resonance_ptrs(seed: pybel.Molecule) -> List[int]:
    return [int(ptr) for ptr in _pipeline.get_radical_resonances_ptr(_get_ptr(seed.OBMol))]


def _cpp_get_tokens_from_ptrs(ptrs: Sequence[int]) -> List[str]:
    return [_pipeline.smiles_token_ptr(ptr) for ptr in ptrs]


def _cpp_process_ptr(mol_ptr: int, charge: int) -> Tuple[int, int]:
    out_ptr, out_charge = _pipeline.process_resonance_ptr(mol_ptr, charge)
    return int(out_ptr), int(out_charge)


_CASES: Sequence[Tuple[str, Sequence[int]]] = [
    ("CC=C", (1,)),
    ("C=CC=C", (2,)),
    ("CC=CC=C", (1,)),
    ("C=CC#N", (2,)),
]


@pytest.mark.parametrize(("smiles", "radical_atom_indices"), _CASES)
def test_get_radical_resonances_cpp_matches_fallback_order(
    smiles: str,
    radical_atom_indices: Sequence[int],
) -> None:
    seed = _make_seed(smiles, radical_atom_indices)

    py_resonances = py_get_radical_resonances(_clone_mol(seed))
    cpp_ptrs = _cpp_get_radical_resonance_ptrs(_clone_mol(seed))
    try:
        py_tokens = [_smiles_token(mol) for mol in py_resonances]
        cpp_tokens = _cpp_get_tokens_from_ptrs(cpp_ptrs)
        assert len(py_tokens) > 1
        assert py_tokens == cpp_tokens
    finally:
        for ptr in cpp_ptrs:
            _core.free_obmol_ptr(ptr)


@pytest.mark.parametrize(("smiles", "radical_atom_indices"), _CASES)
@pytest.mark.parametrize("charge", [0, 1, -1])
def test_process_resonance_cpp_matches_fallback(
    smiles: str,
    radical_atom_indices: Sequence[int],
    charge: int,
) -> None:
    seed = _make_seed(smiles, radical_atom_indices)

    py_resonances = py_get_radical_resonances(_clone_mol(seed))
    cpp_ptrs = _cpp_get_radical_resonance_ptrs(_clone_mol(seed))

    try:
        assert len(py_resonances) == len(cpp_ptrs)
        for py_resonance, cpp_ptr in zip(py_resonances, cpp_ptrs):
            py_processed, py_charge, _py_hit = py_process_resonance(
                _clone_mol(py_resonance), charge
            )
            cpp_processed_ptr, cpp_charge = _cpp_process_ptr(cpp_ptr, charge)
            try:
                assert _smiles_token(py_processed) == _pipeline.smiles_token_ptr(cpp_processed_ptr)
                assert py_charge == cpp_charge
            finally:
                _core.free_obmol_ptr(cpp_processed_ptr)
    finally:
        for ptr in cpp_ptrs:
            _core.free_obmol_ptr(ptr)
