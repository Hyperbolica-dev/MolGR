from __future__ import annotations

from pathlib import Path

import pytest
from rdkit import Chem


pytest.importorskip("openbabel")

from molgr.interface import xyz_to_rdmol
from molgr.utils.equivalence import check_equivalence


AKUFIV_REFERENCE_SMILES = "O=C1[O-]->[Pt+2]2(<-NC3CCCCC3N->2)<-[O-]C1=O"


def _akufiv_xyz_block() -> str:
    return (Path(__file__).parent / "data" / "xyz" / "AKUFIV.xyz").read_text(encoding="utf-8")


@pytest.mark.parametrize("backend", ["python", "cpp"])
def test_akufiv_selects_platinum_ii_without_pentavalent_carbon(backend: str) -> None:
    reconstructed = xyz_to_rdmol(
        _akufiv_xyz_block(),
        total_charge=0,
        spin_multiplicity=1,
        backend=backend,
        make_dative_bonds=True,
    )
    reference = Chem.MolFromSmiles(AKUFIV_REFERENCE_SMILES)
    assert reference is not None
    platinum = next(atom for atom in reconstructed.GetAtoms() if atom.GetAtomicNum() == 78)

    assert platinum.GetFormalCharge() == 2
    assert all(
        sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) <= 4
        for atom in reconstructed.GetAtoms()
        if atom.GetAtomicNum() == 6
    )
    equivalent, info = check_equivalence(reference, reconstructed, use_chirality=False)
    assert equivalent is True, info.reason
