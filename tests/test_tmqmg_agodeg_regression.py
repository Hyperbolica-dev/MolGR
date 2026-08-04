from __future__ import annotations

from pathlib import Path

import pytest
from rdkit import Chem


pytest.importorskip("openbabel")

from molgr.interface import xyz_to_rdmol
from molgr.utils.equivalence import check_equivalence


AGODEG_REFERENCE_SMILES = (
    "CP(C)(C)->[Co+3]1(<-[CH3-])(<-[S-]c2cccc3ccc[c-]->1c23)(<-P(C)(C)C)<-P(C)(C)C"
)


def _agodeg_xyz_block() -> str:
    return (Path(__file__).parent / "data" / "xyz" / "AGODEG.xyz").read_text(encoding="utf-8")


def test_agodeg_cpp_preserves_aromatic_ligand_and_selects_cobalt_iii() -> None:
    reconstructed = xyz_to_rdmol(
        _agodeg_xyz_block(),
        total_charge=0,
        spin_multiplicity=1,
        backend="cpp",
        make_dative_bonds=True,
    )
    cobalt = next(atom for atom in reconstructed.GetAtoms() if atom.GetAtomicNum() == 27)
    reference = Chem.MolFromSmiles(AGODEG_REFERENCE_SMILES)
    assert reference is not None

    equivalent, info = check_equivalence(reference, reconstructed, use_chirality=False)

    assert cobalt.GetFormalCharge() == 3
    assert equivalent is True, info.reason
