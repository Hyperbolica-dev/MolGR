# pyright: reportMissingImports=false

from rdkit import Chem

from molgr.utils import equivalence


def test_equivalence_fallback_uses_inchi_connectivity_when_primary_path_raises(monkeypatch) -> None:
    mol1 = Chem.MolFromSmiles("O=c1cccc[nH]1")
    mol2 = Chem.MolFromSmiles("Oc1ccccn1")

    assert mol1 is not None
    assert mol2 is not None

    def _raise_canon(*_args, **_kwargs):
        raise ValueError("Can't kekulize mol. Unkekulized atoms: 1 3 7")

    monkeypatch.setattr(equivalence, "_canon_smiles", _raise_canon)

    equivalent, info = equivalence.check_equivalence(
        mol1, mol2, use_chirality=True, max_resonance=100
    )

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.INCHI_CONNECTIVITY
    assert "fallback matched" in info.reason
