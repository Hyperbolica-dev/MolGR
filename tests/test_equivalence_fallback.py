# pyright: reportMissingImports=false

import pytest
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


def test_equivalence_accepts_nhc_carbene_and_zwitterion_forms() -> None:
    dataset_mol = Chem.MolFromSmiles("C=CCN1C=CN(C)[C]1->[Pt+2](<-[I-])(<-[I-])<-n1ccccc1")
    molgr_mol = Chem.MolFromSmiles("C=CCn1cc[n+](C)[c-]1->[Pt@SP3+2](<-[I-])(<-[I-])<-n1ccccc1")

    assert dataset_mol is not None
    assert molgr_mol is not None

    equivalent, info = equivalence.check_equivalence(
        dataset_mol,
        molgr_mol,
        use_chirality=False,
        max_resonance=100,
    )

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.CARBENE_ZWITTERION
    assert info.carbene_zwitterion is not None
    assert info.carbene_zwitterion.mol1_normalized == info.carbene_zwitterion.mol2_normalized


@pytest.mark.parametrize(
    ("reference_smiles", "predicted_smiles"),
    [
        (
            "Cc1cc(C)c(N2C=CN(C3=NC(C)(C)CO3)[C]2->[Rh+]234(<-[Br-])<-C5=C->2CCC->3=C->4CC5)c(C)c1",
            "C1=CCCC=CCC1.Cc1cc(C)c(-n2cc[n+](C3=NC(C)(C)CO3)[c-]2->[Rh+]<-[Br-])c(C)c1",
        ),
        (
            "Cc1cc(C)c(N2C=CN(C3=NC(C)(C)CO3)[C]2->[Rh+]234(<-[Br-])<-C5=C->2C2CC5C->3=C->42)c(C)c1",
            "C1=C[C@H]2C=C[C@@H]1C2.Cc1cc(C)c(-n2cc[n+](C3=NC(C)(C)CO3)[c-]2->[Rh+]<-[Br-])c(C)c1",
        ),
        (
            "Cc1cc(C)c(N2C=CN3C4=N(->[Rh+]567(<-[C]32)<-C2=C->5CCC->6=C->7CC2)C(C)(C)CO4)c(C)c1",
            "C1=CCCC=CCC1.Cc1cc(C)c(-n2cc[n+]3[c-]2->[Rh+]<-N2=C3OCC2(C)C)c(C)c1",
        ),
        (
            "Cc1cc(C)c(N2C=CN3Cc4ccccn4->[Ni+2]4567(<-[C]32)<-c2c->4c->5[cH-]->6c->72)c(C)c1",
            "Cc1cc(C)c(-n2cc[n+]3[c-]2->[Ni@SP2+2](<-[cH-]2cccc2)<-n2ccccc2C3)c(C)c1",
        ),
    ],
)
def test_equivalence_accepts_coordination_stripped_nhc_carbene_and_zwitterion_forms(
    reference_smiles: str,
    predicted_smiles: str,
) -> None:
    reference_mol = Chem.MolFromSmiles(reference_smiles)
    predicted_mol = Chem.MolFromSmiles(predicted_smiles)

    assert reference_mol is not None
    assert predicted_mol is not None

    equivalent, info = equivalence.check_equivalence(
        reference_mol,
        predicted_mol,
        use_chirality=False,
        max_resonance=100,
    )

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.COORDINATION_STRIPPED
    assert info.coordination_stripped is not None


def test_equivalence_does_not_overmatch_unrelated_charge_separated_heterocycle() -> None:
    reference_mol = Chem.MolFromSmiles("CC1=[N+]=C(OC=[N-])[N]N1")
    predicted_mol = Chem.MolFromSmiles("Cc1nc(OC=[N])n[nH]1")

    assert reference_mol is not None
    assert predicted_mol is not None

    equivalent, info = equivalence.check_equivalence(
        reference_mol,
        predicted_mol,
        use_chirality=False,
        max_resonance=100,
    )

    assert equivalent is False
    assert info.method is None


def test_equivalence_accepts_metal_coordination_assignment_mismatch() -> None:
    dataset_mol = Chem.MolFromSmiles("C1CCNC1.[Zn+2].N1CCCC1")
    molgr_mol = Chem.MolFromSmiles("C1CCNC1->[Zn+2]<-N1CCCC1")

    assert dataset_mol is not None
    assert molgr_mol is not None

    equivalent, info = equivalence.check_equivalence(
        dataset_mol,
        molgr_mol,
        use_chirality=False,
        max_resonance=100,
    )

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.COORDINATION_STRIPPED
    assert info.coordination_stripped is not None
    assert info.coordination_stripped.mol1_stripped == info.coordination_stripped.mol2_stripped


def test_equivalence_accepts_sulfimide_resonance_forms() -> None:
    reference_mol = Chem.MolFromSmiles(
        "Cc1cc(C)n(->[Ag+](<-n2c(C)cc(C)cc2C)<-[N-](S(=O)(=O)c2ccc(F)cc2)S(=O)(=O)c2ccc(F)cc2)c(C)c1"
    )
    predicted_mol = Chem.MolFromSmiles(
        "Cc1cc(C)n(->[Ag@SP1+](<-N(S(=O)(=O)c2ccc(F)cc2)=[S@](=O)([O-])c2ccc(F)cc2)<-n2c(C)cc(C)cc2C)c(C)c1"
    )

    assert reference_mol is not None
    assert predicted_mol is not None

    equivalent, info = equivalence.check_equivalence(
        reference_mol,
        predicted_mol,
        use_chirality=False,
        max_resonance=100,
    )

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.RESONANCE


def test_equivalence_handles_none_resonance_supplier_entries(monkeypatch) -> None:
    reference_mol = Chem.MolFromSmiles("CC(N=NC(=S)Nc1ccccc1)=C(C)[N-]N=C([S-])Nc1ccccc1.[Ni+2]")
    predicted_mol = Chem.MolFromSmiles("CC(=NN=C([S-])Nc1ccccc1)C(C)=NN=C([S-])Nc1ccccc1.[Ni+2]")

    assert reference_mol is not None
    assert predicted_mol is not None

    monkeypatch.setattr(equivalence, "ResonanceMolSupplier", lambda *_args, **_kwargs: [None, None])

    equivalent, info = equivalence.check_equivalence(
        reference_mol,
        predicted_mol,
        use_chirality=False,
        max_resonance=100,
        _coordination_bonds_already_stripped=True,
    )

    assert equivalent is False
    assert info.method is None
