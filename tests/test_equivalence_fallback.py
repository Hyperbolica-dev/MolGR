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
        ("[C]=[Fe]", "[C-2]=[Fe+2]"),
        ("C[C]=[Fe]", "C[C-2]->[Fe+2]"),
    ],
)
def test_equivalence_marks_metal_carbene_valence_assignment_mismatch(
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
    assert info.method == equivalence.EquivalenceMethod.METAL_CARBENE_VALENCE
    assert info.metal_carbene_valence is not None
    assert info.metal_carbene_valence.mol1_normalized == info.metal_carbene_valence.mol2_normalized


def test_equivalence_marks_aqeriy_metal_carbene_valence_assignment_mismatch() -> None:
    reference_mol = Chem.MolFromSmiles(
        "Cc1cc(C)c(N2CC3CCc4ccccc4N3[C]2->[Ru+2]2(<-[Cl-])(<-[Cl-])<-[CH]c3ccccc3O->2C(C)C)c(C)c1"
    )
    predicted_mol = Chem.MolFromSmiles(
        "Cc1cc(C)c([N+]2=[C-](->[Ru@OH23+4]3(<-[Cl-])(<-[Cl-])<-[CH-2]c4ccccc4[O]->3C(C)C)N3c4ccccc4CC[C@H]3C2)c(C)c1"
    )

    assert reference_mol is not None
    assert predicted_mol is not None

    equivalent, info = equivalence.check_equivalence(
        predicted_mol,
        reference_mol,
        use_chirality=False,
        max_resonance=100,
    )

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.METAL_CARBENE_VALENCE
    assert info.metal_carbene_valence is not None
    assert info.metal_carbene_valence.mol1_transformed_count == 1
    assert info.metal_carbene_valence.mol2_transformed_count == 0


def test_equivalence_marks_desroj_metal_carbene_valence_assignment_mismatch() -> None:
    reference_mol = Chem.MolFromSmiles(
        "CC(C)(C)P1(CC2CCCC3CP(C(C)(C)C)(C(C)(C)C)->[RuH+]<-1(<-[Cl-])<-[C]23)C(C)(C)C"
    )
    predicted_mol = Chem.MolFromSmiles(
        "[H][Ru@OH28+3]12(<-[Cl-])<-[C-2]3[C@@H](C[P]->1(C(C)(C)C)C(C)(C)C)CCC[C@H]3C[P]->2(C(C)(C)C)C(C)(C)C"
    )

    assert reference_mol is not None
    assert predicted_mol is not None

    equivalent, info = equivalence.check_equivalence(
        predicted_mol,
        reference_mol,
        use_chirality=False,
        max_resonance=100,
    )

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.METAL_CARBENE_VALENCE
    assert info.metal_carbene_valence is not None
    assert info.metal_carbene_valence.mol1_transformed_count == 1
    assert info.metal_carbene_valence.mol2_transformed_count == 0


def test_equivalence_accepts_deskuk_explicit_and_implicit_metal_hydride_forms() -> None:
    reference_mol = Chem.MolFromSmiles(
        "CC(C)(C)P1(C=CN2=CCP(C(C)(C)C)(C(C)(C)C)->[RuH+]<-2<-1<-N#N)C(C)(C)C"
    )
    predicted_mol = Chem.MolFromSmiles(
        "[H][Ru@OH28+]12(<-[N]#N)<-[N](C=C[P]->1(C(C)(C)C)C(C)(C)C)=CC[P]->2(C(C)(C)C)C(C)(C)C"
    )

    assert reference_mol is not None
    assert predicted_mol is not None

    equivalent, info = equivalence.check_equivalence(
        predicted_mol,
        reference_mol,
        use_chirality=False,
        max_resonance=100,
    )

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.IDEAL


def test_equivalence_does_not_mark_unbound_dianionic_carbon_as_metal_carbene() -> None:
    reference_mol = Chem.MolFromSmiles("[C]")
    predicted_mol = Chem.MolFromSmiles("[C-2]")

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


def test_equivalence_accepts_abartut_azide_stripped_fragment_topology() -> None:
    reference_mol = Chem.MolFromSmiles(
        "CC1c2cccc[n]2->[Cd+2]2(<-[N-2][N+]#N)(<-[N-2][N+]#N)<-[N]=1CC[N]->21CCOCC1.O"
    )
    predicted_mol = Chem.MolFromSmiles(
        "CC1c2cccc[n]2->[Cd+2]2(<-[N-]=[N+]=[N-])(<-[N-]=[N+]=[N-])<-[N]=1CC[N]->21CCOCC1.O"
    )

    assert reference_mol is not None
    assert predicted_mol is not None

    equivalent, info = equivalence.check_equivalence(
        predicted_mol,
        reference_mol,
        use_chirality=False,
        max_resonance=50,
    )

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.COORDINATION_STRIPPED
    assert info.coordination_stripped is not None


def test_equivalence_accepts_abasaa_disconnected_azide_stripped_fragment_topology() -> None:
    reference_mol = Chem.MolFromSmiles(
        "CC1c2cccc[n]2->[Zn+2]2(<-[N-2][N+]#N)(<-[N-2][N+]#N)<-[N]=1CC[N]->21CCOCC1"
    )
    predicted_mol = Chem.MolFromSmiles(
        "CC1c2cccc[n]2->[Zn+2](<-[N-]=[N+]=[N-])(<-[N-]=[N+]=[N-])<-[N]=1CCN1CCOCC1"
    )

    assert reference_mol is not None
    assert predicted_mol is not None

    equivalent, info = equivalence.check_equivalence(
        predicted_mol,
        reference_mol,
        use_chirality=False,
        max_resonance=50,
    )

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.COORDINATION_STRIPPED
    assert info.coordination_stripped is not None


def test_equivalence_accepts_acaloi_radical_mismatch_stripped_fragment_topology() -> None:
    reference_mol = Chem.MolFromSmiles(
        "CC(C)c1cccc(C(C)C)c1N1C=CN(c2c(C(C)C)cccc2C(C)C)[C]1->[Pd+2]12(<-[Cl-])<-[CH2]=[CH]->1[CH-]->2c1ccccc1"
    )
    predicted_mol = Chem.MolFromSmiles(
        "CC(C)c1cccc(C(C)C)c1-n1cc[n+](-c2c(C(C)C)cccc2C(C)C)[c-]1->[Pd@SP2+2](<-[Cl-])<-[CH2-]C=Cc1ccccc1"
    )

    assert reference_mol is not None
    assert predicted_mol is not None

    equivalent, info = equivalence.check_equivalence(
        predicted_mol,
        reference_mol,
        use_chirality=False,
        max_resonance=50,
    )

    assert equivalent is True
    assert info.checks is not None
    assert info.checks.radical_electrons.passed is False
    assert info.method == equivalence.EquivalenceMethod.COORDINATION_STRIPPED
    assert info.coordination_stripped is not None


def test_equivalence_rejects_stripped_fragment_composition_mismatch() -> None:
    reference_mol = Chem.MolFromSmiles("NN.NN.[Zn+2]")
    predicted_mol = Chem.MolFromSmiles("NNN.N.[Zn+2]")

    assert reference_mol is not None
    assert predicted_mol is not None

    equivalent, info = equivalence.check_equivalence(
        reference_mol,
        predicted_mol,
        use_chirality=False,
        max_resonance=50,
    )

    assert equivalent is False
    assert info.method is None
    assert "fragment element composition differs" in info.reason


def test_equivalence_accepts_denjoz_stripped_fragment_topology() -> None:
    reference_mol = Chem.MolFromSmiles("CC1=C2C=CC=C[N-]2->[Ni+2]2(<-[Cl-])<-S=C(N=N->21)Nc1ccccc1")
    predicted_mol = Chem.MolFromSmiles(
        "CC1c2cccc[n]2->[Ni@SP3+2]2(<-[Cl-])<-[S-]C(=N[N]->2=1)Nc1ccccc1"
    )

    assert reference_mol is not None
    assert predicted_mol is not None

    equivalent, info = equivalence.check_equivalence(
        predicted_mol,
        reference_mol,
        use_chirality=False,
        max_resonance=50,
    )

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.COORDINATION_STRIPPED
    assert info.coordination_stripped is not None


def test_equivalence_accepts_abazek_coordination_stripped_without_full_mapping(
    monkeypatch,
) -> None:
    reference_mol = Chem.MolFromSmiles(
        "CN(C)c1ccc(P(C(C)(C)C)(C(C)(C)C)->[Pd+2]23(<-[Cl-])<-[CH2-]C->2=C->3c2ccccc2)cc1"
    )
    predicted_mol = Chem.MolFromSmiles(
        "CN(C)c1ccc([P](->[Pd+2](<-[Cl-])<-[CH2-]C=Cc2ccccc2)(C(C)(C)C)C(C)(C)C)cc1"
    )

    assert reference_mol is not None
    assert predicted_mol is not None

    def _raise_full_mapping(*_args, **_kwargs):
        raise AssertionError("full topology mapping should be skipped")

    monkeypatch.setattr(equivalence, "_full_simplified_topology_matches", _raise_full_mapping)

    equivalent, info = equivalence.check_equivalence(
        predicted_mol,
        reference_mol,
        use_chirality=False,
        max_resonance=50,
    )

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.COORDINATION_STRIPPED


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
    assert info.method == equivalence.EquivalenceMethod.COORDINATION_STRIPPED
    assert info.coordination_stripped is not None


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
