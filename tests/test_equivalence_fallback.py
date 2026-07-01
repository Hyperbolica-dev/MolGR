from __future__ import annotations

from rdkit import Chem

from molgr.utils import equivalence


def test_equivalence_accepts_resonance_forms_after_standardization() -> None:
    mol1 = Chem.MolFromSmiles("CC1=[N+]=C(OC=[N-])[N]N1")
    mol2 = Chem.MolFromSmiles("Cc1nc(OC=[N])n[nH]1")
    assert mol1 is not None
    assert mol2 is not None

    equivalent, info = equivalence.check_equivalence(mol1, mol2, use_chirality=False)

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.RESONANCE


def test_equivalence_accepts_carbene_zwitterion_forms() -> None:
    mol1 = Chem.MolFromSmiles("[NH3+]CC1=N[N-][N]C1=CO")
    mol2 = Chem.MolFromSmiles("[NH3+]Cc1[n-]nnc1[CH]O")
    assert mol1 is not None
    assert mol2 is not None

    equivalent, info = equivalence.check_equivalence(mol1, mol2, use_chirality=False)

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.RESONANCE


def test_equivalence_rejects_metal_valence_mismatch() -> None:
    mol1 = Chem.MolFromSmiles("C[Cu]N")
    mol2 = Chem.MolFromSmiles("C->[Cu+]<-N")
    assert mol1 is not None
    assert mol2 is not None

    equivalent, info = equivalence.check_equivalence(mol1, mol2, use_chirality=False)

    assert equivalent is False
    assert info.method is None
    assert info.reason == "Not equivalent: metal valence assignment differs."


def test_equivalence_ignores_redundant_metal_radicals() -> None:
    mol1 = Chem.MolFromSmiles("C->[Ru+]<-N", sanitize=False)
    mol2 = Chem.MolFromSmiles("C->[Ru+]<-N", sanitize=False)
    assert mol1 is not None
    assert mol2 is not None
    ru_atom = next(atom for atom in mol1.GetAtoms() if atom.GetSymbol() == "Ru")
    ru_atom.SetNumRadicalElectrons(1)
    mol1.UpdatePropertyCache(strict=False)
    mol2.UpdatePropertyCache(strict=False)

    equivalent, info = equivalence.check_equivalence(mol1, mol2, use_chirality=False)

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.IDEAL


def test_equivalence_treats_metals_as_isolated_ions() -> None:
    mol1 = Chem.MolFromSmiles("C->[Ru+]<-N", sanitize=False)
    mol2 = Chem.MolFromSmiles("C.[Ru+].N", sanitize=False)
    assert mol1 is not None
    assert mol2 is not None
    mol1.UpdatePropertyCache(strict=False)
    mol2.UpdatePropertyCache(strict=False)

    equivalent, info = equivalence.check_equivalence(mol1, mol2, use_chirality=False)

    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.IDEAL


def test_equivalence_rejects_hydrogen_count_mismatch_after_expansion() -> None:
    mol1 = Chem.MolFromSmiles("C=C")
    mol2 = Chem.MolFromSmiles("CC")
    assert mol1 is not None
    assert mol2 is not None

    equivalent, info = equivalence.check_equivalence(mol1, mol2, use_chirality=False)

    assert equivalent is False
    assert info.method is None
    assert info.checks is not None
    assert info.checks.heavy_atom_formula.passed is True
    assert info.checks.explicit_h_formula.passed is False
    assert info.reason == "Not equivalent: explicit-hydrogen element counts differ."


def test_equivalence_rejects_stereoisomers_with_same_inchi_connectivity() -> None:
    mol1 = Chem.MolFromSmiles("C/C=C/C")
    mol2 = Chem.MolFromSmiles("C/C=C\\C")
    assert mol1 is not None
    assert mol2 is not None

    equivalent, info = equivalence.check_equivalence(mol1, mol2)

    assert equivalent is False
    assert info.method is None
    assert info.reason == "Not equivalent: stereochemistry differs."


def test_equivalence_rejects_local_charge_mismatch_after_standardization() -> None:
    mol1 = Chem.MolFromSmiles("CC(C)(C)[N+]#[C-]->[Cu+]1<-[O-]C=C2C=CC=N->12")
    mol2 = Chem.MolFromSmiles("CC(C)(C)[N][C][Cu]1[O][CH][C]2[CH][CH][CH][N]21")
    assert mol1 is not None
    assert mol2 is not None

    equivalent, info = equivalence.check_equivalence(mol1, mol2, use_chirality=False)

    assert equivalent is False
    assert info.method is None
    assert info.reason
