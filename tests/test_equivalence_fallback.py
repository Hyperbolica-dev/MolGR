from __future__ import annotations

import pytest
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


def test_equivalence_rejects_constitutional_isomers_with_same_formula() -> None:
    ethanol = Chem.MolFromSmiles("CCO")
    dimethyl_ether = Chem.MolFromSmiles("COC")
    assert ethanol is not None
    assert dimethyl_ether is not None

    equivalent, info = equivalence.check_equivalence(
        ethanol,
        dimethyl_ether,
        use_chirality=False,
    )

    assert equivalent is False
    assert info.method is None
    assert info.reason == "Not equivalent: non-metal connectivity differs."


def test_equivalence_rejects_charge_transfer_between_components() -> None:
    mol1 = Chem.MolFromSmiles("[OH+].[SH-]")
    mol2 = Chem.MolFromSmiles("[OH-].[SH+]")
    assert mol1 is not None
    assert mol2 is not None

    equivalent, info = equivalence.check_equivalence(mol1, mol2, use_chirality=False)

    assert equivalent is False
    assert info.method is None
    assert info.reason == "Not equivalent: charge or radical count differs within a component."


def test_resonance_match_compares_the_two_generated_sets() -> None:
    mol1 = Chem.MolFromSmiles("CC1=[N+]=C(OC=[N-])[N]N1")
    mol2 = Chem.MolFromSmiles("Cc1nc(OC=[N])n[nH]1")
    assert mol1 is not None
    assert mol2 is not None
    standardized_1 = equivalence._standardize_metal_bonds(mol1)
    standardized_2 = equivalence._standardize_metal_bonds(mol2)
    organic_1 = equivalence._prepare_organic_mol(standardized_1, already_standardized=True)
    organic_2 = equivalence._prepare_organic_mol(standardized_2, already_standardized=True)

    matched, mol1_count, mol2_count, hit_smiles = equivalence._resonance_match(
        organic_1,
        organic_2,
        use_chirality=False,
        max_resonance=50,
        resonance_flags=Chem.ResonanceFlags.UNCONSTRAINED_CATIONS
        | Chem.ResonanceFlags.UNCONSTRAINED_ANIONS,
    )

    assert matched is True
    assert mol1_count > 0
    assert mol2_count > 0
    assert hit_smiles is not None


def test_resonance_timeout_is_not_swallowed(monkeypatch: pytest.MonkeyPatch) -> None:
    class TimeoutSupplier:
        def __init__(self, *args, **kwargs) -> None:  # noqa: ANN002, ANN003
            pass

        def __iter__(self):
            raise TimeoutError("resonance timeout")

    mol1 = Chem.MolFromSmiles("CCO")
    mol2 = Chem.MolFromSmiles("COC")
    assert mol1 is not None
    assert mol2 is not None
    monkeypatch.setattr(equivalence, "ResonanceMolSupplier", TimeoutSupplier)

    with pytest.raises(TimeoutError, match="resonance timeout"):
        equivalence._resonance_match(
            mol1,
            mol2,
            use_chirality=False,
            max_resonance=50,
            resonance_flags=Chem.ResonanceFlags.UNCONSTRAINED_CATIONS
            | Chem.ResonanceFlags.UNCONSTRAINED_ANIONS,
        )


def test_resonance_sets_are_cached(monkeypatch: pytest.MonkeyPatch) -> None:
    real_supplier = equivalence.ResonanceMolSupplier
    supplier_calls = 0

    def counting_supplier(*args, **kwargs):  # noqa: ANN002, ANN003, ANN202
        nonlocal supplier_calls
        supplier_calls += 1
        return real_supplier(*args, **kwargs)

    mol1 = Chem.MolFromSmiles("CCO")
    mol2 = Chem.MolFromSmiles("COC")
    assert mol1 is not None
    assert mol2 is not None
    equivalence._cached_resonance_smiles.cache_clear()
    monkeypatch.setattr(equivalence, "ResonanceMolSupplier", counting_supplier)

    for _ in range(2):
        equivalence._resonance_match(
            mol1,
            mol2,
            use_chirality=False,
            max_resonance=50,
            resonance_flags=Chem.ResonanceFlags.UNCONSTRAINED_CATIONS
            | Chem.ResonanceFlags.UNCONSTRAINED_ANIONS,
        )

    assert supplier_calls == 2
    equivalence._cached_resonance_smiles.cache_clear()


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


def test_equivalence_suppresses_rdkit_warnings_for_metal_hydrides(capfd) -> None:
    mol1 = Chem.MolFromSmiles("C->[IrH+2]<-N", sanitize=False)
    mol2 = Chem.MolFromSmiles("C->[Ir+2]([H])<-N", sanitize=False)
    assert mol1 is not None
    assert mol2 is not None
    mol1.UpdatePropertyCache(strict=False)
    mol2.UpdatePropertyCache(strict=False)

    equivalent, info = equivalence.check_equivalence(mol1, mol2, use_chirality=False)

    captured = capfd.readouterr()
    assert equivalent is True
    assert info.method == equivalence.EquivalenceMethod.IDEAL
    assert captured.err == ""


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
