from __future__ import annotations

import hashlib
import json
import math
import re
import struct
from collections import Counter
from pathlib import Path

import pytest
from rdkit import Chem

from molgr.interface import xyz_to_rdmol
from molgr.utils.equivalence import EquivalenceDecision, evaluate_equivalence


ROOT = Path(__file__).parent / "data/borylxat"
MANIFEST = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))
CASES = MANIFEST["fixtures"]


def _reference(case):
    mol = next(Chem.SDMolSupplier(str(ROOT / case["reference_file"]), removeHs=False))
    assert mol is not None
    return mol


def _atoms(case):
    mol = _reference(case)
    conf = mol.GetConformer()
    return mol.GetNumAtoms(), [
        [atom.GetSymbol(), *conf.GetAtomPosition(atom.GetIdx())] for atom in mol.GetAtoms()
    ]


def _xyz_input(case):
    # Pass only elements/coordinates: reference bonds must not enter reconstruction.
    count, atoms = _atoms(case)
    return (
        "\n".join(
            [str(count), case["case_id"]]
            + [symbol + " " + " ".join(format(x, ".17g") for x in xyz) for symbol, *xyz in atoms]
        )
        + "\n"
    )


def test_selection_and_source_provenance():
    selection = json.loads((ROOT / "selection.json").read_text(encoding="utf-8"))
    assert MANIFEST["schema_version"] == 2
    assert MANIFEST["coordinate_units"] == "angstrom"
    assert MANIFEST["atom_index_base"] == 0
    assert MANIFEST["source_dataset"]["doi"] == "10.5281/zenodo.21330198"
    assert MANIFEST["source_dataset"]["published_md5"] == "d486fbcc36f856bbebf303b5cda4950f"
    assert MANIFEST["source_dataset"]["license"] == "CC-BY-4.0"
    assert len({case["case_id"] for case in CASES}) == len(CASES) == 30
    assert Counter(case["category"] for case in CASES) == {"B": 20, "complex_r": 5, "complex_p": 5}
    assert [
        {key: case[key] for key in selected} for case, selected in zip(CASES, selection)
    ] == selection
    assert len(selection) == len(CASES)
    assert {p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*.sdf")} == {
        case["reference_file"] for case in CASES
    }
    assert not list(ROOT.rglob("*.xyz"))


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["case_id"])
def test_coordinates_and_electronic_state_evidence(case):
    count, atoms = _atoms(case)
    assert count == len(atoms) == case["natoms"]
    assert all(len(atom) == 4 and all(math.isfinite(float(x)) for x in atom[1:]) for atom in atoms)
    symbols = [atom[0] for atom in atoms]
    numbers = [Chem.GetPeriodicTable().GetAtomicNumber(symbol) for symbol in symbols]
    # Compare SDF elements, atom order and full-precision coordinates with the
    # digest of the original SQLite number/position blobs, independently of SDF.
    geometry = struct.pack(f"<{count}i", *numbers) + struct.pack(
        f"<{3 * count}d", *(x for atom in atoms for x in atom[1:])
    )
    assert hashlib.sha256(geometry).hexdigest() == case["source_geometry_sha256"]
    formula = Counter(
        {
            symbol: int(n or 1)
            for symbol, n in re.findall(r"([A-Z][a-z]?)(\d*)", case["source_formula"])
        }
    )
    assert Counter(symbols) == formula
    assert [i for i, symbol in enumerate(symbols) if symbol == "B"] == case["boron_atom_indices"]
    assert len(case["boron_atom_indices"]) == 1
    electrons = (
        sum(Chem.GetPeriodicTable().GetAtomicNumber(symbol) for symbol in symbols) - case["charge"]
    )
    assert electrons % 2 == (case["spin_multiplicity"] - 1) % 2
    assert case["reference_status"] == "corrected_source_graph_geometry_checked"
    assert case["source_smiles"]
    assert case["reference_smiles"]
    spins = case["source_spin_densities"]
    if case["category"] == "complex_p":
        assert case["spin_multiplicity"] == 1
        assert case["electronic_state_basis"].startswith("author workflow: B_N_p singlet")
        assert spins is None
    else:
        assert case["spin_multiplicity"] == 2
        expected_basis = (
            "inferred doublet" if case["category"] == "B" else "author workflow: B_N_r doublet"
        )
        assert case["electronic_state_basis"].startswith(expected_basis)
        assert len(spins) == count
        assert all(math.isfinite(x) for x in spins)
        assert sum(spins) == pytest.approx(1.0, abs=1e-4)
        assert spins[case["boron_atom_indices"][0]] > 0.3


def test_product_controls_match_radical_components_plus_chlorine():
    by_id = {case["case_id"]: case for case in CASES}
    for radical in CASES:
        if radical["category"] != "complex_r":
            continue
        product = by_id[radical["case_id"][:-1] + "p"]
        assert product["category"] == "complex_p"
        assert product["charge"] == radical["charge"]
        _, r_atoms = _atoms(radical)
        _, p_atoms = _atoms(product)
        assert Counter(atom[0] for atom in p_atoms) == Counter(
            atom[0] for atom in r_atoms
        ) + Counter(Cl=1)


def test_source_smiles_hydrogen_mismatch_is_not_silently_used_as_an_answer():
    case = next(case for case in CASES if case["case_id"] == "B_00409")
    reference = Chem.MolFromSmiles(case["source_smiles"])
    assert reference is not None
    with_hydrogens = Chem.AddHs(reference)
    assert with_hydrogens.GetNumAtoms() == case["natoms"] + 1
    _, atoms = _atoms(case)
    assert Counter(atom[0] for atom in atoms) == Counter(B=1, F=2)


def test_delocalized_radical_evidence_is_preserved():
    case = next(case for case in CASES if case["case_id"] == "B_00440")
    _, atoms = _atoms(case)
    spins = case["source_spin_densities"]
    boron_spin = spins[case["boron_atom_indices"][0]]
    assert boron_spin == pytest.approx(0.530277)
    assert max(spin for atom, spin in zip(atoms, spins) if atom[0] == "C") > boron_spin


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["case_id"])
def test_corrected_answer_graph_roundtrip_and_electronic_state(case):
    raw = (ROOT / case["reference_file"]).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == case["reference_sha256"]
    mol = next(Chem.SDMolSupplier(str(ROOT / case["reference_file"]), removeHs=False))
    assert mol is not None
    assert mol.GetNumAtoms() == case["natoms"]
    assert len(Chem.GetMolFrags(mol)) == 1
    assert Chem.GetFormalCharge(mol) == case["charge"]
    assert mol.GetIntProp("TOTAL_CHARGE") == case["charge"]
    assert mol.GetIntProp("SPIN_MULTIPLICITY") == case["spin_multiplicity"]
    assert mol.GetProp("SOURCE_SMILES") == case["source_smiles"]
    assert mol.GetProp("ELECTRONIC_STATE_BASIS") == case["electronic_state_basis"]
    assert json.loads(mol.GetProp("ELECTRONIC_STATE_SOURCES")) == case["electronic_state_sources"]
    for atom in mol.GetAtoms():
        assert atom.GetNumImplicitHs() == atom.GetNumExplicitHs() == 0
        expected_radicals = case["spin_multiplicity"] - 1 if atom.GetSymbol() == "B" else 0
        assert atom.GetNumRadicalElectrons() == expected_radicals
    boron = mol.GetAtomWithIdx(case["boron_atom_indices"][0])
    if case["category"] == "B":
        assert boron.GetFormalCharge() == 0
        assert boron.GetDegree() == 2
    else:
        assert boron.GetFormalCharge() == -1
        assert boron.GetDegree() == (3 if case["category"] == "complex_r" else 4)
        positive = [
            mol.GetAtomWithIdx(i)
            for i in range(mol.GetNumAtoms())
            if mol.GetAtomWithIdx(i).GetFormalCharge() > 0
        ]
        assert len(positive) == 1
        assert positive[0].GetFormalCharge() == 1
        assert mol.GetBondBetweenAtoms(positive[0].GetIdx(), boron.GetIdx()) is not None
    assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == case["reference_smiles"]
    smiles_roundtrip = Chem.AddHs(Chem.MolFromSmiles(case["reference_smiles"]))
    assert smiles_roundtrip.GetNumAtoms() == case["natoms"]
    assert (
        sum(
            smiles_roundtrip.GetAtomWithIdx(i).GetNumRadicalElectrons()
            for i in range(smiles_roundtrip.GetNumAtoms())
        )
        == case["spin_multiplicity"] - 1
    )
    if case["source_spin_densities"] is not None:
        assert json.loads(mol.GetProp("SOURCE_SPIN_DENSITIES")) == case["source_spin_densities"]


def test_corrected_answers_keep_pf3_bonds_and_remove_phantom_boron_hydrogen():
    by_id = {case["case_id"]: case for case in CASES}
    bf2 = Chem.MolFromSmiles(by_id["B_00409"]["reference_smiles"])
    assert Chem.MolToSmiles(bf2) == "F[B]F"
    case = by_id["B_00388_LB_00141_r"]
    mol = next(Chem.SDMolSupplier(str(ROOT / case["reference_file"]), removeHs=False))
    phosphorus = next(
        mol.GetAtomWithIdx(i)
        for i in range(mol.GetNumAtoms())
        if mol.GetAtomWithIdx(i).GetAtomicNum() == 15
    )
    assert Counter(a.GetSymbol() for a in phosphorus.GetNeighbors()) == Counter(F=3, B=1)


def test_corrected_bicyclic_answer_keeps_source_stereochemistry_consistent_with_atoms():
    case = next(case for case in CASES if case["case_id"] == "B_00413")
    mol = next(Chem.SDMolSupplier(str(ROOT / case["reference_file"]), removeHs=False))
    before = Chem.MolToSmiles(mol)
    assert "@" in before
    Chem.AssignStereochemistryFrom3D(mol)
    assert Chem.MolToSmiles(mol) == before


@pytest.mark.parametrize("backend", ["cpp", "python"])
@pytest.mark.parametrize("case_id", ["B_00388", "B_00409"])
def test_simple_boryl_radical_reconstruction_preserves_electronic_state(case_id, backend):
    case = next(case for case in CASES if case["case_id"] == case_id)
    mol = xyz_to_rdmol(
        _xyz_input(case),
        total_charge=case["charge"],
        spin_multiplicity=case["spin_multiplicity"],
        backend=backend,
    )
    assert mol.GetNumAtoms() == case["natoms"]
    assert len(Chem.GetMolFrags(mol)) == 1
    boron = mol.GetAtomWithIdx(case["boron_atom_indices"][0])
    assert boron.GetAtomicNum() == 5
    assert boron.GetDegree() == 2
    assert boron.GetFormalCharge() == 0
    assert boron.GetNumRadicalElectrons() == 1
    assert (
        sum(mol.GetAtomWithIdx(i).GetNumRadicalElectrons() for i in range(mol.GetNumAtoms())) == 1
    )


@pytest.mark.parametrize("backend", ["cpp", "python"])
@pytest.mark.parametrize(
    "case_id",
    [
        "B_00388_LB_00001_r",
        "B_00388_LB_00069_r",
        "B_00388_LB_00141_r",
        "B_00388_LB_00146_r",
        "B_00388_LB_00209_r",
    ],
)
def test_coordinated_boryl_radical_keeps_bonds_and_boron_electronic_state(case_id, backend):
    case = next(case for case in CASES if case["case_id"] == case_id)
    mol = xyz_to_rdmol(
        _xyz_input(case),
        total_charge=case["charge"],
        spin_multiplicity=case["spin_multiplicity"],
        backend=backend,
    )
    reference = _reference(case)
    atoms = [mol.GetAtomWithIdx(i) for i in range(mol.GetNumAtoms())]
    bonds = [mol.GetBondWithIdx(i) for i in range(mol.GetNumBonds())]
    assert [a.GetAtomicNum() for a in atoms] == [a.GetAtomicNum() for a in reference.GetAtoms()]
    assert {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in bonds} == {
        tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in reference.GetBonds()
    }
    assert len(Chem.GetMolFrags(mol)) == 1, Chem.MolToSmiles(mol)
    assert Chem.GetFormalCharge(mol) == case["charge"]
    boron = mol.GetAtomWithIdx(case["boron_atom_indices"][0])
    assert boron.GetFormalCharge() == -1
    assert boron.GetNumRadicalElectrons() == 1
    assert sum(a.GetNumRadicalElectrons() for a in atoms) == 1
    assert evaluate_equivalence(mol, reference).decision == EquivalenceDecision.EQUIVALENT


@pytest.mark.parametrize("backend", ["cpp", "python"])
def test_b00401_accepts_two_step_radical_resonance(backend):
    case = next(case for case in CASES if case["case_id"] == "B_00401")
    reference = _reference(case)
    mol = xyz_to_rdmol(
        _xyz_input(case),
        total_charge=case["charge"],
        spin_multiplicity=case["spin_multiplicity"],
        backend=backend,
    )
    # A fixture-specific witness, not a general relaxation of equivalence:
    # B0*-C4=C3-C2=C1 -> B0=C4-C3*-C2=C1 -> B0=C4-C3=C2-C1*.
    intermediate = Chem.Mol(reference)
    intermediate.GetBondBetweenAtoms(0, 4).SetBondType(Chem.BondType.DOUBLE)
    intermediate.GetBondBetweenAtoms(4, 3).SetBondType(Chem.BondType.SINGLE)
    intermediate.GetAtomWithIdx(0).SetNumRadicalElectrons(0)
    intermediate.GetAtomWithIdx(3).SetNumRadicalElectrons(1)
    Chem.SanitizeMol(intermediate)
    assert evaluate_equivalence(reference, intermediate).decision == EquivalenceDecision.EQUIVALENT
    assert any(
        evaluate_equivalence(mol, witness).decision == EquivalenceDecision.EQUIVALENT
        for witness in (reference, intermediate)
    )


@pytest.mark.parametrize("backend", ["cpp", "python"])
@pytest.mark.parametrize("element,distance", [("Si", 1.85), ("Ge", 1.95)])
def test_synthetic_boron_noncarbon_double_bond_reduction(backend, element, distance):
    # Synthetic probes: preparation produces [BH2-]=XH2 with a +1 budget.
    xyz = (
        f"6\nBH2-XH2 model probe\nB 0 0 0\n{element} {distance} 0 0\n"
        "H -0.6 1.039 0\nH -0.6 -1.039 0\n"
        f"H {distance + 0.7} 1.2 0\nH {distance + 0.7} -1.2 0\n"
    )
    mol = xyz_to_rdmol(xyz, total_charge=0, spin_multiplicity=2, backend=backend)
    atoms = [mol.GetAtomWithIdx(i) for i in range(mol.GetNumAtoms())]
    bonds = [mol.GetBondWithIdx(i) for i in range(mol.GetNumBonds())]
    assert mol.GetNumAtoms() == 6
    assert {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in bonds} == {
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 4),
        (1, 5),
    }
    assert all(b.GetBondType() == Chem.BondType.SINGLE for b in bonds)
    assert [a.GetFormalCharge() for a in atoms] == [-1, 1, 0, 0, 0, 0]
    assert [a.GetNumRadicalElectrons() for a in atoms] == [1, 0, 0, 0, 0, 0]
    assert mol.GetAtomWithIdx(1).GetSymbol() == element


@pytest.mark.parametrize("backend", ["cpp", "python"])
@pytest.mark.parametrize("ligand,distance", [("H", 1.2), ("F", 1.4)])
@pytest.mark.parametrize("charge,multiplicity", [(0, 1), (-1, 2)])
def test_synthetic_bv3_reduction_preserves_all_bonds(
    backend, ligand, distance, charge, multiplicity
):
    # Idealized geometry probes, not additional BorylXAT database records.
    xyz = (
        "4\nsynthetic trigonal borane\nB 0 0 0\n"
        + "\n".join(
            f"{ligand} {distance * math.cos(angle):.17g} {distance * math.sin(angle):.17g} 0"
            for angle in (0, 2 * math.pi / 3, 4 * math.pi / 3)
        )
        + "\n"
    )
    mol = xyz_to_rdmol(xyz, total_charge=charge, spin_multiplicity=multiplicity, backend=backend)
    assert mol.GetNumAtoms() == 4
    assert mol.GetNumBonds() == 3
    assert len(Chem.GetMolFrags(mol)) == 1
    boron = mol.GetAtomWithIdx(0)
    assert boron.GetAtomicNum() == 5
    assert boron.GetDegree() == 3
    assert boron.GetFormalCharge() == charge
    assert boron.GetNumRadicalElectrons() == multiplicity - 1
    assert Chem.GetFormalCharge(mol) == charge
    assert (
        sum(mol.GetAtomWithIdx(i).GetNumRadicalElectrons() for i in range(mol.GetNumAtoms()))
        == multiplicity - 1
    )
