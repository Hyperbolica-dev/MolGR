# pyright: reportMissingImports=false

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import pytest


pytest.importorskip("openbabel")

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom

from molgr.config import make_default_config
from molgr.fallback.pipeline import reconstruct_without_metals as no_metal_module
from molgr.fallback.pipeline.reconstruct_without_metals import (
    xyz_to_omol_no_metal,
    xyz_to_omol_no_metal_state,
)
from molgr.fallback.stages.clean import clean_neighbor_radicals, clean_resonances
from molgr.fallback.stages.eliminate import (
    eliminate_carbene_neighbor_heteroatom,
    eliminate_charge_spliting,
    eliminate_CN_in_doubt,
    eliminate_high_positive_charge_atoms,
    eliminate_negative_charges,
    eliminate_NNN,
)
from molgr.fallback.stages.fresh import fresh_omol_charge_radical
from molgr.fallback.stages.preprocess import make_connections
from molgr.fallback.state import ReconstructionState
from molgr.fallback.utils.no_metals import preparation as no_metal_preparation_module
from molgr.fallback.utils.no_metals import resonance as no_metal_resonance_module
from molgr.fallback.utils.no_metals import selection as no_metal_selection_module
from molgr.fallback.utils.tools import typed_lru_cache
from molgr.utils.converter import pybel_to_rdmol
from molgr.utils.equivalence import check_equivalence
from molgr.utils.post_process import make_stereochemistry


RDLogger.DisableLog("rdApp.*")  # type: ignore


_SEED_LABEL_NORMALIZATION_XYZ = """15

O          1.09630       -0.87290       -0.00190
N          0.73010       -2.07620       -0.00290
C          0.68990       -2.76820       -1.34610
C         -0.16110       -4.04780       -1.24720
C          0.13080       -4.88630       -0.00600
C         -0.15820       -4.05110        1.23820
C          0.69300       -2.77180        1.33850
C          0.08290       -1.81050       -2.38120
C          2.15530       -3.06860       -1.73930
C          2.15940       -3.07320        1.72740
C          0.08850       -1.81700        2.37780
H          0.67210       -0.88080       -2.48180
H         -0.95860       -1.55790       -2.11210
H          0.07560       -2.32530       -3.36140
H          2.16980       -3.47910       -2.76790
"""


def _total_charge_and_radicals(mol: Chem.Mol) -> tuple[int, int]:
    charge = 0
    radicals = 0
    for atom in mol.GetAtoms():  # pyright: ignore[reportCallIssue]
        charge += int(atom.GetFormalCharge())
        radicals += int(atom.GetNumRadicalElectrons())
    return charge, radicals


def _get_ptr(obmol: ob.OBMol) -> int:
    this = getattr(obmol, "this", None)
    if this is not None:
        return int(this)  # type: ignore[arg-type]
    return int(obmol)  # type: ignore[arg-type]


def _pybel_stage_signature(
    omol: pybel.Molecule,
) -> tuple[tuple[tuple[int, int, int], ...], tuple[tuple[int, int, int], ...]]:
    obmol = omol.OBMol
    atoms = tuple(
        (atom.idx, atom.OBAtom.GetFormalCharge(), atom.OBAtom.GetSpinMultiplicity())
        for atom in omol
    )
    bonds = tuple(
        sorted(
            (
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                bond.GetBondOrder(),
            )
            for bond in ob.OBMolBondIter(obmol)
        )
    )
    return atoms, bonds


@typed_lru_cache(maxsize=128, typed=True)
def _seed_case(smiles: str) -> tuple[str, int, int]:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    mol_h = Chem.AddHs(mol)
    embed_code = rdDistGeom.EmbedMolecule(mol_h)  # pyright: ignore[reportCallIssue]
    assert int(embed_code) == 0
    charge, radicals = _total_charge_and_radicals(mol_h)
    return Chem.MolToXYZBlock(mol_h), charge, radicals


def _load_curated_smiles() -> list[str]:
    csv_path = Path(__file__).with_name("test_cases.csv")
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    curated_rows = [1, 2, 5, 10, 17]
    return [rows[idx - 1]["smiles"] for idx in curated_rows]


@pytest.mark.parametrize("smiles", _load_curated_smiles())
def test_fallback_no_metal_reconstructs_curated_cases(smiles: str) -> None:
    xyz_block, total_charge, total_radical_electrons = _seed_case(smiles)

    result = xyz_to_omol_no_metal(xyz_block, total_charge, total_radical_electrons)

    assert result is not None

    expected = Chem.MolFromSmiles(smiles)
    assert expected is not None
    equivalent, info = check_equivalence(
        expected,
        make_stereochemistry(pybel_to_rdmol(result)),
        use_chirality=True,
        max_resonance=100,
    )
    assert equivalent, info.reason


def test_fallback_no_metal_exposes_staged_history() -> None:
    xyz_block = """2
H2
H 0.0 0.0 0.0
H 0.0 0.0 0.74
"""
    state = xyz_to_omol_no_metal_state(xyz_block, 0, 0)

    assert state is not None
    assert state.phase_history[:7] == (
        "read_xyz",
        "normalize_seed_electronic_labels",
        "make_connections",
        "pre_clean",
        "fresh_omol_charge_radical_initial",
        "initialize_charge_budget",
        "eliminate_NNN_negative",
    )
    assert "break_one_bond" in state.phase_history
    assert state.phase_history[-1] in {"clean_resonances", "select_best_resonance_candidate"}


def test_seed_state_normalizes_openbabel_charge_and_spin_labels_for_python_and_cpp() -> None:
    raw_omol = pybel.readstring("xyz", _SEED_LABEL_NORMALIZATION_XYZ)
    raw_nitrogen = raw_omol.OBMol.GetAtom(2)
    raw_oxygen = raw_omol.OBMol.GetAtom(1)
    assert raw_nitrogen.GetFormalCharge() == 1
    assert raw_oxygen.GetFormalCharge() == -1

    state = no_metal_preparation_module._seed_state(_SEED_LABEL_NORMALIZATION_XYZ, 0, 0)
    assert state.phase_history[:2] == ("read_xyz", "normalize_seed_electronic_labels")
    assert all(atom.OBAtom.GetFormalCharge() == 0 for atom in state.omol)
    assert all(atom.OBAtom.GetSpinMultiplicity() == 0 for atom in state.omol)

    from molgr import _core  # type: ignore

    cpp_trace = _core.dev.pipeline.reconstruct_without_metals.debug_linear_pipeline_trace(
        _SEED_LABEL_NORMALIZATION_XYZ,
        0,
        0,
    )
    assert cpp_trace is not None
    assert cpp_trace[0]["phase"] == "read_xyz"
    assert "[N+]" not in cpp_trace[0]["smiles"]
    assert "[O-]" not in cpp_trace[0]["smiles"]


def test_make_connections_reconnects_current_donor_for_python_and_cpp() -> None:
    def make_omol() -> pybel.Molecule:
        obmol = ob.OBMol()
        obmol.BeginModify()

        atom_specs = (
            (7, (0.00, 0.00, 0.00)),
            (1, (-0.60, 0.94, 0.00)),
            (1, (-0.60, -0.47, 0.82)),
            (1, (-0.60, -0.47, -0.82)),
            (5, (1.65, 0.00, 0.00)),
            (1, (2.45, 0.94, 0.00)),
            (1, (2.45, -0.47, 0.82)),
            (1, (2.45, -0.47, -0.82)),
        )
        for atomic_num, (x, y, z) in atom_specs:
            atom = obmol.NewAtom()
            atom.SetAtomicNum(atomic_num)
            atom.SetVector(x, y, z)

        for begin_idx, end_idx in ((1, 2), (1, 3), (1, 4), (5, 6), (5, 7), (5, 8)):
            obmol.AddBond(begin_idx, end_idx, 1)
        obmol.EndModify()
        return pybel.Molecule(obmol)

    omol = make_omol()
    assert omol.OBMol.GetBond(1, 5) is None

    omol, hit = make_connections(omol)

    assert hit
    assert omol.OBMol.GetBond(1, 5) is not None
    assert omol.OBMol.GetBond(1, 5).GetBondOrder() == 1

    from molgr import _core  # type: ignore

    cpp_omol = make_omol()
    assert cpp_omol.OBMol.GetBond(1, 5) is None

    cpp_hit = _core.dev.stages.preprocess.make_connections_ptr(_get_ptr(cpp_omol.OBMol), 0.15)

    assert cpp_hit
    assert _pybel_stage_signature(cpp_omol) == _pybel_stage_signature(omol)


def test_eliminate_nnn_negative_produces_closed_shell_azide() -> None:
    omol = pybel.readstring("smi", "[N][N][N]")

    omol, hit = fresh_omol_charge_radical(omol)

    assert hit
    assert [
        (atom.idx, atom.OBAtom.GetFormalCharge(), atom.OBAtom.GetSpinMultiplicity())
        for atom in omol
    ] == [
        (1, 0, 2),
        (2, 0, 1),
        (3, 0, 2),
    ]

    omol, given_charge, hit = eliminate_NNN(omol, 0, False)

    assert hit
    assert given_charge == 1
    assert [
        (atom.idx, atom.OBAtom.GetFormalCharge(), atom.OBAtom.GetSpinMultiplicity())
        for atom in omol
    ] == [
        (1, -1, 0),
        (2, 1, 0),
        (3, -1, 0),
    ]


def test_eliminate_cn_in_doubt_requires_disjoint_pairs_for_python_and_cpp() -> None:
    smiles = "C[N+](C)=C=[N+](C)C"

    omol = pybel.readstring("smi", smiles)
    before = _pybel_stage_signature(omol)

    omol, given_charge, hit = eliminate_CN_in_doubt(omol, 0)

    assert not hit
    assert given_charge == 0
    assert _pybel_stage_signature(omol) == before

    from molgr import _core  # type: ignore

    cpp_omol = pybel.readstring("smi", smiles)
    before = _pybel_stage_signature(cpp_omol)

    given_charge, hit = _core.dev.stages.eliminate.eliminate_cn_in_doubt_ptr(
        _get_ptr(cpp_omol.OBMol),
        0,
    )

    assert not hit
    assert given_charge == 0
    assert _pybel_stage_signature(cpp_omol) == before


def test_eliminate_negative_charges_converts_cyclopentadienyl_radical_for_python_and_cpp() -> None:
    def make_omol() -> pybel.Molecule:
        obmol = ob.OBMol()
        obmol.BeginModify()
        for _ in range(5):
            atom = obmol.NewAtom()
            atom.SetAtomicNum(6)
            atom.SetFormalCharge(0)
            atom.SetSpinMultiplicity(0)
        obmol.GetAtom(5).SetSpinMultiplicity(1)
        obmol.AddBond(1, 2, 2)
        obmol.AddBond(2, 3, 1)
        obmol.AddBond(3, 4, 2)
        obmol.AddBond(4, 5, 1)
        obmol.AddBond(5, 1, 1)
        obmol.EndModify()
        return pybel.Molecule(obmol)

    omol = make_omol()
    omol, given_charge, hit = eliminate_negative_charges(omol, 0)

    assert hit
    assert given_charge == 1
    assert _pybel_stage_signature(omol) == (
        ((1, 0, 0), (2, 0, 0), (3, 0, 0), (4, 0, 0), (5, -1, 0)),
        ((1, 2, 2), (2, 3, 1), (3, 4, 2), (4, 5, 1), (5, 1, 1)),
    )

    from molgr import _core  # type: ignore

    cpp_omol = make_omol()
    cpp_given_charge, cpp_hit = _core.dev.stages.eliminate.eliminate_negative_charges_ptr(
        _get_ptr(cpp_omol.OBMol),
        0,
    )
    if not cpp_hit:
        pytest.skip("loaded C++ extension does not include ELIM_NEGATIVE_CP; rebuild _core")

    assert cpp_hit
    assert cpp_given_charge == given_charge
    assert _pybel_stage_signature(cpp_omol) == _pybel_stage_signature(omol)


def test_eliminate_negative_charges_prefers_single_radical_heteroatom_even_at_zero_charge_for_python_and_cpp() -> None:
    def make_omol() -> pybel.Molecule:
        obmol = ob.OBMol()
        obmol.BeginModify()
        atom = obmol.NewAtom()
        atom.SetAtomicNum(8)
        atom.SetFormalCharge(0)
        atom.SetSpinMultiplicity(1)
        obmol.EndModify()
        return pybel.Molecule(obmol)

    omol = make_omol()
    omol, given_charge, hit = eliminate_negative_charges(omol, 0)

    assert hit
    assert given_charge == 1
    assert _pybel_stage_signature(omol) == (((1, -1, 0),), ())

    from molgr import _core  # type: ignore

    cpp_omol = make_omol()
    cpp_given_charge, cpp_hit = _core.dev.stages.eliminate.eliminate_negative_charges_ptr(
        _get_ptr(cpp_omol.OBMol),
        0,
    )

    assert cpp_hit
    assert cpp_given_charge == given_charge
    assert _pybel_stage_signature(cpp_omol) == _pybel_stage_signature(omol)


def test_eliminate_high_positive_charge_atoms_skips_overstabilized_match_and_keeps_later_match_for_python_and_cpp() -> None:
    def make_omol() -> pybel.Molecule:
        obmol = ob.OBMol()
        obmol.BeginModify()
        a1 = obmol.NewAtom()
        a1.SetAtomicNum(15)
        a1.SetFormalCharge(1)
        a2 = obmol.NewAtom()
        a2.SetAtomicNum(8)
        a2.SetFormalCharge(0)
        a2.SetSpinMultiplicity(1)
        a3 = obmol.NewAtom()
        a3.SetAtomicNum(8)
        a3.SetFormalCharge(-1)
        a4 = obmol.NewAtom()
        a4.SetAtomicNum(6)
        a4.SetFormalCharge(-1)
        a5 = obmol.NewAtom()
        a5.SetAtomicNum(15)
        a5.SetFormalCharge(1)
        a6 = obmol.NewAtom()
        a6.SetAtomicNum(8)
        a6.SetFormalCharge(0)
        a6.SetSpinMultiplicity(1)
        a7 = obmol.NewAtom()
        a7.SetAtomicNum(6)
        a7.SetFormalCharge(0)
        obmol.AddBond(1, 2, 1)
        obmol.AddBond(1, 3, 1)
        obmol.AddBond(1, 4, 1)
        obmol.AddBond(5, 6, 1)
        obmol.AddBond(5, 7, 1)
        obmol.EndModify()
        return pybel.Molecule(obmol)

    omol = make_omol()
    omol, given_charge, hit = eliminate_high_positive_charge_atoms(omol, 0)

    assert hit
    assert given_charge == 1
    assert _pybel_stage_signature(omol) == (
        (
            (1, 1, 0),
            (2, 0, 1),
            (3, -1, 0),
            (4, -1, 0),
            (5, 1, 0),
            (6, -1, 0),
            (7, 0, 0),
        ),
        ((1, 2, 1), (1, 3, 1), (1, 4, 1), (5, 6, 1), (5, 7, 1)),
    )

    from molgr import _core  # type: ignore

    cpp_omol = make_omol()
    cpp_given_charge, cpp_hit = _core.dev.stages.eliminate.eliminate_high_positive_charge_atoms_ptr(
        _get_ptr(cpp_omol.OBMol),
        0,
    )

    assert cpp_hit
    assert cpp_given_charge == given_charge
    assert _pybel_stage_signature(cpp_omol) == _pybel_stage_signature(omol)


def test_eliminate_high_positive_charge_atoms_skips_non_singlet_donor_for_python_and_cpp() -> None:
    def make_omol() -> pybel.Molecule:
        obmol = ob.OBMol()
        obmol.BeginModify()
        a1 = obmol.NewAtom()
        a1.SetAtomicNum(15)
        a1.SetFormalCharge(1)
        a2 = obmol.NewAtom()
        a2.SetAtomicNum(8)
        a2.SetFormalCharge(0)
        a2.SetSpinMultiplicity(0)
        a3 = obmol.NewAtom()
        a3.SetAtomicNum(6)
        a3.SetFormalCharge(0)
        obmol.AddBond(1, 2, 1)
        obmol.AddBond(1, 3, 1)
        obmol.EndModify()
        return pybel.Molecule(obmol)

    omol = make_omol()
    omol, given_charge, hit = eliminate_high_positive_charge_atoms(omol, 0)

    assert not hit
    assert given_charge == 0
    assert _pybel_stage_signature(omol) == (
        ((1, 1, 0), (2, 0, 0), (3, 0, 0)),
        ((1, 2, 1), (1, 3, 1)),
    )

    from molgr import _core  # type: ignore

    cpp_omol = make_omol()
    cpp_given_charge, cpp_hit = _core.dev.stages.eliminate.eliminate_high_positive_charge_atoms_ptr(
        _get_ptr(cpp_omol.OBMol),
        0,
    )

    assert not cpp_hit
    assert cpp_given_charge == given_charge
    assert _pybel_stage_signature(cpp_omol) == _pybel_stage_signature(omol)


def test_clean_resonances_14_matches_python_and_cpp() -> None:
    def make_omol() -> pybel.Molecule:
        obmol = ob.OBMol()
        obmol.BeginModify()
        phosphorus = obmol.NewAtom()
        phosphorus.SetAtomicNum(15)
        phosphorus.SetFormalCharge(-1)
        nitrogen = obmol.NewAtom()
        nitrogen.SetAtomicNum(7)
        nitrogen.SetFormalCharge(1)
        obmol.AddBond(1, 2, 3)
        obmol.EndModify()
        return pybel.Molecule(obmol)

    omol = make_omol()
    omol, hit = clean_resonances(omol)

    assert hit
    assert _pybel_stage_signature(omol) == (
        ((1, 0, 0), (2, 0, 0)),
        ((1, 2, 2),),
    )

    from molgr import _core  # type: ignore

    cpp_omol = make_omol()
    cpp_hit = _core.dev.stages.clean.clean_resonances_ptr(_get_ptr(cpp_omol.OBMol))

    assert cpp_hit
    assert _pybel_stage_signature(cpp_omol) == _pybel_stage_signature(omol)


def test_clean_resonances_15_matches_python_and_cpp() -> None:
    def make_omol() -> pybel.Molecule:
        obmol = ob.OBMol()
        obmol.BeginModify()
        carbon = obmol.NewAtom()
        carbon.SetAtomicNum(6)
        carbon.SetFormalCharge(0)
        carbon.SetSpinMultiplicity(2)
        oxygen = obmol.NewAtom()
        oxygen.SetAtomicNum(8)
        oxygen.SetFormalCharge(0)
        oxygen.SetSpinMultiplicity(0)
        obmol.AddBond(1, 2, 2)
        obmol.EndModify()
        return pybel.Molecule(obmol)

    omol = make_omol()
    omol, hit = clean_resonances(omol)

    assert hit
    assert _pybel_stage_signature(omol) == (
        ((1, -1, 0), (2, 1, 0)),
        ((1, 2, 3),),
    )

    from molgr import _core  # type: ignore

    cpp_omol = make_omol()
    cpp_hit = _core.dev.stages.clean.clean_resonances_ptr(_get_ptr(cpp_omol.OBMol))

    assert cpp_hit
    assert _pybel_stage_signature(cpp_omol) == _pybel_stage_signature(omol)


def test_clean_resonances_16_matches_python_and_cpp() -> None:
    def make_omol() -> pybel.Molecule:
        obmol = ob.OBMol()
        obmol.BeginModify()
        atom1 = obmol.NewAtom()
        atom1.SetAtomicNum(8)
        atom1.SetFormalCharge(-1)
        atom1.SetSpinMultiplicity(0)
        atom2 = obmol.NewAtom()
        atom2.SetAtomicNum(6)
        atom2.SetFormalCharge(0)
        atom2.SetSpinMultiplicity(0)
        atom3 = obmol.NewAtom()
        atom3.SetAtomicNum(6)
        atom3.SetFormalCharge(0)
        atom3.SetSpinMultiplicity(0)
        atom4 = obmol.NewAtom()
        atom4.SetAtomicNum(6)
        atom4.SetFormalCharge(0)
        atom4.SetSpinMultiplicity(0)
        atom5 = obmol.NewAtom()
        atom5.SetAtomicNum(7)
        atom5.SetFormalCharge(1)
        atom5.SetSpinMultiplicity(0)
        obmol.AddBond(1, 2, 1)
        obmol.AddBond(2, 3, 2)
        obmol.AddBond(3, 4, 1)
        obmol.AddBond(4, 5, 2)
        obmol.EndModify()
        return pybel.Molecule(obmol)

    omol = make_omol()
    omol, hit = clean_resonances(omol)

    assert hit
    assert _pybel_stage_signature(omol) == (
        ((1, 0, 0), (2, 0, 1), (3, 0, 1), (4, 0, 1), (5, 0, 2)),
        ((1, 2, 2), (2, 3, 1), (3, 4, 2), (4, 5, 1)),
    )

    from molgr import _core  # type: ignore

    cpp_omol = make_omol()
    cpp_hit = _core.dev.stages.clean.clean_resonances_ptr(_get_ptr(cpp_omol.OBMol))

    assert cpp_hit
    assert _pybel_stage_signature(cpp_omol) == _pybel_stage_signature(omol)


def test_clean_neighbor_radicals_refreshes_atom_state_for_python_and_cpp() -> None:
    def make_omol() -> pybel.Molecule:
        obmol = ob.OBMol()
        obmol.BeginModify()
        for _ in range(2):
            atom = obmol.NewAtom()
            atom.SetAtomicNum(5)
            atom.SetFormalCharge(-1)
            atom.SetSpinMultiplicity(1)
        obmol.AddBond(1, 2, 1)
        obmol.EndModify()
        return pybel.Molecule(obmol)

    omol = make_omol()
    omol, hit = clean_neighbor_radicals(omol)

    assert hit
    assert _pybel_stage_signature(omol) == (
        ((1, -1, 2), (2, -1, 2)),
        ((1, 2, 2),),
    )

    from molgr import _core  # type: ignore

    cpp_omol = make_omol()
    cpp_hit = _core.dev.stages.clean.clean_neighbor_radicals_ptr(_get_ptr(cpp_omol.OBMol))

    assert cpp_hit
    assert _pybel_stage_signature(cpp_omol) == _pybel_stage_signature(omol)


def test_eliminate_carbene_neighbor_heteroatom_skips_radical_neighbors_for_python_and_cpp() -> None:
    def make_omol() -> pybel.Molecule:
        obmol = ob.OBMol()
        obmol.BeginModify()
        blocked_carbon = obmol.NewAtom()
        blocked_carbon.SetAtomicNum(6)
        blocked_carbon.SetFormalCharge(0)
        blocked_carbon.SetSpinMultiplicity(2)
        radical_neighbor = obmol.NewAtom()
        radical_neighbor.SetAtomicNum(6)
        radical_neighbor.SetFormalCharge(0)
        radical_neighbor.SetSpinMultiplicity(1)
        blocked_oxygen = obmol.NewAtom()
        blocked_oxygen.SetAtomicNum(8)
        blocked_oxygen.SetFormalCharge(0)
        blocked_oxygen.SetSpinMultiplicity(0)
        carbonyl_carbon = obmol.NewAtom()
        carbonyl_carbon.SetAtomicNum(6)
        carbonyl_carbon.SetFormalCharge(0)
        carbonyl_carbon.SetSpinMultiplicity(2)
        carbonyl_oxygen = obmol.NewAtom()
        carbonyl_oxygen.SetAtomicNum(8)
        carbonyl_oxygen.SetFormalCharge(0)
        carbonyl_oxygen.SetSpinMultiplicity(0)
        obmol.AddBond(1, 2, 1)
        obmol.AddBond(1, 3, 1)
        obmol.AddBond(4, 5, 2)
        obmol.EndModify()
        return pybel.Molecule(obmol)

    omol = make_omol()
    omol, given_charge, hit = eliminate_carbene_neighbor_heteroatom(omol, 0)

    assert hit
    assert given_charge == 0
    assert _pybel_stage_signature(omol) == (
        ((1, 0, 2), (2, 0, 1), (3, 0, 0), (4, -1, 0), (5, 1, 0)),
        ((1, 2, 1), (1, 3, 1), (4, 5, 3)),
    )

    from molgr import _core  # type: ignore

    cpp_omol = make_omol()
    given_charge, hit = _core.dev.stages.eliminate.eliminate_carbene_neighbor_heteroatom_ptr(
        _get_ptr(cpp_omol.OBMol),
        0,
    )

    assert hit
    assert given_charge == 0
    assert _pybel_stage_signature(cpp_omol) == _pybel_stage_signature(omol)


@pytest.mark.parametrize("atomic_num", [8, 16, 7])
def test_eliminate_charge_spliting_ignores_spin2_radical_candidates_for_python_and_cpp(
    atomic_num: int,
) -> None:
    def make_omol() -> pybel.Molecule:
        obmol = ob.OBMol()
        obmol.BeginModify()
        target = obmol.NewAtom()
        target.SetAtomicNum(atomic_num)
        target.SetFormalCharge(0)
        target.SetSpinMultiplicity(2)
        carbon = obmol.NewAtom()
        carbon.SetAtomicNum(6)
        carbon.SetFormalCharge(0)
        carbon.SetSpinMultiplicity(1)
        obmol.EndModify()
        return pybel.Molecule(obmol)

    omol = make_omol()

    omol, given_charge, hit = eliminate_charge_spliting(omol, 0)

    assert hit
    assert given_charge == 1
    assert [
        (atom.atomicnum, atom.OBAtom.GetFormalCharge(), atom.OBAtom.GetSpinMultiplicity())
        for atom in omol
    ] == [(atomic_num, 0, 2), (6, -1, 0)]

    from molgr import _core  # type: ignore

    cpp_omol = make_omol()

    given_charge, hit = _core.dev.stages.eliminate.eliminate_charge_spliting_ptr(
        _get_ptr(cpp_omol.OBMol),
        0,
    )

    assert hit
    assert given_charge == 1
    assert [
        (atom.atomicnum, atom.OBAtom.GetFormalCharge(), atom.OBAtom.GetSpinMultiplicity())
        for atom in cpp_omol
    ] == [(atomic_num, 0, 2), (6, -1, 0)]


def test_run_linear_pipeline_passes_current_charge_into_break_stages(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    xyz_block = """2
H2
H 0.0 0.0 0.0
H 0.0 0.0 0.74
"""
    state = ReconstructionState(
        omol=pybel.readstring("xyz", xyz_block),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=3,
        phase_history=("read_xyz",),
    )
    recorded: dict[str, tuple[int, ...]] = {}

    monkeypatch.setattr(no_metal_preparation_module, "make_connections", lambda omol: (omol, False))
    monkeypatch.setattr(no_metal_preparation_module, "pre_clean", lambda omol: (omol, False))
    monkeypatch.setattr(
        no_metal_preparation_module,
        "fresh_omol_charge_radical",
        lambda omol: (omol, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "eliminate_NNN",
        lambda omol, given_charge, positive: (omol, 7 if not positive else given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "eliminate_high_positive_charge_atoms",
        lambda omol, given_charge: (omol, given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "eliminate_CN_in_doubt",
        lambda omol, given_charge: (omol, given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "eliminate_carboxyl",
        lambda omol, given_charge: (omol, given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "eliminate_carbene_neighbor_heteroatom",
        lambda omol, given_charge: (omol, given_charge, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "clean_carbene_neighbor_unsaturated",
        lambda omol: (omol, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "clean_neighbor_radicals",
        lambda omol: (omol, False),
    )
    monkeypatch.setattr(
        no_metal_preparation_module,
        "eliminate_charge_spliting",
        lambda omol, given_charge: (omol, given_charge, False),
    )

    def record_break_deformed_ene(omol, given_charge, total_radical_electrons, tolerance):
        recorded["break_deformed_ene"] = (given_charge, total_radical_electrons, tolerance)
        return omol, False

    def record_break_one_bond(omol, given_charge, total_radical_electrons):
        recorded["break_one_bond"] = (given_charge, total_radical_electrons)
        return omol, given_charge, False

    monkeypatch.setattr(
        no_metal_preparation_module, "break_deformed_ene", record_break_deformed_ene
    )
    monkeypatch.setattr(no_metal_preparation_module, "break_one_bond", record_break_one_bond)

    next_state = no_metal_preparation_module._run_linear_pipeline(state)

    assert recorded["break_deformed_ene"] == (7, 3, 5.0)
    assert recorded["break_one_bond"] == (7, 3)
    assert next_state.given_charge == 7


def test_fallback_no_metal_reuses_resonance_score_metadata(monkeypatch: pytest.MonkeyPatch) -> None:
    no_metal_module._run_no_metal_pipeline_cached.cache_clear()

    base_state = ReconstructionState(
        omol=pybel.readstring("smi", "CC"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz",),
    )
    weaker_candidate = ReconstructionState(
        omol=pybel.readstring("smi", "C=CC=C"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz", "validate_resonance_candidate"),
        metadata={"score": 5.0},
    )
    stronger_candidate = ReconstructionState(
        omol=pybel.readstring("smi", "C=CC=C"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz", "validate_resonance_candidate"),
        metadata={"score": 2.0},
    )

    monkeypatch.setattr(
        no_metal_preparation_module, "_seed_state", lambda *args, **kwargs: base_state
    )
    monkeypatch.setattr(no_metal_preparation_module, "_run_linear_pipeline", lambda state: state)
    monkeypatch.setattr(
        no_metal_preparation_module,
        "validate_omol",
        lambda *args, **kwargs: False,
    )
    monkeypatch.setattr(
        no_metal_resonance_module,
        "_recover_resonance_candidates",
        lambda state, **kwargs: [weaker_candidate, stronger_candidate],
    )
    monkeypatch.setattr(
        no_metal_selection_module,
        "_score_reconstruction_candidate",
        lambda candidate, **kwargs: float(candidate.metadata["score"]),
    )

    result = xyz_to_omol_no_metal_state("cached-resonance", 0, 0)

    assert result is not None
    assert result.metadata["score"] == 2.0
    assert result.phase_history[-1] == "select_best_resonance_candidate"


def test_no_metal_resonance_selection_prefers_aromatic_topology_before_force_field(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    no_metal_module._run_no_metal_pipeline_cached.cache_clear()

    base_state = ReconstructionState(
        omol=pybel.readstring("smi", "CC"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz",),
    )
    lower_force_field_candidate = ReconstructionState(
        omol=pybel.readstring("smi", "C=CC=C"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz", "validate_resonance_candidate"),
        metadata={"score": 1.0},
    )
    aromatic_candidate = ReconstructionState(
        omol=pybel.readstring("smi", "c1ccccc1"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz", "validate_resonance_candidate"),
        metadata={"score": 10.0},
    )

    monkeypatch.setattr(
        no_metal_preparation_module, "_seed_state", lambda *args, **kwargs: base_state
    )
    monkeypatch.setattr(no_metal_preparation_module, "_run_linear_pipeline", lambda state: state)
    monkeypatch.setattr(
        no_metal_preparation_module,
        "validate_omol",
        lambda *args, **kwargs: False,
    )
    monkeypatch.setattr(
        no_metal_resonance_module,
        "_recover_resonance_candidates",
        lambda state, **kwargs: [lower_force_field_candidate, aromatic_candidate],
    )
    monkeypatch.setattr(
        no_metal_selection_module,
        "_score_reconstruction_candidate",
        lambda candidate, **kwargs: float(candidate.metadata["score"]),
    )

    result = xyz_to_omol_no_metal_state("topology-first", 0, 0)

    assert result is not None
    assert result.omol.write("can").strip() == aromatic_candidate.omol.write("can").strip()
    assert result.metadata["organic_aromatic_atom_count"] == 6
    assert result.metadata["organic_aromatic_stability_score"] == pytest.approx(1.0)
    assert result.metadata["organic_topology_selection_key"][:5] == (-1.0, -6, -6, -6, -6)
    assert result.phase_history[-1] == "select_best_resonance_candidate"


def test_no_metal_resonance_selection_penalizes_charge_separated_extra_conjugation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    no_metal_module._run_no_metal_pipeline_cached.cache_clear()

    base_state = ReconstructionState(
        omol=pybel.readstring("smi", "CC"),
        given_charge=0,
        total_charge=-2,
        total_radical_electrons=0,
        phase_history=("read_xyz",),
    )
    lower_force_field_candidate = ReconstructionState(
        omol=pybel.readstring(
            "smi",
            "c1ccc(nc1)c1ccccn1.Oc1cc(O)c2c(c1)oc(c(c2=O)[O-])c1ccc(c(c1)O)O.[Cl-]",
        ),
        given_charge=0,
        total_charge=-2,
        total_radical_electrons=0,
        phase_history=("read_xyz", "validate_resonance_candidate"),
        metadata={"score": 1168.9043528762083},
    )
    charge_separated_candidate = ReconstructionState(
        omol=pybel.readstring(
            "smi",
            "c1ccc(nc1)c1ccccn1.[OH+]=c1cc2-c(c(c1)O)c(c(c(o2)c1ccc(c(c1)O)O)[O-])[O-].[Cl-]",
        ),
        given_charge=0,
        total_charge=-2,
        total_radical_electrons=0,
        phase_history=("read_xyz", "validate_resonance_candidate"),
        metadata={"score": 1226.1886293507894},
    )

    monkeypatch.setattr(
        no_metal_preparation_module, "_seed_state", lambda *args, **kwargs: base_state
    )
    monkeypatch.setattr(no_metal_preparation_module, "_run_linear_pipeline", lambda state: state)
    monkeypatch.setattr(
        no_metal_preparation_module,
        "validate_omol",
        lambda *args, **kwargs: False,
    )
    monkeypatch.setattr(
        no_metal_resonance_module,
        "_recover_resonance_candidates",
        lambda state, **kwargs: [lower_force_field_candidate, charge_separated_candidate],
    )
    monkeypatch.setattr(
        no_metal_selection_module,
        "_score_reconstruction_candidate",
        lambda candidate, **kwargs: float(candidate.metadata["score"]),
    )

    result = xyz_to_omol_no_metal_state("charge-penalized", -2, 0)

    assert result is not None
    assert result.metadata["score"] == pytest.approx(1168.9043528762083)
    assert result.metadata["organic_formal_charge_absolute_sum"] == 2
    assert result.metadata["organic_conjugation_charge_penalty"] == pytest.approx(1.0)
    assert result.omol.write("can").strip() == lower_force_field_candidate.omol.write("can").strip()
    assert (
        lower_force_field_candidate.metadata["organic_adjusted_max_conjugated_component_size"]
        > charge_separated_candidate.metadata["organic_adjusted_max_conjugated_component_size"]
    )
    assert (
        lower_force_field_candidate.metadata["organic_adjusted_conjugated_atom_count"]
        > charge_separated_candidate.metadata["organic_adjusted_conjugated_atom_count"]
    )
    assert (
        lower_force_field_candidate.metadata["organic_adjusted_conjugated_bond_count"]
        > charge_separated_candidate.metadata["organic_adjusted_conjugated_bond_count"]
    )
    assert (
        lower_force_field_candidate.metadata["organic_topology_selection_key"][:5]
        < charge_separated_candidate.metadata["organic_topology_selection_key"][:5]
    )
    assert result.phase_history[-1] == "select_best_resonance_candidate"


def test_no_metal_cached_pipeline_uses_config_in_key_and_execution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    no_metal_module._run_no_metal_pipeline_cached.cache_clear()

    base_state = ReconstructionState(
        omol=pybel.readstring("smi", "CC"),
        given_charge=0,
        total_charge=0,
        total_radical_electrons=0,
        phase_history=("read_xyz",),
    )
    default_config = make_default_config()
    config_a = replace(
        default_config,
        resonance=replace(default_config.resonance, max_depth=2),
    )
    config_b = replace(
        default_config,
        resonance=replace(default_config.resonance, max_depth=3),
    )
    seen_configs = []

    monkeypatch.setattr(
        no_metal_preparation_module, "_seed_state", lambda *args, **kwargs: base_state
    )

    def fake_run_from_state(seed_state: ReconstructionState, *, config=None):
        assert seed_state is base_state
        seen_configs.append(config)
        return base_state

    monkeypatch.setattr(no_metal_module, "_run_no_metal_pipeline_from_state", fake_run_from_state)

    first = xyz_to_omol_no_metal_state("cfg-key", 0, 0, config=config_a)
    second = xyz_to_omol_no_metal_state("cfg-key", 0, 0, config=config_a)
    third = xyz_to_omol_no_metal_state("cfg-key", 0, 0, config=config_b)

    assert first is base_state
    assert second is base_state
    assert third is base_state
    assert seen_configs == [config_a, config_b]
