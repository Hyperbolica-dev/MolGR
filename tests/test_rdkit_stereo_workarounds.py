from __future__ import annotations

from itertools import combinations, permutations

import numpy as np
import pytest
from rdkit import Chem

from molgr.utils import rdkit_stereo
from molgr.utils.post_process import make_stereochemistry


_AXES = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
_TETRAHEDRON = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, -1, -1))
_SP = (_AXES[0], _AXES[2], _AXES[1], _AXES[3])
_TB = ((0, 0, 1), (1, 0, 0), (-0.5, np.sqrt(3) / 2, 0), (-0.5, -np.sqrt(3) / 2, 0), (0, 0, -1))


@pytest.fixture(autouse=True)
def stereo_settings():
    allow = Chem.GetAllowNontetrahedralChirality()
    legacy = Chem.GetUseLegacyStereoPerception()
    Chem.SetAllowNontetrahedralChirality(True)
    Chem.SetUseLegacyStereoPerception(True)
    yield
    Chem.SetAllowNontetrahedralChirality(allow)
    Chem.SetUseLegacyStereoPerception(legacy)


def _molecule(points, *, symbol="Rh", center=(0, 0, 0), order=None, dative=False):
    rw = Chem.RWMol()
    for element in (symbol, "F", "Cl", "Br", "I", "N", "O")[: len(points) + 1]:
        atom = Chem.Atom(element)
        atom.SetNoImplicit(True)
        rw.AddAtom(atom)
    for idx in range(1, len(points) + 1) if order is None else order:
        if dative:
            rw.AddBond(idx, 0, Chem.BondType.DATIVE)
        else:
            rw.AddBond(0, idx, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    if symbol == "C":
        Chem.SanitizeMol(mol)
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for idx, point in enumerate((center, *points)):
        conf.SetAtomPosition(idx, tuple(float(x) for x in point))
    mol.AddConformer(conf)
    return mol


def _native(mol):
    Chem.AssignAtomChiralTagsFromStructure(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    return mol


def _permutation(mol):
    atom = mol.GetAtomWithIdx(0)
    assert atom.GetChiralTag() == Chem.ChiralType.CHI_OCTAHEDRAL
    return atom.GetIntProp("_chiralPermutation")


@pytest.mark.parametrize(
    ("version", "enabled"),
    [
        ("2023.09.5", False),
        ("2023.9.6", True),
        ("2024.03.5", True),
        ("2026.09.1", True),
        ("2026.9.1", True),
        ("2026.09.2", False),
        ("2027.03.1", False),
        ("2026.09.1rc1", False),
        ("2026.09.1.dev1", False),
        ("unknown", False),
    ],
)
def test_version_bounds(version, enabled):
    assert rdkit_stereo._needs_stereo_workaround(version) is enabled


@pytest.mark.parametrize("version", ["2026.09.2", "2027.03.1", "unknown"])
def test_out_of_range_keeps_native_results(monkeypatch, version):
    monkeypatch.setattr(rdkit_stereo.rdBase, "rdkitVersion", version)
    for mol in (
        _molecule(_TETRAHEDRON, symbol="C", center=(3, 0, 0)),
        _molecule((_AXES[0], _AXES[2], _AXES[3], _AXES[4])),
    ):
        expected = _native(Chem.Mol(mol))
        actual = make_stereochemistry(mol)
        actual_atom, expected_atom = actual.GetAtomWithIdx(0), expected.GetAtomWithIdx(0)
        for prop in ("_CIPCode", "_chiralPermutation", "_NonExplicit3DChirality"):
            assert actual_atom.HasProp(prop) == expected_atom.HasProp(prop)
            if expected_atom.HasProp(prop):
                assert actual_atom.GetProp(prop) == expected_atom.GetProp(prop)
        assert actual_atom.GetChiralTag() == expected_atom.GetChiralTag()


@pytest.mark.parametrize("center", [(0, 0, 0), (3, 0, 0)])
@pytest.mark.parametrize("mirror", [False, True])
def test_tetrahedral_all_neighbor_orders_and_cip(center, mirror):
    points = np.array(_TETRAHEDRON, dtype=float)
    if mirror:
        points[:, 0] *= -1
    expected = _native(_molecule(points, symbol="C"))
    expected_smiles = Chem.MolToSmiles(expected)
    expected_cip = expected.GetAtomWithIdx(0).GetProp("_CIPCode")
    for order in permutations(range(1, 5)):
        mol = _molecule(points, symbol="C", center=center, order=order)
        make_stereochemistry(mol)
        assert Chem.MolToSmiles(mol) == expected_smiles
        assert mol.GetAtomWithIdx(0).GetProp("_CIPCode") == expected_cip
        make_stereochemistry(mol)
        assert Chem.MolToSmiles(mol) == expected_smiles
        assert mol.GetAtomWithIdx(0).GetProp("_CIPCode") == expected_cip


@pytest.mark.parametrize("scale", [0.01, 1, 100])
@pytest.mark.parametrize("height", [0.0, 1e-5])
def test_degenerate_ligands_clear_tags_independent_of_order_and_scale(scale, height):
    points = scale * np.array(((1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, height)))
    for order in permutations(range(1, 5)):
        mol = _molecule(points, symbol="C", center=(0, 0, scale), order=order)
        make_stereochemistry(mol)
        atom = mol.GetAtomWithIdx(0)
        assert atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
        assert not atom.HasProp("_CIPCode")


def test_tetrahedral_scale_rotation_translation_and_atom_renumbering():
    rotation = np.array(((0, -1, 0), (1, 0, 0), (0, 0, 1)))
    expected = Chem.MolToSmiles(_native(_molecule(_TETRAHEDRON, symbol="C")))
    for scale in (0.01, 1, 100):
        points = scale * np.array(_TETRAHEDRON) @ rotation + (4, 5, 6)
        center = scale * np.array((3, 0, 0)) @ rotation + (4, 5, 6)
        mol = _molecule(points, symbol="C", center=center)
        for order in ((0, 1, 2, 3, 4), (4, 3, 0, 2, 1)):
            permuted = Chem.RenumberAtoms(mol, order)
            make_stereochemistry(permuted)
            assert Chem.MolToSmiles(permuted) == expected


@pytest.mark.parametrize("scale", [0.01, 1.0, 100.0])
def test_tetrahedral_geometry_matches_determinant_near_cutoff_and_with_distortion(scale):
    rng = np.random.default_rng(9127)
    geometries = list(rng.normal(size=(128, 4, 3)))
    # On both sides of the dimensionless 1e-3 volume cutoff, including mirrors.
    geometries += [
        np.array(((1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, height)))
        for height in (-0.0041, -0.0039, 0.0, 0.0039, 0.0041)
    ]
    geometries.append(np.zeros((4, 3)))
    for points in geometries:
        points = scale * points
        diameter = max(np.linalg.norm(a - b) for a, b in combinations(points, 2))
        volume = np.linalg.det((points[:3] - points[3]) / diameter) if diameter else 0.0
        expected = Chem.ChiralType.CHI_UNSPECIFIED
        if volume > 1e-3:
            expected = Chem.ChiralType.CHI_TETRAHEDRAL_CCW
        elif volume < -1e-3:
            expected = Chem.ChiralType.CHI_TETRAHEDRAL_CW
        mol = _molecule(points, symbol="C", center=scale * rng.normal(size=3))
        make_stereochemistry(mol)
        assert mol.GetAtomWithIdx(0).GetChiralTag() == expected


def test_coordinate_updates_and_requested_conformer_are_observed():
    mol = _molecule(_TETRAHEDRON, symbol="C", center=(3, 0, 0))
    mirrored = Chem.Conformer(mol.GetConformer())
    mirrored.SetId(7)
    for idx in range(mol.GetNumAtoms()):
        point = mirrored.GetAtomPosition(idx)
        mirrored.SetAtomPosition(idx, (-point.x, point.y, point.z))
    mol.AddConformer(mirrored, assignId=False)
    for conf_id, expected_cip in ((0, "S"), (7, "R"), (0, "S")):
        Chem.AssignStereochemistryFrom3D(mol, confId=conf_id)
        rdkit_stereo.apply_3d_stereo_workarounds(mol, conf_id=conf_id)
        assert mol.GetAtomWithIdx(0).GetProp("_CIPCode") == expected_cip
    conf = mol.GetConformer(0)
    for idx in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(idx, mirrored.GetAtomPosition(idx))
    make_stereochemistry(mol)
    assert mol.GetAtomWithIdx(0).GetProp("_CIPCode") == "R"


@pytest.mark.parametrize("legacy", [False, True])
def test_duplicate_ligands_do_not_create_tetrahedral_stereo(legacy):
    Chem.SetUseLegacyStereoPerception(legacy)
    mol = _molecule(_TETRAHEDRON, symbol="C", center=(3, 0, 0))
    mol.GetAtomWithIdx(2).SetAtomicNum(9)
    make_stereochemistry(mol)
    assert mol.GetAtomWithIdx(0).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
    assert not mol.GetAtomWithIdx(0).HasProp("_CIPCode")


@pytest.mark.parametrize("isotopes", [(0, 0), (0, 2), (2, 2), (2, 3)])
@pytest.mark.parametrize("legacy", [False, True])
def test_terminal_hydrogens_respect_isotopic_stereo_and_coordinate_scale(isotopes, legacy):
    Chem.SetUseLegacyStereoPerception(legacy)

    def build(points, center, order=None):
        mol = _molecule(points, symbol="C", center=center, order=order)
        for idx, isotope in enumerate(isotopes, 1):
            mol.GetAtomWithIdx(idx).SetAtomicNum(1)
            mol.GetAtomWithIdx(idx).SetIsotope(isotope)
        mol.UpdatePropertyCache(strict=False)
        return mol

    expected = Chem.MolToSmiles(_native(build(_TETRAHEDRON, (0, 0, 0))))
    for scale in (0.01, 1.0):
        for order in permutations(range(1, 5)):
            mol = build(scale * np.array(_TETRAHEDRON), (3 * scale, 0, 0), order)
            make_stereochemistry(mol)
            assert Chem.MolToSmiles(mol) == expected
            # Only achiral centers are eligible for the shortcut. Older RDKit
            # versions may omit this marker on already-correct chiral atoms.
            if isotopes[0] == isotopes[1]:
                assert mol.GetAtomWithIdx(0).GetIntProp("_NonExplicit3DChirality") == 1


@pytest.mark.parametrize("metadata", ["permutation", "cip", "missing_3d_marker"])
def test_achiral_carbon_still_refreshes_stale_metadata(metadata):
    mol = _molecule(_TETRAHEDRON, symbol="C", center=(3, 0, 0))
    mol.GetAtomWithIdx(1).SetAtomicNum(1)
    mol.GetAtomWithIdx(2).SetAtomicNum(1)
    mol.UpdatePropertyCache(strict=False)
    _native(mol)
    atom = mol.GetAtomWithIdx(0)
    if metadata == "permutation":
        atom.SetIntProp("_chiralPermutation", 19)
    elif metadata == "cip":
        atom.SetProp("_CIPCode", "R")
    else:
        atom.ClearProp("_NonExplicit3DChirality")
    rdkit_stereo.apply_3d_stereo_workarounds(mol)
    assert atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
    assert not atom.HasProp("_chiralPermutation")
    assert not atom.HasProp("_CIPCode")
    assert atom.GetIntProp("_NonExplicit3DChirality") == 1


@pytest.mark.parametrize("dative", [False, True])
def test_all_288_sparse_octahedral_arrangements_match_completed_geometry(dative):
    count = 0
    for occupied in permutations(_AXES, 4):
        pairs = [
            (i, j) for i, j in combinations(range(4), 2) if np.dot(occupied[i], occupied[j]) == -1
        ]
        if len(pairs) != 1:
            continue
        vacant = tuple(p for p in _AXES if p not in occupied)
        expected = {
            _permutation(_native(_molecule(occupied + vacant, dative=dative))),
            _permutation(_native(_molecule(occupied + vacant[::-1], dative=dative))),
        }
        sparse = _molecule(occupied, dative=dative)
        make_stereochemistry(sparse)
        assert _permutation(sparse) in expected, occupied
        make_stereochemistry(sparse)
        assert _permutation(sparse) in expected
        count += 1
    assert count == 288


def test_distorted_octahedral_geometry_with_bond_reordering():
    occupied = np.array((_AXES[0], _AXES[2], _AXES[3], _AXES[4]), dtype=float)
    occupied += ((0.01, 0.02, 0), (0.01, 0, -0.02), (0.01, 0.01, 0), (-0.02, 0.01, 0))
    vacant = np.array((_AXES[1], _AXES[5]), dtype=float)
    for order in permutations(range(1, 5)):
        expected = {
            _permutation(
                _native(_molecule(np.concatenate((occupied, missing)), order=order + (5, 6)))
            )
            for missing in (vacant, vacant[::-1])
        }
        mol = _molecule(occupied, order=order)
        # Renumbering preserves the bond insertion order but changes atom IDs.
        for renumbering in ((0, 1, 2, 3, 4), (0, 4, 2, 1, 3)):
            permuted = Chem.RenumberAtoms(mol, renumbering)
            make_stereochemistry(permuted)
            assert _permutation(permuted) in expected


@pytest.mark.parametrize(
    "points",
    [
        (_AXES[0], _AXES[1], _AXES[2], _AXES[3]),  # square planar
        (_AXES[4], _AXES[5], _AXES[0], (-0.5, 0.8660254, 0)),  # TB seesaw
        (_AXES[4], _AXES[5], _AXES[0], (-0.5, 0.8660254, 0), (-0.5, -0.8660254, 0)),
        _AXES[:5],  # square pyramidal
        _AXES,  # complete octahedral
    ],
)
def test_other_coordination_geometries_keep_native_assignment(points):
    mol = _molecule(points)
    expected = _native(Chem.Mol(mol)).GetAtomWithIdx(0)
    make_stereochemistry(mol)
    atom = mol.GetAtomWithIdx(0)
    assert atom.GetChiralTag() == expected.GetChiralTag()
    assert atom.GetIntProp("_chiralPermutation") == expected.GetIntProp("_chiralPermutation")


def test_outgoing_dative_bond_is_not_a_tetrahedral_ligand():
    mol = Chem.RWMol(_molecule(_TETRAHEDRON, symbol="C", center=(3, 0, 0)))
    metal = mol.AddAtom(Chem.Atom("Zn"))
    mol.AddBond(0, metal, Chem.BondType.DATIVE)
    mol.GetConformer().SetAtomPosition(metal, (5, 0, 0))
    make_stereochemistry(mol)
    expected = _native(_molecule(_TETRAHEDRON, symbol="C"))
    assert mol.GetAtomWithIdx(0).GetChiralTag() == expected.GetAtomWithIdx(0).GetChiralTag()


@pytest.mark.parametrize(
    "kind", ["implicit_h", "unknown_bond", "unknown_property", "no_3d", "no_conf"]
)
def test_ineligible_inputs_are_not_repaired(kind):
    mol = _molecule(_TETRAHEDRON, symbol="C", center=(3, 0, 0))
    if kind == "implicit_h":
        rw = Chem.RWMol(mol)
        rw.RemoveAtom(4)
        mol = rw.GetMol()
        mol.GetAtomWithIdx(0).SetNoImplicit(False)
        mol.UpdatePropertyCache()
    elif kind == "unknown_bond":
        mol.GetBondWithIdx(0).SetBondDir(Chem.BondDir.UNKNOWN)
    elif kind == "unknown_property":
        mol.GetBondWithIdx(0).SetIntProp("_UnknownStereo", 1)
    elif kind == "no_3d":
        mol.GetConformer().Set3D(False)
    else:
        mol.RemoveAllConformers()
    _native(mol)
    tag = mol.GetAtomWithIdx(0).GetChiralTag()
    rdkit_stereo.apply_3d_stereo_workarounds(mol)
    assert mol.GetAtomWithIdx(0).GetChiralTag() == tag


def test_refresh_preserves_double_bond_stereo():
    from rdkit.Chem import rdDistGeom

    alkene = Chem.AddHs(Chem.MolFromSmiles("F/C=C/Cl"))
    assert rdDistGeom.EmbedMolecule(alkene, randomSeed=42) == 0
    _native(alkene)
    expected = Chem.MolToSmiles(alkene)
    tetrahedral = _molecule(_TETRAHEDRON, symbol="C", center=(3, 0, 0))
    mol = Chem.CombineMols(alkene, tetrahedral)
    _native(mol)
    rdkit_stereo.apply_3d_stereo_workarounds(mol)
    assert expected in Chem.MolToSmiles(mol).split(".")
    # Also check the public post-processing entry point's retained bond state.
    make_stereochemistry(mol)
    assert expected in Chem.MolToSmiles(mol).split(".")
    assert mol.GetBondWithIdx(1).GetStereo() == alkene.GetBondWithIdx(1).GetStereo()
    assert list(mol.GetBondWithIdx(1).GetStereoAtoms()) == list(
        alkene.GetBondWithIdx(1).GetStereoAtoms()
    )
    assert mol.GetAtomWithIdx(alkene.GetNumAtoms()).GetProp("_CIPCode") == "S"


def _coordination_state(mol):
    atom = mol.GetAtomWithIdx(0)
    return atom.GetChiralTag(), atom.GetIntProp("_chiralPermutation")


def _geometry_signature(points):
    points = np.asarray(points)
    # Gram matrix distinguishes ligand relationships; signed triple products
    # additionally distinguish reflections. Both are invariant under rotation.
    return tuple(np.round(points @ points.T, 6).ravel()) + tuple(
        round(float(np.linalg.det(points[list(triple)])), 6)
        for triple in combinations(range(len(points)), 3)
    )


@pytest.mark.parametrize(("points", "classes"), [(_SP, 3), (_TB, 20), (_AXES, 30)])
def test_all_complete_coordination_permutations_preserve_native_geometry(points, classes):
    signatures = {}
    for occupied in permutations(points):
        native = _native(_molecule(occupied))
        state = _coordination_state(native)
        signature = _geometry_signature(occupied)
        assert signatures.setdefault(state, signature) == signature
        patched = make_stereochemistry(_molecule(occupied))
        assert _coordination_state(patched) == state
    assert len(signatures) == classes


@pytest.mark.parametrize(
    ("points", "count", "eligible"), [(_SP, 3, 24), (_TB, 4, 72), (_AXES, 5, 720)]
)
@pytest.mark.parametrize("dative", [False, True])
def test_sparse_coordination_matches_all_missing_slot_labelings(points, count, eligible, dative):
    checked = 0
    for occupied_indices in permutations(range(len(points)), count):
        occupied = tuple(points[i] for i in occupied_indices)
        missing = tuple(points[i] for i in range(len(points)) if i not in occupied_indices)
        expected = {
            _coordination_state(_native(_molecule(occupied + vacant, dative=dative)))
            for vacant in permutations(missing)
        }
        sparse = make_stereochemistry(_molecule(occupied, dative=dative))
        # TB subsets missing an axial ligand have no trans pair and are not
        # identified as TB by RDKit. Geometry classification remains upstream.
        if sparse.GetAtomWithIdx(0).GetChiralTag() != next(iter(expected))[0]:
            continue
        assert _coordination_state(sparse) in expected
        checked += 1
    assert checked == eligible


@pytest.mark.parametrize(
    "points",
    [
        # A complete TB whose equatorial triangle does not enclose the axis:
        # RDKit's first-two-equatorial-ligand triples can disagree in sign.
        ((0, 0, 1), (1, 0, 0), (0.8, 0.6, 0), (0.8, -0.6, 0), (0, 0, -1)),
        # Coplanar OH seesaw, square pyramid, and six-ligand environments.
        ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0.7, 0.7, 0)),
        ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0.7, 0.7, 0)),
        (
            (1, 0, 0),
            (-1, 0, 0),
            (0.5, 0.866, 0),
            (-0.5, -0.866, 0),
            (-0.5, 0.866, 0),
            (0.5, -0.866, 0),
        ),
    ],
)
def test_conflicting_or_degenerate_coordination_is_unspecified_for_every_order(points):
    for order in permutations(range(1, len(points) + 1)):
        mol = _molecule(points, order=order)
        for _ in range(2):
            make_stereochemistry(mol)
            atom = mol.GetAtomWithIdx(0)
            assert atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
            assert not atom.HasProp("_chiralPermutation")


@pytest.mark.parametrize(("points", "seed"), [(_SP, 3), (_TB, 4), (_AXES[:5], 4)])
def test_distorted_coordination_is_invariant_under_neighbor_order(points, seed):
    # An independent rotation/reflection signature decodes RDKit's full class
    # numbers. Sparse, physically equivalent missing-slot permutations compare
    # equal without relying on SMILES canonicalization or embedding.
    templates = {}
    for template in (_SP, _TB, _AXES):
        for occupied in permutations(template):
            templates[_coordination_state(_native(_molecule(occupied)))] = np.asarray(occupied)
    distorted = np.asarray(points) + np.random.default_rng(seed).normal(0, 0.5, (len(points), 3))
    fingerprints = set()
    for order in permutations(range(len(points))):
        mol = make_stereochemistry(_molecule(distorted, order=tuple(i + 1 for i in order)))
        atom = mol.GetAtomWithIdx(0)
        if atom.GetChiralTag() in rdkit_stereo._NONTETRAHEDRAL_TAGS:
            ideal = templates[_coordination_state(mol)][: len(points)][np.argsort(order)]
            fingerprints.add((str(atom.GetChiralTag()), _geometry_signature(ideal)))
        else:
            # Raw CW/CCW must flip for odd permutations, unlike the actual
            # molecular configuration represented by canonical SMILES.
            fingerprints.add(Chem.MolToSmiles(mol))
    assert len(fingerprints) == 1


def _alkene(order, angle):
    rw = Chem.RWMol()
    for symbol in ("C", "C", "F", "Cl", "Br", "I"):
        rw.AddAtom(Chem.Atom(symbol))
    bonds = ((0, 1, Chem.BondType.DOUBLE),) + tuple(
        (a, b, Chem.BondType.SINGLE) for a, b in ((0, 2), (0, 3), (1, 4), (1, 5))
    )
    for idx in order:
        rw.AddBond(*bonds[idx])
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    conf = Chem.Conformer(6)
    conf.Set3D(True)
    points = (
        (-0.7, 0, 0),
        (0.7, 0, 0),
        (-1, 1, 0),
        (-1, -1, 0),
        (1, np.cos(angle), np.sin(angle)),
        (1, -np.cos(angle), -np.sin(angle)),
    )
    for idx, point in enumerate(points):
        conf.SetAtomPosition(idx, tuple(float(x) for x in point))
    mol.AddConformer(conf)
    return mol


@pytest.mark.parametrize("angle", [0.0, 0.4, np.pi / 2, np.pi / 2 + 1e-5, 2.8, np.pi])
@pytest.mark.parametrize("legacy", [False, True])
def test_double_bond_order_invariance_and_roundtrip(angle, legacy):
    Chem.SetUseLegacyStereoPerception(legacy)
    outputs = set()
    for order in permutations(range(5)):
        mol = _alkene(order, angle)
        for _ in range(2):
            make_stereochemistry(mol)
            smiles = Chem.MolToSmiles(mol)
            assert Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) == smiles
            outputs.add(smiles)
            if abs(np.cos(angle)) < 1e-3:
                assert mol.GetBondBetweenAtoms(0, 1).GetStereo() == Chem.BondStereo.STEREOANY
                assert "/" not in smiles and "\\" not in smiles
            else:
                assert "/" in smiles or "\\" in smiles
    assert len(outputs) == 1


def test_ambiguous_double_bond_preserves_adjacent_valid_double_bond():
    mol = Chem.MolFromSmiles("C(=CC=CCl)F")
    conf = Chem.Conformer(6)
    conf.Set3D(True)
    # Atom order: C=C-C=C-Cl, with F attached to the first carbon.
    for idx, point in enumerate(
        ((0, 0, 0), (1, 0, 0), (2, 1, 0), (3, 1, 0), (4, 2, 0), (-1, 0, 1))
    ):
        conf.SetAtomPosition(idx, point)
    mol.AddConformer(conf)
    make_stereochemistry(mol)
    assert mol.GetBondBetweenAtoms(0, 1).GetStereo() == Chem.BondStereo.STEREOANY
    assert mol.GetBondBetweenAtoms(2, 3).GetStereo() == Chem.BondStereo.STEREOE
    restored = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    assert sum(b.GetStereo() == Chem.BondStereo.STEREOE for b in restored.GetBonds()) == 1


@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize("update", ["coordinates", "conformer", "copy"])
def test_derived_unknown_double_bond_is_recomputed_after_geometry_changes(legacy, update):
    Chem.SetUseLegacyStereoPerception(legacy)
    for order in permutations(range(5)):
        mol = _alkene(order, 0.0)
        for angle in (0.0, np.pi / 2, np.pi, np.pi / 2, 0.0):
            reference = _alkene(order, angle)
            conf_id = -1
            if update == "conformer":
                conf = Chem.Conformer(reference.GetConformer())
                conf.SetId(7)
                if mol.GetNumConformers() > 1:
                    mol.RemoveConformer(7)
                mol.AddConformer(conf, assignId=False)
                conf_id = 7
            else:
                if update == "copy":
                    mol = Chem.Mol(mol)
                for idx in range(mol.GetNumAtoms()):
                    mol.GetConformer().SetAtomPosition(
                        idx, reference.GetConformer().GetAtomPosition(idx)
                    )
            make_stereochemistry(reference)
            if update == "conformer":
                Chem.AssignStereochemistryFrom3D(mol, confId=conf_id)
                rdkit_stereo.apply_3d_stereo_workarounds(mol, conf_id=conf_id)
            else:
                make_stereochemistry(mol)
            assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(reference)
            bond = mol.GetBondBetweenAtoms(0, 1)
            assert bond.GetStereo() == reference.GetBondBetweenAtoms(0, 1).GetStereo()
            assert bool(bond.HasProp(rdkit_stereo._DERIVED_UNKNOWN_STEREO)) == (
                bond.GetStereo() == Chem.BondStereo.STEREOANY
            )


@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize("annotation", ["any", "crossed", "wiggly", "property"])
def test_reperception_preserves_explicit_double_bond_uncertainty(legacy, annotation):
    Chem.SetUseLegacyStereoPerception(legacy)
    mol = _alkene(range(5), np.pi / 2)
    # Test input STEREOANY separately from explicit annotations applied after
    # the workaround has already marked a geometry as ambiguous.
    if annotation != "any":
        make_stereochemistry(mol)
    bond = mol.GetBondBetweenAtoms(0, 1)
    bond.SetStereo(Chem.BondStereo.STEREOANY)
    if annotation == "crossed":
        bond.SetBondDir(Chem.BondDir.EITHERDOUBLE)
    elif annotation == "wiggly":
        mol.GetBondBetweenAtoms(0, 2).SetBondDir(Chem.BondDir.UNKNOWN)
    elif annotation == "property":
        bond.SetIntProp("_UnknownStereo", 1)
    reference = _alkene(range(5), 0.0)
    for idx in range(mol.GetNumAtoms()):
        mol.GetConformer().SetAtomPosition(idx, reference.GetConformer().GetAtomPosition(idx))
    for _ in range(2):
        make_stereochemistry(mol)
        assert bond.GetStereo() == Chem.BondStereo.STEREOANY
        assert "/" not in Chem.MolToSmiles(mol)


def test_reperception_of_derived_unknown_keeps_other_stereo_repairs():
    alkene = _alkene(range(5), np.pi / 2)
    tetrahedron = _molecule(_TETRAHEDRON, symbol="C", center=(3, 0, 0))
    mol = Chem.CombineMols(alkene, tetrahedron)
    make_stereochemistry(mol)
    valid = _alkene(range(5), np.pi)
    for idx in range(valid.GetNumAtoms()):
        mol.GetConformer().SetAtomPosition(idx, valid.GetConformer().GetAtomPosition(idx))
    # Exercise the repair-only entry as well as the integrated perception path.
    Chem.AssignStereochemistryFrom3D(mol)
    rdkit_stereo.apply_3d_stereo_workarounds(mol)
    reference = Chem.CombineMols(valid, tetrahedron)
    make_stereochemistry(reference)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(reference)
    assert mol.GetAtomWithIdx(6).HasProp("_CIPCode")


@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize("geometry", ["collinear", "collapsed"])
def test_native_geometric_unknown_is_recomputed(legacy, geometry):
    Chem.SetUseLegacyStereoPerception(legacy)
    mol = _alkene(range(5), 0.0)
    conf = mol.GetConformer()
    if geometry == "collinear":
        conf.SetAtomPosition(2, (-2, 0, 0))
        conf.SetAtomPosition(3, (-3, 0, 0))
    else:
        conf.SetAtomPosition(1, conf.GetAtomPosition(0))
    for _ in range(2):
        make_stereochemistry(mol)
        assert mol.GetBondWithIdx(0).GetStereo() == Chem.BondStereo.STEREOANY
    reference = _alkene(range(5), 0.0)
    for idx in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(idx, reference.GetConformer().GetAtomPosition(idx))
    make_stereochemistry(mol)
    make_stereochemistry(reference)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(reference)


@pytest.mark.parametrize("annotation", ["crossed", "wiggly", "property"])
def test_input_unknown_directions_are_not_marked_as_geometric_abstention(annotation):
    mol = _alkene(range(5), 0.0)
    bond = mol.GetBondWithIdx(0)
    if annotation == "crossed":
        bond.SetBondDir(Chem.BondDir.EITHERDOUBLE)
    elif annotation == "wiggly":
        mol.GetBondBetweenAtoms(0, 2).SetBondDir(Chem.BondDir.UNKNOWN)
    else:
        bond.SetIntProp("_UnknownStereo", 1)
    for _ in range(2):
        make_stereochemistry(mol)
        assert not bond.HasProp(rdkit_stereo._DERIVED_UNKNOWN_STEREO)


def test_unknown_tetrahedral_stereo_survives_repeated_postprocessing():
    mol = _molecule(_TETRAHEDRON, symbol="C")
    mol.GetBondWithIdx(0).SetBondDir(Chem.BondDir.UNKNOWN)
    for _ in range(3):
        make_stereochemistry(mol)
        assert mol.GetAtomWithIdx(0).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
        assert mol.GetBondWithIdx(0).GetBondDir() == Chem.BondDir.UNKNOWN


def test_obsolete_coordination_permutation_is_removed():
    mol = _molecule(_TETRAHEDRON, symbol="C")
    mol.GetAtomWithIdx(0).SetIntProp("_chiralPermutation", 19)
    make_stereochemistry(mol)
    assert not mol.GetAtomWithIdx(0).HasProp("_chiralPermutation")


def test_candidate_filter_cleans_all_stale_properties_and_preserves_explicit_unknown():
    mol = Chem.MolFromSmiles("CC." + ".".join(["[He]"] * 1001))
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    mol.AddConformer(conf)
    mol.GetBondWithIdx(0).SetBondDir(Chem.BondDir.UNKNOWN)
    for idx in range(mol.GetNumAtoms()):
        mol.GetAtomWithIdx(idx).SetIntProp("_chiralPermutation", 19)
    rdkit_stereo.apply_3d_stereo_workarounds(mol)
    # Preserve the original unknown-stereo handling at the annotation's begin
    # atom. Clearing properties must not skip subsequent query matches.
    assert mol.GetAtomWithIdx(0).HasProp("_chiralPermutation")
    assert all(
        not mol.GetAtomWithIdx(idx).HasProp("_chiralPermutation")
        for idx in range(1, mol.GetNumAtoms())
    )


def test_double_bond_candidate_search_does_not_truncate_at_1000_matches():
    mol = Chem.MolFromSmiles(".".join(["FC=CCl"] * 1001))
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    points = ((-1, 1, 0), (0, 0, 0), (1, 0, 0), (2, 0, 1))
    for idx in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(idx, points[idx % 4])
    mol.AddConformer(conf)
    make_stereochemistry(mol)
    assert all(
        mol.GetBondBetweenAtoms(idx + 1, idx + 2).GetStereo() == Chem.BondStereo.STEREOANY
        for idx in range(0, mol.GetNumAtoms(), 4)
    )


@pytest.mark.parametrize("version", ["2026.09.1", "2026.09.2"])
def test_direction_preservation_is_not_an_rdkit_version_workaround(monkeypatch, version):
    monkeypatch.setattr(rdkit_stereo.rdBase, "rdkitVersion", version)
    mol = _alkene(range(5), 0)
    make_stereochemistry(mol)
    assert "/" in Chem.MolToSmiles(mol)
