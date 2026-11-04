"""Local 3D stereo workarounds for RDKit 2023.09.6 through 2026.09.1.

See https://github.com/rdkit/rdkit/issues/9760 (four explicit tetrahedral
ligands) and https://github.com/rdkit/rdkit/issues/9761 (four-ligand OH).
Also covers the analogous four-ligand TB sign error, order-dependent reference
vectors in distorted TB/OH environments, and ambiguous/twisted double bonds.
The upper bound is the last reviewed release, not an assumed upstream fix.
Newer, prerelease, and unrecognized versions use RDKit's native behavior.
"""

from __future__ import annotations

import re
from functools import lru_cache
from itertools import combinations, product
from math import isfinite, sqrt
from typing import Optional, Sequence, Tuple

from rdkit import Chem, rdBase
from rdkit.Chem import rdqueries


# Dimensionless affine volume, normalized by the largest ligand-ligand distance
# cubed. This is a MolGR degeneracy policy, not RDKit's center-based 0.1 A^3 cutoff.
_TETRAHEDRAL_VOLUME_TOL = 1e-3
_DIRECTION_TOL = 1e-3
_DERIVED_UNKNOWN_STEREO = "_MolGRDerivedUnknownStereo"
# Only the confirmed incorrect sparse branches are replaced. Complete
# environments retain RDKit's permutation tables.
_OH_CORRECTIONS = {
    (1, 2): ((0, 2, 3), (10, 8)),
    (1, 3): ((0, 1, 2), (1, 2)),
    (2, 3): ((0, 1, 2), (4, 14)),
}
_NONTETRAHEDRAL_TAGS = (
    Chem.ChiralType.CHI_SQUAREPLANAR,
    Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL,
    Chem.ChiralType.CHI_OCTAHEDRAL,
)
_TETRAHEDRAL_TAGS = (
    Chem.ChiralType.CHI_UNSPECIFIED,
    Chem.ChiralType.CHI_TETRAHEDRAL_CW,
    Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
)
_Vector = Tuple[float, float, float]


def _candidate_atom_query() -> Chem.Mol:
    # Run the degree/property scan inside RDKit. This only selects candidates;
    # bond filtering and geometric perception still use the rules below.
    atom = rdqueries.ExplicitDegreeGreaterQueryAtom(3)
    atom.ExpandQuery(
        rdqueries.HasPropQueryAtom("_chiralPermutation"), Chem.CompositeQueryType.COMPOSITE_OR
    )
    query = Chem.RWMol()
    query.AddAtom(atom)
    return query.GetMol()


def _achiral_carbon_query() -> Chem.Mol:
    # A carbon with two identical ordinary terminal H ligands cannot be a
    # tetrahedral stereocenter. Only skip centers already perceived/cleaned by
    # RDKit; leave stale metadata and unusual hydrogen states on the slow path.
    center = rdqueries.AtomNumEqualsQueryAtom(6)
    for condition in (
        rdqueries.ExplicitDegreeEqualsQueryAtom(4),
        rdqueries.HasChiralTagQueryAtom(negate=True),
        rdqueries.HasIntPropWithValueQueryAtom("_NonExplicit3DChirality", 1),
        rdqueries.HasPropQueryAtom("_chiralPermutation", negate=True),
        rdqueries.HasPropQueryAtom("_CIPCode", negate=True),
    ):
        center.ExpandQuery(condition)
    hydrogen = rdqueries.AtomNumEqualsQueryAtom(1)
    for condition in (
        rdqueries.ExplicitDegreeEqualsQueryAtom(1),
        rdqueries.IsotopeEqualsQueryAtom(0),
        rdqueries.FormalChargeEqualsQueryAtom(0),
        rdqueries.NumRadicalElectronsEqualsQueryAtom(0),
        rdqueries.HasChiralTagQueryAtom(negate=True),
        rdqueries.HasPropQueryAtom("molAtomMapNumber", negate=True),
    ):
        hydrogen.ExpandQuery(condition)
    query = Chem.RWMol()
    query.AddAtom(center)
    query.AddAtom(hydrogen)
    query.AddAtom(hydrogen)
    query.AddBond(0, 1, Chem.BondType.SINGLE)
    query.AddBond(0, 2, Chem.BondType.SINGLE)
    return query.GetMol()


_CANDIDATE_ATOM_QUERY = _candidate_atom_query()
_ACHIRAL_CARBON_QUERY = _achiral_carbon_query()
_DOUBLE_BOND_QUERY = Chem.MolFromSmarts("*=*")
assert _DOUBLE_BOND_QUERY is not None


def _subtract(a: Sequence[float], b: Sequence[float]) -> _Vector:
    return a[0] - b[0], a[1] - b[1], a[2] - b[2]


def _scale(a: Sequence[float], factor: float) -> _Vector:
    return a[0] * factor, a[1] * factor, a[2] * factor


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _triple(a: Sequence[float], b: Sequence[float], c: Sequence[float]) -> float:
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        + a[1] * (b[2] * c[0] - b[0] * c[2])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


@lru_cache(maxsize=32)
def _needs_stereo_workaround(version: str) -> bool:
    if re.fullmatch(r"\d+\.\d+\.\d+", version) is None:
        return False
    return (2023, 9, 6) <= tuple(map(int, version.split("."))) <= (2026, 9, 1)


def _has_unknown_bond(atom: Chem.Atom) -> bool:
    # Match RDKit's isWigglyBond(), including the direction of the annotation.
    return any(
        bond.GetBeginAtomIdx() == atom.GetIdx()
        and bond.GetBondType() == Chem.BondType.SINGLE
        and (
            bond.GetBondDir() == Chem.BondDir.UNKNOWN
            or (bond.HasProp("_UnknownStereo") and bond.GetIntProp("_UnknownStereo"))
        )
        for bond in atom.GetBonds()
    )


def _tetrahedral_neighbors(atom: Chem.Atom) -> Optional[list[int]]:
    # Match bondAffectsAtomChirality(): a donor's outgoing dative bond is not
    # a tetrahedral ligand; an acceptor's incoming dative bond is.
    # Detect explicit uncertainty in the same traversal; None preserves it.
    atom_idx = atom.GetIdx()
    neighbors = []
    for bond in atom.GetBonds():
        bond_type = bond.GetBondType()
        if bond_type == Chem.BondType.SINGLE:
            if bond.GetBeginAtomIdx() == atom_idx and (
                bond.GetBondDir() == Chem.BondDir.UNKNOWN
                or (bond.HasProp("_UnknownStereo") and bond.GetIntProp("_UnknownStereo"))
            ):
                return None
        elif bond_type in (Chem.BondType.UNSPECIFIED, Chem.BondType.ZERO) or (
            bond_type == Chem.BondType.DATIVE and bond.GetBeginAtomIdx() == atom_idx
        ):
            continue
        neighbors.append(bond.GetOtherAtomIdx(atom_idx))
    return neighbors


def _repair_tetrahedral(
    atom: Chem.Atom, positions: list[list[float]], neighbors: list[int]
) -> bool:
    if len(neighbors) != 4 or atom.GetTotalNumHs() != 0:
        return False
    points = [positions[idx] for idx in neighbors]
    if not all(isfinite(x) for point in points for x in point):
        return False
    differences = [_subtract(points[i], points[j]) for i, j in combinations(range(4), 2)]
    scale = sqrt(max(_dot(v, v) for v in differences))
    volume = (
        _triple(*(_scale(_subtract(point, points[3]), 1 / scale) for point in points[:3]))
        if scale > 0
        else 0.0
    )
    tag = Chem.ChiralType.CHI_UNSPECIFIED
    if volume > _TETRAHEDRAL_VOLUME_TOL:
        tag = Chem.ChiralType.CHI_TETRAHEDRAL_CCW
    elif volume < -_TETRAHEDRAL_VOLUME_TOL:
        tag = Chem.ChiralType.CHI_TETRAHEDRAL_CW
    if atom.GetChiralTag() == tag:
        return False
    atom.SetChiralTag(tag)
    if tag == Chem.ChiralType.CHI_UNSPECIFIED:
        if atom.HasProp("_NonExplicit3DChirality"):
            atom.ClearProp("_NonExplicit3DChirality")
    else:
        atom.SetIntProp("_NonExplicit3DChirality", 1)
    return True


def _repair_nontetrahedral(atom: Chem.Atom, positions: list[list[float]]) -> bool:
    # Unlike tetrahedral perception, RDKit's non-tetrahedral branch uses all
    # graph neighbors in their internal order. The geometry class itself is
    # still selected by RDKit; this fixes its orientation, not the classifier.
    count = atom.GetDegree()
    if count not in (4, 5, 6):
        return False
    center = positions[atom.GetIdx()]
    vectors = []
    for neighbor in atom.GetNeighbors():
        vector = _subtract(positions[neighbor.GetIdx()], center)
        length = sqrt(_dot(vector, vector))
        if not all(isfinite(x) for x in vector) or length == 0:
            return False
        vectors.append(_scale(vector, 1 / length))
    trans_pairs = [
        (i, j) for i, j in combinations(range(count), 2) if _dot(vectors[i], vectors[j]) < -0.9
    ]
    paired = {idx for pair in trans_pairs for idx in pair}
    if len(paired) != 2 * len(trans_pairs):
        return False
    unpaired = [idx for idx in range(count) if idx not in paired]
    if atom.GetChiralTag() == Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL:
        if count not in (4, 5) or len(trans_pairs) != 1:
            return False
        i, j = trans_pairs[0]
        equatorial_pairs = [(unpaired[0], unpaired[1])]
        if count == 5:
            equatorial_pairs += [(unpaired[1], unpaired[2]), (unpaired[2], unpaired[0])]
        volumes = [
            _triple(axis, vectors[a], vectors[b])
            for axis in (vectors[i], _scale(vectors[j], -1))
            for a, b in equatorial_pairs
        ]
    else:
        if len(trans_pairs) != count - 3:
            return False
        alternatives = [[vectors[i], _scale(vectors[j], -1)] for i, j in trans_pairs]
        alternatives += [[vectors[idx]] for idx in unpaired]
        volumes = [_triple(a, b, c) for a, b, c in product(*alternatives)]
    # If interchangeable references disagree, this distorted environment has
    # no reliable RDKit orientation. Abstain consistently instead of inventing
    # a geometry fit or keeping whichever reference was inserted first.
    if not (all(v > _DIRECTION_TOL for v in volumes) or all(v < -_DIRECTION_TOL for v in volumes)):
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        if atom.HasProp("_NonExplicit3DChirality"):
            atom.ClearProp("_NonExplicit3DChirality")
        return True
    if count != 4:
        return False
    if atom.GetChiralTag() == Chem.ChiralType.CHI_OCTAHEDRAL:
        if trans_pairs[0] not in _OH_CORRECTIONS:
            return False
        triple, permutations = _OH_CORRECTIONS[trans_pairs[0]]
    elif trans_pairs[0] == (2, 3):
        triple, permutations = (2, 0, 1), (16, 19)
    else:
        return False
    i, j, k = triple
    volume = _triple(vectors[i], vectors[j], vectors[k])
    permutation = permutations[0 if volume >= 0 else 1]
    if atom.HasProp("_chiralPermutation") and atom.GetIntProp("_chiralPermutation") == permutation:
        return False
    atom.SetIntProp("_chiralPermutation", permutation)
    return True


def _repair_double_bonds(
    mol: Chem.Mol, positions: list[list[float]], candidates: Sequence[Sequence[int]]
) -> bool:
    """Use both explicit substituents at each end; reject ambiguous projections."""
    changed = False
    for begin_idx, end_idx in candidates:
        bond = mol.GetBondBetweenAtoms(begin_idx, end_idx)
        if bond.GetStereo() not in (
            Chem.BondStereo.STEREOCIS,
            Chem.BondStereo.STEREOTRANS,
            Chem.BondStereo.STEREOE,
            Chem.BondStereo.STEREOZ,
        ):
            continue
        stereo_atoms = list(bond.GetStereoAtoms())
        if len(stereo_atoms) != 2:
            continue
        begin, end = bond.GetBeginAtom(), bond.GetEndAtom()
        axis = _subtract(positions[end.GetIdx()], positions[begin.GetIdx()])
        axis_length = sqrt(_dot(axis, axis))
        if not all(isfinite(x) for x in axis) or axis_length == 0:
            continue
        axis = _scale(axis, 1 / axis_length)
        projections = []
        for atom, selected in zip((begin, end), stereo_atoms):
            neighbors = [
                b.GetOtherAtomIdx(atom.GetIdx())
                for b in atom.GetBonds()
                if b.GetIdx() != bond.GetIdx()
                and b.GetBondType() in (Chem.BondType.SINGLE, Chem.BondType.AROMATIC)
            ]
            if selected not in neighbors or len(neighbors) not in (1, 2):
                break
            side = []
            for neighbor in neighbors:
                vector = _subtract(positions[neighbor], positions[atom.GetIdx()])
                if not all(isfinite(x) for x in vector):
                    break
                projected = _subtract(vector, _scale(axis, _dot(vector, axis)))
                norm = sqrt(_dot(projected, projected))
                length = sqrt(_dot(vector, vector))
                direction = (
                    _scale(projected, 1 / norm)
                    if norm > _DIRECTION_TOL * length and norm > 0
                    else (0.0, 0.0, 0.0)
                )
                side.append(direction if neighbor == selected else _scale(direction, -1))
            if len(side) != len(neighbors):
                break
            projections.append(side)
        if len(projections) != 2:
            continue
        cosines = [_dot(a, b) for a, b in product(*projections)]
        if all(c > _DIRECTION_TOL for c in cosines) or all(c < -_DIRECTION_TOL for c in cosines):
            continue
        bond.SetStereo(Chem.BondStereo.STEREOANY)
        bond.SetBoolProp(_DERIVED_UNKNOWN_STEREO, True)
        changed = True
    if changed:
        # A single directional bond can serve two adjacent double bonds.
        # Regenerate all slash directions together, preserving explicit unknowns
        # and wedges. This avoids leaving E/Z information for a cleared bond.
        for bond_idx in range(mol.GetNumBonds()):
            bond = mol.GetBondWithIdx(bond_idx)
            if bond.GetBondDir() in (Chem.BondDir.ENDUPRIGHT, Chem.BondDir.ENDDOWNRIGHT):
                bond.SetBondDir(Chem.BondDir.NONE)
        Chem.SetDoubleBondNeighborDirections(mol)
    return changed


def apply_3d_stereo_workarounds(
    mol: Chem.Mol, conf_id: int = -1, *, perceive: bool = False
) -> None:
    """Repair freshly perceived 3D tags in place, then refresh stereo metadata.

    With perceive=True, also run RDKit's native 3D perception, after resetting
    uncertainty generated by a previous repair. This preserves input crossed
    bonds before RDKit clears their direction flags. Otherwise call after the
    last AssignStereochemistryFrom3D(): running that function again would
    overwrite the repairs. Four explicit tetrahedral ligands use
    their affine orientation even when the center lies outside their hull;
    nearly coplanar ligands are left unspecified. Chemical stereo eligibility
    (e.g. duplicate ligands) is still cleaned by RDKit, without re-reading 3D.
    """
    if not _needs_stereo_workaround(rdBase.rdkitVersion) or not mol.GetNumConformers():
        if perceive:
            Chem.AssignStereochemistryFrom3D(mol, confId=conf_id)
        return
    conf = mol.GetConformer(conf_id)
    if not conf.Is3D():
        if perceive:
            Chem.AssignStereochemistryFrom3D(mol, confId=conf_id)
        return
    # Override the default match limit: large or disconnected molecules can
    # have more than 1000 candidates. These queries do not inspect chirality.
    double_bonds = mol.GetSubstructMatches(_DOUBLE_BOND_QUERY, maxMatches=mol.GetNumBonds())
    # RDKit deliberately preserves STEREOANY even when replacing 3D tags.
    # Geometric abstention must instead be reconsidered after coordinate edits
    # or a conformer switch. Only reset uncertainty whose geometric origin is
    # known; input STEREOANY and explicit crossed/wiggly annotations survive.
    refresh = False
    perceived_bonds = []
    for begin_idx, end_idx in double_bonds:
        bond = mol.GetBondBetweenAtoms(begin_idx, end_idx)
        if (
            perceive
            and bond.GetStereo() != Chem.BondStereo.STEREOANY
            and bond.GetBondDir() == Chem.BondDir.NONE
            and not (bond.HasProp("_UnknownStereo") and bond.GetIntProp("_UnknownStereo"))
        ):
            perceived_bonds.append(bond)
        if not bond.HasProp(_DERIVED_UNKNOWN_STEREO):
            continue
        bond.ClearProp(_DERIVED_UNKNOWN_STEREO)
        if (
            bond.GetStereo() == Chem.BondStereo.STEREOANY
            and bond.GetBondDir() == Chem.BondDir.NONE
            and not (bond.HasProp("_UnknownStereo") and bond.GetIntProp("_UnknownStereo"))
        ):
            bond.SetStereo(Chem.BondStereo.STEREONONE)
            perceived_bonds.append(bond)
            refresh = True
    if perceive or refresh:
        # Rerun before repairing atoms: a later native pass would overwrite
        # the tetrahedral and coordination corrections below.
        Chem.AssignStereochemistryFrom3D(mol, confId=conf_id)
        for bond in perceived_bonds:
            if bond.GetStereo() != Chem.BondStereo.STEREOANY:
                continue
            # A wiggly neighboring bond also produces STEREOANY. It is input
            # uncertainty, not an abstention caused by collinear coordinates.
            related = [bond, *bond.GetBeginAtom().GetBonds(), *bond.GetEndAtom().GetBonds()]
            if any(
                b.GetBondDir() in (Chem.BondDir.UNKNOWN, Chem.BondDir.EITHERDOUBLE)
                or (b.HasProp("_UnknownStereo") and b.GetIntProp("_UnknownStereo"))
                for b in related
            ):
                continue
            bond.SetBoolProp(_DERIVED_UNKNOWN_STEREO, True)
    atoms = mol.GetSubstructMatches(_CANDIDATE_ATOM_QUERY, maxMatches=mol.GetNumAtoms())
    if atoms and not mol.HasQuery():
        # At most six H pairs per degree-four carbon; two times the total atom
        # count also covers methane. This is only an optimization: a missed
        # match would still receive the full geometric check below.
        achiral = {
            match[0]
            for match in mol.GetSubstructMatches(
                _ACHIRAL_CARBON_QUERY, maxMatches=2 * mol.GetNumAtoms()
            )
        }
        atoms = tuple(match for match in atoms if match[0] not in achiral)
    if not atoms and not double_bonds:
        return
    # One bulk read avoids repeated Point3D construction and Python/C++ calls.
    # Coordinates are local to this invocation, so conformer edits stay visible.
    positions: list[list[float]] = conf.GetPositions().tolist()
    changed = False
    for (atom_idx,) in atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Filtering can only remove neighbors, so degree < 4 cannot reach any
        # repaired geometry. Keep cleanup of pre-existing stale metadata.
        if atom.GetDegree() < 4:
            if (
                atom.HasProp("_chiralPermutation")
                and atom.GetChiralTag() not in _NONTETRAHEDRAL_TAGS
                and not _has_unknown_bond(atom)
            ):
                atom.ClearProp("_chiralPermutation")
            continue
        tag = atom.GetChiralTag()
        if tag in _TETRAHEDRAL_TAGS:
            neighbors = _tetrahedral_neighbors(atom)
            if neighbors is None:
                continue
            changed |= _repair_tetrahedral(atom, positions, neighbors)
        else:
            if _has_unknown_bond(atom):
                continue
            if tag in (Chem.ChiralType.CHI_OCTAHEDRAL, Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL):
                changed |= _repair_nontetrahedral(atom, positions)
        if atom.GetChiralTag() not in _NONTETRAHEDRAL_TAGS and atom.HasProp("_chiralPermutation"):
            atom.ClearProp("_chiralPermutation")
    changed |= _repair_double_bonds(mol, positions, double_bonds)
    if changed:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
