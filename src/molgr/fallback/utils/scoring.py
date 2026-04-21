from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Callable, Dict, List, Tuple, TypeVar, Union, cast

from openbabel import openbabel as ob
from openbabel import pybel
from typing_extensions import Protocol

from molgr.fallback.state import ReconstructionState
from molgr.fallback.utils import consts, smarts
from molgr.fallback.utils.dataclasses import MetalAtomPosition
from molgr.fallback.utils.tools import typed_lru_cache


ScoreKey = str
MetalStateKey = Tuple[Tuple[int, int, int, int, int, int, int], ...]
PostReinsertionKey = Tuple[ScoreKey, MetalStateKey]
ChargedAtomSnapshot = Tuple[int, float, float, float]
_OMOL_SCORE_CACHE_MAXSIZE = 4096
_OMOL_SCORE_INPUTS: Dict[ScoreKey, pybel.Molecule] = {}
_ORGANIC_SCORE_INPUTS: Dict[ScoreKey, pybel.Molecule] = {}
_POST_REINSERTION_SCORE_INPUTS: Dict[PostReinsertionKey, pybel.Molecule] = {}
_POST_REINSERTION_METAL_STATE_INPUTS: Dict[
    PostReinsertionKey,
    Tuple[float, Tuple[ChargedAtomSnapshot, ...], Tuple[MetalAtomPosition, ...]],
] = {}
_COORDINATE_SCALE = 1_000_000
_EMPTY_METAL_STATE_KEY: MetalStateKey = ()
T = TypeVar("T")
ScoringInput = Union[pybel.Molecule, ReconstructionState]


class CacheInfoLike(Protocol):
    hits: int
    misses: int
    currsize: int


def _resolve_cached_omol_value(
    omol_or_state: ScoringInput,
    cache_name: str,
    builder: Callable[[pybel.Molecule], T],
) -> Tuple[pybel.Molecule, T]:
    if isinstance(omol_or_state, ReconstructionState):
        return omol_or_state.omol, omol_or_state.get_cached_omol_value(cache_name, builder)
    return omol_or_state, builder(omol_or_state)


def _quantized_coordinate(value: float) -> int:
    return int(round(float(value) * _COORDINATE_SCALE))


def _build_score_key(omol: pybel.Molecule) -> ScoreKey:
    return omol.write("molreport") or ""


def build_post_reinsertion_base_key(omol_or_state: ScoringInput) -> ScoreKey:
    if isinstance(omol_or_state, ReconstructionState):
        return omol_or_state.post_reinsertion_base_key()
    _omol, score_key = _resolve_cached_omol_value(
        omol_or_state,
        "organic_score_key",
        _build_score_key,
    )
    return score_key


def build_metal_state_key(metal_states: Sequence[MetalAtomPosition]) -> MetalStateKey:
    return tuple(
        (
            metal_state.idx,
            metal_state.element_idx,
            metal_state.valence,
            metal_state.radical_num,
            _quantized_coordinate(metal_state.position_x),
            _quantized_coordinate(metal_state.position_y),
            _quantized_coordinate(metal_state.position_z),
        )
        for metal_state in metal_states
    )


def _compute_omol_score(omol: pybel.Molecule) -> float:
    obmol = cast(ob.OBMol, omol.OBMol)
    obmol.SetAromaticPerceived(False)
    score = _compute_organic_core_score(omol)
    score += _compute_post_reinsertion_score(omol)
    return score


def _compute_organic_core_score(omol: pybel.Molecule) -> float:
    obmol = cast(ob.OBMol, omol.OBMol)
    score = calc_remain_bond_order_reward(omol)

    for atom in ob.OBMolAtomIter(obmol):
        atom = cast(ob.OBAtom, atom)
        if atom.IsMetal():
            continue
        if atom.IsAromatic():
            score -= 5 - abs(cast(int, atom.GetFormalCharge())) * 3

        deviation_score = None
        if atom.GetSpinMultiplicity() > 0 or atom.GetFormalCharge() != 0:
            deviation_score = get_deviation_score(omol, atom.GetIdx())

        if atom.GetSpinMultiplicity() > 0 and deviation_score is not None:
            score += deviation_score * 10.0
        if atom.GetFormalCharge() > 0 and deviation_score is not None:
            score += deviation_score * 10.0
        if atom.GetFormalCharge() < 0 and deviation_score is not None:
            score += (1 - deviation_score) * 10.0

        if atom.GetAtomicNum() == 6 and all(
            cast(int, bond.GetBondOrder()) == 2 for bond in ob.OBAtomBondIter(atom)
        ):
            score += 5

    score += calculate_physchem_penalty(omol)
    score += calculate_heteroatom_penalty(omol)
    score -= calculate_conjugation_reward(omol)
    return score


def _compute_post_reinsertion_score(omol: pybel.Molecule) -> float:
    score = calc_symmetry_penalty(omol)
    score += calculate_metal_penalty(omol)
    return score


@typed_lru_cache(maxsize=_OMOL_SCORE_CACHE_MAXSIZE, typed=True)
def _omol_score_cached(score_key: ScoreKey) -> float:
    omol = _OMOL_SCORE_INPUTS[score_key]
    return _compute_omol_score(omol)


@typed_lru_cache(maxsize=_OMOL_SCORE_CACHE_MAXSIZE, typed=True)
def _organic_core_score_cached(score_key: ScoreKey) -> float:
    omol = _ORGANIC_SCORE_INPUTS[score_key]
    return _compute_organic_core_score(omol)


@typed_lru_cache(maxsize=_OMOL_SCORE_CACHE_MAXSIZE, typed=True)
def _post_reinsertion_score_cached(score_key: PostReinsertionKey) -> float:
    omol = _POST_REINSERTION_SCORE_INPUTS[score_key]
    return _compute_post_reinsertion_score(omol)


@typed_lru_cache(maxsize=_OMOL_SCORE_CACHE_MAXSIZE, typed=True)
def _post_reinsertion_metal_state_score_cached(score_key: PostReinsertionKey) -> float:
    base_symmetry_penalty, charged_atom_snapshots, metal_states = (
        _POST_REINSERTION_METAL_STATE_INPUTS[score_key]
    )
    return _compute_post_reinsertion_score_from_metal_states(
        base_symmetry_penalty,
        charged_atom_snapshots,
        metal_states,
    )


def omol_score_cache_info() -> Tuple[int, int, int]:
    info = cast(CacheInfoLike, _omol_score_cached.cache_info())
    hits = info.hits
    misses = info.misses
    currsize = info.currsize
    return hits, misses, currsize


def omol_score_cache_clear() -> None:
    _OMOL_SCORE_INPUTS.clear()
    _ORGANIC_SCORE_INPUTS.clear()
    _POST_REINSERTION_SCORE_INPUTS.clear()
    _POST_REINSERTION_METAL_STATE_INPUTS.clear()
    _omol_score_cached.cache_clear()
    _organic_core_score_cached.cache_clear()
    _post_reinsertion_score_cached.cache_clear()
    _post_reinsertion_metal_state_score_cached.cache_clear()


def organic_core_score_cache_info() -> Tuple[int, int, int]:
    info = cast(CacheInfoLike, _organic_core_score_cached.cache_info())
    return info.hits, info.misses, info.currsize


def post_reinsertion_score_cache_info() -> Tuple[int, int, int]:
    direct_info = cast(CacheInfoLike, _post_reinsertion_score_cached.cache_info())
    metal_state_info = cast(CacheInfoLike, _post_reinsertion_metal_state_score_cached.cache_info())
    return (
        direct_info.hits + metal_state_info.hits,
        direct_info.misses + metal_state_info.misses,
        direct_info.currsize + metal_state_info.currsize,
    )


def organic_core_score(omol_or_state: ScoringInput) -> float:
    if isinstance(omol_or_state, ReconstructionState):
        return omol_or_state.organic_core_score()
    omol, score_key = _resolve_cached_omol_value(
        omol_or_state,
        "organic_score_key",
        _build_score_key,
    )
    return organic_core_score_from_key(score_key, omol)


def organic_core_score_from_key(score_key: ScoreKey, omol: pybel.Molecule) -> float:
    _ORGANIC_SCORE_INPUTS[score_key] = omol
    try:
        return _organic_core_score_cached(score_key)
    finally:
        _ORGANIC_SCORE_INPUTS.pop(score_key, None)


def post_reinsertion_score_from_base_key(
    post_reinsertion_base_key: ScoreKey,
    omol: pybel.Molecule,
) -> float:
    score_key = (post_reinsertion_base_key, _EMPTY_METAL_STATE_KEY)
    _POST_REINSERTION_SCORE_INPUTS[score_key] = omol
    try:
        return _post_reinsertion_score_cached(score_key)
    finally:
        _POST_REINSERTION_SCORE_INPUTS.pop(score_key, None)


def build_post_reinsertion_base_components(
    omol: pybel.Molecule,
) -> Tuple[float, Tuple[ChargedAtomSnapshot, ...]]:
    obmol = cast(ob.OBMol, omol.OBMol)
    charged_atom_snapshots: List[ChargedAtomSnapshot] = []
    for atom in ob.OBMolAtomIter(obmol):
        obatom = cast(ob.OBAtom, atom)
        charge = cast(int, obatom.GetFormalCharge())
        if charge == 0:
            continue
        charged_atom_snapshots.append((charge, obatom.GetX(), obatom.GetY(), obatom.GetZ()))
    return calc_symmetry_penalty(omol), tuple(charged_atom_snapshots)


def post_reinsertion_score_from_metal_states(
    post_reinsertion_base_key: ScoreKey,
    base_symmetry_penalty: float,
    charged_atom_snapshots: Sequence[ChargedAtomSnapshot],
    metal_states: Sequence[MetalAtomPosition],
) -> float:
    score_key = (post_reinsertion_base_key, build_metal_state_key(metal_states))
    _POST_REINSERTION_METAL_STATE_INPUTS[score_key] = (
        base_symmetry_penalty,
        tuple(charged_atom_snapshots),
        tuple(metal_states),
    )
    try:
        return _post_reinsertion_metal_state_score_cached(score_key)
    finally:
        _POST_REINSERTION_METAL_STATE_INPUTS.pop(score_key, None)


def combined_candidate_score_from_metal_states(
    organic_score: float,
    post_reinsertion_base_key: ScoreKey,
    base_symmetry_penalty: float,
    charged_atom_snapshots: Sequence[ChargedAtomSnapshot],
    metal_states: Sequence[MetalAtomPosition],
) -> float:
    return organic_score + post_reinsertion_score_from_metal_states(
        post_reinsertion_base_key,
        base_symmetry_penalty,
        charged_atom_snapshots,
        metal_states,
    )


def omol_score_from_parts(
    organic_score: float,
    post_reinsertion_base_key: ScoreKey,
    omol: pybel.Molecule,
) -> float:
    return organic_score + post_reinsertion_score_from_base_key(post_reinsertion_base_key, omol)


def calculate_tetrahedron_volume(
    p1: Sequence[float], p2: Sequence[float], p3: Sequence[float], p4: Sequence[float]
) -> float:
    a = [p1[0] - p4[0], p1[1] - p4[1], p1[2] - p4[2]]
    b = [p2[0] - p4[0], p2[1] - p4[1], p2[2] - p4[2]]
    c = [p3[0] - p4[0], p3[1] - p4[1], p3[2] - p4[2]]
    det = (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )
    return abs(det) / 6.0


def calculate_shape_quality(
    p1: Sequence[float], p2: Sequence[float], p3: Sequence[float], p4: Sequence[float]
) -> float:
    volume = calculate_tetrahedron_volume(p1, p2, p3, p4)
    if math.isclose(volume, 0.0):
        return 0.0

    def d2(a: Sequence[float], b: Sequence[float]) -> float:
        return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2

    edges_sq = [d2(p1, p2), d2(p1, p3), d2(p1, p4), d2(p2, p3), d2(p2, p4), d2(p3, p4)]
    sum_edges_sq = sum(edges_sq)
    l_rms_cubed = (sum_edges_sq / 6.0) ** 1.5
    if math.isclose(l_rms_cubed, 0.0):
        return 0.0
    quality = (6.0 * math.sqrt(2.0)) * (volume / l_rms_cubed)
    return max(0.0, min(1.0, quality))


def get_deviation_score(omol: pybel.Molecule, atom_idx: int) -> float:
    def _get_coords(a: ob.OBAtom) -> Sequence[float]:
        vec = cast(ob.vector3, a.GetVector())
        return [vec.GetX(), vec.GetY(), vec.GetZ()]

    obmol = cast(ob.OBMol, omol.OBMol)
    atom = cast(ob.OBAtom, obmol.GetAtom(atom_idx))
    neighbor_atoms: List[ob.OBAtom] = list(ob.OBAtomAtomIter(atom))
    if len(neighbor_atoms) == 2:
        angle = obmol.GetAngle(neighbor_atoms[0], atom, neighbor_atoms[1])
        return abs(angle - 108) / 108.0
    if len(neighbor_atoms) == 3:
        quality = calculate_shape_quality(
            _get_coords(neighbor_atoms[0]),
            _get_coords(neighbor_atoms[1]),
            _get_coords(neighbor_atoms[2]),
            _get_coords(atom),
        )
        return 1.0 - quality
    return 0.0


def calc_remain_bond_order_reward(omol: pybel.Molecule) -> float:
    obmol = cast(ob.OBMol, omol.OBMol)
    return (
        sum(
            consts.NON_METAL_DICT[atom.GetAtomicNum()].default_valence
            for atom in ob.OBMolAtomIter(obmol)
            if not atom.IsMetal()
        )
        - sum(cast(int, bond.GetBondOrder()) for bond in ob.OBMolBondIter(obmol)) * 2
    ) * 5


def calc_symmetry_penalty(omol: pybel.Molecule) -> float:
    obmol = cast(ob.OBMol, omol.OBMol)
    gs = ob.OBGraphSym(obmol)
    symmetry_ids_vec = ob.vectorUnsignedInt()
    gs.GetSymmetry(symmetry_ids_vec)
    return (len(set(symmetry_ids_vec)) - obmol.NumAtoms()) * 2.0


def _metal_state_symmetry_penalty(metal_states: Sequence[MetalAtomPosition]) -> float:
    return (
        len(
            {
                (metal_state.element_idx, metal_state.valence, metal_state.radical_num)
                for metal_state in metal_states
            }
        )
        - len(metal_states)
    ) * 2.0


def calculate_charge_penalty(atom: ob.OBAtom) -> float:
    return _calculate_charge_penalty_from_data(
        cast(int, atom.GetAtomicNum()),
        cast(int, atom.GetFormalCharge()),
        cast(int, atom.GetTotalValence()),
        cast(int, atom.GetSpinMultiplicity()),
    )


def _calculate_charge_penalty_from_data(
    atomic_num: int,
    charge: int,
    total_valence: int,
    spin_multiplicity: int,
) -> float:
    if charge == 0:
        return 0.0
    total_electrons = cast(
        int,
        consts.NON_METAL_DICT[atomic_num].num_outer_electrons
        + total_valence
        - charge
        + spin_multiplicity,
    )
    if total_electrons == 8 or total_electrons == 2:
        return 0.0
    en = cast(float, ob.GetElectroNeg(atomic_num))
    if charge > 0:
        return abs(charge) * max(0.0, en - 2) * 3.0
    return abs(charge) * max(0.0, 4 - en) * 3.0


def calculate_radical_penalty(atom: ob.OBAtom) -> float:
    return _calculate_radical_penalty_from_data(
        cast(int, atom.GetAtomicNum()),
        cast(int, atom.GetSpinMultiplicity()),
        cast(int, atom.GetHvyDegree()),
    )


def _calculate_radical_penalty_from_data(
    atomic_num: int,
    radical_num: int,
    heavy_degree: int,
) -> float:
    if radical_num == 0:
        return 0.0
    if atomic_num in consts.HETEROATOM:
        return radical_num * 10.0
    return (3 - heavy_degree) * 1.5


def calculate_coulombic_penalty(bond: ob.OBBond) -> float:
    q1 = cast(int, cast(ob.OBAtom, bond.GetBeginAtom()).GetFormalCharge())
    q2 = cast(int, cast(ob.OBAtom, bond.GetEndAtom()).GetFormalCharge())
    if q1 == 0 or q2 == 0:
        return 0.0
    if q1 * q2 > 0:
        return 15.0
    return -0.5


def calculate_physchem_penalty(omol: pybel.Molecule) -> float:
    total_penalty = 0.0
    obmol = cast(ob.OBMol, omol.OBMol)
    for atom in ob.OBMolAtomIter(obmol):
        obatom = cast(ob.OBAtom, atom)
        if obatom.IsMetal():
            continue
        total_penalty += calculate_charge_penalty(obatom)
        total_penalty += calculate_radical_penalty(obatom)
    for bond in ob.OBMolBondIter(obmol):
        total_penalty += calculate_coulombic_penalty(bond)
    return total_penalty


def get_metal_coordination_sphere(
    omol: pybel.Molecule, metal_atom: ob.OBAtom, cutoff: float = 3
) -> List[Tuple[ob.OBAtom, float]]:
    neighbors: List[Tuple[ob.OBAtom, float]] = []
    metal_idx = metal_atom.GetIdx()
    obmol = cast(ob.OBMol, omol.OBMol)
    for atom in ob.OBMolAtomIter(obmol):
        obatom = cast(ob.OBAtom, atom)
        if obatom.GetIdx() == metal_idx:
            continue
        dist = metal_atom.GetDistance(obatom)
        if dist <= cutoff:
            neighbors.append((obatom, dist))
    return neighbors


def calculate_metal_penalty(omol: pybel.Molecule) -> float:
    penalty = 0.0
    obmol = cast(ob.OBMol, omol.OBMol)
    for atom in ob.OBMolAtomIter(obmol):
        metal_atom = cast(ob.OBAtom, atom)
        if not metal_atom.IsMetal():
            continue
        symbol = cast(str, ob.GetSymbol(metal_atom.GetAtomicNum()))
        valence = cast(int, metal_atom.GetFormalCharge())
        if valence <= 0:
            penalty += 10 * max(abs(valence), 1)
        prior_list = consts.METAL_VALENCE_AVAILABLE_PRIOR.get(symbol, [])
        minor_list = consts.METAL_VALENCE_AVAILABLE_MINOR.get(symbol, [])
        if valence not in prior_list:
            if valence in minor_list:
                penalty += 10.0
            else:
                penalty += 20.0
        neighbors = get_metal_coordination_sphere(omol, metal_atom, cutoff=2.6)
        for ligand_atom, dist in neighbors:
            ligand_charge = cast(int, ligand_atom.GetFormalCharge())
            if valence > 0:
                if ligand_charge > 0:
                    penalty += 100.0 * (ligand_charge * valence) / (dist**2)
                elif ligand_charge < 0:
                    penalty -= 5 * (abs(ligand_charge) * valence) / (dist**2)
    return penalty


def _calculate_charge_interaction_penalty(
    metal_valence: int,
    ligand_charge: int,
    dist_sq: float,
) -> float:
    if ligand_charge > 0:
        return 100.0 * (ligand_charge * metal_valence) / dist_sq
    if ligand_charge < 0:
        return -5.0 * (abs(ligand_charge) * metal_valence) / dist_sq
    return 0.0


def calculate_metal_penalty_from_metal_states(
    charged_atom_snapshots: Sequence[ChargedAtomSnapshot],
    metal_states: Sequence[MetalAtomPosition],
    cutoff: float = 2.6,
) -> float:
    penalty = 0.0
    cutoff_sq = cutoff * cutoff

    for metal_state in metal_states:
        valence = metal_state.valence
        if valence <= 0:
            penalty += 10 * max(abs(valence), 1)

        prior_list = consts.METAL_VALENCE_AVAILABLE_PRIOR.get(metal_state.symbol, [])
        minor_list = consts.METAL_VALENCE_AVAILABLE_MINOR.get(metal_state.symbol, [])
        if valence not in prior_list:
            if valence in minor_list:
                penalty += 10.0
            else:
                penalty += 20.0

        if valence <= 0:
            continue

        mx = metal_state.position_x
        my = metal_state.position_y
        mz = metal_state.position_z
        for ligand_charge, x, y, z in charged_atom_snapshots:
            dx = mx - x
            dy = my - y
            dz = mz - z
            dist_sq = dx * dx + dy * dy + dz * dz
            if dist_sq > cutoff_sq:
                continue
            penalty += _calculate_charge_interaction_penalty(valence, ligand_charge, dist_sq)

        for other_metal_state in metal_states:
            if other_metal_state is metal_state:
                continue
            ligand_charge = other_metal_state.valence
            if ligand_charge == 0:
                continue
            dx = mx - other_metal_state.position_x
            dy = my - other_metal_state.position_y
            dz = mz - other_metal_state.position_z
            dist_sq = dx * dx + dy * dy + dz * dz
            if dist_sq > cutoff_sq:
                continue
            penalty += _calculate_charge_interaction_penalty(valence, ligand_charge, dist_sq)

    return penalty


def _compute_post_reinsertion_score_from_metal_states(
    base_symmetry_penalty: float,
    charged_atom_snapshots: Sequence[ChargedAtomSnapshot],
    metal_states: Sequence[MetalAtomPosition],
) -> float:
    return (
        base_symmetry_penalty
        + _metal_state_symmetry_penalty(metal_states)
        + calculate_metal_penalty_from_metal_states(charged_atom_snapshots, metal_states)
    )


def calculate_heteroatom_penalty(omol: pybel.Molecule) -> float:
    penalty = 0.0
    obmol = cast(ob.OBMol, omol.OBMol)
    for atom in ob.OBMolAtomIter(obmol):
        obatom = cast(ob.OBAtom, atom)
        if obatom.GetAtomicNum() not in consts.HETEROATOM:
            continue
        penalty += 10 * (obatom.GetFormalCharge() - obatom.GetTotalValence())
    return penalty


def calculate_conjugation_reward(omol: pybel.Molecule) -> float:
    res: List[Tuple[int, int, int, int]] = list(smarts.SCORING_CONJUGATION.findall(omol))
    return len(res) * 2.0


def omol_score(omol_or_state: ScoringInput) -> float:
    if isinstance(omol_or_state, ReconstructionState):
        return omol_or_state.full_score()
    omol, score_key = _resolve_cached_omol_value(
        omol_or_state,
        "omol_score_key",
        _build_score_key,
    )
    _OMOL_SCORE_INPUTS[score_key] = omol
    try:
        return _omol_score_cached(score_key)
    finally:
        _OMOL_SCORE_INPUTS.pop(score_key, None)


__all__ = [
    "combined_candidate_score_from_metal_states",
    "build_post_reinsertion_base_key",
    "build_post_reinsertion_base_components",
    "omol_score",
    "omol_score_from_parts",
    "omol_score_cache_clear",
    "omol_score_cache_info",
    "organic_core_score",
    "organic_core_score_from_key",
    "organic_core_score_cache_info",
    "post_reinsertion_score_from_metal_states",
    "post_reinsertion_score_from_base_key",
    "post_reinsertion_score_cache_info",
]
