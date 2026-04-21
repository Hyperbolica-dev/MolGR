"""Metal-aware reconstruction for fallback.

The production path is:
1. Strip metals and enumerate a small set of metal valence/radical assignments.
2. Group assignments by the no-metal target they induce via meet-in-the-middle DP.
3. Reconstruct each no-metal target once, reuse that state across the bucket, and
   only materialize the winning combined molecule at the end.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from typing import DefaultDict, Dict, List, Optional, Set, Tuple, cast

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.pipeline.reconstruct_without_metals import xyz_to_omol_no_metal_state
from molgr.fallback.state import (
    MetalCandidateState,
    MetalCandidateStateMachine,
    MetalPreparationState,
    ReconstructionState,
)
from molgr.fallback.utils import consts, dataclasses


_ValenceBoundsKey = Tuple[Tuple[str, int, int], ...]
_ChargeGroupedAssignments = Dict[
    int,
    List[Tuple[_ValenceBoundsKey, List["_PartialMetalAssignment"]]],
]
_RadicalBucketLevel = Tuple[int, _ChargeGroupedAssignments]
_RadicalPrefixReachability = Tuple[Tuple[int, ...], Tuple[_RadicalBucketLevel, ...]]


@dataclass(frozen=True)
class _PartialMetalAssignment:
    metal_states: Tuple[dataclasses.MetalAtomPosition, ...]
    total_metal_charge: int
    total_metal_radicals: int
    metal_assignment_rank: float
    valence_bounds: _ValenceBoundsKey
    order: int


def get_possible_metal_radicals(metal: str, valence: int) -> Set[int]:
    """Return allowed radical counts for a metal under a candidate valence."""

    f_d_s_p = consts.METAL_F_D_S_P_ELECTRONS.get(metal)
    if f_d_s_p is None:
        return set()

    f, d, s, p = f_d_s_p.f, f_d_s_p.d, f_d_s_p.s, f_d_s_p.p

    if valence <= s + p:
        if 0 <= d < len(consts.D_ELECTRONS_SPIN):
            base = (f + s + p - valence) % 2
            return {base + dd for dd in consts.D_ELECTRONS_SPIN[d]}
        return set()

    if valence <= s + p + d:
        idx = d - valence + s + p
        if 0 <= idx < len(consts.D_ELECTRONS_SPIN):
            return {f % 2 + dd for dd in consts.D_ELECTRONS_SPIN[idx]}
        return set()

    if valence <= s + p + d + f:
        return {f % 2}

    return set()


def _build_metal_states(obatom: ob.OBAtom) -> List[dataclasses.MetalAtomPosition]:
    symbol = ob.GetSymbol(obatom.GetAtomicNum())

    def _default_state() -> dataclasses.MetalAtomPosition:
        return dataclasses.MetalAtomPosition(
            idx=obatom.GetIdx(),
            symbol=symbol,
            element_idx=obatom.GetAtomicNum(),
            valence=0,
            radical_num=0,
            position_x=obatom.GetX(),
            position_y=obatom.GetY(),
            position_z=obatom.GetZ(),
        )

    prior = consts.METAL_VALENCE_AVAILABLE_PRIOR.get(symbol, [])
    minor = consts.METAL_VALENCE_AVAILABLE_MINOR.get(symbol, [])
    seen_valences = set()
    valences: List[int] = []
    for valence in prior + minor:
        if valence in seen_valences:
            continue
        seen_valences.add(valence)
        valences.append(valence)
    if not valences:
        valences = [0]

    if symbol not in consts.METAL_F_D_S_P_ELECTRONS:
        return [_default_state()]

    states: List[dataclasses.MetalAtomPosition] = []
    for valence in valences:
        try:
            radicals = get_possible_metal_radicals(symbol, valence)
        except ValueError:
            continue
        for radical_num in sorted(radicals):
            states.append(
                dataclasses.MetalAtomPosition(
                    idx=obatom.GetIdx(),
                    symbol=symbol,
                    element_idx=obatom.GetAtomicNum(),
                    valence=valence,
                    radical_num=radical_num,
                    position_x=obatom.GetX(),
                    position_y=obatom.GetY(),
                    position_z=obatom.GetZ(),
                )
            )

    if not states:
        return [_default_state()]
    return states


def combine_metal_with_omol(
    omol: pybel.Molecule, metal_list: Sequence[dataclasses.MetalAtomPosition]
) -> pybel.Molecule:
    """Insert the selected metal states back into the no-metal winner."""

    obmol = cast(ob.OBMol, omol.clone.OBMol)
    obmol.BeginModify()
    try:
        num_organic = obmol.NumAtoms()
        num_metals = len(metal_list)
        total_atoms = num_organic + num_metals

        for metal in metal_list:
            atom = cast(ob.OBAtom, obmol.NewAtom())
            atom.SetAtomicNum(metal.element_idx)
            atom.SetFormalCharge(metal.valence)
            atom.SetSpinMultiplicity(metal.radical_num)
            atom.SetVector(metal.position_x, metal.position_y, metal.position_z)

        new_order = [0] * total_atoms
        has_error = False
        for i, metal in enumerate(metal_list):
            current_idx = num_organic + 1 + i
            target_slot = metal.idx - 1
            if target_slot < 0 or target_slot >= total_atoms:
                has_error = True
                continue
            if new_order[target_slot] != 0:
                has_error = True
                continue
            new_order[target_slot] = current_idx

        if not has_error:
            current_organic_idx = 1
            for i in range(total_atoms):
                if new_order[i] != 0:
                    continue
                if current_organic_idx > num_organic:
                    has_error = True
                    break
                new_order[i] = current_organic_idx
                current_organic_idx += 1

        if not has_error and all(idx > 0 for idx in new_order):
            obmol.RenumberAtoms(new_order)
    finally:
        obmol.EndModify()
    return pybel.Molecule(obmol)


def prepare_metal_state(
    xyz_block: str,
    total_charge: int,
    total_radical_electrons: int,
) -> MetalPreparationState:
    """Split the input into a no-metal XYZ block plus per-metal state options."""

    omol = pybel.readstring("xyz", xyz_block)
    removable_metal_atoms = [cast(ob.OBAtom, atom.OBAtom) for atom in omol.atoms if atom.OBAtom.IsMetal()]
    available_valence_radical_states = tuple(
        tuple(_build_metal_states(obatom)) for obatom in removable_metal_atoms
    )
    for obatom in removable_metal_atoms:
        omol.OBMol.DeleteAtom(obatom)
    state = MetalPreparationState(
        no_metal_xyz_block=cast(str, omol.write("xyz")),
        available_valence_radical_states=available_valence_radical_states,
        total_charge=total_charge,
        total_radical_electrons=total_radical_electrons,
        phase_history=(
            "read_xyz",
            "build_metal_state_options",
            "remove_metal_atoms",
            "serialize_no_metal_xyz",
        ),
        metadata={"metal_atom_count": len(removable_metal_atoms)},
    )
    return state


def _score_candidate_with_no_metal_state(
    candidate: MetalCandidateState,
    no_metal_state: ReconstructionState,
) -> MetalCandidateState:
    """Attach the shared no-metal reconstruction and score the metal candidate."""

    candidate_machine = MetalCandidateStateMachine.from_candidate_state(candidate)
    candidate_machine.set_no_metal_state("reconstruct_no_metal", no_metal_state)
    candidate_machine.annotate("score_candidate")
    scored_candidate = candidate_machine.freeze()
    scored_candidate.combined_score()
    return scored_candidate


def _metal_state_assignment_penalty(metal_state: dataclasses.MetalAtomPosition) -> float:
    penalty = 0.0
    if metal_state.valence <= 0:
        penalty += 10.0 * max(abs(metal_state.valence), 1)

    prior_list = consts.METAL_VALENCE_AVAILABLE_PRIOR.get(metal_state.symbol, [])
    minor_list = consts.METAL_VALENCE_AVAILABLE_MINOR.get(metal_state.symbol, [])
    if metal_state.valence not in prior_list:
        if metal_state.valence in minor_list:
            penalty += 10.0
        else:
            penalty += 20.0
    return penalty


def _update_valence_bounds(
    bounds: _ValenceBoundsKey,
    metal_state: dataclasses.MetalAtomPosition,
    max_mixed_valence_spread: Optional[int],
) -> Optional[_ValenceBoundsKey]:
    if max_mixed_valence_spread is None or max_mixed_valence_spread < 0:
        return ()

    updated_bounds = list(bounds)
    for idx, (symbol, lower, upper) in enumerate(updated_bounds):
        if symbol != metal_state.symbol:
            continue
        next_lower = min(lower, metal_state.valence)
        next_upper = max(upper, metal_state.valence)
        if next_upper - next_lower > max_mixed_valence_spread:
            return None
        updated_bounds[idx] = (symbol, next_lower, next_upper)
        break
    else:
        updated_bounds.append((metal_state.symbol, metal_state.valence, metal_state.valence))
        updated_bounds.sort(key=lambda item: item[0])

    return tuple(updated_bounds)


def _trim_partial_assignments(
    entries: List[_PartialMetalAssignment],
    max_assignments_per_target: int,
) -> List[_PartialMetalAssignment]:
    limit = max(1, max_assignments_per_target)
    if len(entries) <= limit:
        entries.sort(key=lambda entry: (entry.metal_assignment_rank, entry.order))
        return entries
    entries.sort(key=lambda entry: (entry.metal_assignment_rank, entry.order))
    return entries[:limit]


def _enumerate_partial_assignment_frontier(
    available_valence_radical_states: Sequence[Sequence[dataclasses.MetalAtomPosition]],
    *,
    max_mixed_valence_spread: Optional[int] = 3,
    max_total_metal_radicals: Optional[int] = None,
    max_assignments_per_state: int = 64,
) -> Dict[Tuple[int, int, _ValenceBoundsKey], List[_PartialMetalAssignment]]:
    """Enumerate one half of the metal search space while pruning dominated prefixes."""

    partial_assignments: Dict[Tuple[int, int, _ValenceBoundsKey], List[_PartialMetalAssignment]] = {
        (0, 0, ()): [
            _PartialMetalAssignment(
                metal_states=(),
                total_metal_charge=0,
                total_metal_radicals=0,
                metal_assignment_rank=0.0,
                valence_bounds=(),
                order=0,
            )
        ]
    }
    next_order = 1

    for metal_state_options in available_valence_radical_states:
        next_partial_assignments: Dict[
            Tuple[int, int, _ValenceBoundsKey], List[_PartialMetalAssignment]
        ] = defaultdict(list)
        for entries in partial_assignments.values():
            for entry in entries:
                for metal_state in metal_state_options:
                    next_total_metal_radicals = entry.total_metal_radicals + metal_state.radical_num
                    if (
                        max_total_metal_radicals is not None
                        and next_total_metal_radicals > max_total_metal_radicals
                    ):
                        continue

                    next_valence_bounds = _update_valence_bounds(
                        entry.valence_bounds,
                        metal_state,
                        max_mixed_valence_spread,
                    )
                    if next_valence_bounds is None:
                        continue

                    next_entry = _PartialMetalAssignment(
                        metal_states=entry.metal_states + (metal_state,),
                        total_metal_charge=entry.total_metal_charge + metal_state.valence,
                        total_metal_radicals=next_total_metal_radicals,
                        metal_assignment_rank=entry.metal_assignment_rank
                        + _metal_state_assignment_penalty(metal_state),
                        valence_bounds=next_valence_bounds,
                        order=next_order,
                    )
                    next_order += 1
                    next_key = (
                        next_entry.total_metal_charge,
                        next_entry.total_metal_radicals,
                        next_entry.valence_bounds,
                    )
                    next_partial_assignments[next_key].append(next_entry)

        partial_assignments = {
            key: _trim_partial_assignments(entries, max_assignments_per_state)
            for key, entries in next_partial_assignments.items()
        }
        if not partial_assignments:
            break

    return partial_assignments


def _merge_valence_bounds(
    left_bounds: _ValenceBoundsKey,
    right_bounds: _ValenceBoundsKey,
    max_mixed_valence_spread: Optional[int],
) -> Optional[_ValenceBoundsKey]:
    if max_mixed_valence_spread is None or max_mixed_valence_spread < 0:
        return ()

    merged_bounds: Dict[str, Tuple[int, int]] = {
        symbol: (lower, upper) for symbol, lower, upper in left_bounds
    }
    for symbol, lower, upper in right_bounds:
        previous_bounds = merged_bounds.get(symbol)
        if previous_bounds is None:
            merged_bounds[symbol] = (lower, upper)
            continue
        next_lower = min(previous_bounds[0], lower)
        next_upper = max(previous_bounds[1], upper)
        if next_upper - next_lower > max_mixed_valence_spread:
            return None
        merged_bounds[symbol] = (next_lower, next_upper)

    return tuple(
        sorted(
            (symbol, lower, upper)
            for symbol, (lower, upper) in merged_bounds.items()
        )
    )


def _bucket_partial_assignments_by_charge_radicals(
    frontier: Dict[Tuple[int, int, _ValenceBoundsKey], List[_PartialMetalAssignment]],
) -> Dict[int, _ChargeGroupedAssignments]:
    bucket_index: DefaultDict[
        int,
        DefaultDict[int, List[Tuple[_ValenceBoundsKey, List[_PartialMetalAssignment]]]],
    ] = defaultdict(lambda: defaultdict(list))
    for (charge, radicals, valence_bounds), entries in frontier.items():
        bucket_index[radicals][charge].append((valence_bounds, entries))

    return {
        radicals: {charge: charge_groups[charge] for charge in sorted(charge_groups)}
        for radicals, charge_groups in sorted(bucket_index.items())
    }


def _build_radical_prefix_reachability(
    bucket_index: Dict[int, _ChargeGroupedAssignments],
) -> _RadicalPrefixReachability:
    radical_levels = tuple(sorted(bucket_index))
    bucket_levels = tuple((radicals, bucket_index[radicals]) for radicals in radical_levels)
    return radical_levels, bucket_levels


def _iter_reachable_radical_buckets(
    prefix_reachability: _RadicalPrefixReachability,
    max_radicals: int,
) -> Sequence[_RadicalBucketLevel]:
    radical_levels, bucket_levels = prefix_reachability
    prefix_end = bisect_right(radical_levels, max_radicals)
    if prefix_end <= 0:
        return ()
    return bucket_levels[:prefix_end]


def _combine_partial_assignment_frontiers(
    left_frontier: Dict[Tuple[int, int, _ValenceBoundsKey], List[_PartialMetalAssignment]],
    right_frontier: Dict[Tuple[int, int, _ValenceBoundsKey], List[_PartialMetalAssignment]],
    total_charge: int,
    total_radical_electrons: int,
    *,
    max_mixed_valence_spread: Optional[int] = 3,
    max_total_metal_radicals: Optional[int] = None,
    max_assignments_per_target: int = 64,
) -> Dict[Tuple[int, int], List[_PartialMetalAssignment]]:
    """Join the left/right DP frontiers by compatible charge and radical buckets."""

    grouped_entries: Dict[Tuple[int, int], List[_PartialMetalAssignment]] = defaultdict(list)
    next_order = 0
    trim_trigger = max(1, max_assignments_per_target) * 4
    left_bucket_index = _bucket_partial_assignments_by_charge_radicals(left_frontier)
    right_bucket_index = _bucket_partial_assignments_by_charge_radicals(right_frontier)
    right_prefix_reachability = _build_radical_prefix_reachability(right_bucket_index)
    max_combined_metal_radicals = total_radical_electrons
    if max_total_metal_radicals is not None:
        max_combined_metal_radicals = min(max_combined_metal_radicals, max_total_metal_radicals)

    for left_radicals, left_charge_groups in left_bucket_index.items():
        max_right_radicals = max_combined_metal_radicals - left_radicals
        if max_right_radicals < 0:
            continue

        for right_radicals, right_charge_groups in _iter_reachable_radical_buckets(
            right_prefix_reachability,
            max_right_radicals,
        ):
            total_metal_radicals = left_radicals + right_radicals
            target_radicals = total_radical_electrons - total_metal_radicals
            for left_charge, left_valence_groups in left_charge_groups.items():
                for right_charge, right_valence_groups in right_charge_groups.items():
                    target = (
                        total_charge - (left_charge + right_charge),
                        target_radicals,
                    )
                    bucket = grouped_entries[target]
                    initial_bucket_size = len(bucket)
                    for left_bounds, left_entries in left_valence_groups:
                        for right_bounds, right_entries in right_valence_groups:
                            merged_bounds = _merge_valence_bounds(
                                left_bounds,
                                right_bounds,
                                max_mixed_valence_spread,
                            )
                            if merged_bounds is None:
                                continue

                            for left_entry in left_entries:
                                for right_entry in right_entries:
                                    bucket.append(
                                        _PartialMetalAssignment(
                                            metal_states=left_entry.metal_states + right_entry.metal_states,
                                            total_metal_charge=left_entry.total_metal_charge
                                            + right_entry.total_metal_charge,
                                            total_metal_radicals=left_entry.total_metal_radicals
                                            + right_entry.total_metal_radicals,
                                            metal_assignment_rank=left_entry.metal_assignment_rank
                                            + right_entry.metal_assignment_rank,
                                            valence_bounds=merged_bounds,
                                            order=next_order,
                                        )
                                    )
                                    next_order += 1

                    if len(bucket) == initial_bucket_size:
                        if initial_bucket_size == 0:
                            grouped_entries.pop(target, None)
                        continue

                    if len(bucket) > trim_trigger:
                        grouped_entries[target] = _trim_partial_assignments(
                            bucket,
                            max_assignments_per_target,
                        )

    return {
        target: _trim_partial_assignments(entries, max_assignments_per_target)
        for target, entries in grouped_entries.items()
    }


def _group_candidates_by_target_dp(
    base_phase_history: Sequence[str],
    available_valence_radical_states: Sequence[Sequence[dataclasses.MetalAtomPosition]],
    total_charge: int,
    total_radical_electrons: int,
    *,
    max_mixed_valence_spread: Optional[int] = 3,
    max_total_metal_radicals: Optional[int] = None,
    max_assignments_per_target: int = 64,
) -> Dict[Tuple[int, int], List[MetalCandidateState]]:
    """Build target buckets directly, without enumerating the full metal Cartesian product."""

    split_index = len(available_valence_radical_states) // 2
    left_frontier = _enumerate_partial_assignment_frontier(
        available_valence_radical_states[:split_index],
        max_mixed_valence_spread=max_mixed_valence_spread,
        max_total_metal_radicals=max_total_metal_radicals,
        max_assignments_per_state=max_assignments_per_target,
    )
    right_frontier = _enumerate_partial_assignment_frontier(
        available_valence_radical_states[split_index:],
        max_mixed_valence_spread=max_mixed_valence_spread,
        max_total_metal_radicals=max_total_metal_radicals,
        max_assignments_per_state=max_assignments_per_target,
    )
    if not left_frontier or not right_frontier:
        return {}

    grouped_entries = _combine_partial_assignment_frontiers(
        left_frontier,
        right_frontier,
        total_charge,
        total_radical_electrons,
        max_mixed_valence_spread=max_mixed_valence_spread,
        max_total_metal_radicals=max_total_metal_radicals,
        max_assignments_per_target=max_assignments_per_target,
    )
    grouped_candidates: Dict[Tuple[int, int], List[MetalCandidateState]] = {}
    combination_index = 0
    for target, entries in grouped_entries.items():
        bucket: List[MetalCandidateState] = []
        for entry in entries:
            candidate_machine = MetalCandidateStateMachine(
                entry.metal_states,
                target[0],
                target[1],
                phase_history=base_phase_history,
                metadata={"combination_index": combination_index},
            )
            candidate_machine.annotate("enumerate_metal_combination")
            candidate_machine.annotate("reconstruct_no_metal_candidate")
            candidate_machine.annotate(
                "rank_metal_assignment_for_target",
                metal_assignment_rank=entry.metal_assignment_rank,
            )
            candidate = candidate_machine.freeze()
            bucket.append(candidate)
            combination_index += 1
        grouped_candidates[target] = bucket
    return grouped_candidates


def xyz2omol_state(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
    *,
    max_mixed_valence_spread: Optional[int] = 3,
    max_total_metal_radicals: Optional[int] = 0,
    max_assignments_per_target: int = 64,
) -> Optional[MetalCandidateState]:
    """Return the best scored metal candidate state for the input XYZ block."""

    base_state = prepare_metal_state(xyz_block, total_charge, total_radical_electrons)
    possible_candidates: List[MetalCandidateState] = []
    grouped_candidates = _group_candidates_by_target_dp(
        base_state.phase_history,
        base_state.available_valence_radical_states,
        total_charge,
        total_radical_electrons,
        max_mixed_valence_spread=max_mixed_valence_spread,
        max_total_metal_radicals=(
            total_radical_electrons if max_total_metal_radicals == 0 else max_total_metal_radicals
        ),
        max_assignments_per_target=max_assignments_per_target,
    )

    for _, candidates in grouped_candidates.items():
        if not candidates:
            continue
        prototype = candidates[0]
        try:
            no_metal_state = xyz_to_omol_no_metal_state(
                base_state.no_metal_xyz_block,
                prototype.no_metal_charge_target,
                prototype.no_metal_radical_target,
            )
        except (OSError, ValueError):
            continue
        if no_metal_state is None:
            continue
        possible_candidates.extend(
            _score_candidate_with_no_metal_state(candidate, no_metal_state)
            for candidate in candidates
        )
    if not possible_candidates:
        return None
    possible_candidates.sort(key=lambda candidate: cast(float, candidate.score))
    best_candidate = possible_candidates[0]
    if best_candidate.combined_omol is None:
        best_candidate.materialize_combined_omol(combine_metal_with_omol)
        winner_machine = MetalCandidateStateMachine.from_candidate_state(best_candidate)
        winner_machine.annotate("combine_metal_with_omol")
        best_candidate = winner_machine.freeze()
    winner_machine = MetalCandidateStateMachine.from_candidate_state(best_candidate)
    winner_machine.annotate("select_best_candidate")
    return winner_machine.freeze()


def xyz2omol(
    xyz_block: str,
    total_charge: int = 0,
    total_radical_electrons: int = 0,
    *,
    max_mixed_valence_spread: Optional[int] = 3,
    max_total_metal_radicals: Optional[int] = 0,
    max_assignments_per_target: int = 64,
) -> Optional[pybel.Molecule]:
    """Materialize the winning metal-aware reconstruction."""

    candidate = xyz2omol_state(
        xyz_block,
        total_charge,
        total_radical_electrons,
        max_mixed_valence_spread=max_mixed_valence_spread,
        max_total_metal_radicals=max_total_metal_radicals,
        max_assignments_per_target=max_assignments_per_target,
    )
    if candidate is None:
        return None
    return candidate.combined_omol


__all__ = [
    "_build_metal_states",
    "_group_candidates_by_target_dp",
    "combine_metal_with_omol",
    "get_possible_metal_radicals",
    "prepare_metal_state",
    "xyz2omol",
    "xyz2omol_state",
]
