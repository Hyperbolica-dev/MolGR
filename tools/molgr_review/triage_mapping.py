# pyright: reportCallIssue=false, reportArgumentType=false
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from typing import Any, Callable, Optional

from rdkit import Chem
from rdkit.Chem import rdFMCS


METAL_ATOMIC_NUMBERS = {
    3,
    4,
    11,
    12,
    13,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    37,
    38,
    39,
    40,
    41,
    42,
    43,
    44,
    45,
    46,
    47,
    48,
    49,
    50,
    55,
    56,
    57,
    58,
    59,
    60,
    61,
    62,
    63,
    64,
    65,
    66,
    67,
    68,
    69,
    70,
    71,
    72,
    73,
    74,
    75,
    76,
    77,
    78,
    79,
    80,
    81,
    82,
    83,
    87,
    88,
    89,
    90,
    91,
    92,
    93,
    94,
    95,
    96,
    97,
    98,
    99,
    100,
    101,
    102,
    103,
}


def is_metal(atom: Chem.Atom) -> bool:
    return atom.GetAtomicNum() in METAL_ATOMIC_NUMBERS


def bond_kind(bond: Chem.Bond) -> str:
    if bond.GetBondType() in {
        Chem.BondType.DATIVE,
        Chem.BondType.DATIVEL,
        Chem.BondType.DATIVER,
        Chem.BondType.DATIVEONE,
    }:
        return "dative"
    return {
        Chem.BondType.SINGLE: "single",
        Chem.BondType.DOUBLE: "double",
        Chem.BondType.TRIPLE: "triple",
        Chem.BondType.AROMATIC: "aromatic",
    }.get(bond.GetBondType(), str(bond.GetBondType()).lower())


def atom_h_count(atom: Chem.Atom) -> int:
    attached = sum(neighbour.GetAtomicNum() == 1 for neighbour in atom.GetNeighbors())
    try:
        implicit = atom.GetNumImplicitHs()
    except RuntimeError:
        implicit = 0
    return attached + atom.GetNumExplicitHs() + implicit


def _prepare_heavy(
    mol: Chem.Mol,
) -> tuple[Chem.Mol, dict[int, int], dict[int, int], dict[int, tuple[int, ...]]]:
    source = Chem.Mol(mol)
    source.UpdatePropertyCache(strict=False)
    h_counts: dict[int, int] = {}
    attached_h: dict[int, tuple[int, ...]] = {}
    for atom in source.GetAtoms():
        atom.SetIntProp("_triage_source_index", atom.GetIdx())
        if atom.GetAtomicNum() == 1:
            continue
        h_counts[atom.GetIdx()] = atom_h_count(atom)
        attached_h[atom.GetIdx()] = tuple(
            sorted(
                neighbour.GetIdx()
                for neighbour in atom.GetNeighbors()
                if neighbour.GetAtomicNum() == 1
            )
        )
    heavy = Chem.RemoveHs(source, sanitize=False)
    heavy.UpdatePropertyCache(strict=False)
    heavy_to_source = {
        atom.GetIdx(): atom.GetIntProp("_triage_source_index") for atom in heavy.GetAtoms()
    }
    return heavy, heavy_to_source, h_counts, attached_h


def _edge_map(mol: Chem.Mol) -> dict[tuple[int, int], str]:
    return {
        tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))): bond_kind(bond)
        for bond in mol.GetBonds()
    }


def _mapping_score(
    candidate: Chem.Mol,
    reference: Chem.Mol,
    mapping: dict[int, int],
    candidate_edges: Optional[dict[tuple[int, int], str]] = None,
    reference_edges: Optional[dict[tuple[int, int], str]] = None,
) -> tuple[int, int, int, int]:
    inverse = {
        reference_index: candidate_index for candidate_index, reference_index in mapping.items()
    }
    candidate_edges = candidate_edges if candidate_edges is not None else _edge_map(candidate)
    reference_edges = reference_edges if reference_edges is not None else _edge_map(reference)
    edge_delta = 0
    kind_delta = 0
    for pair, kind in candidate_edges.items():
        if not all(index in mapping for index in pair):
            continue
        reference_pair = tuple(sorted(mapping[index] for index in pair))
        reference_kind = reference_edges.get(reference_pair)
        edge_delta += reference_kind is None
        kind_delta += reference_kind is not None and reference_kind != kind
    for pair in reference_edges:
        if not all(index in inverse for index in pair):
            continue
        candidate_pair = tuple(sorted(inverse[index] for index in pair))
        edge_delta += candidate_pair not in candidate_edges
    charge_delta = sum(
        candidate.GetAtomWithIdx(candidate_index).GetFormalCharge()
        != reference.GetAtomWithIdx(reference_index).GetFormalCharge()
        for candidate_index, reference_index in mapping.items()
    )
    environment_delta = sum(
        candidate.GetAtomWithIdx(candidate_index).IsInRing()
        != reference.GetAtomWithIdx(reference_index).IsInRing()
        for candidate_index, reference_index in mapping.items()
    )
    return edge_delta, kind_delta, charge_delta, environment_delta


def _candidate_xyz_check(
    candidate: Chem.Mol,
    xyz_atoms: Optional[list[tuple[str, tuple[float, float, float]]]],
    tolerance: float,
) -> tuple[bool, dict[int, int], str]:
    if xyz_atoms is None:
        return False, {}, "xyz_atoms_missing"
    if candidate.GetNumAtoms() != len(xyz_atoms) or candidate.GetNumConformers() == 0:
        return False, {}, "atom_count_or_conformer_mismatch"
    conformer = candidate.GetConformer()
    for atom in candidate.GetAtoms():
        index = atom.GetIdx()
        position = conformer.GetAtomPosition(index)
        xyz_element, xyz_position = xyz_atoms[index]
        delta = math.dist((position.x, position.y, position.z), xyz_position)
        if atom.GetSymbol() != xyz_element or delta > tolerance:
            return False, {}, f"candidate_xyz_mismatch_at_{index}"
    return True, {atom.GetIdx(): atom.GetIdx() for atom in candidate.GetAtoms()}, ""


def parse_xyz_atoms(xyz_block: str) -> list[tuple[str, tuple[float, float, float]]]:
    lines = xyz_block.splitlines()
    if not lines:
        return []
    count = int(lines[0].strip())
    atoms = []
    for line in lines[2 : 2 + count]:
        fields = line.split()
        atoms.append((fields[0], tuple(float(value) for value in fields[1:4])))
    return atoms


def _signature(
    candidate: Chem.Mol,
    reference: Chem.Mol,
    mapping: dict[int, int],
    candidate_to_source: dict[int, int],
    reference_to_source: dict[int, int],
    candidate_to_xyz: dict[int, int],
    candidate_h_counts: dict[int, int],
    reference_h_counts: dict[int, int],
    candidate_attached_h: dict[int, tuple[int, ...]],
) -> tuple[Any, ...]:
    inverse = {
        reference_index: candidate_index for candidate_index, reference_index in mapping.items()
    }
    candidate_edges = _edge_map(candidate)
    reference_edges = _edge_map(reference)
    differences = []

    def xyz_index(candidate_heavy_index: int) -> int:
        return candidate_to_xyz[candidate_to_source[candidate_heavy_index]]

    for pair, candidate_kind in candidate_edges.items():
        if not all(index in mapping for index in pair):
            continue
        reference_pair = tuple(sorted(mapping[index] for index in pair))
        reference_kind = reference_edges.get(reference_pair)
        involves_metal = any(is_metal(candidate.GetAtomWithIdx(index)) for index in pair)
        if reference_kind is None or reference_kind != candidate_kind:
            kind = "metal_bond" if involves_metal else "organic_bond"
            differences.append(
                (
                    kind,
                    tuple(sorted(xyz_index(index) for index in pair)),
                    True,
                    reference_kind is not None,
                    candidate_kind,
                    reference_kind or "none",
                )
            )
    for pair, reference_kind in reference_edges.items():
        if not all(index in inverse for index in pair):
            continue
        candidate_pair = tuple(sorted(inverse[index] for index in pair))
        if candidate_pair in candidate_edges:
            continue
        involves_metal = any(is_metal(reference.GetAtomWithIdx(index)) for index in pair)
        kind = "metal_bond" if involves_metal else "organic_bond"
        differences.append(
            (
                kind,
                tuple(sorted(xyz_index(index) for index in candidate_pair)),
                False,
                True,
                "none",
                reference_kind,
            )
        )

    h_differences = []
    for candidate_index, reference_index in mapping.items():
        candidate_source = candidate_to_source[candidate_index]
        reference_source = reference_to_source[reference_index]
        candidate_h = candidate_h_counts.get(candidate_source, 0)
        reference_h = reference_h_counts.get(reference_source, 0)
        if candidate_h == reference_h:
            continue
        h_differences.append(
            {
                "center": xyz_index(candidate_index),
                "delta": candidate_h - reference_h,
                "attached_h": tuple(
                    candidate_to_xyz[index]
                    for index in candidate_attached_h.get(candidate_source, ())
                    if index in candidate_to_xyz
                ),
            }
        )
    surplus = [item for item in h_differences if item["delta"] > 0]
    deficit = [item for item in h_differences if item["delta"] < 0]
    for item in surplus:
        differences.append(
            (
                "hydrogen_assignment",
                item["attached_h"],
                item["center"],
                tuple(sorted(other["center"] for other in deficit)),
                "candidate",
                "reference",
            )
        )
    if not surplus and deficit:
        differences.extend(
            (
                "hydrogen_assignment",
                (),
                None,
                (item["center"],),
                "candidate",
                "reference",
            )
            for item in deficit
        )
    return tuple(sorted(differences, key=lambda item: json.dumps(item, sort_keys=True)))


@dataclass(frozen=True)
class MappingResult:
    chosen_mapping: dict[int, int]
    candidate_to_xyz: dict[int, int]
    reference_to_candidate: dict[int, int]
    decision_relevant_signatures: tuple[tuple[Any, ...], ...]
    mapping_count_examined: int
    enumeration_truncated: bool
    timeout: bool
    confidence: str
    error: str = ""
    equal_best_mapping_count: int = 0

    @property
    def mapping_signature_count(self) -> int:
        return len(self.decision_relevant_signatures)


def map_candidate_reference_xyz(
    candidate: Chem.Mol,
    reference: Chem.Mol,
    xyz_atoms: list[tuple[str, tuple[float, float, float]]],
    *,
    max_matches: int = 256,
    max_mapping_combinations: int = 16384,
    timeout_seconds: int = 2,
    source_reference_to_xyz: Optional[dict[int, int]] = None,
    mcs_finder: Callable[..., Any] = rdFMCS.FindMCS,
) -> MappingResult:
    candidate_exact, candidate_to_xyz, xyz_error = _candidate_xyz_check(candidate, xyz_atoms, 1e-8)
    if not candidate_exact:
        return MappingResult({}, {}, {}, (), 0, False, False, "failed", xyz_error)
    try:
        candidate_heavy, candidate_to_source, candidate_h, candidate_attached_h = _prepare_heavy(
            candidate
        )
        reference_heavy, reference_to_source, reference_h, _ = _prepare_heavy(reference)
        mcs = mcs_finder(
            [candidate_heavy, reference_heavy],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringMatchesRingOnly=False,
            completeRingsOnly=False,
            timeout=timeout_seconds,
        )
        timed_out = bool(mcs.canceled)
        if mcs.queryMol is None or mcs.numAtoms == 0:
            return MappingResult(
                {}, candidate_to_xyz, {}, (), 0, False, timed_out, "failed", "mcs_failed"
            )
        if (
            mcs.numAtoms != candidate_heavy.GetNumAtoms()
            or mcs.numAtoms != reference_heavy.GetNumAtoms()
        ):
            return MappingResult(
                {},
                candidate_to_xyz,
                {},
                (),
                0,
                False,
                timed_out,
                "failed",
                "incomplete_heavy_atom_mapping",
            )
        candidate_matches = candidate_heavy.GetSubstructMatches(
            mcs.queryMol, uniquify=False, useChirality=False, maxMatches=max_matches
        )
        reference_matches = reference_heavy.GetSubstructMatches(
            mcs.queryMol, uniquify=False, useChirality=False, maxMatches=max_matches
        )
        possible_combinations = len(candidate_matches) * len(reference_matches)
        truncated = (
            len(candidate_matches) >= max_matches
            or len(reference_matches) >= max_matches
            or possible_combinations > max_mapping_combinations
        )
        best_score: Optional[tuple[int, int, int, int]] = None
        best_mappings: dict[tuple[tuple[int, int], ...], dict[int, int]] = {}
        examined = 0
        candidate_edges = _edge_map(candidate_heavy)
        reference_edges = _edge_map(reference_heavy)
        for candidate_match in candidate_matches:
            for reference_match in reference_matches:
                if examined >= max_mapping_combinations:
                    break
                examined += 1
                mapping = dict(zip(candidate_match, reference_match))
                score = _mapping_score(
                    candidate_heavy,
                    reference_heavy,
                    mapping,
                    candidate_edges,
                    reference_edges,
                )
                key = tuple(sorted(mapping.items()))
                if best_score is None or score < best_score:
                    best_score = score
                    best_mappings = {key: mapping}
                elif score == best_score:
                    best_mappings[key] = mapping
            if examined >= max_mapping_combinations:
                break
        if not best_mappings:
            return MappingResult(
                {}, candidate_to_xyz, {}, (), examined, truncated, timed_out, "failed", "no_mapping"
            )
        signatures = {
            _signature(
                candidate_heavy,
                reference_heavy,
                mapping,
                candidate_to_source,
                reference_to_source,
                candidate_to_xyz,
                candidate_h,
                reference_h,
                candidate_attached_h,
            )
            for mapping in best_mappings.values()
        }
        chosen = next(iter(best_mappings.values()))
        reference_to_candidate = {
            reference_to_source[reference_index]: candidate_to_source[candidate_index]
            for candidate_index, reference_index in chosen.items()
        }
        exact = source_reference_to_xyz is not None and all(
            source_reference_to_xyz.get(reference_index) == candidate_to_xyz[candidate_index]
            for reference_index, candidate_index in reference_to_candidate.items()
        )
        if exact:
            confidence = "exact"
        elif timed_out or truncated or len(signatures) != 1:
            confidence = "ambiguous"
        else:
            confidence = "unique_graph_mapping"
        return MappingResult(
            chosen,
            candidate_to_xyz,
            reference_to_candidate,
            tuple(sorted(signatures, key=repr)),
            examined,
            truncated,
            timed_out,
            confidence,
            equal_best_mapping_count=len(best_mappings),
        )
    except Exception as exc:
        return MappingResult(
            {},
            candidate_to_xyz,
            {},
            (),
            0,
            False,
            False,
            "failed",
            f"{type(exc).__name__}:{exc}",
        )
