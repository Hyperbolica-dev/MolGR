from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Deque, Dict, List, Optional, Sequence, Tuple, cast


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from openbabel import openbabel as ob
from openbabel import pybel

from molgr.fallback.pipeline import resonance as resonance_module
from molgr.fallback.pipeline.reconstruct_with_metals import prepare_metal_state
from molgr.fallback.pipeline.reconstruct_without_metals import (
    _DEFAULT_RESONANCE_TRAVERSAL_POLICY,
    _run_linear_pipeline,
    _seed_state,
)
from molgr.fallback.pipeline.resonance import (
    build_processed_resonance_key,
    process_resonance,
)
from molgr.fallback.stages.preprocess import validate_omol
from molgr.fallback.state import OmolStateMachine, ReconstructionState
from scripts.molgr_cases_molfile import load_molfile_cases


@dataclass(frozen=True)
class RawResonanceRecord:
    index: int
    parent_index: Optional[int]
    depth: int
    move_path: Optional[Tuple[int, int, int]]
    search_key: resonance_module.ResonanceStateKey
    omol: pybel.Molecule


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Trace the no-metal radical resonance traversal for a molfile/SDF case and "
            "dump the created resonances, processed candidates, and frontier decisions."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("tests/data/sdf/MoNNMo.sdf"),
        help="Molfile/SDF path to inspect.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=2,
        help="Maximum resonance traversal depth.",
    )
    parser.add_argument(
        "--policy",
        choices=("default", "none"),
        default="default",
        help=(
            "Traversal policy to apply. 'default' matches the current fallback pipeline; "
            "'none' disables traversal-policy pruning."
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("/tmp/molgr_fallback_resonance_trace"),
        help="Directory to write CSV summaries and per-state serialization artifacts.",
    )
    return parser.parse_args()


def _short_digest(value: object) -> str:
    if isinstance(value, str):
        payload = value.encode("utf-8", errors="replace")
    else:
        payload = repr(value).encode("utf-8", errors="replace")
    return hashlib.sha1(payload).hexdigest()[:12]


def _compact_paths(paths: Sequence[Tuple[int, int, int]]) -> str:
    if not paths:
        return ""
    return ";".join(f"{a}-{b}-{c}" for a, b, c in paths)


def _atom_annotations(omol: pybel.Molecule) -> Tuple[str, str, str]:
    obmol = cast(ob.OBMol, omol.OBMol)
    charges: List[str] = []
    radicals: List[str] = []
    aromatic_atoms: List[str] = []
    for atom in ob.OBMolAtomIter(obmol):
        obatom = cast(ob.OBAtom, atom)
        atom_idx = obatom.GetIdx()
        charge = cast(int, obatom.GetFormalCharge())
        spin = cast(int, obatom.GetSpinMultiplicity())
        if charge != 0:
            charges.append(f"{atom_idx}:{charge}")
        if spin > 0:
            radicals.append(f"{atom_idx}:{spin}")
        if bool(obatom.IsAromatic()):
            aromatic_atoms.append(str(atom_idx))
    return "|".join(charges), "|".join(radicals), "|".join(aromatic_atoms)


def _preview_text(omol: pybel.Molecule) -> str:
    try:
        return cast(str, omol.write("can")).strip()
    except Exception:
        return ""


def _write_artifact(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _build_raw_resonance_records(
    state: ReconstructionState,
    *,
    max_depth: int,
    policy_name: str,
) -> Tuple[List[RawResonanceRecord], List[dict[str, object]]]:
    traversal_policy = _DEFAULT_RESONANCE_TRAVERSAL_POLICY if policy_name == "default" else None
    root_omol = state.omol
    root_key, bond_index_map = resonance_module._build_resonance_search_context(root_omol)
    raw_records: List[RawResonanceRecord] = [
        RawResonanceRecord(
            index=0,
            parent_index=None,
            depth=0,
            move_path=None,
            search_key=root_key,
            omol=root_omol,
        )
    ]
    frontier_rows: List[dict[str, object]] = []
    seen: Dict[resonance_module.ResonanceStateKey, int] = {root_key: 0}
    frontier: Deque[Tuple[int, pybel.Molecule, resonance_module.ResonanceStateKey, int]] = deque(
        [(0, root_omol, root_key, 0)]
    )

    while frontier:
        current_index, current_omol, current_key, depth = frontier.popleft()
        if depth >= max_depth:
            frontier_rows.append(
                {
                    "current_index": current_index,
                    "depth": depth,
                    "raw_move_count": 0,
                    "selected_move_count": 0,
                    "selected_new_child_indices": "",
                    "selected_seen_child_indices": "",
                    "raw_move_paths": "",
                    "selected_move_paths": "",
                    "pruned_move_paths": "",
                }
            )
            continue

        indexed_moves = resonance_module._enumerate_one_step_resonance_moves(
            current_omol,
            current_key,
            bond_index_map,
        )
        selected_moves = resonance_module._apply_resonance_traversal_policy(
            resonance_module.ResonanceTraversalContext(
                root_omol=root_omol,
                current_omol=current_omol,
                current_state_key=current_key,
                depth=depth,
                max_depth=max_depth,
            ),
            indexed_moves,
            traversal_policy,
        )

        selected_move_keys = {(move.idxs, move.next_state_key) for move in selected_moves}
        raw_paths = [move.idxs for move in indexed_moves]
        selected_paths = [move.idxs for move in selected_moves]
        pruned_paths = [
            move.idxs
            for move in indexed_moves
            if (move.idxs, move.next_state_key) not in selected_move_keys
        ]
        new_child_indices: List[str] = []
        seen_child_indices: List[str] = []

        for move in selected_moves:
            existing_index = seen.get(move.next_state_key)
            if existing_index is not None:
                seen_child_indices.append(str(existing_index))
                continue

            child_omol = resonance_module._materialize_one_step_resonance(current_omol, move.idxs)
            child_index = len(raw_records)
            seen[move.next_state_key] = child_index
            raw_records.append(
                RawResonanceRecord(
                    index=child_index,
                    parent_index=current_index,
                    depth=depth + 1,
                    move_path=move.idxs,
                    search_key=move.next_state_key,
                    omol=child_omol,
                )
            )
            frontier.append((child_index, child_omol, move.next_state_key, depth + 1))
            new_child_indices.append(str(child_index))

        frontier_rows.append(
            {
                "current_index": current_index,
                "depth": depth,
                "raw_move_count": len(indexed_moves),
                "selected_move_count": len(selected_moves),
                "selected_new_child_indices": "|".join(new_child_indices),
                "selected_seen_child_indices": "|".join(seen_child_indices),
                "raw_move_paths": _compact_paths(raw_paths),
                "selected_move_paths": _compact_paths(selected_paths),
                "pruned_move_paths": _compact_paths(pruned_paths),
            }
        )

    return raw_records, frontier_rows


def _process_resonance_records(
    linear_state: ReconstructionState,
    raw_records: Sequence[RawResonanceRecord],
    out_dir: Path,
) -> Tuple[List[dict[str, object]], Optional[int]]:
    summary_rows: List[dict[str, object]] = []
    seen_processed_keys: Dict[str, int] = {}
    best_row_index: Optional[int] = None
    best_score = float("inf")
    base_machine = OmolStateMachine.from_reconstruction_state(linear_state)

    for raw_record in raw_records:
        raw_charges, raw_radicals, raw_aromatic_atoms = _atom_annotations(raw_record.omol)
        raw_artifact_path = (
            out_dir / "raw" / f"{raw_record.index:02d}_depth{raw_record.depth}.molreport.txt"
        )
        _write_artifact(raw_artifact_path, cast(str, raw_record.omol.write("molreport")))

        candidate_machine = base_machine.branch(
            "branch_resonance_candidate",
            omol=raw_record.omol,
        )
        candidate_machine.run_omol_charge_stage("process_resonance", process_resonance)
        processed_key = candidate_machine.get_cached_omol_value(
            "resonance_state_key",
            build_processed_resonance_key,
        )
        processed_key_digest = _short_digest(processed_key)
        first_seen_index = seen_processed_keys.get(processed_key)
        is_duplicate = first_seen_index is not None
        if first_seen_index is None:
            seen_processed_keys[processed_key] = raw_record.index

        processed_omol = candidate_machine.omol
        processed_charges, processed_radicals, processed_aromatic_atoms = _atom_annotations(
            processed_omol
        )
        processed_artifact_path = (
            out_dir / "processed" / f"{raw_record.index:02d}_depth{raw_record.depth}.molreport.txt"
        )
        _write_artifact(processed_artifact_path, cast(str, processed_omol.write("molreport")))

        is_valid = False
        organic_core_score: Optional[float] = None
        full_score: Optional[float] = None
        if not is_duplicate:
            is_valid = validate_omol(
                processed_omol,
                linear_state.total_charge,
                linear_state.total_radical_electrons,
            )
            if is_valid:
                candidate_machine.annotate(
                    "validate_resonance_candidate",
                    resonance_index=raw_record.index,
                )
                candidate_state = candidate_machine.freeze_like(linear_state)
                organic_core_score = candidate_state.score("organic_core")
                candidate_state.post_reinsertion_base_components()
                full_score = candidate_state.full_score()
                if full_score < best_score:
                    best_score = full_score
                    best_row_index = raw_record.index

        summary_rows.append(
            {
                "raw_index": raw_record.index,
                "depth": raw_record.depth,
                "parent_index": "" if raw_record.parent_index is None else raw_record.parent_index,
                "move_path": ""
                if raw_record.move_path is None
                else f"{raw_record.move_path[0]}-{raw_record.move_path[1]}-{raw_record.move_path[2]}",
                "raw_search_key_digest": _short_digest(raw_record.search_key),
                "raw_preview": _preview_text(raw_record.omol),
                "raw_charges": raw_charges,
                "raw_radicals": raw_radicals,
                "raw_aromatic_atoms": raw_aromatic_atoms,
                "raw_artifact": str(raw_artifact_path),
                "processed_key_digest": processed_key_digest,
                "dedup_duplicate": is_duplicate,
                "dedup_first_raw_index": "" if first_seen_index is None else first_seen_index,
                "processed_preview": _preview_text(processed_omol),
                "processed_charges": processed_charges,
                "processed_radicals": processed_radicals,
                "processed_aromatic_atoms": processed_aromatic_atoms,
                "processed_artifact": str(processed_artifact_path),
                "validate_passed": is_valid,
                "organic_core_score": "" if organic_core_score is None else organic_core_score,
                "full_score": "" if full_score is None else full_score,
                "selected_best": False,
            }
        )

    if best_row_index is not None:
        for row in summary_rows:
            if row["raw_index"] == best_row_index:
                row["selected_best"] = True
                break

    return summary_rows, best_row_index


def _write_csv(path: Path, rows: Sequence[dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    args = _parse_args()
    cases = load_molfile_cases(args.input, limit=1)
    if not cases:
        raise ValueError(f"No cases loaded from {args.input}")
    case = cases[0]
    if case.get("provider_error"):
        raise ValueError(f"Case provider failed: {case['provider_error']}")

    xyz_block = case.get("xyz_block")
    total_charge = case.get("total_charge")
    total_radical_electrons = case.get("total_radical_electrons")
    if not isinstance(xyz_block, str) or not isinstance(total_charge, int):
        raise TypeError("Loaded case is missing xyz_block or total_charge")
    if not isinstance(total_radical_electrons, int):
        raise TypeError("Loaded case is missing total_radical_electrons")

    prepared = prepare_metal_state(xyz_block, total_charge, total_radical_electrons)
    linear_state = _run_linear_pipeline(
        _seed_state(
            prepared.no_metal_xyz_block,
            prepared.total_charge,
            prepared.total_radical_electrons,
        )
    )

    out_dir = args.out_dir / args.input.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_artifact(
        out_dir / "linear_state.molreport.txt", cast(str, linear_state.omol.write("molreport"))
    )

    direct_valid = validate_omol(
        linear_state.omol,
        linear_state.total_charge,
        linear_state.total_radical_electrons,
    )
    raw_records, frontier_rows = _build_raw_resonance_records(
        linear_state,
        max_depth=args.max_depth,
        policy_name=args.policy,
    )
    summary_rows, best_raw_index = _process_resonance_records(linear_state, raw_records, out_dir)

    _write_csv(
        out_dir / "frontier.csv",
        frontier_rows,
        (
            "current_index",
            "depth",
            "raw_move_count",
            "selected_move_count",
            "selected_new_child_indices",
            "selected_seen_child_indices",
            "raw_move_paths",
            "selected_move_paths",
            "pruned_move_paths",
        ),
    )
    _write_csv(
        out_dir / "summary.csv",
        summary_rows,
        (
            "raw_index",
            "depth",
            "parent_index",
            "move_path",
            "raw_search_key_digest",
            "raw_preview",
            "raw_charges",
            "raw_radicals",
            "raw_aromatic_atoms",
            "raw_artifact",
            "processed_key_digest",
            "dedup_duplicate",
            "dedup_first_raw_index",
            "processed_preview",
            "processed_charges",
            "processed_radicals",
            "processed_aromatic_atoms",
            "processed_artifact",
            "validate_passed",
            "organic_core_score",
            "full_score",
            "selected_best",
        ),
    )

    unique_processed = sum(1 for row in summary_rows if not row["dedup_duplicate"])
    valid_candidates = sum(1 for row in summary_rows if row["validate_passed"])
    print(f"input={args.input}")
    print(f"policy={args.policy} max_depth={args.max_depth}")
    print(f"direct_valid={direct_valid}")
    print(f"created_raw_resonances={len(raw_records)}")
    print(f"unique_processed_resonances={unique_processed}")
    print(f"valid_processed_candidates={valid_candidates}")
    print(f"best_raw_index={best_raw_index}")
    print(f"out_dir={out_dir}")
    print(f"summary_csv={out_dir / 'summary.csv'}")
    print(f"frontier_csv={out_dir / 'frontier.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
