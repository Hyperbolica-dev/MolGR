from __future__ import annotations

import argparse
import ast
import csv
import gzip
import hashlib
import heapq
import json
from collections import Counter
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


SAMPLE_QUOTAS = {
    "closed_shell": 4,
    "c_radical": 8,
    "n_radical": 3,
    "o_radical": 3,
    "charged": 1,
    "other_radical": 1,
}


def _score(seed: int, record_index: int, case_id: str) -> int:
    value = f"{seed}:{record_index}:{case_id}".encode()
    return int.from_bytes(hashlib.sha256(value).digest()[:8], "big")


def _bonds(mol: Chem.Mol) -> str:
    return json.dumps(
        [
            [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), str(bond.GetBondType())]
            for bond in mol.GetBonds()
        ],
        separators=(",", ":"),
    )


def _sample_stratum(charge: int, radical_atoms: list[Chem.Atom]) -> str:
    if charge != 0:
        return "charged"
    electrons = sum(atom.GetNumRadicalElectrons() for atom in radical_atoms)
    if electrons == 0:
        return "closed_shell"
    if electrons == 1 and len(radical_atoms) == 1:
        element = radical_atoms[0].GetSymbol().lower()
        if element in {"c", "n", "o"}:
            return f"{element}_radical"
    return "other_radical"


def _property_list_length(mol: Chem.Mol, name: str) -> int | None:
    if not mol.HasProp(name):
        return None
    value = ast.literal_eval(mol.GetProp(name))
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"{name} is not a list")
    return len(value)


def _without_hs_for_graph(mol: Chem.Mol) -> Chem.Mol:
    if all(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()):
        return Chem.Mol(mol)
    return Chem.RemoveAllHs(mol)


def inspect_dataset(input_path: Path, out_dir: Path, seed: int = 0) -> dict[str, Any]:
    counters: dict[str, Counter[Any]] = {
        "properties": Counter(),
        "formal_charge": Counter(),
        "radical_electrons": Counter(),
        "radical_centers": Counter(),
        "radical_elements": Counter(),
    }
    samples: dict[str, list[tuple[int, int, dict[str, Any]]]] = {}
    failures: list[dict[str, Any]] = []
    total_records = 0
    parsed_records = 0
    graph_disagreements = 0
    stereo_disagreements = 0
    atom_charges_length_mismatches = 0
    atom_spins_length_mismatches = 0
    explicit_multiplicity_fields: Counter[str] = Counter()
    graph_disagreement_examples: list[dict[str, Any]] = []
    atom_spins_on_closed_shell: list[dict[str, Any]] = []
    charged_records: list[dict[str, Any]] = []
    other_radical_records: list[dict[str, Any]] = []

    with gzip.open(input_path, "rb") as stream:
        supplier = Chem.ForwardSDMolSupplier(
            stream,
            sanitize=True,
            removeHs=False,
            strictParsing=True,
        )
        for record_index, mol in enumerate(supplier, start=1):
            total_records += 1
            if mol is None:
                failures.append(
                    {"record_index": record_index, "error": "RDKit failed to parse SDF record"}
                )
                continue
            parsed_records += 1
            property_names = list(mol.GetPropNames(includePrivate=False, includeComputed=False))
            counters["properties"].update(property_names)
            explicit_multiplicity_fields.update(
                name for name in property_names if "multiplic" in name.lower()
            )
            charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
            radical_atoms = [atom for atom in mol.GetAtoms() if atom.GetNumRadicalElectrons()]
            radical_electrons = sum(atom.GetNumRadicalElectrons() for atom in radical_atoms)
            counters["formal_charge"][charge] += 1
            counters["radical_electrons"][radical_electrons] += 1
            counters["radical_centers"][len(radical_atoms)] += 1
            counters["radical_elements"].update(atom.GetSymbol() for atom in radical_atoms)
            try:
                atom_charges_length = _property_list_length(mol, "AtomCharges")
                atom_spins_length = _property_list_length(mol, "AtomSpins")
            except Exception as exc:
                failures.append(
                    {
                        "record_index": record_index,
                        "error": f"property parse failed: {type(exc).__name__}: {exc}",
                    }
                )
                continue
            atom_count = mol.GetNumAtoms()
            if atom_charges_length is not None and atom_charges_length != atom_count:
                atom_charges_length_mismatches += 1
            if atom_spins_length is not None and atom_spins_length != atom_count:
                atom_spins_length_mismatches += 1

            source_smiles = mol.GetProp("SMILES") if mol.HasProp("SMILES") else ""
            source_mol = Chem.MolFromSmiles(source_smiles) if source_smiles else None
            sdf_no_h = _without_hs_for_graph(mol)
            sdf_isomeric = Chem.MolToSmiles(sdf_no_h, canonical=True, isomericSmiles=True)
            sdf_graph = Chem.MolToSmiles(sdf_no_h, canonical=True, isomericSmiles=False)
            source_isomeric = None
            source_graph = None
            if source_mol is None:
                graph_disagreements += 1
            else:
                source_isomeric = Chem.MolToSmiles(source_mol, canonical=True, isomericSmiles=True)
                source_graph = Chem.MolToSmiles(source_mol, canonical=True, isomericSmiles=False)
                graph_disagreements += source_graph != sdf_graph
                stereo_disagreements += source_isomeric != sdf_isomeric

            case_id = mol.GetProp("_Name") if mol.HasProp("_Name") else str(record_index)
            stratum = _sample_stratum(charge, radical_atoms)
            row = {
                "identifier": case_id,
                "record_index": record_index,
                "stratum": stratum,
                "smiles": source_smiles,
                "sdf_smiles": sdf_isomeric,
                "formula": rdMolDescriptors.CalcMolFormula(mol),
                "total_charge": charge,
                "radical_electrons": radical_electrons,
                "radical_sites": json.dumps(
                    [atom.GetIdx() for atom in radical_atoms], separators=(",", ":")
                ),
                "radical_elements": json.dumps(
                    [atom.GetSymbol() for atom in radical_atoms], separators=(",", ":")
                ),
                "spin_multiplicity": 1 + radical_electrons
                if radical_electrons in {0, 1}
                else "UNKNOWN",
                "atom_count": atom_count,
                "sdf_reference_bonds": _bonds(mol),
                "atom_charges_present": atom_charges_length is not None,
                "atom_spins_present": atom_spins_length is not None,
                "graph_agrees": source_graph == sdf_graph,
                "stereo_agrees": source_isomeric == sdf_isomeric,
            }
            if source_graph != sdf_graph and len(graph_disagreement_examples) < 20:
                graph_disagreement_examples.append(row)
            if radical_electrons == 0 and atom_spins_length is not None:
                atom_spins_on_closed_shell.append(row)
            if charge != 0:
                charged_records.append(row)
            if stratum == "other_radical":
                other_radical_records.append(row)
            quota = SAMPLE_QUOTAS.get(stratum, 0)
            if quota:
                heap = samples.setdefault(stratum, [])
                score = _score(seed, record_index, case_id)
                item = (-score, record_index, row)
                if len(heap) < quota:
                    heapq.heappush(heap, item)
                elif score < -heap[0][0]:
                    heapq.heapreplace(heap, item)

    sample_rows = [
        item[2]
        for stratum in SAMPLE_QUOTAS
        for item in sorted(samples.get(stratum, []), key=lambda value: (-value[0], value[1]))
    ]
    out_dir.mkdir(parents=True, exist_ok=True)
    if sample_rows:
        with (out_dir / "inspected_cases.csv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(sample_rows[0]))
            writer.writeheader()
            writer.writerows(sample_rows)
    summary = {
        "input": str(input_path),
        "total_records": total_records,
        "parsed_records": parsed_records,
        "failures": failures,
        "property_counts": dict(counters["properties"]),
        "formal_charge_counts": dict(counters["formal_charge"]),
        "radical_electron_counts": dict(counters["radical_electrons"]),
        "radical_center_counts": dict(counters["radical_centers"]),
        "radical_element_counts": dict(counters["radical_elements"]),
        "explicit_multiplicity_fields": dict(explicit_multiplicity_fields),
        "graph_disagreements": graph_disagreements,
        "stereo_disagreements": stereo_disagreements,
        "atom_charges_length_mismatches": atom_charges_length_mismatches,
        "atom_spins_length_mismatches": atom_spins_length_mismatches,
        "sample_count": len(sample_rows),
        "graph_disagreement_examples": graph_disagreement_examples,
        "atom_spins_on_closed_shell": atom_spins_on_closed_shell,
        "charged_records": charged_records,
        "other_radical_records": other_radical_records,
    }
    (out_dir / "inspection_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Inspect the official BDE-db SDF.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    print(json.dumps(inspect_dataset(args.input, args.out, args.seed), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
