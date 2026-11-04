#!/usr/bin/env python3
"""Freeze corrected BorylXAT-DB SDF fixtures and electronic-state evidence.

Requires a local copy of the ASE SQLite file from Zenodo record 21330198.
No ASE or Parquet dependency is needed; this reads the selected SQLite rows.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sqlite3
import struct
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdqueries


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests/data/borylxat"
SOURCE = {
    "name": "BorylXAT-DB",
    "record_url": "https://zenodo.org/records/21330198",
    "doi": "10.5281/zenodo.21330198",
    "record_version": "v2",
    "file": "BorylXAT-DB.db",
    "download_url": "https://zenodo.org/records/21330198/files/BorylXAT-DB.db?download=1",
    "size_bytes": 291348480,
    "published_md5": "d486fbcc36f856bbebf303b5cda4950f",
    "license": "CC-BY-4.0",
    "license_url": "https://creativecommons.org/licenses/by/4.0/",
}


REFERENCE_STATUS = "corrected_source_graph_geometry_checked"
UPSTREAM = "https://github.com/jackie-illiilli/BorylXAT-DB/blob/main/"
STRUCTURE_DOC = UPSTREAM + "Database_Structure.md"
WORKFLOW_DOC = UPSTREAM + "1_Calc_Reactant.ipynb"


def electronic_state_basis(category: str) -> str:
    if category == "B":
        return (
            "inferred doublet: author open-shell category, odd electron count, spin-density sum ~1"
        )
    if category == "complex_r":
        return "author workflow: B_N_r doublet (multiplicity 2); spin-density sum ~1"
    return "author workflow: B_N_p singlet (multiplicity 1); author closed-shell category"


def freeze(connection, selections: list[dict], output: Path, *, checksum_verified=False) -> dict:
    """Preserve source geometry and correct source graphs independently of MolGR."""
    records = []
    geometries = []
    for selection in selections:
        row = connection.execute(
            "SELECT id, natoms, charge, key_value_pairs, numbers, positions, data "
            "FROM systems WHERE id = ?",
            (selection["ase_row_id"],),
        ).fetchone()
        if row is None:
            raise ValueError(f"Missing source row: {selection}")
        row_id, natoms, charge, metadata, numbers_blob, positions_blob, data_blob = row
        metadata = json.loads(metadata)
        case_id = selection["case_id"]
        if metadata["key"] != case_id or metadata["category"] != selection["category"]:
            raise ValueError(f"Source identity changed: {selection}")
        numbers = struct.unpack(f"<{natoms}i", numbers_blob)
        flat_positions = struct.unpack(f"<{3 * natoms}d", positions_blob)
        if not all(math.isfinite(x) for x in flat_positions):
            raise ValueError(f"Nonfinite positions: {case_id}")
        json_offset = int.from_bytes(data_blob[:8], "little")
        data = json.loads(data_blob[json_offset:])
        spins = data.get("spin_densities")
        radical = metadata["category"] != "complex_p"
        if radical and (
            not isinstance(spins, list)
            or len(spins) != natoms
            or not all(math.isfinite(x) for x in spins)
            or not math.isclose(sum(spins), 1.0, abs_tol=1e-4)
        ):
            raise ValueError(f"Missing or inconsistent open-shell evidence: {case_id}")
        if int(charge) != charge or charge != metadata["charge"]:
            raise ValueError(f"Inconsistent charge: {case_id}")
        multiplicity = 2 if radical else 1
        if (sum(numbers) - int(charge)) % 2 != (multiplicity - 1) % 2:
            raise ValueError(f"Inconsistent electron parity: {case_id}")
        geometry = Chem.RWMol()
        conf = Chem.Conformer(natoms)
        conf.Set3D(True)
        for i, number in enumerate(numbers):
            atom = Chem.Atom(number)
            atom.SetNoImplicit(True)
            geometry.AddAtom(atom)
            conf.SetAtomPosition(i, flat_positions[3 * i : 3 * i + 3])
        geometry.AddConformer(conf)
        geometries.append(geometry.GetMol())
        records.append(
            {
                **selection,
                "ase_row_id": row_id,
                "source_geometry_sha256": hashlib.sha256(numbers_blob + positions_blob).hexdigest(),
                "natoms": natoms,
                "source_formula": data["formula"],
                "charge": int(charge),
                "spin_multiplicity": multiplicity,
                "electronic_state_basis": electronic_state_basis(metadata["category"]),
                "electronic_state_sources": [STRUCTURE_DOC]
                + ([WORKFLOW_DOC] if metadata["category"] != "B" else []),
                "boron_atom_indices": [i for i, z in enumerate(numbers) if z == 5],
                "source_spin_densities": spins,
                "source_smiles": metadata["smiles"],
                "reference_status": "unreviewed_source_smiles_not_an_answer_graph",
            }
        )
    manifest = {
        "schema_version": 2,
        "source_dataset": {**SOURCE, "full_file_checksum_verified": checksum_verified},
        "coordinate_units": "angstrom",
        "atom_index_base": 0,
        "fixtures": records,
    }
    write_references(output, manifest, geometries)
    return manifest


def corrected_reference(geometry: Chem.Mol, case: dict) -> tuple[Chem.Mol, list[int], int]:
    params = Chem.SmilesParserParams()
    params.removeHs = False
    params.sanitize = False
    source = Chem.MolFromSmiles(case["source_smiles"], params)
    if source is None or source.GetNumAtoms() != case["natoms"]:
        raise ValueError(f"Source graph atom count mismatch: {case['case_id']}")
    removed_boron_h = 0
    for atom in source.GetAtoms():
        # The geometry already contains every H as an atom. Bracket B-H counts in
        # the source are additional phantom H, not missing coordinate atoms.
        atom.SetNoImplicit(True)
        if atom.GetAtomicNum() == 5:
            removed_boron_h += atom.GetNumExplicitHs()
            atom.SetNumExplicitHs(0)
            atom.SetNumRadicalElectrons(case["spin_multiplicity"] - 1)
    Chem.SanitizeMol(source)
    if any(a.GetNumExplicitHs() or a.GetNumImplicitHs() for a in source.GetAtoms()):
        raise ValueError(f"Unresolved virtual hydrogens: {case['case_id']}")
    if (
        Chem.GetFormalCharge(source) != case["charge"]
        or sum(a.GetNumRadicalElectrons() for a in source.GetAtoms())
        != case["spin_multiplicity"] - 1
    ):
        raise ValueError(f"Reference electronic-state mismatch: {case['case_id']}")

    geometry = Chem.Mol(geometry)
    # Check adjacency from coordinates independently of source bond orders.
    rdDetermineBonds.DetermineConnectivity(geometry, useVdw=True, covFactor=1.2)
    if case["category"] == "B":
        # These source atom-map numbers are the original zero-based coordinate indices.
        mapping = [a.GetAtomMapNum() for a in source.GetAtoms()]
        if sorted(mapping) != list(range(case["natoms"])):
            raise ValueError(f"Incomplete source atom mapping: {case['case_id']}")
    else:
        # Complex SMILES lack atom maps. Match element-labelled adjacency only;
        # equivalent hydrogens may be exchanged, but coordinates are never moved.
        query = Chem.RWMol()
        for atom in source.GetAtoms():
            query.AddAtom(rdqueries.AtomNumEqualsQueryAtom(atom.GetAtomicNum()))
        for bond in source.GetBonds():
            query.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), Chem.BondType.SINGLE)
        mapping = list(geometry.GetSubstructMatch(query))
        if len(mapping) != case["natoms"]:
            raise ValueError(f"Source/geometry connectivity mismatch: {case['case_id']}")
    order = sorted(range(case["natoms"]), key=lambda idx: mapping[idx])
    answer = Chem.RenumberAtoms(source, order)
    if (
        answer.GetNumBonds() != geometry.GetNumBonds()
        or any(
            geometry.GetBondBetweenAtoms(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) is None
            for b in answer.GetBonds()
        )
        or any(
            answer.GetAtomWithIdx(i).GetAtomicNum() != geometry.GetAtomWithIdx(i).GetAtomicNum()
            for i in range(case["natoms"])
        )
    ):
        raise ValueError(f"Mapped source/geometry mismatch: {case['case_id']}")
    for atom in answer.GetAtoms():
        atom.SetAtomMapNum(0)
    answer.AddConformer(Chem.Conformer(geometry.GetConformer()))
    # Validate supplied stereochemistry against the coordinates without
    # replacing the reference with a reconstruction or newly perceived graph.
    perceived = Chem.Mol(answer)
    Chem.AssignStereochemistryFrom3D(perceived)
    if Chem.MolToSmiles(answer) != Chem.MolToSmiles(perceived):
        raise ValueError(f"Source stereochemistry disagrees with coordinates: {case['case_id']}")
    answer.SetProp("_Name", case["case_id"])
    return answer, mapping, removed_boron_h


def write_references(root: Path, manifest: dict, geometries: list[Chem.Mol]) -> None:
    manifest_path = root / "manifest.json"
    results = []
    for case, geometry in zip(manifest["fixtures"], geometries):
        mol, mapping, removed_h = corrected_reference(geometry, case)
        params = Chem.MolWriterParams()
        params.forceV3000 = True
        params.precision = 17
        block = Chem.MolToV3KMolBlock(mol, params)
        properties = {
            "SOURCE_DOI": manifest["source_dataset"]["doi"],
            "SOURCE_KEY": case["case_id"],
            "SOURCE_SMILES": case["source_smiles"],
            "TOTAL_CHARGE": str(case["charge"]),
            "SPIN_MULTIPLICITY": str(case["spin_multiplicity"]),
            "ELECTRONIC_STATE_BASIS": case["electronic_state_basis"],
            "ELECTRONIC_STATE_SOURCES": json.dumps(case["electronic_state_sources"]),
            "REFERENCE_STATUS": REFERENCE_STATUS,
        }
        if case["source_spin_densities"] is not None:
            properties["SOURCE_SPIN_DENSITIES"] = json.dumps(case["source_spin_densities"])
        sdf = (
            block
            + "".join(f">  <{key}>\n{value}\n\n" for key, value in properties.items())
            + "$$$$\n"
        )
        sdf_bytes = sdf.encode("utf-8")
        reference_file = f"reference_graph/{case['case_id']}.sdf"
        case.update(
            reference_status=REFERENCE_STATUS,
            reference_file=reference_file,
            reference_smiles=Chem.MolToSmiles(Chem.RemoveHs(mol)),
            reference_source_atom_to_geometry=mapping,
            reference_corrections={
                "implicit_hydrogens_disabled": True,
                "extra_boron_bracket_hydrogens_removed": removed_h,
                "boron_radical_electrons": case["spin_multiplicity"] - 1,
                "source_formal_charges_preserved": True,
                "source_bond_orders_preserved_before_aromatic_normalization": True,
                "geometry_connectivity_check": "RDKit covalent radii, factor 1.2",
                "representation": "one B-centred Lewis form; source covalent/charge-separated convention",
            },
        )
        results.append((root / reference_file, sdf_bytes))
    # Complete every validation before replacing any of the reference files.
    (root / "reference_graph").mkdir(parents=True, exist_ok=True)
    for path, contents in results:
        path.write_bytes(contents)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=FIXTURES)
    args = parser.parse_args()
    digest = hashlib.md5()
    with args.db.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    if (
        args.db.stat().st_size != SOURCE["size_bytes"]
        or digest.hexdigest() != SOURCE["published_md5"]
    ):
        raise ValueError("Database does not match the pinned Zenodo file")
    selections = json.loads((FIXTURES / "selection.json").read_text(encoding="utf-8"))
    connection = sqlite3.connect(args.db.resolve().as_uri() + "?mode=ro", uri=True)
    try:
        manifest = freeze(connection, selections, args.output, checksum_verified=True)
    finally:
        connection.close()
    print(f"Wrote {len(manifest['fixtures'])} fixtures to {args.output}")


if __name__ == "__main__":
    main()
