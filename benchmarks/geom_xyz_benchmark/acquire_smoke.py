"""Derive the frozen 100-molecule smoke fixture from official GEOM QM9 crude data."""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import tarfile
from pathlib import Path
from typing import Any, Iterator

import msgpack
from rdkit import Chem, RDLogger
from rdkit.Chem import Lipinski

from benchmarks.geom_xyz_benchmark.adapter import stable_molecule_id


QUOTAS = {"stereo": 20, "aromatic": 20, "hetero_rich": 20, "flexible": 15, "small": 10, "other": 15}
RDLogger.DisableLog("rdApp.*")


def _entries(archive: Path) -> Iterator[tuple[str, dict[str, Any]]]:
    with tarfile.open(archive, "r|gz") as tar:  # type: ignore[call-overload]
        member = next(item for item in tar if item.isfile() and item.name.endswith(".msgpack"))
        stream = tar.extractfile(member)
        if stream is None:
            raise RuntimeError(f"cannot read {member.name}")
        unpacker = msgpack.Unpacker(stream, raw=False, max_buffer_size=1024**3)
        while True:
            try:
                count = unpacker.read_map_header()
            except msgpack.OutOfData:
                return
            for _ in range(count):
                yield unpacker.unpack(), unpacker.unpack()


def _stratum(mol: Chem.Mol) -> str:
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    if "@" in smiles:
        return "stereo"
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return "aromatic"
    hetero = sum(atom.GetAtomicNum() not in (1, 6) for atom in mol.GetAtoms())
    if hetero >= 3:
        return "hetero_rich"
    if Lipinski.NumRotatableBonds(mol) >= 2:
        return "flexible"
    if mol.GetNumHeavyAtoms() <= 5:
        return "small"
    return "other"


def _fixture_record(smiles: str, payload: dict[str, Any]) -> tuple[dict[str, Any] | None, str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "reference_parse_failure"
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    radicals = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    if charge != 0 or radicals != 0:
        return None, "electronic_state_excluded"
    conformers = payload.get("conformers") or []
    viable = [(i, conf) for i, conf in enumerate(conformers) if conf.get("xyz")]
    if not viable:
        return None, "missing_xyz"
    conf_idx, conf = min(
        viable,
        key=lambda item: (float(item[1].get("relativeenergy", float("inf"))), item[0]),
    )
    table = Chem.GetPeriodicTable()
    xyz_rows = conf["xyz"]
    xyz = [str(len(xyz_rows)), f"GEOM QM9 {smiles} conformer {conf_idx}"]
    xyz.extend(
        f"{table.GetElementSymbol(int(row[0]))} {float(row[1]):.10f} {float(row[2]):.10f} {float(row[3]):.10f}"
        for row in xyz_rows
    )
    molecule_id = stable_molecule_id(smiles)
    return {
        "case_id": f"geom-qm9-{molecule_id}-c{conf_idx}",
        "molecule_id": molecule_id,
        "conformer_id": str(conf_idx),
        "xyz": "\n".join(xyz) + "\n",
        "total_charge": 0,
        "spin_multiplicity": 1,
        "reference_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        "source_metadata": {
            "dataset": "GEOM-QM9",
            "release": "Harvard Dataverse v4 (2022-02-11)",
            "source_smiles": smiles,
            "selection": "minimum_relativeenergy_then_source_index",
            "relativeenergy_kcal_mol": conf.get("relativeenergy"),
            "totalenergy_hartree": conf.get("totalenergy"),
            "source_conformer_count": len(conformers),
        },
    }, "eligible"


def acquire(archive: Path, output: Path, *, seed: int) -> None:
    heaps: dict[str, list[tuple[int, str, dict[str, Any]]]] = {key: [] for key in QUOTAS}
    audit: dict[str, int] = {"source_molecules": 0}
    for smiles, payload in _entries(archive):
        audit["source_molecules"] += 1
        record, reason = _fixture_record(smiles, payload)
        audit[reason] = audit.get(reason, 0) + 1
        if record is None:
            continue
        mol = Chem.MolFromSmiles(record["reference_smiles"])
        assert mol is not None
        stratum = _stratum(mol)
        score = int.from_bytes(hashlib.sha256(f"{seed}:{smiles}".encode()).digest()[:8], "big")
        item = (-score, smiles, record)
        heap = heaps[stratum]
        if len(heap) < QUOTAS[stratum]:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)
    records = [item[2] for heap in heaps.values() for item in heap]
    if len(records) != sum(QUOTAS.values()):
        raise RuntimeError(f"insufficient eligible records: selected {len(records)}")
    records.sort(key=lambda row: row["case_id"])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in records), encoding="utf-8"
    )
    audit.update(selected_molecules=len(records), seed=seed, quota_total=sum(QUOTAS.values()))
    output.with_suffix(".acquisition.json").write_text(
        json.dumps({"counts": audit, "stratum_quotas": QUOTAS}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260825)
    args = parser.parse_args()
    acquire(args.archive, args.output, seed=args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
