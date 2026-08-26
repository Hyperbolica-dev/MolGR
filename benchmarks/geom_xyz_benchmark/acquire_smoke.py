"""Derive frozen 100-molecule smoke fixtures from official GEOM crude data."""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import math
import tarfile
from pathlib import Path
from typing import Any, Iterator

from rdkit import Chem, RDLogger
from rdkit.Chem import Lipinski

from benchmarks.geom_xyz_benchmark.adapter import stable_molecule_id


QM9_QUOTAS = {
    "stereo": 20,
    "aromatic": 20,
    "hetero_rich": 20,
    "flexible": 15,
    "small": 10,
    "other": 15,
}
DRUGS_QUOTAS = {
    "heavy_01_15": 15,
    "heavy_16_25": 25,
    "heavy_26_35": 25,
    "heavy_36_50": 20,
    "heavy_51_plus": 15,
}
OFFICIAL_ARTIFACTS = {
    "qm9": {
        "filename": "qm9_crude.msgpack.tar.gz",
        "dataverse_file_id": 4327190,
        "md5": "aad0081ed5d9b8c93c2bd0235987573b",
    },
    "drugs": {
        "filename": "drugs_crude.msgpack.tar.gz",
        "dataverse_file_id": 4360331,
        "md5": "7778e84c50b7cde755cca670d1f75091",
    },
}
RDLogger.DisableLog("rdApp.*")


def _entries(archive: Path) -> Iterator[tuple[str, dict[str, Any]]]:
    import msgpack

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


def _qm9_stratum(mol: Chem.Mol) -> str:
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


def _drugs_stratum(mol: Chem.Mol) -> str:
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms <= 15:
        return "heavy_01_15"
    if heavy_atoms <= 25:
        return "heavy_16_25"
    if heavy_atoms <= 35:
        return "heavy_26_35"
    if heavy_atoms <= 50:
        return "heavy_36_50"
    return "heavy_51_plus"


def _fixture_record(
    smiles: str, payload: dict[str, Any], *, dataset: str
) -> tuple[dict[str, Any] | None, str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "reference_parse_failure"
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    radicals = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    if charge != 0 or radicals != 0:
        return None, "electronic_state_excluded"
    if len(Chem.GetMolFrags(mol)) != 1:
        return None, "fragmented_reference"
    conformers = payload.get("conformers") or []
    with_xyz = [(i, conf) for i, conf in enumerate(conformers) if conf.get("xyz")]
    if not with_xyz:
        return None, "missing_xyz"
    viable: list[tuple[int, dict[str, Any], float]] = []
    for conf_idx, conf in with_xyz:
        try:
            relative_energy = float(conf["relativeenergy"])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(relative_energy):
            viable.append((conf_idx, conf, relative_energy))
    if not viable:
        return None, "missing_relativeenergy"
    conf_idx, conf, relative_energy = min(
        viable,
        key=lambda item: (item[2], item[0]),
    )
    table = Chem.GetPeriodicTable()
    xyz_rows = conf["xyz"]
    if Chem.AddHs(mol).GetNumAtoms() != len(xyz_rows):
        return None, "reference_xyz_atom_count_mismatch"
    dataset_label = dataset.upper()
    xyz = [str(len(xyz_rows)), f"GEOM {dataset_label} {smiles} conformer {conf_idx}"]
    xyz.extend(
        f"{table.GetElementSymbol(int(row[0]))} {float(row[1]):.10f} {float(row[2]):.10f} {float(row[3]):.10f}"
        for row in xyz_rows
    )
    molecule_id = stable_molecule_id(smiles)
    return {
        "case_id": f"geom-{dataset}-{molecule_id}-c{conf_idx}",
        "molecule_id": molecule_id,
        "conformer_id": str(conf_idx),
        "xyz": "\n".join(xyz) + "\n",
        "total_charge": 0,
        "spin_multiplicity": 1,
        "reference_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        "source_metadata": {
            "dataset": f"GEOM-{dataset_label}",
            "release": "Harvard Dataverse v4 (2022-02-11)",
            "source_artifact": f"{dataset}_crude.msgpack.tar.gz",
            "source_smiles": smiles,
            "selection": "minimum_relativeenergy_then_source_index",
            "relativeenergy_kcal_mol": relative_energy,
            "totalenergy_hartree": conf.get("totalenergy"),
            "source_conformer_count": len(conformers),
        },
    }, "eligible"


def acquire(
    archive: Path,
    output: Path,
    *,
    seed: int,
    dataset: str = "qm9",
    sample_size: int = 100,
    sampling: str = "representative_smoke",
) -> None:
    if sample_size <= 0:
        raise ValueError("sample_size must be positive")
    if sampling not in ("representative_smoke", "population"):
        raise ValueError(f"unsupported sampling mode: {sampling}")
    if sampling == "representative_smoke" and sample_size != 100:
        raise ValueError("representative_smoke uses the fixed 100-molecule quotas")
    quotas = DRUGS_QUOTAS if dataset == "drugs" else QM9_QUOTAS
    heap_keys = ("population",) if sampling == "population" else tuple(quotas)
    heap_limits = {"population": sample_size} if sampling == "population" else quotas
    heaps: dict[str, list[tuple[int, str, dict[str, Any]]]] = {key: [] for key in heap_keys}
    audit: dict[str, int] = {
        "source_molecules": 0,
        "eligible": 0,
        "reference_parse_failure": 0,
        "electronic_state_excluded": 0,
        "fragmented_reference": 0,
        "missing_xyz": 0,
        "missing_relativeenergy": 0,
        "reference_xyz_atom_count_mismatch": 0,
        "duplicate_canonical_reference": 0,
    }
    seen_reference_smiles: set[str] = set()
    for smiles, payload in _entries(archive):
        audit["source_molecules"] += 1
        record, reason = _fixture_record(smiles, payload, dataset=dataset)
        audit[reason] += 1
        if record is None:
            continue
        reference_smiles = record["reference_smiles"]
        if reference_smiles in seen_reference_smiles:
            audit["duplicate_canonical_reference"] += 1
            audit["eligible"] -= 1
            continue
        seen_reference_smiles.add(reference_smiles)
        mol = Chem.MolFromSmiles(record["reference_smiles"])
        assert mol is not None
        stratum = (
            "population"
            if sampling == "population"
            else (_drugs_stratum(mol) if dataset == "drugs" else _qm9_stratum(mol))
        )
        score = int.from_bytes(hashlib.sha256(f"{seed}:{smiles}".encode()).digest()[:8], "big")
        item = (-score, smiles, record)
        heap = heaps[stratum]
        if len(heap) < heap_limits[stratum]:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)
    records = [item[2] for heap in heaps.values() for item in heap]
    if len(records) != sample_size:
        raise RuntimeError(f"insufficient eligible records: selected {len(records)}")
    records.sort(key=lambda row: row["case_id"])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in records), encoding="utf-8"
    )
    audit.update(selected_molecules=len(records), seed=seed, sample_size=sample_size)
    output.with_suffix(".acquisition.json").write_text(
        json.dumps(
            {
                "counts": audit,
                "dataset": dataset,
                "source": {
                    **OFFICIAL_ARTIFACTS[dataset],
                    "release": "Harvard Dataverse v4 (2022-02-11)",
                    "doi": "10.7910/DVN/JNGTDF",
                },
                "stratum_quotas": quotas if sampling == "representative_smoke" else None,
                "selection": {
                    "conformer": "minimum_relativeenergy_then_source_index",
                    "molecule": (
                        "lowest_sha256(seed:source_smiles)_over_all_eligible_unique_molecules"
                        if sampling == "population"
                        else "lowest_sha256(seed:source_smiles)_within_size_stratum"
                    ),
                    "primary_unit": "molecule",
                    "sampling": sampling,
                },
                "xyz2mol_provenance": "unavailable_in_crude_artifact",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260825)
    parser.add_argument("--dataset", choices=("qm9", "drugs"), default="qm9")
    parser.add_argument("--sample-size", type=int, default=100)
    parser.add_argument(
        "--sampling", choices=("representative_smoke", "population"), default="representative_smoke"
    )
    args = parser.parse_args()
    acquire(
        args.archive,
        args.output,
        seed=args.seed,
        dataset=args.dataset,
        sample_size=args.sample_size,
        sampling=args.sampling,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
