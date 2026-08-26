"""Freeze an already evaluated GEOM-Drugs formal run into a paper export."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path
from typing import Any


RUN_COMMAND = (
    ".venv/bin/python benchmarks/geom_xyz_benchmark/formal_run.py "
    "--archive benchmarks/_data/geom/drugs_crude.msgpack.tar.gz "
    "--out benchmarks/geom_xyz_benchmark/_runs/drugs_formal_full "
    "--expected-eligible 291709 --case-timeout-seconds 10"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def freeze(run_dir: Path, export_dir: Path, *, molgr_git_sha: str) -> None:
    summary: dict[str, Any] = json.loads((run_dir / "summary.json").read_text(encoding="utf-8"))
    expected = {"results.csv.gz", "failures.csv", "review_cases.csv", "summary.json"}
    missing = sorted(name for name in expected if not (run_dir / name).is_file())
    if missing:
        raise FileNotFoundError(f"formal run is incomplete: {missing}")
    export_dir.mkdir(parents=True, exist_ok=True)
    for name in sorted(expected):
        shutil.copyfile(run_dir / name, export_dir / name)

    provenance = f"""# GEOM-Drugs formal benchmark provenance

## Dataset

- Release: GEOM Harvard Dataverse v4, 2022-02-11
- DOI: `10.7910/DVN/JNGTDF`
- Dataverse file ID: `4360331`
- Source archive: `drugs_crude.msgpack.tar.gz`
- Source archive MD5: `7778e84c50b7cde755cca670d1f75091` (verified)
- Source records: 292,035
- Eligible unique molecules: 291,709
- Exclusions: 15 non-neutral/open-shell references; 311 duplicate canonical references; zero
  parse, fragment, XYZ, finite-relative-energy, or reference/XYZ atom-count exclusions
- `xyz2mol_smiles` was not used as ground truth and the featurized release was not acquired.

## Frozen protocol

- Primary unit: one molecule, one deterministic conformer
- Conformer: minimum finite `relativeenergy`, original source index tie-break
- Reference: source/canonical SMILES molecular graph
- Electronic-state scope: neutral, zero-radical, single fragment; charge 0, multiplicity 1
- Reconstruction: `molgr_cpp`, git SHA `{molgr_git_sha}`
- Evaluator: evaluator v1; Candidate = predicted, Reference = reference,
  `use_chirality=False`
- Decisions: `equivalent`, `not_equivalent`, `inconclusive`
- Exact SMILES and chirality are diagnostics only.
- No permissive sulfoxide `S=O` / `[S+][O-]` rule was added.

Exact command:

```bash
{RUN_COMMAND}
```

## Environment and execution

- Python: {summary["environment"]["python"]}
- RDKit: {summary["environment"]["rdkit"]}
- Platform: {summary["environment"]["platform"]}
- Wall time: {summary["wall_time_seconds"]:.6f} seconds

`results.csv.gz`, `failures.csv`, `review_cases.csv`, and `summary.json` are byte-for-byte copies
of the frozen evaluated run. Creating this export did not invoke reconstruction or the evaluator.
"""
    (export_dir / "provenance.md").write_text(provenance, encoding="utf-8")
    checksums = {
        name: {"sha256": _sha256(export_dir / name), "bytes": (export_dir / name).stat().st_size}
        for name in (
            "summary.json",
            "results.csv.gz",
            "failures.csv",
            "review_cases.csv",
            "provenance.md",
        )
    }
    manifest = {
        "schema_version": 1,
        "dataset": "GEOM-Drugs",
        "release": "Harvard Dataverse v4 (2022-02-11)",
        "doi": "10.7910/DVN/JNGTDF",
        "dataverse_file_id": 4360331,
        "source_archive": {
            "filename": "drugs_crude.msgpack.tar.gz",
            "md5": "7778e84c50b7cde755cca670d1f75091",
            "verified": True,
        },
        "source_molecules": 292035,
        "eligible_unique_molecules": 291709,
        "molgr_git_sha": molgr_git_sha,
        "evaluator": "molgr evaluator v1, Candidate=predicted, Reference=reference, use_chirality=False",
        "run_command": RUN_COMMAND,
        "files": checksums,
    }
    (export_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument("--molgr-git-sha", required=True)
    args = parser.parse_args()
    freeze(args.run_dir, args.export_dir, molgr_git_sha=args.molgr_git_sha)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
