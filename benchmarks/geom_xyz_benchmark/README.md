# GEOM molecular-graph benchmark (READY protocol)

This adapter evaluates graph reconstruction from source GEOM coordinates.  It deliberately uses
the molecule's original/canonical input SMILES as the reference graph, not GEOM's
`xyz2mol_smiles`, because the latter was re-inferred from the same coordinates and is therefore
not independent ground truth.

## Unit and electronic-state policy

The primary unit is a **molecule**.  One conformer is selected deterministically per molecule
(minimum reported conformer energy, then source order as a tie-breaker).  Multiple conformers
must not be counted as independent chemical evidence.  A later robustness study may report
per-molecule all/mean/worst-conformer recovery separately.

GEOM does not publish a trustworthy spin-multiplicity field in the audited MessagePack release.
The v1 adapter therefore admits only references whose RDKit graph is neutral and has zero radical
electrons, and records charge `0` and multiplicity `1`.  Charged or open-shell records are loader
failures, not silently coerced.  This is an inclusion policy, not a claim that all GEOM records
have known singlet states.

## Evaluation

Primary equivalence is `molgr.utils.equivalence.evaluate_equivalence(..., use_chirality=False)`.
Its invariant gates and bounded resonance handling remain unchanged.  Exact non-isomeric
canonical SMILES is diagnostic only.  Chirality is a secondary diagnostic because the GEOM paper
states that source SMILES can omit stereo while each 3D conformer necessarily realizes it; the
paper's own graph recovery validation removed stereochemical indicators.

Reference preprocessing is limited to parsing the published canonical source SMILES.  No
prediction-side sanitization, kekulization, charge repair, fragment selection, or valence repair is
performed by this adapter.  Aromatic SMILES and Kekule forms are representation choices; primary
equivalence handles them without mutating only the prediction.

Run the frozen 100-molecule smoke fixture:

```bash
uv run python benchmarks/geom_xyz_benchmark/run.py \
  --input benchmarks/geom_xyz_benchmark/data/qm9_smoke100.jsonl \
  --limit 100 --seed 20260825 --out benchmarks/geom_xyz_benchmark/_runs/smoke100
```

Outputs are `summary.json`, `results.csv`, and `failures.csv`.

Rebuild the fixture only from the checksummed official archive (the source archive remains
read-only and is gitignored):

```bash
uv run --with msgpack python -m benchmarks.geom_xyz_benchmark.acquire_smoke \
  --archive benchmarks/_data/geom/qm9_crude.msgpack.tar.gz \
  --output benchmarks/geom_xyz_benchmark/data/qm9_smoke100.jsonl --seed 20260825
```
