# BDE-db XYZ Benchmark

This benchmark evaluates molecular-graph reconstruction from the optimized coordinates in
BDE-db. It is not a bond-dissociation-energy prediction benchmark.

The expected source is the official Springer Nature Figshare release:

- DOI: `10.6084/m9.figshare.12158646`
- file: `20200415_radical_database.sdf.gz`
- Figshare file id: `22357962`

Keep the source file outside Git, for example at
`.local/bde_db/raw/20200415_radical_database.sdf.gz`.

The pilot uses only `molgr_cpp`. It derives total charge from SDF atom formal charges. Spin
multiplicity is accepted only when the reference graph encodes zero radical electrons
(singlet) or exactly one radical electron (doublet); other records fail explicitly.

The primary metric is resonance-aware molecular-graph equivalence using
`check_equivalence(..., use_chirality=False, max_resonance=100)`. Exact isomeric SMILES is
diagnostic only. Exact formal-radical localization is reported as
`formal_radical_atom_index_match` only after an atom-identity guard verifies that predicted
atom order, elements, and coordinates correspond to the original XYZ/reference atoms.

Sampling is deterministic and proportional to the observed closed-shell, C-radical,
N-radical, and O-radical populations. The isolated `[H]` and `[O-]` records are always
included when they are inside the requested record range. There is no artificial charged
quota.

Run a deterministic, stratified 100-case pilot with:

```bash
bash scripts/benchmark_env.sh run python benchmarks/bde_db_benchmark/run.py \
  --input .local/bde_db/raw/20200415_radical_database.sdf.gz \
  --limit 100 \
  --seed 0 \
  --out .local/bde_db/runs/pilot-100
```

Use `--start` and `--end` for an inclusive one-based SDF record range. Outputs are
`results.csv`, `failures.csv`, `review_cases.csv`, `run_metadata.json`, and `summary.md`.
The review CSV prioritizes reconstruction failures, non-equivalent cases,
resonance-equivalent cases, formal-radical atom-index mismatches, charge/radical-electron
mismatches, and exact-SMILES mismatches.

The current loader and result writer are suitable for bounded pilots. Before running all
289,639 records, replace in-memory case/result collection with chunked loading and streaming
CSV output, and move large XYZ, bond, and source-metadata payloads into opt-in review
artifacts.
