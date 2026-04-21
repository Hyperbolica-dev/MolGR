# Benchmarks

This directory contains benchmark entrypoints and benchmark-specific docs.

## Available benchmarks

- `smiles_xyz_benchmark`: compares molecule reconstruction methods from XYZ-like inputs.
- `molfile_xyz_benchmark`: compares the same reconstruction methods using `.mol` / `.molfile` / `.sdf` inputs converted to XYZ cases.

Methods currently wired in `smiles_xyz_benchmark`:

- `rdkit_determine_bonds`
- `openbabel_read_xyz`
- `cell2mol_v2`
- `molgr_fallback`
- `xyzgraph_cheminf_full`

`molfile_xyz_benchmark` currently reuses the same method registry and result format as `smiles_xyz_benchmark`, but loads cases from molfile/SDF inputs through `scripts/molgr_cases_molfile.py`. The recommended fixture root is `tests/data/sdf/`, which can contain nested category directories such as cations, anions, or metal complexes.

## Requirements

- Project runtime: Python `>=3.8`
- Benchmark runtime: Python `>=3.10` (for `xyzgraph_cheminf_full`)
- Recommended: `uv`

## Quickstart (benchmark env)

From repo root:

```bash
bash scripts/benchmark_env.sh create
bash scripts/benchmark_env.sh run python benchmarks/smiles_xyz_benchmark/run.py \
  --input tests/test_cases.csv \
  --limit 10 \
  --out benchmarks/_runs/demo
```

Molfile / SDF benchmark example:

```bash
bash scripts/benchmark_env.sh run python benchmarks/molfile_xyz_benchmark/run.py \
  --input tests/data/sdf/MoNNMo.sdf \
  --limit 1 \
  --out benchmarks/_runs/molfile-demo
```

Optional shell switch for repeated benchmark commands:

```bash
eval "$(bash scripts/benchmark_env.sh env)"
```

Expected output directory:

- `benchmarks/_runs/demo`

## Run commands

Run the benchmark with the dedicated env script:

```bash
bash scripts/benchmark_env.sh run python benchmarks/smiles_xyz_benchmark/run.py --input tests/test_cases.csv --out benchmarks/_runs/run1
```

Run the molfile / SDF benchmark with the same environment:

```bash
bash scripts/benchmark_env.sh run python benchmarks/molfile_xyz_benchmark/run.py --input tests/data/sdf/MoNNMo.sdf --out benchmarks/_runs/molfile-run1
```

Useful flags:

- `--input`: input cases file or molfile/SDF fixture directory
- `--limit`: cap number of cases for quick checks
- `--out`: output run directory

## Outputs

Each run writes two main files under your `--out` directory:

- `results.csv`: one row per `(case, method)` attempt with status, prediction, equivalence result, and timing fields.
- `summary.csv`: aggregated metrics by method (counts and latency stats like average, p50, p95).

Timing columns in `results.csv` are flattened; `timing_ms_breakdown_json` preserves the full timing breakdown dict.

How to interpret quickly:

- `results.csv`: use for per-case debugging and failure analysis.
- `summary.csv`: use for method-level comparisons across success/failure rates and runtime.

## Optional `cell2mol_v2` setup and GPL note

`cell2mol_v2` is optional. Benchmarks still run without vendoring `cell2mol`.

We pin `numpy<2` in the benchmark dependency set for compatibility with optional `cell2mol`/`cosymlib` stacks.

If you enable `cell2mol_v2`, install/configure it separately in your environment and ensure your usage/redistribution complies with its license terms, including GPL obligations where applicable.

See benchmark-specific details in `benchmarks/smiles_xyz_benchmark/README.md`.
