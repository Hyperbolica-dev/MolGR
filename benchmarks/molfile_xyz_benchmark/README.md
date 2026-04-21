# Molfile/SDF XYZ Benchmark

## What it does

This benchmark compares molecule reconstruction approaches from `.mol`, `.molfile`, and `.sdf` inputs.

The recommended fixture root is `tests/data/sdf/`. You can place files directly there or organize them into nested categories such as `tests/data/sdf/cations/`, `tests/data/sdf/anions/`, and `tests/data/sdf/metal_complexes/`.

Each input structure is first converted into an XYZ case with:

- `xyz_block`
- `total_charge`
- `total_radical_electrons`
- `ground_truth_rdmol`

using `scripts/molgr_cases_molfile.py`, then benchmarked with the shared method registry from `benchmarks/smiles_xyz_benchmark`.

## Methods

This benchmark currently reuses the same methods as `smiles_xyz_benchmark`:

- `rdkit_determine_bonds`
- `openbabel_read_xyz`
- `cell2mol_v2`
- `molgr_fallback`
- `molgr_fallback`
- `molgr_cpp`
- `xyzgraph_cheminf_full`

## Benchmark Environment (Python >=3.10)

Reproduce the benchmark environment exactly as documented in `benchmarks/README.md`.

From repo root:

```bash
bash scripts/benchmark_env.sh create
```

Run directly inside the dedicated benchmark environment:

```bash
bash scripts/benchmark_env.sh run python benchmarks/molfile_xyz_benchmark/run.py \
  --input tests/data/sdf/MoNNMo.sdf \
  --limit 1 \
  --out benchmarks/_runs/molfile-demo
```

Optional shell switch for repeated commands:

```bash
eval "$(bash scripts/benchmark_env.sh env)"
```

Then run:

```bash
uv run python benchmarks/molfile_xyz_benchmark/run.py \
  --input tests/data/sdf/MoNNMo.sdf \
  --limit 1 \
  --out benchmarks/_runs/molfile-demo
```

Switch back to the default project environment:

```bash
unset UV_PROJECT_ENVIRONMENT UV_PYTHON
```

## Inputs

- `--input` accepts either:
  - a single `.mol` / `.molfile` / `.sdf` file, or
  - a directory, recursively scanned for those suffixes.
- Recommended fixture root: `tests/data/sdf/`.
- Nested directories are supported, so you can benchmark a single class with a path like `tests/data/sdf/metal_complexes/` or the whole fixture tree with `tests/data/sdf/`.
- `--limit` optionally caps the number of discovered files.

## Outputs

The output directory contains:

- `results.csv`: one row per `(case, method)` run.
- `summary.csv`: aggregated metrics by method.

The CSV schema is intentionally shared with `smiles_xyz_benchmark` so comparisons stay straightforward.
