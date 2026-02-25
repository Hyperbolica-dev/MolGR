# Benchmarks

This directory contains benchmark entrypoints and benchmark-specific docs.

## Available benchmarks

- `smiles_xyz_benchmark`: compares molecule reconstruction methods from XYZ-like inputs.

Methods currently wired in `smiles_xyz_benchmark`:

- `rdkit_determine_bonds`
- `openbabel_read_xyz`
- `cell2mol_v2`
- `molgr_fallback`

## Requirements

- Python `>=3.8`
- Recommended: `uv`
- Fallback: `pip` (`python -m venv` + `python -m pip`)

## Quickstart (uv)

From repo root:

```bash
uv venv
source .venv/bin/activate
uv pip install -e .
uv run python benchmarks/smiles_xyz_benchmark/run.py \
  --input tests/test_cases.csv \
  --limit 10 \
  --out benchmarks/_runs/demo
```

Expected output directory:

- `benchmarks/_runs/demo`

## Environment install (pip fallback)

From repo root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
python benchmarks/smiles_xyz_benchmark/run.py \
  --input tests/test_cases.csv \
  --limit 10 \
  --out benchmarks/_runs/demo
```

## Run commands

Run the benchmark with `uv`:

```bash
uv run python benchmarks/smiles_xyz_benchmark/run.py --input tests/test_cases.csv --out benchmarks/_runs/run1
```

Useful flags:

- `--input`: input cases file
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

If you enable `cell2mol_v2`, install/configure it separately in your environment and ensure your usage/redistribution complies with its license terms, including GPL obligations where applicable.

See benchmark-specific details in `benchmarks/smiles_xyz_benchmark/README.md`.
