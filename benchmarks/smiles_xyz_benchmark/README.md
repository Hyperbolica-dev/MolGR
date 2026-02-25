<!--
 * @Author: TMJ
 * @Date: 2026-02-25 14:46:29
 * @LastEditors: TMJ
 * @LastEditTime: 2026-02-25 19:22:16
 * @Description: 请填写简介
-->
# SMILES/XYZ Benchmark

## What it does

This benchmark compares molecule reconstruction approaches from XYZ-like inputs and reports per-case outputs plus aggregate metrics.

## Methods

- `rdkit_determine_bonds`
- `openbabel_read_xyz`
- `cell2mol_v2`
- `molgr_fallback`

## Run

```bash
uv run python benchmarks/smiles_xyz_benchmark/run.py --input tests/test_cases.csv --limit 10 --out benchmarks/_runs/demo
```

## Outputs

The run directory contains:

- `results.csv`: one row per test case and method.
- `summary.csv`: aggregated metrics by method.

Timing columns in `results.csv` are flattened; `timing_ms_breakdown_json` preserves the full timing breakdown dict.

### `results.csv` columns

- `case_idx`
- `method_id`
- `input_smiles`
- `ground_truth_smiles`
- `status`
- `error`
- `predicted_smiles`
- `equivalent`
- `equivalence_method`
- `timing_ms_total`
- `timing_ms_breakdown_json`

### `summary.csv` columns

- `method_id`
- `count`
- `success_count`
- `fail_count`
- `skip_count`
- `avg_ms_total`
- `p50_ms_total`
- `p95_ms_total`

## Licensing note (`cell2mol_v2`)

`cell2mol_v2` is optional. If enabled, ensure your use and redistribution comply with its license terms (including GPL obligations where applicable).

## Reproducibility

- Pin the environment and dependency versions (for example via `uv.lock`).
- Run with the same input file and flags (`--input`, `--limit`, `--out`).
- Keep method set fixed when comparing runs.
- Record the commit hash and run timestamp with each benchmark output.
