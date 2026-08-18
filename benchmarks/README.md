# Benchmarks

[English](README.md) | [中文](README.zh-CN.md)

This directory contains benchmark entrypoints, benchmark-specific dependencies, and
benchmark documentation. Benchmark dependencies are intentionally isolated from the
root project dependency graph so they do not affect package releases.

## Available Benchmarks

- `smiles_xyz_benchmark`: builds XYZ cases from SMILES inputs and compares reconstruction methods.
- `molfile_xyz_benchmark`: loads `.mol`, `.molfile`, and `.sdf` inputs, converts them to XYZ cases,
  and runs the same method registry.
- `tmqmg_xyz_benchmark`: runs the shared method registry on tmQMg CSV/XYZ pairs and supports
  row/id subset selection.

All benchmark entrypoints use the same shared method registry:

- `rdkit_determine_bonds`
- `openbabel_read_xyz`
- `cell2mol_v2`
- `molgr_fallback`
- `molgr_cpp`
- `xyzgraph_cheminf_full`

`cell2mol_v2` and `xyzgraph_cheminf_full` are competitor baselines. `molgr_fallback`
uses the Python reference backend, and `molgr_cpp` uses the default C++ backend.

## Environment

- Project runtime: Python `>=3.8`
- Benchmark runtime: Python `>=3.10,<3.12`
- Recommended package manager: `uv`
- Benchmark dependency file: [`pyproject.toml`](pyproject.toml)
- Benchmark lock file: [`uv.lock`](uv.lock)
- Default benchmark virtualenv: `.venv-benchmark`
- Default benchmark Python executable: `python3.10`

The benchmark runtime is narrower than the package runtime because optional comparison
stacks such as `xyzgraph`, `cell2mol`, and `cosymlib` have older dependency constraints.
The benchmark dependency set also keeps `numpy<2` for compatibility with those stacks.

## Quickstart

Create or refresh the dedicated benchmark environment from the repository root:

```bash
bash scripts/benchmark_env.sh create
```

Run a small SMILES/XYZ benchmark:

```bash
bash scripts/benchmark_env.sh run python benchmarks/smiles_xyz_benchmark/run.py \
  --input tests/test_cases.csv \
  --limit 10 \
  --out benchmarks/_runs/smiles-demo
```

Run a small Molfile/SDF benchmark:

```bash
bash scripts/benchmark_env.sh run python benchmarks/molfile_xyz_benchmark/run.py \
  --input tests/data/sdf/MoNNMo.sdf \
  --limit 1 \
  --out benchmarks/_runs/molfile-demo
```

Run a small tmQMg benchmark on a selected subset:

```bash
bash scripts/benchmark_env.sh run python benchmarks/tmqmg_xyz_benchmark/run.py \
  --csv /path/to/tmqmg.csv \
  --xyz-dir /path/to/xyz \
  --start-row 1 \
  --end-row 100 \
  --ids ABC123,DEF456 \
  --limit 10 \
  --out benchmarks/_runs/tmqmg-demo
```

For backend parity checks, restrict tmQMg to MolGR's two backends:

```bash
bash scripts/benchmark_env.sh run python benchmarks/tmqmg_xyz_benchmark/run.py \
  --csv /path/to/tmqmg.csv \
  --xyz-dir /path/to/xyz \
  --limit 1000 \
  --out benchmarks/_runs/tmqmg-parity \
  --case-timeout-seconds 1.0 \
  --cpp-accelerations all \
  --methods molgr_cpp,molgr_fallback
```

`tmqmg_xyz_benchmark` defaults to `--process-workers 1` and runs rows serially; the
C++ single-molecule path still keeps internal target-bucket parallelism. For
`molgr_cpp` with a larger value, one subprocess passes the worker count to the
native batch pool instead of stacking external processes. Measure batch worker
counts on the target machine rather than assuming the largest value is fastest.

For repeated commands, switch the current shell to the benchmark project:

```bash
eval "$(bash scripts/benchmark_env.sh env)"
```

Then use normal `uv run` commands:

```bash
uv run python benchmarks/smiles_xyz_benchmark/run.py \
  --input tests/test_cases.csv \
  --limit 10 \
  --out benchmarks/_runs/smiles-demo
```

Switch back to the default project environment:

```bash
unset UV_PROJECT UV_PROJECT_ENVIRONMENT UV_PYTHON
```

## Inputs

SMILES/XYZ benchmark:

- `--input` points to a text or CSV-style file containing SMILES.
- Header rows named `smiles`, `canonicalsmiles`, or `canonicalsmi` are accepted.
- Each SMILES is embedded with RDKit, optimized with UFF, converted to XYZ, and compared
  against the original RDKit molecule.

Molfile/SDF benchmark:

- `--input` accepts one `.mol`, `.molfile`, or `.sdf` file, or a directory scanned recursively.
- Recommended fixture root: [`../tests/data/sdf/`](../tests/data/sdf/).
- Nested fixture categories under [`../tests/data/sdf/`](../tests/data/sdf/) are supported.
- Each input structure is loaded with RDKit, converted to an XYZ case, and compared against
  the source molecule.

Shared flags:

- `--limit`: cap the number of cases for quick checks.
- `--out`: output run directory.
- `--methods`: for tmQMg, restrict the shared method registry by method id.
- `--process-workers`: tmQMg serial/native-batch worker budget; `1` is the serial
  single-molecule path, and larger C++ values use one native batch pool.
- `--cpp-accelerations`: for tmQMg, choose the C++ acceleration preset.

## Outputs

Each run writes two main files under the `--out` directory:

- `results.csv`: one row per `(case, method)` attempt with status, prediction,
  equivalence result, and timing fields.
- `summary.csv`: aggregated metrics by method, including counts and latency statistics.

The important `results.csv` columns are:

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
- `method_ms`
- `equivalence_ms`
- `timing_ms_breakdown_json`

Timing columns are flattened when present. `timing_ms_breakdown_json` preserves the full
method-specific timing dictionary for profiling.

## Interpreting Results

- Use `summary.csv` for method-level comparison.
- Use `results.csv` for per-case debugging and failure analysis.
- A row is counted as successful when equivalence is `True`.
- A row is counted as failed when the method errors or equivalence is `False`.
- A row is counted as skipped when the input case could not be prepared.

## Optional `cell2mol_v2` Note

`cell2mol_v2` is optional as a benchmark baseline. The benchmark dependency file points
to the upstream `cell2mol` `v2` Git revision. If you enable or redistribute this baseline,
ensure your usage complies with its license terms, including GPL obligations where applicable.

## Benchmark-Specific Docs

- [`smiles_xyz_benchmark/README.md`](smiles_xyz_benchmark/README.md)
- [`smiles_xyz_benchmark/README.zh-CN.md`](smiles_xyz_benchmark/README.zh-CN.md)
- [`molfile_xyz_benchmark/README.md`](molfile_xyz_benchmark/README.md)
- [`molfile_xyz_benchmark/README.zh-CN.md`](molfile_xyz_benchmark/README.zh-CN.md)
- [`tmqmg_xyz_benchmark/README.md`](tmqmg_xyz_benchmark/README.md)
- [`tmqmg_xyz_benchmark/README.zh-CN.md`](tmqmg_xyz_benchmark/README.zh-CN.md)
