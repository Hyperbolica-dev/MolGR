# tmQMg XYZ Benchmark

This benchmark compares molecule reconstruction methods on tmQMg CSV/XYZ pairs.
It reuses the shared method registry from the other benchmark entrypoints.

## Inputs

`--csv` points to a tmQMg metadata CSV with at least `id`, `smiles`, and `charge` columns.
`--xyz-dir` points to a directory containing `<id>.xyz` files.
The official source, pinned revision, checksums, and download commands are in the
[molecule review development guide](../../docs/development/MOLECULE_REVIEW_TOOL.zh-CN.md#获取-tmqmg-数据).

Subset selection flags:

- `--start-row` and `--end-row` filter by 1-based CSV row index.
- `--ids` filters by tmQMg id. Repeat the flag or pass comma-separated values.
- `--limit` caps the number of selected rows after filtering.

Execution flags:

- `--methods` filters the shared method registry. Use `molgr_cpp,molgr_fallback`
  for C++/Python backend parity checks.
- `--cpp-accelerations default|all` selects the C++ backend acceleration preset.
  Both presets keep target-bucket parallelism enabled; `all` additionally enables
  optional C++ accelerations such as the vendor UFF atom-typing cache.
- `--enable-uff-atom-typing-cache` enables the vendor UFF atom-typing cache when
  using the default preset.
- `--process-workers N` splits each method across `N` worker subprocesses. This can
  stack with C++ internal target-bucket threading, but high values may compete for CPU.
- `--case-timeout-seconds` sets the per-method per-case wall-time limit; use `0` to
  disable it.

Backend parity check:

```bash
bash scripts/benchmark_env.sh run python benchmarks/tmqmg_xyz_benchmark/run.py \
  --csv /path/to/tmqmg.csv \
  --xyz-dir /path/to/xyz \
  --limit 1000 \
  --out benchmarks/_runs/tmqmg-parity \
  --progress-every 50 \
  --case-timeout-seconds 1.0 \
  --cpp-accelerations all \
  --methods molgr_cpp,molgr_fallback
```

## Outputs

The benchmark writes `results.csv` and `summary.csv` to `--out`.
