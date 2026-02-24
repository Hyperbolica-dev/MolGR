# MolGR

MolGR provides molecule graph reasoning utilities plus debug/backtest tooling.

## Project Structure

- `src/`: MolGR library code.
- `scripts/`: backtest/debug scripts and HTML report tooling.
- `tests/`: test data and test cases.
- `pyproject.toml`: project metadata and Python tooling config.

## Debug Trace / Backtest

Run from repo root with `uv run python3 ...`.

### SMILES backtest trace

```bash
uv run python3 scripts/smiles_backtest_debug_trace.py \
  --input tests/test_cases.csv \
  --out-root .molgr_backtest_smiles_snapshots_full \
  --limit 50 \
  --no-chirality \
  --max-resonance 100
```

- Script: `scripts/smiles_backtest_debug_trace.py`
- Flags:
  - `--input` (default: `tests/test_cases.csv`)
  - `--out-root` (default: `.molgr_backtest_smiles_snapshots_full/`)
  - `--limit` (default: all)
  - `--no-chirality` (optional)
  - `--max-resonance` (default: `100`)

### Molfile backtest trace

```bash
uv run python3 scripts/molfile_backtest_debug_trace.py \
  --input tests/path/to_cases.mol \
  --out-root .molgr_backtest_molfile_snapshots \
  --limit 50 \
  --no-chirality \
  --max-resonance 100
```

- Script: `scripts/molfile_backtest_debug_trace.py`
- Flags:
  - `--input` (required)
  - `--out-root` (default: `.molgr_backtest_molfile_snapshots/`)
  - `--limit` (default: all)
  - `--no-chirality` (optional)
  - `--max-resonance` (default: `100`)

### Render HTML report from trace

```bash
uv run python3 scripts/molgr_debug_html.py \
  --trace-dir .molgr_backtest_smiles_snapshots_full/run-case001 \
  --out .molgr_backtest_smiles_snapshots_full/run-case001/report.html \
  --max-events-rendered 2000 \
  --max-atoms-svg 200
```

- Script: `scripts/molgr_debug_html.py`
- Flags:
  - `--trace-dir` (required; directory containing `trace.jsonl`)
  - `--out` (optional; default: `<trace-dir>/report.html`)
  - `--max-events-rendered` (default: `2000`)
  - `--max-atoms-svg` (default: `200`)

### Output layout

- Backtest scripts write per-case directories under `--out-root` as `run-caseNNN/`.
- Each case directory contains `trace.jsonl` and `report.html`.
- Default output roots are `.molgr_backtest_smiles_snapshots_full/` and `.molgr_backtest_molfile_snapshots/`.

### HTML viewer notes

- Candidate pool lists `omol_score` and `#node_id` for each scoring candidate.
- Clicking a candidate highlights its reasoning chain and jumps to the corresponding event details.
