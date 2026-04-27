# MolGR

## Environment

- Python: `>=3.8`
- Package manager: `uv`
- Core runtime dependencies are declared in [pyproject.toml](/home/tmj/proj/MolGR/pyproject.toml)

## Development Setup

Create the default development environment:

```bash
uv sync --dev
```

Install the package in editable mode and build the C++ extension:

```bash
uv pip install -e . -v --no-build-isolation
```

If you only need to rebuild the extension after C++ changes:

```bash
make cpp-build
```

## Verification

Lint:

```bash
uv run ruff check .
```

Format check:

```bash
uv run ruff format --check .
```

Tests:

```bash
uv run pytest
```

Type check:

```bash
uv run mypy src
```

Build distributions:

```bash
uv build
```

## Benchmark Environment

Benchmarks use a dedicated environment instead of the default `.venv`:

```bash
bash scripts/benchmark_env.sh create
```

Run a benchmark inside that environment:

```bash
bash scripts/benchmark_env.sh run python benchmarks/smiles_xyz_benchmark/run.py --input tests/test_cases.csv --out benchmarks/_runs/run1
```

Optional shell switch for repeated benchmark commands:

```bash
eval "$(bash scripts/benchmark_env.sh env)"
```

## Optional IDE Setup

Generate VSCode / clangd C++ config from the active environment:

```bash
uv run python scripts/gen_vscode_config_with_ob.py --build-dir build
```
