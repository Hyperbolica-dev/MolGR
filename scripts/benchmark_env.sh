#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(CDPATH='' cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_DIR="${BENCH_VENV_DIR:-"${ROOT_DIR}/.venv-benchmark"}"
PYTHON="${BENCH_PYTHON:-python3.10}"

usage() {
  cat <<'EOF'
Usage: bash scripts/benchmark_env.sh <command> [args...]

Commands:
  env
    Print export lines for uv project environment switching.

  create
    Create/sync a dedicated benchmark environment and build the C++ extension.

  run <cmd> [args...]
    Run a command inside the benchmark environment.

Environment variables:
  BENCH_VENV_DIR   Path to the benchmark venv directory (default: .venv-benchmark)
  BENCH_PYTHON     Python interpreter for the benchmark env (default: python3.10)
EOF
}

cmd="${1:-}"
shift || true

case "$cmd" in
  env)
    printf 'export UV_PROJECT_ENVIRONMENT=%q\n' "$VENV_DIR"
    printf 'export UV_PYTHON=%q\n' "$PYTHON"
    ;;
  create)
    UV_PROJECT_ENVIRONMENT="$VENV_DIR" UV_PYTHON="$PYTHON" uv venv --allow-existing "$VENV_DIR"
    UV_PROJECT_ENVIRONMENT="$VENV_DIR" UV_PYTHON="$PYTHON" SETUPTOOLS_USE_DISTUTILS=1 uv sync --group benchmark  --no-build-isolation
    UV_PROJECT_ENVIRONMENT="$VENV_DIR" UV_PYTHON="$PYTHON" uv pip install --python "$VENV_DIR/bin/python" -e .
    ;;
  run)
    if [ "$#" -lt 1 ]; then
      usage
      exit 2
    fi
    UV_PROJECT_ENVIRONMENT="$VENV_DIR" UV_PYTHON="$PYTHON" uv run --frozen --group benchmark --group dev "$@"
    ;;
  -h|--help|help|"")
    usage
    ;;
  *)
    echo "Unknown command: $cmd" >&2
    usage >&2
    exit 2
    ;;
esac
