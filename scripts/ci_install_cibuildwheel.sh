#!/usr/bin/env sh
set -eu

: "${PYPI_INDEX_URL:?PYPI_INDEX_URL is required}"

UV_BIN="${UV_BIN:-}"
USER_UV_BIN="$(python3 -c 'import site; print(site.getuserbase())')/bin/uv"
if [ -z "$UV_BIN" ] && [ -x "$USER_UV_BIN" ]; then
  UV_BIN="$USER_UV_BIN"
fi
if [ -z "$UV_BIN" ] && command -v uv >/dev/null 2>&1; then
  UV_BIN="$(command -v uv)"
fi
if [ -z "$UV_BIN" ]; then
  UV_BIN="$USER_UV_BIN"
fi

"$UV_BIN" python install 3.11

CIBUILDWHEEL_ENV_DIR="${CIBUILDWHEEL_ENV_DIR:-$(pwd)/.tmp/cibuildwheel-venv}"
"$UV_BIN" venv --clear --python 3.11 "$CIBUILDWHEEL_ENV_DIR"
"$UV_BIN" pip install --python "$CIBUILDWHEEL_ENV_DIR/bin/python" --refresh-package cibuildwheel --reinstall-package cibuildwheel --index "$PYPI_INDEX_URL" --default-index "$PYPI_INDEX_URL" cibuildwheel==3.4.1

CIBUILDWHEEL_BIN_DIR="$CIBUILDWHEEL_ENV_DIR/bin"
CIBUILDWHEEL="$CIBUILDWHEEL_BIN_DIR/cibuildwheel"
export CIBUILDWHEEL
export CIBUILDWHEEL_ENV_DIR
if [ -z "${CIBW_CACHE_PATH:-}" ]; then
  CIBW_CACHE_PATH="$HOME/.cache/cibuildwheel"
fi
mkdir -p "$CIBW_CACHE_PATH"
export CIBW_CACHE_PATH
if [ ! -x "$CIBUILDWHEEL" ]; then
  echo "cibuildwheel executable not found at $CIBUILDWHEEL" >&2
  exit 1
fi

if [ -n "${GITHUB_PATH:-}" ]; then
  echo "$CIBUILDWHEEL_BIN_DIR" >> "$GITHUB_PATH"
fi
if [ -n "${GITHUB_ENV:-}" ]; then
  echo "CIBUILDWHEEL=$CIBUILDWHEEL" >> "$GITHUB_ENV"
  echo "CIBUILDWHEEL_ENV_DIR=$CIBUILDWHEEL_ENV_DIR" >> "$GITHUB_ENV"
  echo "CIBW_CACHE_PATH=$CIBW_CACHE_PATH" >> "$GITHUB_ENV"
fi

"$CIBUILDWHEEL" --help >/dev/null
