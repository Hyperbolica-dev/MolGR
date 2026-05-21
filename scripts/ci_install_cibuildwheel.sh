#!/usr/bin/env sh
set -eu

: "${PYPI_INDEX_URL:?PYPI_INDEX_URL is required}"

UV_BIN="${UV_BIN:-}"
USER_UV_BIN="$(python3 -m site --user-base)/bin/uv"
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
"$UV_BIN" tool install --force --python 3.11 --index "$PYPI_INDEX_URL" --default-index "$PYPI_INDEX_URL" cibuildwheel==3.4.1

CIBUILDWHEEL_BIN_DIR="$("$UV_BIN" tool dir --bin)"
CIBUILDWHEEL="$CIBUILDWHEEL_BIN_DIR/cibuildwheel"
export CIBUILDWHEEL
if [ ! -x "$CIBUILDWHEEL" ]; then
  echo "cibuildwheel executable not found at $CIBUILDWHEEL" >&2
  exit 1
fi

if [ -n "${GITHUB_PATH:-}" ]; then
  echo "$CIBUILDWHEEL_BIN_DIR" >> "$GITHUB_PATH"
fi
if [ -n "${GITHUB_ENV:-}" ]; then
  echo "CIBUILDWHEEL=$CIBUILDWHEEL" >> "$GITHUB_ENV"
fi

"$CIBUILDWHEEL" --help >/dev/null
