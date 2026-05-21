#!/usr/bin/env sh
set -eu

: "${PYPI_INDEX_URL:?PYPI_INDEX_URL is required}"

PYTHON_BIN="${PYTHON_BIN:-}"
for candidate in python3 /usr/bin/python3 /usr/local/bin/python3; do
  if [ -z "$PYTHON_BIN" ] && command -v "$candidate" >/dev/null 2>&1 && "$candidate" -m pip --version >/dev/null 2>&1; then
    PYTHON_BIN="$(command -v "$candidate")"
  fi
done
if [ -z "$PYTHON_BIN" ]; then
  echo "No python3 with pip found" >&2
  exit 1
fi

USER_BIN="$("$PYTHON_BIN" -c 'import site; print(site.getuserbase())')/bin"
UV_BIN="$USER_BIN/uv"
PIP_BREAK_SYSTEM_PACKAGES=""

if "$PYTHON_BIN" -m pip install --help | grep -q -- '--break-system-packages'; then
  PIP_BREAK_SYSTEM_PACKAGES="--break-system-packages"
fi

"$PYTHON_BIN" -m pip install --user ${PIP_BREAK_SYSTEM_PACKAGES} -i "$PYPI_INDEX_URL" uv
if [ ! -x "$UV_BIN" ]; then
  if command -v uv >/dev/null 2>&1; then
    UV_BIN="$(command -v uv)"
    USER_BIN="$(dirname "$UV_BIN")"
  else
    echo "uv executable not found at $UV_BIN or on PATH" >&2
    exit 1
  fi
fi

if [ -n "${GITHUB_PATH:-}" ]; then
  echo "$USER_BIN" >> "$GITHUB_PATH"
fi
if [ -n "${GITHUB_ENV:-}" ]; then
  echo "UV_BIN=$UV_BIN" >> "$GITHUB_ENV"
fi

"$UV_BIN" --version
