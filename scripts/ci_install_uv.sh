#!/usr/bin/env sh
set -eu

: "${PYPI_INDEX_URL:?PYPI_INDEX_URL is required}"

USER_BIN="$(python3 -c 'import site; print(site.getuserbase())')/bin"
PIP_BREAK_SYSTEM_PACKAGES=""

if python3 -m pip install --help | grep -q -- '--break-system-packages'; then
  PIP_BREAK_SYSTEM_PACKAGES="--break-system-packages"
fi

python3 -m pip install --user ${PIP_BREAK_SYSTEM_PACKAGES} -i "$PYPI_INDEX_URL" uv
echo "$USER_BIN" >> "$GITHUB_PATH"
"$USER_BIN/uv" --version
