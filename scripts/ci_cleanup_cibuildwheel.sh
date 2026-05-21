#!/usr/bin/env sh
set -eu

cleanup_tmp_path() {
  path="$1"
  case "$path" in
    "$(pwd)/.tmp/"* | .tmp/*)
      if [ -d "$path" ]; then
        rm -rf "$path"
      fi
      ;;
    *)
      echo "Refusing to remove unexpected path: $path" >&2
      exit 1
      ;;
  esac
}

CIBUILDWHEEL_ENV_DIR="${CIBUILDWHEEL_ENV_DIR:-$(pwd)/.tmp/cibuildwheel-venv}"
cleanup_tmp_path "$CIBUILDWHEEL_ENV_DIR"
