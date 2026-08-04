#!/bin/sh
set -eu

review_root="${MOLGR_REVIEW_ROOT:-/workspace}"
review_db="${MOLGR_REVIEW_DB:-/var/lib/molgr-review/review.sqlite}"
fixtures_dir="${MOLGR_REVIEW_FIXTURES_DIR:-${review_root}/tests/data/reviewed/tmqmg}"
host="${MOLGR_REVIEW_HOST:-0.0.0.0}"
port="${MOLGR_REVIEW_PORT:-8765}"
xyz_dir="${MOLGR_XYZ_DIR:-}"

mkdir -p "$(dirname "$review_db")" "$fixtures_dir"

echo "Installing current MolGR checkout: ${review_root}"
uv pip install --system --no-cache -e "$review_root"

for benchmark_venv in "${MOLGR_PY38_VENV:-}" "${MOLGR_PY310_VENV:-}"; do
    benchmark_python="${benchmark_venv}/bin/python"
    if [ -n "$benchmark_venv" ] && [ -x "$benchmark_python" ]; then
        echo "Installing current MolGR checkout into ${benchmark_venv}"
        SETUPTOOLS_SCM_PRETEND_VERSION="${MOLGR_BUILD_VERSION:-0.0.0}" \
        SETUPTOOLS_USE_DISTUTILS=stdlib uv pip install \
            --python "$benchmark_python" \
            --no-build-isolation \
            --no-cache \
            --reinstall-package molgr \
            -e "$review_root"
    fi
done

if [ ! -f "$review_db" ]; then
    echo "Initializing review database: ${review_db}"
    python -c 'import sqlite3, sys; db, schema = sys.argv[1:]; conn = sqlite3.connect(db); conn.executescript(open(schema, encoding="utf-8").read()); conn.commit(); conn.close()' \
        "$review_db" "$review_root/tools/molgr_review/schema.sql"
fi

set -- python "$review_root/tools/molgr_review/server.py" \
    --db "$review_db" \
    --fixtures-dir "$fixtures_dir" \
    --host "$host" \
    --port "$port"
if [ -n "$xyz_dir" ]; then
    set -- "$@" --xyz-dir "$xyz_dir"
fi

echo "Starting MolGR review server on ${host}:${port}"
exec "$@"
