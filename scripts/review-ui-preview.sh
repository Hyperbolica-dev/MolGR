#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Start an isolated Docker preview of the MolGR Review UI.

Usage:
  scripts/review-ui-preview.sh [options]

Options:
  --source-db PATH      Source review.sqlite (read-only; copied with SQLite backup)
  --xyz-dir PATH        Source XYZ directory (mounted read-only)
  --fixtures-dir PATH   Source reviewed fixtures directory (copied)
  --port PORT           Host port (default: 8766)
  --work-dir PATH       Preview state directory
  --container NAME      Docker container name
  --image IMAGE         Docker image name/tag
  --fresh               Remove an existing preview and copy sources again
  --clean               Remove the preview container and preview state, then exit
  -h, --help            Show this help

Environment equivalents:
  MOLGR_REVIEW_PREVIEW_SOURCE_DB
  MOLGR_REVIEW_PREVIEW_XYZ_DIR
  MOLGR_REVIEW_PREVIEW_SOURCE_FIXTURES
  MOLGR_REVIEW_PREVIEW_PORT
  MOLGR_REVIEW_PREVIEW_DIR
  MOLGR_REVIEW_PREVIEW_CONTAINER
  MOLGR_REVIEW_PREVIEW_IMAGE

Without --fresh, existing preview data is never overwritten.
EOF
}

die() {
  printf 'review-ui-preview: %s\n' "$*" >&2
  exit 1
}

command -v git >/dev/null 2>&1 || die "git is required"

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(git -C "$script_dir/.." rev-parse --show-toplevel)"

source_db="${MOLGR_REVIEW_PREVIEW_SOURCE_DB:-}"
xyz_dir="${MOLGR_REVIEW_PREVIEW_XYZ_DIR:-}"
source_fixtures="${MOLGR_REVIEW_PREVIEW_SOURCE_FIXTURES:-}"
port="${MOLGR_REVIEW_PREVIEW_PORT:-8766}"
preview_dir="${MOLGR_REVIEW_PREVIEW_DIR:-$repo_root/.local/molgr_review_preview}"
container_name="${MOLGR_REVIEW_PREVIEW_CONTAINER:-molgr-review-ui-preview}"
image_name="${MOLGR_REVIEW_PREVIEW_IMAGE:-}"
fresh=0
clean=0

while (($#)); do
  case "$1" in
    --source-db)
      (($# >= 2)) || die "--source-db requires a path"
      source_db="$2"
      shift 2
      ;;
    --xyz-dir)
      (($# >= 2)) || die "--xyz-dir requires a path"
      xyz_dir="$2"
      shift 2
      ;;
    --fixtures-dir)
      (($# >= 2)) || die "--fixtures-dir requires a path"
      source_fixtures="$2"
      shift 2
      ;;
    --port)
      (($# >= 2)) || die "--port requires a value"
      port="$2"
      shift 2
      ;;
    --work-dir)
      (($# >= 2)) || die "--work-dir requires a path"
      preview_dir="$2"
      shift 2
      ;;
    --container)
      (($# >= 2)) || die "--container requires a name"
      container_name="$2"
      shift 2
      ;;
    --image)
      (($# >= 2)) || die "--image requires a name"
      image_name="$2"
      shift 2
      ;;
    --fresh)
      fresh=1
      shift
      ;;
    --clean)
      clean=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown argument: $1 (use --help)"
      ;;
  esac
done

if command -v python3 >/dev/null 2>&1; then
  host_python="$(command -v python3)"
elif command -v python >/dev/null 2>&1; then
  host_python="$(command -v python)"
else
  die "Python is required to create a consistent SQLite preview copy"
fi

canonical_path() {
  "$host_python" -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).expanduser().resolve())' "$1"
}

preview_dir="$(canonical_path "$preview_dir")"
repo_root="$(canonical_path "$repo_root")"

validate_preview_target() {
  [[ -n "$preview_dir" && "$preview_dir" != "/" && "$preview_dir" != "$repo_root" ]] || \
    die "unsafe preview working directory: $preview_dir"
}

remove_preview_dir() {
  validate_preview_target
  if [[ -e "$preview_dir" ]]; then
    [[ -f "$preview_dir/.molgr-review-preview" ]] || die \
      "refusing to remove an unrecognized directory without the preview marker: $preview_dir"
    rm -rf -- "$preview_dir"
  fi
}

container_exists() {
  docker container inspect "$container_name" >/dev/null 2>&1
}

container_is_preview() {
  [[ "$(docker container inspect --format '{{ index .Config.Labels "org.molgr.review-preview" }}' "$container_name" 2>/dev/null || true)" == "true" ]]
}

if ((clean)); then
  validate_preview_target
  if [[ -e "$preview_dir" && ! -f "$preview_dir/.molgr-review-preview" ]]; then
    die "refusing to clean an unrecognized directory without the preview marker: $preview_dir"
  fi
  if command -v docker >/dev/null 2>&1 && container_exists; then
    container_is_preview || die "refusing to remove an existing non-preview container: $container_name"
    docker rm -f "$container_name" >/dev/null
  fi
  remove_preview_dir
  printf 'Removed preview container: %s\n' "$container_name"
  printf 'Removed preview data: %s\n' "$preview_dir"
  exit 0
fi

command -v docker >/dev/null 2>&1 || die \
  "Docker is required. Install Docker Engine/Desktop and ensure the docker command is on PATH."
docker info >/dev/null 2>&1 || die \
  "Docker daemon is unavailable or not accessible to the current user"

if [[ -z "$source_db" ]]; then
  default_db="$repo_root/.local/molgr_review/review.sqlite"
  [[ -f "$default_db" ]] || die \
    "source DB was not provided and the default does not exist: $default_db"
  source_db="$default_db"
fi
if [[ -z "$xyz_dir" ]]; then
  default_xyz="$repo_root/.local/tmQMg/data/xyz"
  [[ -d "$default_xyz" ]] || die \
    "XYZ directory was not provided and the default does not exist: $default_xyz"
  xyz_dir="$default_xyz"
fi
if [[ -z "$source_fixtures" ]]; then
  default_fixtures="$repo_root/tests/data/reviewed/tmqmg"
  [[ -d "$default_fixtures" ]] || die \
    "fixtures directory was not provided and the default does not exist: $default_fixtures"
  source_fixtures="$default_fixtures"
fi

source_db="$(canonical_path "$source_db")"
xyz_dir="$(canonical_path "$xyz_dir")"
source_fixtures="$(canonical_path "$source_fixtures")"

[[ -f "$source_db" ]] || die "source DB does not exist: $source_db"
[[ -r "$source_db" ]] || die "source DB is not readable: $source_db"
[[ -d "$xyz_dir" ]] || die "XYZ directory does not exist: $xyz_dir"
[[ -d "$source_fixtures" ]] || die "fixtures directory does not exist: $source_fixtures"
[[ "$port" =~ ^[0-9]+$ ]] && ((port >= 1 && port <= 65535)) || die "invalid port: $port"
[[ -n "$container_name" ]] || die "container name must not be empty"
validate_preview_target

case "$source_db" in
  "$preview_dir"/*) die "source DB must be outside the preview working directory" ;;
esac
case "$source_fixtures"/ in
  "$preview_dir"/*) die "source fixtures must be outside the preview working directory" ;;
esac
case "$preview_dir"/ in
  "$source_fixtures"/*) die "preview working directory must not be inside source fixtures" ;;
esac

if container_exists || [[ -e "$preview_dir" ]]; then
  if ((!fresh)); then
    die "preview already exists. Use 'docker start $container_name' to resume it, or rerun with --fresh to replace its isolated data."
  fi
  if container_exists; then
    container_is_preview || die "refusing to remove an existing non-preview container: $container_name"
    docker rm -f "$container_name" >/dev/null
  fi
  remove_preview_dir
fi

revision="$(git -C "$repo_root" rev-parse --short HEAD)"
revision_label="$revision"
if [[ -n "$(git -C "$repo_root" status --porcelain --untracked-files=no)" ]]; then
  revision_label="${revision}-dirty"
fi
if [[ -z "$image_name" ]]; then
  image_name="molgr-review-ui-preview:${revision_label}"
fi

preview_parent="$(dirname -- "$preview_dir")"
mkdir -p "$preview_parent"
staging_dir="$(mktemp -d "$preview_parent/.molgr-review-preview.XXXXXX")"
cleanup_staging() {
  if [[ -n "${staging_dir:-}" && -d "$staging_dir" ]]; then
    rm -rf -- "$staging_dir"
  fi
}
trap cleanup_staging EXIT
mkdir -p "$staging_dir/state" "$staging_dir/fixtures"
printf '%s\n' "MolGR Review UI preview state" >"$staging_dir/.molgr-review-preview"

"$host_python" -c '
import pathlib
import sqlite3
import sys

source = pathlib.Path(sys.argv[1]).resolve()
destination = pathlib.Path(sys.argv[2]).resolve()
uri = source.as_uri() + "?mode=ro"
with sqlite3.connect(uri, uri=True) as source_db:
    with sqlite3.connect(str(destination)) as preview_db:
        source_db.backup(preview_db)
        result = preview_db.execute("PRAGMA integrity_check").fetchone()
        if not result or result[0] != "ok":
            raise SystemExit("preview SQLite integrity_check failed")
' "$source_db" "$staging_dir/state/review.sqlite"

cp -a "$source_fixtures/." "$staging_dir/fixtures/"
mv -- "$staging_dir" "$preview_dir"
staging_dir=""
trap - EXIT

preview_db="$preview_dir/state/review.sqlite"
preview_fixtures="$preview_dir/fixtures"

printf 'Building preview image: %s\n' "$image_name"
docker build \
  --file "$repo_root/Dockerfile.review" \
  --tag "$image_name" \
  --build-arg "MOLGR_BUILD_VERSION=0.0.0+preview.${revision}" \
  "$repo_root"

container_id="$(docker run --detach \
  --name "$container_name" \
  --label org.molgr.review-preview=true \
  --publish "127.0.0.1:${port}:8765" \
  --env MOLGR_REVIEW_HOST=0.0.0.0 \
  --env MOLGR_REVIEW_PORT=8765 \
  --env MOLGR_REVIEW_DB=/var/lib/molgr-review/review.sqlite \
  --env MOLGR_REVIEW_FIXTURES_DIR=/var/lib/molgr-review-fixtures \
  --env MOLGR_XYZ_DIR=/data/xyz \
  --mount "type=bind,src=$preview_dir/state,dst=/var/lib/molgr-review" \
  --mount "type=bind,src=$preview_fixtures,dst=/var/lib/molgr-review-fixtures" \
  --mount "type=bind,src=$xyz_dir,dst=/data/xyz,readonly" \
  "$image_name")"

printf '\nMolGR Review UI preview started.\n'
printf 'Preview URL:       http://127.0.0.1:%s\n' "$port"
printf 'UI git revision:   %s\n' "$revision_label"
printf 'Source DB:         %s\n' "$source_db"
printf 'Preview DB:        %s\n' "$preview_db"
printf 'Source fixtures:   %s\n' "$source_fixtures"
printf 'Preview fixtures:  %s\n' "$preview_fixtures"
printf 'XYZ (read-only):   %s\n' "$xyz_dir"
printf 'Container:         %s (%s)\n' "$container_name" "${container_id:0:12}"
printf 'Stop preview:      docker stop %q\n' "$container_name"
printf 'Resume preview:    docker start %q\n' "$container_name"
printf 'Reset preview:     %q --fresh\n' "$repo_root/scripts/review-ui-preview.sh"
printf 'Clean preview:     %q --clean\n' "$repo_root/scripts/review-ui-preview.sh"
