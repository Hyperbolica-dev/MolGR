#!/usr/bin/env sh
set -eu

if [ -z "${GITHUB_ENV:-}" ]; then
  echo "GITHUB_ENV is not set; skip proxy cleanup."
  exit 0
fi

{
  echo "HTTP_PROXY="
  echo "HTTPS_PROXY="
  echo "NO_PROXY="
  echo "http_proxy="
  echo "https_proxy="
  echo "no_proxy="
} >> "$GITHUB_ENV"
