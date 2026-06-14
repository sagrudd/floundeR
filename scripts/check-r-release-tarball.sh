#!/usr/bin/env sh
set -eu

ROOT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
CHECK_ROOT="${TMPDIR:-/tmp}/flounder-release-check"

rm -rf "$CHECK_ROOT"
mkdir -p "$CHECK_ROOT"

cd "$CHECK_ROOT"
R CMD build "$ROOT_DIR"

TARBALL=$(find "$CHECK_ROOT" -maxdepth 1 -type f -name 'floundeR_*.tar.gz' | sort | tail -n 1)

if [ -z "$TARBALL" ]; then
  echo "No floundeR source tarball was produced." >&2
  exit 1
fi

R CMD check --no-manual "$TARBALL"
