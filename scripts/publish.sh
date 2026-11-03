#!/usr/bin/env bash
# Builds the site and uploads web/ to a server over rsync.
# The target (host:path) comes from GLASSWORKS_TARGET or .publish-target, which git ignores.
set -euo pipefail
cd "$(dirname "$0")/.."

TARGET="${GLASSWORKS_TARGET:-$(cat .publish-target 2>/dev/null || true)}"
if [ -z "$TARGET" ]; then
  echo "Set GLASSWORKS_TARGET or write host:/path/to/site/ into .publish-target" >&2
  exit 1
fi

# Publish only what's committed, so the site always matches a commit.
if [ -n "$(git status --porcelain)" ]; then
  echo "Uncommitted changes; commit or stash them first" >&2
  exit 1
fi

./scripts/build-web.sh

# --delay-updates renames every file into place at the end, so visitors never get
# new JS glue with an old .wasm. --delete drops pages the build no longer makes.
# No --chmod: macOS ships openrsync, which rejects most --chmod forms.
rsync -rlt --delete --delay-updates --exclude '.DS_Store' web/ "${TARGET%/}/"

echo "Published $(git rev-parse --short HEAD) to $TARGET"
