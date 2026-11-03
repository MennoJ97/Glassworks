#!/usr/bin/env bash
# Builds the site and uploads web/ to a server over rsync.
# The target (host:path) comes from GLASSWORKS_TARGET or .publish-target, which git ignores.
# If GLASSWORKS_FAMILY or .publish-family (also ignored) names the home of a family of sites
# (https://example.org), the published pages link that site's shared family.css and
# family.js in place of the local copy, and their top bar links back to it.
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

FAMILY="${GLASSWORKS_FAMILY:-$(cat .publish-family 2>/dev/null || true)}"
FAMILY="${FAMILY%/}"

./scripts/build-web.sh

# Publish from a copy, so web/ keeps the self-contained pages that git holds.
SITE=target/site
mkdir -p "$SITE"
rsync -a --delete --exclude '.DS_Store' web/ "$SITE/"
if [ -n "$FAMILY" ]; then
  find "$SITE" -name '*.html' -print0 | FAMILY="$FAMILY" xargs -0 perl -pi -e '
    BEGIN { ($host = $ENV{FAMILY}) =~ s{^https?://}{} }
    s{<link rel="stylesheet" href="(?:\.\./)*family\.css">}{<link rel="stylesheet" href="$ENV{FAMILY}/family.css">\n<script src="$ENV{FAMILY}/family.js"></script>};
    s{<!--apex-->}{<a class="apex" href="$ENV{FAMILY}">$host</a> / };
  '
  grep -rL "$FAMILY/family.css" --include='*.html' "$SITE" | grep . && { echo "Pages above were not pointed at $FAMILY" >&2; exit 1; }
fi

# --delay-updates renames every file into place at the end, so visitors never get
# new JS glue with an old .wasm. --delete drops pages the build no longer makes.
# No --chmod: macOS ships openrsync, which rejects most --chmod forms.
rsync -rlt --delete --delay-updates "$SITE/" "${TARGET%/}/"

echo "Published $(git rev-parse --short HEAD) to $TARGET"
