#!/bin/bash
set -euo pipefail

out="git_hash.h"

if [[ ! -d ../.git ]]; then
  echo "Repo not found, git hash set to zero"
  hash="0"
else
  PATH="$PATH:/usr/bin"

  hash="$(git rev-parse --short HEAD)"

  # Mark dirty only if tracked files differ (ignores untracked)
  if ! git diff --quiet --no-ext-diff --; then
    hash="${hash}-dirty"
  fi
fi

new="#define GIT_HASH \"${hash}\""

# Only rewrite if content changed
if [[ ! -f "$out" ]] || [[ "$(cat "$out")" != "$new" ]]; then
  printf '%s\n' "$new" > "$out"
fi

cat "$out"
