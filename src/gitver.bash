#!/bin/bash
set -euo pipefail

# repo root (works no matter where called from)
repo_root="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
out_rel="src/git_hash.h"
out="$repo_root/$out_rel"

hash="0"
if git -C "$repo_root" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  hash="$(git -C "$repo_root" rev-parse --short HEAD)"
  if ! git -C "$repo_root" diff --quiet --no-ext-diff -- \
     || ! git -C "$repo_root" diff --cached --quiet --no-ext-diff --
  then
    hash="${hash}-dirty"
  fi
fi

new="#define GIT_HASH \"${hash}\""

# only update if content differs
if [[ ! -f "$out" ]] || [[ "$(cat "$out")" != "$new" ]]; then
  mkdir -p "$(dirname "$out")"
  printf '%s\n' "$new" > "$out"
fi
