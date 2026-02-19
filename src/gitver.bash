#!/bin/bash

set -euo pipefail

PATH=$PATH:/usr/bin

out_rel="src/git_hash.h"

repo_root="$(git rev-parse --show-toplevel 2>/dev/null || true)"
if [[ -z "$repo_root" ]]; then
  echo "Repo not found, git hash set to zero"
  hash="0"
else
  hash="$(git -C "$repo_root" rev-parse --short HEAD)"

  # Dirty if any tracked changes EXCEPT the generated header
  if ! git -C "$repo_root" diff --quiet --no-ext-diff -- . ":(exclude)$out_rel" \
     || ! git -C "$repo_root" diff --cached --quiet --no-ext-diff -- . ":(exclude)$out_rel"
  then
    hash="${hash}-dirty"
  fi
fi

new="#define GIT_HASH \"${hash}\""
out_abs="$repo_root/$out_rel"

# Only rewrite if content actually changes (prevents churn)
if [[ ! -f "$out_abs" ]] || [[ "$(cat "$out_abs")" != "$new" ]]; then
  printf '%s\n' "$new" > "$out_abs"
fi
