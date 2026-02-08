#!/usr/bin/env bash
set -euo pipefail

commit=$(cat vcxproj_make_commit.txt)

curl -fsSL "https://raw.githubusercontent.com/rcedgar/vcxproj_make/$commit/vcxproj_make.py" \
  > vcxproj_make.py

mkdir -p ../bin

portable_opt=""
if [ "${PORTABLE:-0}" = "1" ]; then
  portable_opt="--nonative"
fi

python3 ./vcxproj_make.py --profile --binary reseekp --ccompiler "gcc -std=gnu17" $portable_opt
