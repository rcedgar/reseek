#!/usr/bin/env bash
set -euo pipefail

commit=$(cat vcxproj_make_commit.txt)

curl -fsSL "https://raw.githubusercontent.com/rcedgar/vcxproj_make/$commit/vcxproj_make.py" \
  > vcxproj_make.py

mkdir -p ../bin

cmd=(python3 ./vcxproj_make.py --debug --binary reseekd --ccompiler "gcc -std=gnu17")
if [ "${PORTABLE:-0}" = "1" ]; then
  cmd+=(--nonative --cppcompiler "g++ -mavx2")
fi

"${cmd[@]}"
