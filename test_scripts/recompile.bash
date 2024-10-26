#!/bin/bash -e

echo
echo recompile.bash
echo

cd ../src
rm -rf o/ ../bin/reseek*
./build_linux_x86.bash
