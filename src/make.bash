mkdir -p ../bin

python3 $src/vcxproj_make/vcxproj_make.py --openmp --bash 2> make.stderr
rc=$?

echo
echo
echo '=== rc=$rc tail make.stderr ==='
tail make.stderr
echo
echo
ls -lh ../bin/reseek
