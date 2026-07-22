mkdir -p ../bin

./gitver.bash
python3 ../py/vcxproj_make_0c8d9d1.py --openmp --x86-64-v3 --bash 2> make.stderr
rc=$?

echo
echo
echo '=== tail make.stderr ==='
cat make.stderr \
	| grep -v "lto-wrapper: warning: using serial compilation" \
	| grep -v "lto-wrapper: note: see the" \
	| grep -v " warning: Using .dlopen. in statically linked applications" \
	| tail
echo
echo
if [ $rc == 0 ] ; then
	echo SUCCESS
else
	echo ERROR
fi
echo

ls -lh ../bin/reseek
