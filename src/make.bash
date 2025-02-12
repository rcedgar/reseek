mkdir -p ../bin

python3 $src/vcxproj_make/vcxproj_make.py --bash 2> make.stderr

tail make.stderr
