curl -fsSL https://raw.githubusercontent.com/rcedgar/vcxproj_make/15594f1/vcxproj_make.py \
  > vcxproj_make.py

mkdir -p ../bin

python3 ./vcxproj_make.py --openmp
