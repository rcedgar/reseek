curl -fsSL https://raw.githubusercontent.com/rcedgar/vcxproj_make/d27afda/vcxproj_make.py \
  > vcxproj_make.py

mkdir -p ../bin

python3 ./vcxproj_make.py --debug
