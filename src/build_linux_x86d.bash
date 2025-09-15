commit=`cat vcxproj_make_commit.txt`

curl -fsSL https://raw.githubusercontent.com/rcedgar/vcxproj_make/$commit/vcxproj_make.py \
  > vcxproj_make.py

mkdir -p ../bin

python3 ./vcxproj_make.py --debug -binary reseekd
