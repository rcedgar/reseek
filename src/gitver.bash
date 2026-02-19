#!/bin/bash

if [ ! -d ../.git ] ; then
  if [ ! -f git_hash.h ] ; then
    echo "0" > git_hash.h
  fi
  echo "Repo not found, git hash set to zero"
  exit 0
fi

PATH=$PATH:/usr/bin

hold=`cat git_hash.h`
h=`git describe --abbrev=7 --dirty --long --always`
h=$h
hold=$hold
if [[ $h == $hold ]] ; then
    echo Same version githash=$h
else
    echo New version githash=$h
    echo "#define GIT_HASH \"$h\"" > git_hash.h
fi
