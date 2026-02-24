#!/bin/bash

if [ ! -d ../.git ] ; then
  if [ ! -f git_hash.h ] ; then
    echo "0" > git_hash.h
  fi
  echo "Repo not found, git hash set to zero"
  exit 0
fi

PATH=$PATH:/usr/bin
h=`git describe --abbrev=7 --dirty --long --always`
echo $h
echo "#define GIT_HASH \"$h\"" > /tmp/git_hash.h

if [[ ! -s git_hash.h || `sum git_hash.h` != `sum /tmp/git_hash.h` ]] ; then
	echo Update git_hash.h
else
	echo No change git_hash.h
fi
