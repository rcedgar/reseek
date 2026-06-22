#!/bin/bash

if [ ! -d ../.git ] ; then
  if [ ! -f git_hash.h ] ; then
    echo "0" > git_hash.h
  fi
#  echo "Repo not found, git hash set to zero"
  exit 0
fi

PATH=$PATH:/usr/bin
h=`git describe --abbrev=7 --dirty --long --always`
echo $h
echo "#define GIT_HASH \"$h\"" > /tmp/git_hash.h

old_sum=`cat git_hash.h | sum`
new_sum=`cat /tmp/git_hash.h | sum`

if [[ ! -s git_hash.h ||  "$old_sum" != "$new_sum" ]] ; then
	/bin/mv /tmp/git_hash.h .
fi
