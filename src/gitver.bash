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

if [[ ! -s git_hash.h ]] ; then
	echo Not found git_hash.h
fi

if [[ "$old_sum" != "$new_sum" ]] ; then
	echo Sum changed old=$old_sum new=$new_sum
fi

if [[ ! -s git_hash.h ||  "$old_sum" != "$new_sum" ]] ; then
#	echo sum1=`sum git_hash.h` sum2=`sum /tmp/git_hash.h`
	echo Update git_hash.h
	/bin/mv /tmp/git_hash.h .
else
	echo No change git_hash.h
fi
