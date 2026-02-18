#!/bin/bash

if [ ! -d ../.git ] ; then
	echo "Repo not found, git hash set to zero"
	hash=0
else
	PATH=$PATH:/usr/bin
	hash=$(git rev-parse --short HEAD)$([ -n "$(git status --porcelain)" ] && echo "-dirty")
fi
echo "#define GIT_HASH \"$hash\"" > git_hash.h
cat git_hash.h
