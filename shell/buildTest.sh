#!/bin/bash

#Assume this is called in the LADEL directory
curdir=`pwd`

#Build direcetories
if [ ! -d "build" ]; then
  mkdir build
fi

if [ ! -d "build/debug" ]; then
  mkdir build/debug
fi

if [ ! -d "build/lib" ]; then
  mkdir build/lib
fi

builddir=$curdir/build/debug
cd $builddir

cmake $curdir -DCMAKE_BUILD_TYPE=debug -DUNITTESTS=ON -DSIMPLE_COL_COUNTS=ON -DAMD=1
make
ctest -VV
# valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose bin/run_all_tests