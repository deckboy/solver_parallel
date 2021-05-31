#!/bin/sh

module load intel

mkdir ../libraries

cd ../files/lis-1.7.17

./configure --enable-f90 --enable-omp --prefix=$(pwd)/../../libraries/

make

make install
