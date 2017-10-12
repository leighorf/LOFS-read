#!/bin/bash
module unload PrgEnv-cray
module unload PrgEnv-gnu
module unload PrgEnv-pgi

module unload darshan
module load visit

make clean

rm -rf CMakeCache.txt CMakeFiles/ cmake_install.cmake Makefile CMakeLists.txt
xml2cmake cm1visit-bw-2.12.0.xml

mv CMakeLists.txt temp
head -n 65 temp > CMakeLists.txt
cat newLine.txt >> CMakeLists.txt
tail -n 33 temp >> CMakeLists.txt
rm temp

cmake -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ .
make
