#!/usr/bin/env bash
echo "Building LOFS from local libraries..."
echo "Loading module files"
spack load netcdf-c
spack load h5z-zfp
cd src
make -f Makefile.spack clean
make -f Makefile.spack all
make -f Makefile.spack install
