#!/usr/bin/env bash
echo "Building LOFS from local libraries..."
echo "Loading module files"
spack load netcdf-c
spack load h5z-zfp
cd src
PREFIX=$PREFIX make -f Makefile.spack clean 
PREFIX=$PREFIX make -f Makefile.spack all
PREFIX=$PREFIX make -f Makefile.spack install
