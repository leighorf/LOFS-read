#!/usr/bin/env bash
echo "Building LOFS from local libraries..."
echo "Loading module files"
spack load hdf5
spack load netcdf-c
spack load netcdf-fortran
spack load netcdf-cxx4
spack load openmpi
spack load h5z-zfp
make -f Makefile.orf.spack all
make -f Makefile.orf.spack install
