#!/bin/bash
module unload darshan

if [ "$1" == "--local" ]; then
    echo "Building LOFS from local libraries..."
    module load cray-hdf5-parallel
    module load cray-netcdf-hdf5parallel
    make -f Makefile.bw.local all
    make -f Makefile.bw.local install
else
    echo "Building LOFS from system libraries..."
    module use /sw/bw/thg/modulefiles
    module load phdf5/1.8.18
    module load zfp
    module load h5zzfp/1.8.18
    module load szip
    #module load cray-netcdf
    module load cray-netcdf-hdf5parallel
    make -f Makefile.bw install
fi
