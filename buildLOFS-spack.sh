#!/usr/bin/env bash

## KELTON NOTE:
## This is an install script for 
## LOFS using Spack to handle the 
## package management. I can't
## seem to figure out how to do this
## natively in Spack just yet but this
## will do the trick for now. I think.

## This assumes the spack  binary is 
## already in your user path. If it isn't,
## add '. /path/to/spack/share/setup-env.sh' 
## to your ~/.bash_profile or ~/.bashrc



if [ "$1" == "--local" ]; then
    echo "Building LOFS from local libraries..."
    echo "Loading module files"
    spack load hdf5
    spack load netcdf
    spack load netcdf-fortran
    spack load netcdf-cxx4
    spack load openmpi
    make -f Makefile.spack.local all
    make -f Makefile.spack.local install
else
## Install HDF5 first
    echo "Installing LOFS Dependencies"
    echo "Installing OpenMPI"
    spack install openmpi
    echo "Installing HDF5"
    spack install hdf5@1.10.3+cxx~debug+fortran+hl+mpi+pic+shared+szip~threadsafe 
    echo "Installing NetCDF"
    spack install netcdf@4.6.1~dap~hdf4+mpi~parallel-netcdf+shared ^hdf5@1.10.3+cxx~debug+fortran+hl+mpi+pic+shared+szip~threadsafe
    echo "Installing NetCDF C++ and Fortran APIs"
    spack install netcdf-cxx4@4.3.0 ^hdf5@1.10.3+cxx~debug+fortran+hl+mpi+pic+shared+szip~threadsafe
    spack install netcdf-fortran@4.4.4 ^hdf5@1.10.3+cxx~debug+fortran+hl+mpi+pic+shared+szip~threadsafe
    echo "Installing ZFP Plugin for HDF5"
    spack install h5z-zfp@develop+fortran ^hdf5@1.10.3+cxx~debug+fortran+hl+mpi+pic+shared+szip~threadsafe


    echo "Loading module files"
    spack load hdf5
    spack load netcdf
    spack load netcdf-fortran
    spack load netcdf-cxx4
    spack load openmpi

    echo "Building LOFS"
    make -f Makefile.spack all
    make -f Makefile.spack install
fi
