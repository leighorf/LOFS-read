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

##ORF NOTE:

## I have updated this script so the latest versions of code are installed, and no
## duplicates occur. Also installing ncview. This works as of 2020-03-10

## module load has an accompanying module list command. spack load does not.
## Fuck you spack load. If you do module load you need to specify the full name that is
## listed with module list. spack load lets you just pick the app.


if [ "$1" != "--local" ]; then
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
    spack install hdf5@1.10.6+cxx~debug+fortran+hl+mpi+pic+shared+szip~threadsafe 
    echo "Installing NetCDF"
    spack install netcdf-c@4.7.3~dap~hdf4+mpi~parallel-netcdf+shared ^hdf5@1.10.6+cxx~debug+fortran+hl+mpi+pic+shared+szip~threadsafe
    echo "Installing NetCDF C++ and Fortran APIs"
    spack install netcdf-cxx4@4.3.1 ^hdf5@1.10.6+cxx~debug+fortran+hl+mpi+pic+shared+szip~threadsafe
    spack install netcdf-fortran@4.5.2 ^hdf5@1.10.6+cxx~debug+fortran+hl+mpi+pic+shared+szip~threadsafe
    echo "Installing ZFP Plugin for HDF5"
    spack install h5z-zfp@0.8.0+fortran ^hdf5@1.10.6+cxx~debug+fortran+hl+mpi+pic+shared+szip~threadsafe
	spack install ncview ^netcdf-c@4.7.3~dap~hdf4+mpi~parallel-netcdf+shared


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
