#!/usr/bin/env bash
echo "Building LOFS from local libraries..."
echo "Loading module files"
spack load netcdf-c
spack load h5z-zfp
cd src


if spack find --format "{compiler}" zfp | grep -q intel; then
    echo "It appears ZFP was compiled with Intel. Using Intel Makefile..."
    PREFIX=$PREFIX make -f Makefile.spack.intel clean 
    PREFIX=$PREFIX make -f Makefile.spack.intel all
    PREFIX=$PREFIX make -f Makefile.spack.intel install
elif spack find --format "{compiler}" zfp | grep -q intel; then
    PREFIX=$PREFIX make -f Makefile.spack.gcc clean 
    PREFIX=$PREFIX make -f Makefile.spack.gcc all
    PREFIX=$PREFIX make -f Makefile.spack.gcc install
else
    echo "ZFP wasn't compiled by either gcc or intel, which are currently the only supported configurations. You may have to manually modify the Makefile to work."
fi

