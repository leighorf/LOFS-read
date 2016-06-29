module add cmake/2.8.8
module add szip
module add netcdf
xml2cmake -v 2.5.2 -clobber cm1visit-nautilus.xml
cmake .
make

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/analysis/visit-thirdparty/1.0/sles11.1_gnu4.3.4/visit/hdf5/1.8.7/linux-x86_64_gcc-4.3/lib:/sw/analysis/szip/2.1/sles11.1_intel11.1/lib
