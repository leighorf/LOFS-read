module unload darshan

module use /sw/bw/thg/modulefiles
module load phdf5
module load h5zzfp
#module swap cce cce/8.3.14
module add cray-netcdf

#make -f Makefile.zfp clean
#make -f Makefile.zfp install
cmake .
make
