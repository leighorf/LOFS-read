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
spack install h5z-zfp@develop+fortran ^hdf5@1.10.6+cxx~debug+fortran+hl+mpi+pic+shared+szip~threadsafe
echo "Installing ncview"
spack install ncview ^netcdf-c@4.7.3~dap~hdf4+mpi~parallel-netcdf+shared ^hdf5@1.10.6+cxx~debug+fortran+hl+mpi+pic+shared+szip~threadsafe
