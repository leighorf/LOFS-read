# v5d.c, binio.c are from the vis5d+ distribution
#
# You will need to set these appropriately for your system

#NETCDF_DIR=/usr/local
#HDF5_DIR=/u/staff/sisneros/visit/trunk/visit/hdf5/1.8.7/linux-x86_64_gcc-4.3
#HDF5_DIR=/home/orf/build/visit/visit/hdf5/1.8.7/linux-x86_64_gcc-4.4

#NETC
#HDF5LIBS = -lhdf5
CC = icc
LN = cc
AR = ar
CFLAGS = -O3 -fPIC
FFLAGS = -O3
LARGE=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE

all: hdf2 readtraj hdf2d2nc libcm.a

hdf2: hdf2.o parsedir.o hdfio.o readmult.o binio.o v5d.o
	$(CC) -o hdf2 hdf2.o parsedir.o hdfio.o readmult.o binio.o v5d.o $(NETCDF_CLIB) $(HDF5_CLIB)
readtraj: readtraj.o hdfio.o 
	$(CC) -o readtraj readtraj.o hdfio.o  $(HDF5_CLIB)
hdf2d2nc: hdf2d2nc.o hdfio.o
	$(CC) $(CFLAGS) -o hdf2d2nc hdf2d2nc.o hdfio.o $(NETCDF_CLIB) $(HDF5_CLIB)
libcm.a: parsedir.o hdfio.o readmult.o
	$(AR) cr libcm.a parsedir.o hdfio.o readmult.o
readmult.o: readmult.c
	$(CC) $(CFLAGS) $(HDF5_CLIB) -c readmult.c
hdfio.o: hdfio.c
	$(CC) $(CFLAGS) $(HDF5_CLIB) -c hdfio.c
parsedir.o: parsedir.c
	$(CC) $(CFLAGS) $(HDF5_CLIB) -c parsedir.c
hdf2.o: hdf2.c
	$(CC) $(CFLAGS) $(NETCDF_CLIB) $(HDF5_CLIB) -c hdf2.c
binio.o: binio.c
	$(CC) $(CFLAGS) $(LARGE) -c binio.c
v5d.o: v5d.c
	$(CC) $(CFLAGS) $(LARGE) -c v5d.c
readtraj.o: readtraj.c
	$(CC) -c $(HDF5_CLIB) readtraj.c
hdf2d2nc.o: hdf2d2nc.c
	$(CC) $(CFLAGS) $(NETCDF_CLIB) $(HDF5_CLIB) -c hdf2d2nc.c
clean:
	rm -rf *.o hdf2 hdf2d2nc libcm.a
install: hdf2 readtraj
	install -m 0755 hdf2 $(HOME)/bin
	install -m 0755 readtraj $(HOME)/bin
	cd $(HOME)/bin && ln -fs hdf2 hdf2nc
	cd $(HOME)/bin && ln -fs hdf2 hdf2v5d
	cd $(HOME)/bin && ln -fs hdf2 makevisit

.PHONY: install
