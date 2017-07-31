# v5d.c, binio.c are from the vis5d+ distribution
#
# You will need to set these appropriately for your system

NETCDF_DIR=/usr/local/ncview/io-netcdf4-hdf5-gun/lib
#HDF5_DIR=/u/staff/sisneros/visit/trunk/visit/hdf5/1.8.7/linux-x86_64_gcc-4.3
#HDF5_DIR=/home/orf/build/visit/visit/hdf5/1.8.7/linux-x86_64_gcc-4.4

NETCDFLIBS = -lnetcdf -lm -lhdf5_hl
#HDF5LIBS = -lhdf5
CC = /usr/local/ncview/io-netcdf4-hdf5-gnu/bin/h5cc
LN = /usr/local/ncview/io-netcdf4-hdf5-gnu/bin/h5cc
AR = ar
CFLAGS = -O3 -fPIC
FFLAGS = -O3
LARGE=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE

all: hdf2 readtraj hdf2d2nc libcm.a

hdf2: hdf2.o parsedir.o hdfio.o readmult.o binio.o v5d.o
	$(CC) -o hdf2 hdf2.o parsedir.o hdfio.o readmult.o binio.o v5d.o -L${NETCDF_DIR}/lib ${NETCDFLIBS}
readtraj: readtraj.o hdfio.o 
	$(CC) -o readtraj readtraj.o hdfio.o 
hdf2d2nc: hdf2d2nc.o hdfio.o
	$(CC) $(CFLAGS) -o hdf2d2nc hdf2d2nc.o hdfio.o -L${NETCDF_DIR}/lib  ${NETCDFLIBS}
libcm.a: parsedir.o hdfio.o readmult.o
	$(AR) cr libcm.a parsedir.o hdfio.o readmult.o
readmult.o: readmult.c
	$(CC) $(CFLAGS) -I -c readmult.c
hdfio.o: hdfio.c
	$(CC) $(CFLAGS) -c hdfio.c
parsedir.o: parsedir.c
	$(CC) $(CFLAGS) -c parsedir.c
hdf2.o: hdf2.c
	$(CC) $(CFLAGS) -I$(NETCDF_DIR)/include -c hdf2.c
binio.o: binio.c
	$(CC) $(LARGE) -c binio.c
v5d.o: v5d.c
	$(CC) $(LARGE) -c v5d.c
readtraj.o: readtraj.c
	$(CC) -c readtraj.c
hdf2d2nc.o: hdf2d2nc.c
	$(CC) $(CFLAGS) -I${NETCDF_DIR}/include -c hdf2d2nc.c
clean:
	rm -rf *.o hdf2 readtraj hdf2d2nc libcm.a
install: hdf2 hdf2d2nc readtraj
	install -m 0755 hdf2 /home/orf/bin
	install -m 0755 hdf2d2nc /home/orf/bin
	install -m 0755 readtraj /home/orf/bin
	cd /usr/local/bin && ln -fs hdf2 hdf2nc
	cd /usr/local/bin && ln -fs hdf2 hdf2v5d

.PHONY: install
