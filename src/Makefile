CC = mpicc
AR = ar
CFLAGS = -fopenmp -O2 -fPIC 

LINKOPTS  = -lh5zzfp -lzfp -lhdf5_hl -lhdf5 -lsz -lz -lm -lnetcdf -fopenmp 

all: hdf2nc 

hdf2nc: hdf2nc.o hdf2nc-util.o parse_cmdline_hdf2nc.o parsedir.o hdfio.o readlofs.o
	$(CC) $(CFLAGS) -o hdf2nc hdf2nc.o hdf2nc-util.o parse_cmdline_hdf2nc.o parsedir.o hdfio.o readlofs.o $(LINKOPTS)
parse_cmdline_hdf2nc.o: parse_cmdline_hdf2nc.c
	$(CC) $(CFLAGS) -c parse_cmdline_hdf2nc.c
parsedir.o: parsedir.c
	$(CC) $(CFLAGS) -c parsedir.c
readlofs.o: readlofs.c
	$(CC) $(CFLAGS) -c readlofs.c
hdfio.o: hdfio.c
	$(CC) $(CFLAGS) -c hdfio.c
hdf2nc.o: hdf2nc.c
	$(CC) $(CFLAGS) -c hdf2nc.c
hdf2nc-util.o: hdf2nc-util.c
	$(CC) $(CFLAGS) -c hdf2nc-util.c
clean:
	rm -rf *.o hdf2nc
install: hdf2nc
	install -m 0755 hdf2nc $(HOME)/bin

.PHONY: install
