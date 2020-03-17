#spack load hdf5-zfp netcdf-c first

CC = mpicc
CFLAGS = -fopenmp -g -fPIC 
LINKOPTS  = -lh5zzfp -lzfp -lhdf5_hl -lhdf5 -lsz -lz -lm -lnetcdf -fopenmp 
SRC=hdf2nc.c hdf2nc-util.c parse_cmdline_hdf2nc.c parsedir.c hdfio.c readlofs.c do_requested_variables.c

OBJS = $(addsuffix .o, $(basename $(SRC)))

all: hdf2nc

hdf2nc: $(OBJS)
	$(CC) $(CFLAGS) -o hdf2nc  $(OBJS) $(LINKOPTS)
.c.o:	
	$(CC) $(CFLAGS) -c $*.c
clean:
	rm -rf *.o hdf2nc
install: hdf2nc
	install -m 0755 hdf2nc $(HOME)/bin
