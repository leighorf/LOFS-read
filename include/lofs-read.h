#ifndef LOFS_READ_H
#define LOFS_READ_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <netcdf.h>
#include <hdf5.h> 
#include <hdf5_hl.h>
#include <H5Zzfp.h>
#include <getopt.h>
#include "macros.h"
#include "limits.h"

void get0dint(hid_t sd_id, char *varname, int *var);
void get0dfloat(hid_t sd_id, char *varname, float *var);
void get1dint(hid_t sd_id, char *varname, int *var, int p0, int np);
void get1dfloat(hid_t sd_id, char *varname, float *var, int p0, int np);
void get1ddouble(hid_t sd_id, char *varname, double *var, int p0, int np);

#endif
