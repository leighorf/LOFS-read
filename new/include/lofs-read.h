/*
    This file is part of LOFS, written by Leigh Orf (http://orf.media)
*/

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
void get_sorted_node_dirs(char *topdir, char *timedir,char **nodedir, int *dn, int nnodedirs, int debug,int regenerate_cahce);
void get_sorted_time_dirs(char *topdir, char **timedir, double *times, int ntimedirs, int debug,int regenerate_cahce);
int get_num_time_dirs(char *topdir, int debug,int regenerate_cahce);
int get_num_node_dirs(char *topdir, char *timedir, int debug,int regenerate_cahce);
double * get_all_available_times (char *topdir, char **timedir, int ntimedirs, char **nodedir,
		int nnodedirs, int *ntottimes,char *firstfilename, int *firsttimedirindex,
		int *saved_snx0, int *saved_sny0, int *saved_snx1, int *saved_sny1, int debug,int regenerate_cahce);
