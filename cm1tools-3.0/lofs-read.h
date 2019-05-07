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

#define ERROR_STOP(string) { fprintf(stderr," *** Fatal error in %s line %i: %s\n",__FILE__,__LINE__,string); exit(0); }
#define ERROR_WARN(string) fprintf(stderr," *** Warning: %s line %i: %s\n",__FILE__,__LINE__,string);

#define P3(x,y,z,mx,my) (((z)*(mx)*(my))+((y)*(mx))+(x))
#define P2(x,y,mx) (((y)*(mx))+(x))
#define TRUE (1)
#define FALSE (0)

void get_hdf_metadata(char *base, int *nx, int *ny, int *nz, int *nodex, int *nodey);
void read_hdf_mult_md(float *gf, char *topdir, char **timedir, char **nodedir, int ntimedirs, int dn,
    double *dirtimes, double *alltimes, int ntottimes, double dtime, char *varname,
    int gx0, int gy0, int gxf, int gyf, int gz0, int gzf,
    int nx, int ny, int nz, int nodex, int nodey);
void open_cm1_hdf_file(hid_t *sd_id, char *base, int itime, int inode);
void close_cm1_hdf_file(hid_t sd_id);
void get0dint(hid_t sd_id, char *varname, int *var);
void get0dfloat(hid_t sd_id, char *varname, float *var);
void get1dint(hid_t sd_id, char *varname, int *var, int p0, int np);
void get1dfloat(hid_t sd_id, char *varname, float *var, int p0, int np);
void get1ddouble(hid_t sd_id, char *varname, double *var, int p0, int np);
void get2dfloat (hid_t file_id, char *varname, float *var, int y0, int ny, int x0, int nx);
void put0dint(hid_t sd_id, char *varname, int *var);
void put0dfloat(hid_t sd_id, char *varname, float *var);
void put1dfloat(hid_t sd_id, char *varname, float *var, int p0, int np);
int isNumeric (const char * s);
void sortarray(char **strarray,int nel,int csize);
void get_sorted_node_dirs(char *topdir, char *timedir,char **nodedir, int *dn, int nnodedirs, int debug,int regenerate_cahce);
void get_sorted_time_dirs(char *topdir, char **timedir, double *times, int ntimedirs, int debug,int regenerate_cahce);
int get_num_time_dirs(char *topdir, int debug,int regenerate_cahce);
int get_num_node_dirs(char *topdir, char *timedir, int debug,int regenerate_cahce);
double * get_all_available_times (char *topdir, char **timedir, int ntimedirs, char **nodedir, int nnodedirs, int *ntottimes,char *firstfilename, int *firsttimedirindex, int *saved_snx0, int *saved_sny0, int *saved_snx1, int *saved_sny1, int debug,int regenerate_cahce);
