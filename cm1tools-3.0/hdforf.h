/*
    This file is part of cm1tools, written by Leigh Orf (http://orf5.com)

    cm1tools is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cm1tools is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cm1tools.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#define EATSHIT
#ifdef EATSHIT
#include <hdf5.h> 
#endif
void get_hdf_metadata(char *base, int node, int itime, int *nx, int *ny, int *nz, int *nodex, int *nodey);
void read_hdf_mult_md(float *gf, char *topdir, char **timedir, char **nodedir, int ntimedirs, int dn,
    int *dirtimes, int *alltimes, int ntottimes, int itime, char *varname,
    int gx0, int gy0, int gxf, int gyf, int gz0, int gzf,
    int nx, int ny, int nz, int nodex, int nodey);
void link_hdf_files (char *topdir, char **timedir, char **nodedir, int ntimedirs, int dn, int *dirtimes, int *alltimes, int ntottimes,
	int gx0, int gy0, int gxf, int gyf,
	int nx, int ny, int nz, int nodex, int nodey);
void open_cm1_hdf_file(hid_t *sd_id, char *base, int itime, int inode);
void close_cm1_hdf_file(hid_t sd_id);
void get0dint(hid_t sd_id, char *varname, int *var);
void get0dfloat(hid_t sd_id, char *varname, float *var);
void get1dint(hid_t sd_id, char *varname, int *var, int p0, int np);
void get1dfloat(hid_t sd_id, char *varname, float *var, int p0, int np);
void put0dint(hid_t sd_id, char *varname, int *var);
void put0dfloat(hid_t sd_id, char *varname, float *var);
void put1dfloat(hid_t sd_id, char *varname, float *var, int p0, int np);

int isNumeric (const char * s);
static int cmpstringp(const void *p1, const void *p2);
void sortarray(char **strarray,int nel,int csize);
void get_sorted_node_dirs(char *topdir, char *timedir,char **nodedir, int *dn, int nnodedirs);
void get_sorted_time_dirs(char *topdir, char **timedir, int *times, int ntimedirs, char *base);
int get_num_time_dirs(char *topdir);
int get_num_node_dirs(char *topdir, char *timedir);
void get_first_hdf_file_name(char *topdir, char *timedir, char *nodedir, char *filename);
int * get_all_available_times (char *topdir, char **timedir, int ntimedirs, char **nodedir, int nnodedirs, int *ntottimes,char *firstfilename, int *firsttimedirindex);
