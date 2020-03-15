#ifndef HDF2NC_H
#define HDF2NC_H

#include <stdlib.h>
//Above is for size_t, not picked up from lofs-read.h

//Trashcan struct for all things netcdf I guess
typedef struct netcdf_struct
{
	int nxh_dimid,nyh_dimid,nzh_dimid;
	int nxf_dimid,nyf_dimid,nzf_dimid,time_dimid,timeid;
	int thsfcid,dbzsfcid;
	int x0id,y0id,z0id,x1id,y1id,z1id;
	int xhid,yhid,zhid;
	int xfid,yfid,zfid;
	int ncid;
	int *varnameid;
	char **varname;
	char *ncfilename;
	int dims[4],start[4],edges[4],d2[3];
	size_t s2[3],e2[3];
	int u0id,v0id,pres0id,pi0id,th0id,qv0id;
	int twodslice;
} ncstruct;

typedef struct sounding
{
	float *u0,*v0,*pres0,*pi0,*th0,*qv0;
} sounding;

typedef struct buffers
{
	float *ustag, *vstag, *wstag;
	float *twodfield,*twodfield0;
	float *twodbuf,*twodbuf0;
	float *buf,*buf0;
	float *dum0,*dum00;
	float *dum1,*dum10;
	float *writeptr;
} buffers;

typedef struct readahead
{
	int u,v,w;
	int vortmag,hvort,streamvort;//Not really readahead, used for mallocs
} readahead;


void parse_cmdline_hdf2nc(int argc, char *argv[], cmdline *cmd, dir_meta *dm, grid *gd);
void get_saved_base(char *timedir, char *saved_base);

/*
void read_hdf_mult_md(float *gf, char *topdir, char **timedir, char **nodedir, int ntimedirs, int dn,
    double *dirtimes, double *alltimes, int ntottimes, double dtime, char *varname,
    int gx0, int gy0, int gxf, int gyf, int gz0, int gzf,
    int nx, int ny, int nz, int nodex, int nodey);

*/

#endif
