#ifndef DIRSTRUCT_H
#define DIRSTRUCT_H

#include <hdf5.h>
//for hid_t below

typedef struct dir_meta
{
	char **timedir;
	char **nodedir;
	char *firstfilename;
	char *saved_base;      //Contains the base name of the run
	char *topdir;          //Full path to '3D' directory or equivalent
	int firsttimedirindex; //Not currently utilized, but set
	double *dirtimes,*alltimes;
	int ntimedirs;
	int ntottimes;
	int nnodedirs;
	int dn;
	int regenerate_cache;
} dir_meta;

typedef struct hdf_meta
{
	int nx,ny,nz,nodex,nodey;
	int nvar_available;
	int n2dswaths;
	char **varname_available;

} hdf_meta;

// grid for user-selected data
// X,Y,Z values with respect to full saved model domain
typedef struct gridstuct
{
	int X0,Y0,X1,Y1,Z0,Z1;
	int NX,NY,NZ;
	int saved_X0,saved_Y0;
	int saved_X1,saved_Y1;
	int saved_Z0,saved_Z1;
} grid;

//Our requested Cartesian grid points for any LOFS read
typedef struct requested_cube
{
	int X0,Y0,Z0,X1,Y1,Z1;
	int NX,NY,NZ;
} requested_cube;

typedef struct meshstruct
{
	float *xhfull,*yhfull,*xffull,*yffull;
	float *xhout,*yhout,*zhout;
	float *xfout,*yfout,*zfout;
	float *zh, *zf;
	float *uh,*uf,*vh,*vf,*mh,*mf;
	float dx,dy,dz;
	float rdx,rdy,rdz;
	float umove,vmove;
	double dt;
} mesh;

typedef struct cmdline
{
	int debug;
	int verbose;
	int do_swaths;
	int do_allvars;
	int gzip;
	int use_interp;
	int use_box_offset;
	int filetype;
	int nthreads;
	int got_base;
	int optcount;
	int nvar,nvar_cmdline;
	char *histpath,*base;
	char **varname_cmdline;
	float time;
	int argc_hdf2nc_min;
} cmdline;

//void get_sorted_time_dirs    (dir_meta *dm, cmdline cmd);
//void get_sorted_node_dirs    (dir_meta *dm, cmdline cmd);
//void get_num_time_dirs       (dir_meta *dm, cmdline cmd);
//void get_num_node_dirs       (dir_meta *dm, cmdline cmd);
//void get_all_available_times (dir_meta *dm, grid *gd, cmdline cmd);
//void get_hdf_metadata (dir_meta dm, hdf_meta *hm, cmdline *cm, char *argv[], hid_t *f_id);

#endif
