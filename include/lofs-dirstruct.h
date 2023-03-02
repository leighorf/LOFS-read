#ifndef DIRSTRUCT_H
#define DIRSTRUCT_H

typedef struct dir_meta
{
	char **timedir;
	char **nodedir;
	char *firstfilename;
	char *saved_base;      //Contains the base name of the run
	char *topdir;          //Full path to '3D' directory or equivalent
	int firsttimedirindex; //Not currently utilized, but set
	char *devshmdir;
	char *cachedir;//Can now live in /dev/shm or in the usual history-directory-relative spot
	double *dirtimes,*alltimes;
	int ntimedirs;
	int ntottimes;
	int nnodedirs;
	int dn;
	int regenerate_cache;
} dir_meta;

typedef struct hdf_meta
{
	int nx,ny,nz,rankx,ranky;
	int nvar_available;
	int n2dswaths;
	char **varname_available;
	char **zfpacc_LOFS_all;
	int nzfplofs;

} hdf_meta;

// grid for user-selected data
// X,Y,Z values with respect to full saved model domain
typedef struct gridstuct
{
	int X0,Y0,X1,Y1,Z0,Z1;
	float XC,YC,ZC;//for grabpoint, Cartesian location input
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
	int zfp;
	int zfplossless;
	int bitgroom1,bitgroom2,bitgroom3,bitgroom_nsd;
	int use_interp;
	int use_box_offset;
	int filetype;
	int nthreads;
	int got_base;
	int got_ncdir;
	int twodwrite;
	int write_cmd_file;
	int optcount;
	int centiseconds;
	int devshmcache;
	int checkcmd;
	int nvar,nvar_cmdline;
	int header;
	char *histpath,*base,*ncdir;
	char **varname_cmdline;
	double time; //ORF 2021-11-10
	int argc_hdf2nc_min;
} cmdline;

#endif
