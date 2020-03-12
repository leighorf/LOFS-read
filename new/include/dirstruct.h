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
	double *dirtimes,*alltimes;
	int ntimedirs;
	int notottimes;
	int dn;
	int regenerate_cache;
} dir_meta;

typedef struct hdf_meta
{
	int nx,ny,nz,nodex,nodey;
} hdf_meta;

// grid for user-selected data
// X,Y,Z values with respect to full saved model domain
typedef struct gridstuct
{
	int X0,Y0,X1,Y1,Z0,Z1;
	int saved_X0,saved_Y0;
	int saved_X1,saved_Y1;
	int saved_Z0,saved_Z1;
	float umove,vmove;
} grid;

#endif
