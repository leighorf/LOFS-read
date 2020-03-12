#ifndef HDF2NC_H
#define HDF2NC_H
typedef struct cmdline
{
	/* MUST INITIALIZE THESE */
	int debug;
	int do_swaths;
	int do_allvars;
	int gzip;
	int use_interp;
	int use_box_offset;
	int filetype;
	int nthreads;
	int got_base;
	int optcount;
	char *histpath,*base;
	float time;

} cmdline;

/*
void parse_cmdline_hdf2nc(int argc, char *argv[],
		char *histpath, char *base, int *got_base,
		double *time, int *X0, int *Y0, int *X1, int *Y1, int *Z0, int *Z1, int *regenerate_cache);
*/
void parse_cmdline_hdf2nc(int argc, char *argv[], cmdline *cmd, dir_meta *dm, grid *gd);

void get_saved_base(char *timedir, char *saved_base);

void read_hdf_mult_md(float *gf, char *topdir, char **timedir, char **nodedir, int ntimedirs, int dn,
    double *dirtimes, double *alltimes, int ntottimes, double dtime, char *varname,
    int gx0, int gy0, int gxf, int gyf, int gz0, int gzf,
    int nx, int ny, int nz, int nodex, int nodey);

void get_hdf_metadata(char *base, int *nx, int *ny, int *nz, int *nodex, int *nodey);

#endif
