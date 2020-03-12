#ifndef HDF2NC_H
#define HDF2NC_H

void parse_cmdline_hdf2nc(int argc, char *argv[], cmdline *cmd, dir_meta *dm, grid *gd);
void get_saved_base(char *timedir, char *saved_base);

/*
void read_hdf_mult_md(float *gf, char *topdir, char **timedir, char **nodedir, int ntimedirs, int dn,
    double *dirtimes, double *alltimes, int ntottimes, double dtime, char *varname,
    int gx0, int gy0, int gxf, int gyf, int gz0, int gzf,
    int nx, int ny, int nz, int nodex, int nodey);

void get_hdf_metadata(char *base, int *nx, int *ny, int *nz, int *nodex, int *nodey);
*/

#endif
