#ifndef HDF2NC_H
#define HDF2NC_H

#include <stdlib.h>
#include <hdf5.h>
#include <dirent.h>
#include <sys/stat.h>

//Above is for size_t, not picked up from lofs-read.h

//Trashcan struct for all things netcdf I guess
typedef struct netcdf_struct
{
	int nxh_dimid,nyh_dimid,nzh_dimid;
	int nxf_dimid,nyf_dimid,nzf_dimid,time_dimid,timeid;
	int thsfcid,dbzsfcid;
	int umoveid,vmoveid;
	int x0id,y0id,z0id,x1id,y1id,z1id;
	int xhid,yhid,zhid;
	int xfid,yfid,zfid;
	int ncid;
	int *varnameid;
	char **varname;
	char *ncfilename;
	int dims[4],d2[3];
	size_t start[4],edges[4];//not used anymore, se do_requested_variables
	int u0id,v0id,pres0id,pi0id,th0id,qv0id,rho0id;
	int twodslice;
} ncstruct;

typedef struct sounding
{
	float *u0,*v0,*pres0,*pi0,*th0,*qv0,*rho0;
} sounding;

typedef struct buffers
{
	float *ustag, *vstag, *wstag;
	float *ppert, *thrhopert;
	float *buf0, *buf, *dum0, *dum1;
	float *threedbuf;
} buffers;

typedef struct readahead
{
	int u,v,w,ppert,thrhopert;
	int vortmag,hvort,streamvort;//Not really readahead, used for mallocs
	int budgets;
} readahead;
void dealloc_structs(cmdline *cmd,dir_meta *dm, grid *gd,ncstruct *nc, readahead *rh);
void parse_cmdline_hdf2nc(int argc, char *argv[], cmdline *cmd, dir_meta *dm, grid *gd);
void get_saved_base(char *timedir, char *saved_base);
void init_structs(cmdline *cmd,dir_meta *dm, grid *gd,ncstruct *nc, readahead *rh);
void get_num_time_dirs (dir_meta *dm,cmdline cmd);
void get_sorted_time_dirs (dir_meta *dm,cmdline cmd);
void get_num_node_dirs (dir_meta *dm,cmdline cmd);
void get_sorted_node_dirs (dir_meta *dm,cmdline cmd);
void get_all_available_times (dir_meta *dm, grid *gd, cmdline cmd);
void get_hdf_metadata(dir_meta dm, hdf_meta *hm, cmdline *cmd, char *argv[], hid_t *f_id);
void set_span(grid *gd,hdf_meta hm,cmdline cmd);
void allocate_1d_arrays(hdf_meta hm, grid gd, mesh *msh, sounding *snd);
void set_1d_arrays(hdf_meta hm, grid gd, mesh *msh, sounding *snd, hid_t *f_id);
void set_netcdf_attributes(ncstruct *nc, grid gd, cmdline *cmd, buffers *b, hdf_meta *hm, hid_t *f_id);
void nc_write_1d_data (ncstruct nc, grid gd, mesh msh, sounding snd, cmdline cmd);
void set_readahead(readahead *rh,ncstruct nc, cmdline cmd);
void malloc_3D_arrays (buffers *b, grid gd, readahead rh,cmdline cmd);
void free_3D_arrays (buffers *b, grid gd, readahead rh,cmdline cmd);
void do_the_swaths(hdf_meta hm, ncstruct nc, dir_meta dm, grid gd, cmdline cmd);
void do_readahead(buffers *b,grid gd,readahead rh,dir_meta dm,hdf_meta hm,cmdline cmd);
void do_requested_variables(buffers *b, ncstruct nc, grid gd, mesh msh, sounding *snd, readahead rh,dir_meta dm,hdf_meta hm,cmdline cmd);
void read_lofs_buffer(float *buf, char *varname, dir_meta dm, hdf_meta hm, requested_cube rc, cmdline cmd);
void copy_grid_to_requested_cube (requested_cube *rc, grid gd);
void sortchararray (char **strarray, int nel);
void compress_with_nccopy(ncstruct nc,cmdline cmd);
void write_hdf2nc_command_txtfile(int argc, char *argv[],ncstruct nc);

#endif
