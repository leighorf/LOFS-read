#ifndef HDF2NC_H
#define HDF2NC_H

#include <stdlib.h>
#include <hdf5.h>
#include <dirent.h>
#include <sys/stat.h>

//Above is for size_t, not picked up from lofs-read.h

typedef struct var3d_struct
{
	char varname[MAXSTR];
	int is_LOFS_var;
	double zfpacc_LOFS;
	double zfpacc_netcdf;
	int varnameid;
} var3dstruct;

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
	var3dstruct *var3d;
	char *ncfilename;
	int dims[4],d2[3];
	size_t start[4],edges[4];//not used anymore, see do_requested_variables
	int u0id,v0id,pres0id,pi0id,th0id,thv0id,qv0id,rho0id;
	int twodslice;
} ncstruct;

/* ORF 2022-08-23 
 * My solution to the whole ZFP thing is to just create zfp accuracy
 * variables that have a specific naming scheme (identical to varname).
 * Could be oh so much more elegant if I could do some object oriented
 * stuff... so everything can be specified at the command line, with
 * sane defaults in the init code.
 * The LOFS parameters are read only - they will be read from the LOFS
 * data, as I save the metadata for those.*/

typedef struct lofs_zfp
{
	float u,v,w,uinterp,vinterp,winterp;
	float prespert,thrhopert,dbz;
	float qc,qi,qr,qg,qs;
	float nci,ncg,ncr,ncs;
	/* next line NSSL microphysics only */
	float qhl,ccn,ccw,crw,cci,csw,chw,chl,vhw,vhl;
	float qvpert,thpert,th,prs;
	float pi,pipert,rho,rhopert;
	float tke_sg,km,kh,qv;
	/* I don't save derived quantities mucn anymore,
	 * calculate them on the fly is the way to go but
	 * here is vorticity anyway */
	float xvort,yvort,zvort,vortmag;
} lofs;

typedef struct netcdf_zfp
{
	/* first, CM1 / LOFS saved variables, identical to above */
	float u,v,w;
	float prespert,thrhopert,dbz;
	float qc,qi,qr,qg,qs;
	float nci,ncg,ncr,ncs;
	/* next line NSSL microphysics only */
	float qhl,ccn,ccw,crw,cci,csw,chw,chl,vhw,vhl;
	float qvpert,thpert,th,prs;
	float pi,pipert,rho,rhopert;
	float tke_sg,km,kh,qv;

	/* Now, derived variables... add at your leisure */
	/* For instance, budget stuff will need to be added here */
	/* For now I only include what I'm using at the moment */
	float uinterp,vinterp,winterp;
	float xvort,yvort,zvort,vortmag;
	float hwin_sr,hwin_gr,windmag_sr;

	float wb_buoy,ub_pgrad,vb_pgrad,wb_pgrad;
	float xvort_stretch,yvort_stretch,zvort_stretch;
	float xvort_baro,yvort_baro;
	float xvort_solenoid,yvort_solenoid,zvort_solenoid;
	float hvort,streamvort,qiqvpert,qtot,tempC;
	float hdiv;
} netcdf;

typedef struct zfp_acc
{
	lofs *lofs;
	netcdf *netcdf;
}zfpacc;

typedef struct sounding
{
	float *u0,*v0,*pres0,*pi0,*th0,*thv0,*qv0,*rho0;
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
	int vortmag,hvort;
	int streamvort,qiqvpert,qtot;
	int tempC;
	int budgets;
} readahead;

typedef struct point
{
	float XC,YC,ZC;//The Cartesian locations that lie on the appropriate mesh (not interpolated)
	float val;
} point;

void dealloc_structs(cmdline *cmd,dir_meta *dm, grid *gd,ncstruct *nc, readahead *rh);
void parse_cmdline_hdf2nc(int argc, char *argv[], cmdline *cmd, dir_meta *dm, grid *gd, zfpacc *zfpacc);
void get_saved_base(char *timedir, char *saved_base);
void init_structs(cmdline *cmd,dir_meta *dm, grid *gd,ncstruct *nc, readahead *rh, hdf_meta *hm, zfpacc *zfpacc);
void get_num_time_dirs (dir_meta *dm,cmdline cmd);
void get_sorted_time_dirs (dir_meta *dm,cmdline cmd);
void get_num_node_dirs (dir_meta *dm,cmdline cmd);
void get_sorted_node_dirs (dir_meta *dm,cmdline cmd);
void get_all_available_times (dir_meta *dm, grid *gd, cmdline cmd);
void get_hdf_metadata(dir_meta dm, hdf_meta *hm, cmdline *cmd, ncstruct *nc, char *argv[], hid_t *f_id, zfpacc *zfpacc);
void set_span(grid *gd,hdf_meta hm,cmdline cmd);
void allocate_1d_arrays(hdf_meta hm, grid gd, mesh *msh, sounding *snd);
void set_1d_arrays(hdf_meta hm, grid gd, mesh *msh, sounding *snd, hid_t *f_id);
void set_netcdf_attributes(ncstruct *nc, grid gd, cmdline *cmd, buffers *b, hdf_meta *hm, hid_t *f_id, zfpacc *zfpacc);
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
int mkdir_p(const char *path);
void add_CM1_LOFS_zfp_metadata_to_netcdf_file (hdf_meta *hm, hid_t *f_id, ncstruct nc);
void parse_cmdline_grabpoint(int argc, char *argv[], cmdline *cmd, dir_meta *dm, grid *gd, zfpacc *zfpacc);
float grabpoint(grid *gd,hdf_meta hm,dir_meta dm, cmdline cmd,mesh msh, char *varname);

#endif
