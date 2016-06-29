#include <stdio.h>
#include <math.h>
#include <string.h>
#include "hdforf.h"
#include <netcdf.h>

#define MAXVARS 20

#define P2(x,y,mx) ((y)*(mx)+(x))
#define P3(x,y,z,mx,my) (((z)*(mx)*(my))+((y)*(mx))+(x))

/* Leigh Orf 8/31/09 */

/* hdf2nc: Convert distributed cm1hdf5 format to single netcdf4 file */

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

FILE *fp_data;

int snx, sny, snz;

int X0, X1;
int Y0, Y1;
int Z0, Z1;


main (int argc, char *argv[])
{
	/* hdf vars */

	char *datafile = "hdf2nc-mult.dat";
	char cpcmd[100];

	char basename[100], ncfilename[100], ncbase[100];
	int X0, X1, Y0, Y1, Z0, Z1;	/* cube parameters */
	float dx, dy, dz;
	float *dzc, *dze, *zcntr, *ubar, *vbar;
	float *xh,*yh,*zh;
	float ubox, vbox;
	float dxkm, dykm;
	int dt, t0, t1;
	int ivar = 0;
	int i, ix, iy, iz, count, ib, nblocks, status;
	int nx,ny,nz,nodex,nodey;

	/* 3d data to be written to v5d file */
	float *datacube, *xyslice;
	float *datacube2;
	float *diffcube,*d0;
	float *u, *v;
	float *qr, *qs, *qg, *pi, *th, *dbz, *tmp;
	float *qrz, *qsz, *qgz, *piz, *thz;
	float dbz_sfc, dbz_max, cnor, cnos, cnoh, rwdn, swdn, hwdn;

	int arsize,nvar;
	int firsthdf;
	char varname[MAXVARS][20];
	//ORF TEST
	hid_t group_id;
//	H5G_info_t group_info;
	char name[100];
	//END ORF TEST
	hid_t sd_id;
// from b2nc.c
	int ncid,nxid,nyid,nzid,x0id,y0id,z0id,x1id,y1id,z1id;
	int xhid,yhid,zhid;
	int varnameid[MAXVARS];
	int anx,any,anz;
	int twod=0;
	size_t start[3],edges[3];
	int dims[3];

	if (argc < 10)
	{
		fprintf(stderr,"Usage: %s [hdfbase] [ncbase] [X0] [Y0] [X1] [Y1] [Z0] [Z1] [time]\n",argv[0]);
		exit(0);
	}

	strcpy(basename,argv[1]);
	strcpy(ncbase,argv[2]);
	X0 = atoi(argv[3]);
	Y0 = atoi(argv[4]);
	X1 = atoi(argv[5]);
	Y1 = atoi(argv[6]);
	Z0 = atoi(argv[7]);
	Z1 = atoi(argv[8]);
	t0 = atoi(argv[9]);

	fprintf(stderr,"%s %s %i %i %i %i %i %i %i\n", basename,ncbase,X0,Y0,X1,Y1,Z0,Z1,t0);

	for (i=0; i < argc-10; i++)
	{
	    strcpy(varname[i],argv[i+10]);
	    fprintf(stderr,"%s \n",varname[i]);
	}
	nvar = argc - 10;
	printf("nvar = %i\n",nvar);

	sprintf(ncfilename,"%s.%05i.nc",ncbase,t0);

	snx = X1 - X0 + 1;
	sny = Y1 - Y0 + 1;
	snz = Z1 - Z0 + 1;

	datacube = (float *) malloc ((size_t) (snx * sny * snz * sizeof (float)));
	if (datacube == NULL) die_here(3000,"hdf2b","Cannot allocate","datacube");


	status = nc_create (ncfilename, NC_CLOBBER, &ncid); if (status != NC_NOERR) printf ("ERROR: nc_create\n");
//	status = nc_create (ncfilename, NC_NETCDF4, &ncid); if (status != NC_NOERR) printf ("ERROR: nc_create\n");
	status = nc_def_dim (ncid, "nx", snx, &nxid); if (status != NC_NOERR) printf ("ERROR: nc_def_dim 1\n");
	status = nc_def_dim (ncid, "ny", sny, &nyid); if (status != NC_NOERR) printf ("ERROR: nc_def_dim 2\n");
	status = nc_def_dim (ncid, "nz", snz, &nzid); if (status != NC_NOERR) printf ("ERROR: nc_def_dim 3\n");
	status = nc_def_var (ncid, "X0", NC_INT, 0, dims, &x0id); if (status != NC_NOERR) printf ("ERROR: nc_def_var\n");
	status = nc_def_var (ncid, "Y0", NC_INT, 0, dims, &y0id); if (status != NC_NOERR) printf ("ERROR: nc_def_var\n");
	status = nc_def_var (ncid, "X1", NC_INT, 0, dims, &x1id); if (status != NC_NOERR) printf ("ERROR: nc_def_var\n");
	status = nc_def_var (ncid, "Y1", NC_INT, 0, dims, &y1id); if (status != NC_NOERR) printf ("ERROR: nc_def_var\n");
	status = nc_def_var (ncid, "Z0", NC_INT, 0, dims, &z0id); if (status != NC_NOERR) printf ("ERROR: nc_def_var\n");
	status = nc_def_var (ncid, "Z1", NC_INT, 0, dims, &z1id); if (status != NC_NOERR) printf ("ERROR: nc_def_var\n");




//	dims[0] = nxid;
//	status = nc_def_var (ncid, "xh", NC_FLOAT, 1, dims, &xhid); if (status != NC_NOERR) printf ("ERROR: nc_def_var\n");
//	dims[0] = nyid;
//	status = nc_def_var (ncid, "yh", NC_FLOAT, 1, dims, &yhid); if (status != NC_NOERR) printf ("ERROR: nc_def_var\n");
//	dims[0] = nzid;
//	status = nc_def_var (ncid, "zh", NC_FLOAT, 1, dims, &zhid); if (status != NC_NOERR) printf ("ERROR: nc_def_var\n");

	dims[0] = nzid;
	dims[1] = nyid;
	dims[2] = nxid;

	for (ivar = 0; ivar < nvar; ivar++)
	{
		status = nc_def_var (ncid, varname[ivar], NC_FLOAT, 3, dims, &(varnameid[ivar]));
		if (status != NC_NOERR)
			printf ("Cannot nc_def_var for var #%i %s\n", ivar, varname[ivar]);
			status=nc_def_var_deflate(ncid, byte_varid, 0, 1, DEFLATE_LEVEL_1)
		if (status != NC_NOERR)
			printf ("Cannot nc_def_var_deflate for var #%i %s\n", ivar, varname[ivar]);
	}
	status = nc_enddef (ncid);
	if (status != NC_NOERR)
	{
		printf("status = %i\n",status);
		printf ("nc_enddef\n");
	}

      status = nc_put_var_int (ncid,x0id,&X0);
      status = nc_put_var_int (ncid,y0id,&Y0);
      status = nc_put_var_int (ncid,x1id,&X1);
      status = nc_put_var_int (ncid,y1id,&Y1);
      status = nc_put_var_int (ncid,z0id,&Z0);
      status = nc_put_var_int (ncid,z1id,&Z1);

	firsthdf = first_hdf_index(basename, t0);
	get_hdf_metadata(basename,firsthdf,t0,&nx,&ny,&nz,&nodex,&nodey);
	open_cm1_hdf_file (&sd_id, basename, t0, firsthdf);
//	get0dfloat (sd_id, "mesh/dx", &dx);
//	get0dfloat (sd_id, "mesh/dy", &dy);

	//xh = (float *) malloc ((size_t) (snx * sizeof (float)));
//	yh = (float *) malloc ((size_t) (sny * sizeof (float)));
	zh = (float *) malloc ((size_t) (snz * sizeof (float)));

	printf("snx = %i\n",snx);
	//ORF 5/12/11
	//This only works when you are operating on a single HDF5 file
	//which spans the full domain
//	get1dfloat (sd_id,"mesh/xh", xh, X0, snx);
//	get1dfloat (sd_id,"mesh/yh", yh, Y0, sny);
//	get1dfloat (sd_id,"mesh/zh", zh, Z0, snz);

//	start[0] = 0; edges[0] = snx;
//	status = nc_put_vara_float (ncid, xhid, start, edges, xh);
//	edges[0] = sny;
//	status = nc_put_vara_float (ncid, yhid, start, edges, yh);
//	edges[0] = snz;
//	status = nc_put_vara_float (ncid, zhid, start, edges, zh);

	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	edges[0] = snz;
	edges[1] = sny;
	edges[2] = snx;

	/**************
	//ORF TEST
	// mesh first
	group_id = H5Gopen(sd_id,"/mesh",H5P_DEFAULT);
	H5Gget_info(group_id,&group_info);
	printf("There are %i links in the mesh group\n",group_info.nlinks);
	for (i=0; i < group_info.nlinks; i++)
	{
		H5Lget_name_by_idx(group_id,".",H5_INDEX_NAME,H5_ITER_INC,i,name,20,H5P_DEFAULT);
		printf("Name = %s\n",name);
	}
	H5Gclose(group_id);

	// 2d vars
	group_id = H5Gopen(sd_id,"/2d",H5P_DEFAULT);
	H5Gget_info(group_id,&group_info);
	printf("There are %i links in the 2d var group\n",group_info.nlinks);
	for (i=0; i < group_info.nlinks; i++)
	{
		H5Lget_name_by_idx(group_id,".",H5_INDEX_NAME,H5_ITER_INC,i,name,20,H5P_DEFAULT);
		printf("Name = %s\n",name);
	}
	H5Gclose(group_id);

	// 3d vars
	group_id = H5Gopen(sd_id,"/3d",H5P_DEFAULT);
	H5Gget_info(group_id,&group_info);
	printf("There are %i links in the 3d group\n",group_info.nlinks);
	for (i=0; i < group_info.nlinks; i++)
	{
		H5Lget_name_by_idx(group_id,".",H5_INDEX_NAME,H5_ITER_INC,i,name,20,H5P_DEFAULT);
		printf("Name = %s\n",name);
	}
	H5Gclose(group_id);
	************/

	H5Fclose(sd_id);

	for (ivar = 0; ivar < nvar; ivar++)
	{
		fprintf(stderr,"Working on %s\n",varname[ivar]);
		if (!strcmp (varname[ivar], "2Dswaths"))
		{
			read_hdf_mult (&(datacube[0]), basename, t0, "sfcrain_trans", X0, Y0, X1, Y1, 0, 0,nx,ny,1,nodex,nodey);
			read_hdf_mult (&(datacube[1*snx*sny]), basename, t0, "maxsws_trans", X0, Y0, X1, Y1, 0, 0,nx,ny,1,nodex,nodey);
			read_hdf_mult (&(datacube[2*snx*sny]), basename, t0, "maxw_trans 5km", X0, Y0, X1, Y1, 0, 0,nx,ny,1,nodex,nodey);
			read_hdf_mult (&(datacube[3*snx*sny]), basename, t0, "maxw_trans 1km", X0, Y0, X1, Y1, 0, 0,nx,ny,1,nodex,nodey);
			read_hdf_mult (&(datacube[4*snx*sny]), basename, t0, "maxzeta_trans 1km", X0, Y0, X1, Y1, 0, 0,nx,ny,1,nodex,nodey);
			read_hdf_mult (&(datacube[5*snx*sny]), basename, t0, "maxzeta_trans sfc", X0, Y0, X1, Y1, 0, 0,nx,ny,1,nodex,nodey);
		}
		else
		{
			fprintf(stderr,"\n Reading %s\n",basename);
			read_hdf_mult (datacube, basename, t0, varname[ivar], X0, Y0, X1, Y1, Z0, Z1,nx,ny,nz,nodex,nodey);
		}
		status = nc_put_vara_float (ncid, varnameid[ivar], start, edges, datacube);
	}
	nc_close(ncid);
	fprintf(stderr,"\nDone.\n");
	exit (0);
}
