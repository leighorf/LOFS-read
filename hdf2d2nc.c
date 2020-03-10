#include "lofs-read.h"

FILE *fp_data;
int n2d;
int i2d;
const char **varname;
int  **varid;

char *replace_str(char *str, char *orig, char *rep)
{
  static char buffer[4096];
  char *p;

  if(!(p = strstr(str, orig)))  // Is 'orig' even in 'str'?
    return str;

  strncpy(buffer, str, p-str); // Copy characters from 'str' start to 'orig' st$
  buffer[p-str] = '\0';

  sprintf(buffer+(p-str), "%s%s", rep, p+strlen(orig));

  return buffer;
}

herr_t first_pass(hid_t loc_id, const char *name, void *opdata)
{
    H5G_stat_t statbuf;
    H5Gget_objinfo(loc_id, name, FALSE, &statbuf);
    n2d++;
    return 0;
}

herr_t second_pass(hid_t loc_id, const char *name, void *opdata)
{
    H5G_stat_t statbuf;
    H5Gget_objinfo(loc_id, name, FALSE, &statbuf);
    printf("Got %s\n",name);
    strcpy((char *)varname[i2d],name);
    i2d++;
    return 0;
}


int main (int argc, char *argv[])
{
	/* hdf vars */

	char hdf5filename[100], *ncfilename, *tmpfilename;
	float dx, dy;
	float *twodvar;
	int status;
	int nx,ny,nz;

	// Should nuke most of these don't need them after iterate //
	

	char fullvarname[100];

	float *xh,*yh,*xf,*yf,*th0,*pi0,*qv0,*u0,*v0;
	int xh_id,yh_id,th0_id,pi0_id,qv0_id,u0_id,v0_id;
	int xf_id,yf_id;
	int dx_id, dy_id;

	hid_t file_id;
	int ncid,nxid,nyid,nzid;
	int unxid,vnyid,wnzid;
	size_t start[3],edges[3];
	int dims[3];

	printf("argc = %i\n",argc);
	if (argc < 2)
	{
		fprintf(stderr,"Usage: %s [2Dcm1hdf5file]",argv[0]);
		exit(0);
	}

	strcpy(hdf5filename,argv[1]);
	ncfilename = (char *)malloc(100*sizeof(char));
	tmpfilename = (char *)malloc(100*sizeof(char));
	tmpfilename = replace_str(hdf5filename,".2Dcm1hdf5","_2D.nc");
	if(argc==3) strcpy(ncfilename,argv[2]); else strcpy(ncfilename,".");
	strcat(ncfilename,"/");
	strcat(ncfilename,tmpfilename);

	file_id = H5Fopen (hdf5filename, H5F_ACC_RDONLY,H5P_DEFAULT);


	get0dint (file_id, "/grid/nx", &nx);
	get0dint (file_id, "/grid/ny", &ny);
	get0dint (file_id, "/grid/nz", &nz);
	get0dfloat (file_id, "/mesh/dx", &dx);
	get0dfloat (file_id, "/mesh/dy", &dy);
//	get0dfloat (file_id, "/shear/angle",&angle);
//	get0dfloat (file_id, "/shear/dl1",&dl1);
//	get0dfloat (file_id, "/shear/dl2",&dl2);
//	get0dfloat (file_id, "/shear/sl1",&sl1);
//	get0dfloat (file_id, "/shear/sl2",&sl2);
//	get0dfloat (file_id, "/time",&time);
	printf("nx = %3i ny = %3i nz = %3i\n",nx,ny,nz);
//	printf("dx = %f dy = %f angle = %f\n",dx,dy,angle);
//	printf("dl1 = %f dl2 = %f sl1 = %f sl2 = %f time = %f\n",dl1,dl2,sl1,sl2,time);

	twodvar = (float *) malloc (nx * ny * sizeof(float));
//	zh = (float *) malloc (nz * sizeof(float));
	yh = (float *) malloc (ny * sizeof(float));
	xh = (float *) malloc (nx * sizeof(float));
//	zf = (float *) malloc ((nz+1) * sizeof(float));
	yf = (float *) malloc ((ny+1) * sizeof(float));
	xf = (float *) malloc ((nx+1) * sizeof(float));
	th0 = (float *) malloc (nz * sizeof(float));
	pi0 = (float *) malloc (nz * sizeof(float));
	qv0 = (float *) malloc (nz * sizeof(float));
	u0 = (float *) malloc (nz * sizeof(float));
	v0 = (float *) malloc (nz * sizeof(float));


//      get1dfloat(file_id,"/mesh/zh",zh,0,nz);
 //     get1dfloat(file_id,"/mesh/zf",zf,0,nz+1);

	get1dfloat(file_id, "/mesh/xh",xh,0,nx);
	get1dfloat(file_id, "/mesh/xf",xf,0,nx); //ORF 6.4.12 was nx+1
	get1dfloat(file_id, "/mesh/yh",yh,0,ny); // ORF BUG IN CM1 10/5/12 CHANGE JUST TO GET CODE TO RUN
	get1dfloat(file_id, "/mesh/yf",yf,0,ny); //ORF 6.4.12 was ny+1
	get1dfloat(file_id, "/basestate/u0",u0,0,nz);
	get1dfloat(file_id, "/basestate/v0",v0,0,nz);
	get1dfloat(file_id, "/basestate/pi0",pi0,0,nz);
	get1dfloat(file_id, "/basestate/th0",th0,0,nz);
	get1dfloat(file_id, "/basestate/qv0",qv0,0,nz);

	status = nc_create (ncfilename, NC_CLOBBER|NC_64BIT_OFFSET , &ncid); if (status != NC_NOERR) printf ("ERROR: nc_create\n");
	status = nc_def_dim (ncid, "nx", nx, &nxid); if (status != NC_NOERR) printf ("ERROR: nc_def_dim 1\n");
	status = nc_def_dim (ncid, "nxu", nx+1, &unxid); if (status != NC_NOERR) printf ("ERROR: nc_def_dim 1\n");
	status = nc_def_dim (ncid, "ny", ny, &nyid); if (status != NC_NOERR) printf ("ERROR: nc_def_dim 2\n");
	status = nc_def_dim (ncid, "nyv", ny+1, &vnyid); if (status != NC_NOERR) printf ("ERROR: nc_def_dim 2\n");
	status = nc_def_dim (ncid, "nz", nz, &nzid); if (status != NC_NOERR) printf ("ERROR: nc_def_dim 3\n");
	status = nc_def_dim (ncid, "nzw", nz+1, &wnzid); if (status != NC_NOERR) printf ("ERROR: nc_def_dim 3\n");

	status = nc_def_var (ncid, "dx", NC_FLOAT, 0, dims, &dx_id);
	status = nc_def_var (ncid, "dy", NC_FLOAT, 0, dims, &dy_id);
//	status = nc_def_var (ncid, "angle", NC_FLOAT, 0, dims, &angle_id);
//	status = nc_def_var (ncid, "dl1", NC_FLOAT, 0, dims, &dl1_id);
//	status = nc_def_var (ncid, "dl2", NC_FLOAT, 0, dims, &dl2_id);
//	status = nc_def_var (ncid, "sl1", NC_FLOAT, 0, dims, &sl1_id);
//	status = nc_def_var (ncid, "sl2", NC_FLOAT, 0, dims, &sl2_id);
//	status = nc_def_var (ncid, "time", NC_FLOAT, 0, dims, &time_id);

//	dims[0] = nzid; status = nc_def_var (ncid, "zh", NC_FLOAT, 1, dims, &zh_id);
	dims[0] = nyid; status = nc_def_var (ncid, "yh", NC_FLOAT, 1, dims, &yh_id);
	dims[0] = nxid; status = nc_def_var (ncid, "xh", NC_FLOAT, 1, dims, &xh_id);
//	dims[0] = wnzid; status = nc_def_var (ncid, "zf", NC_FLOAT, 1, dims, &zf_id);
	dims[0] = vnyid; status = nc_def_var (ncid, "yf", NC_FLOAT, 1, dims, &yf_id);
	dims[0] = unxid; status = nc_def_var (ncid, "xf", NC_FLOAT, 1, dims, &xf_id);
	dims[0] = nzid; status = nc_def_var (ncid, "th0", NC_FLOAT, 1, dims, &th0_id);
	dims[0] = nzid; status = nc_def_var (ncid, "pi0", NC_FLOAT, 1, dims, &pi0_id);
	dims[0] = nzid; status = nc_def_var (ncid, "qv0", NC_FLOAT, 1, dims, &qv0_id);
	dims[0] = nzid; status = nc_def_var (ncid, "u0", NC_FLOAT, 1, dims, &u0_id);
	dims[0] = nzid; status = nc_def_var (ncid, "v0", NC_FLOAT, 1, dims, &v0_id);

	n2d=0;
#ifdef OLD_FORMAT
	H5Giterate(file_id, "/2dfull",NULL,first_pass,NULL);
#else
	H5Giterate(file_id, "/2Dfull",NULL,first_pass,NULL);
#endif

	varname = (const char **)malloc(n2d*sizeof(char *));
	varid =   (int **)       malloc(n2d*sizeof(int  *));

	for (i2d=0; i2d<n2d; i2d++)
	{
		varname[i2d] = (char *)malloc(50*sizeof(char)); // 50 characters per variable
		varid  [i2d] = (int *) malloc(sizeof(int));
	}

	i2d=0;
#ifdef OLD_FORMAT
	H5Giterate(file_id, "/2dfull",NULL,second_pass,NULL);
#else
	H5Giterate(file_id, "/2Dfull",NULL,second_pass,NULL);
#endif
	dims[0] = nyid; dims[1] = nxid;
	for (i2d=0; i2d<n2d; i2d++)
	{
		status = nc_def_var (ncid, varname[i2d], NC_FLOAT, 2, dims, varid[i2d]);
		//status = nc_def_var (ncid, "rain", NC_FLOAT, 2, dims, varid[i2d]);
		if (status != NC_NOERR)
		{
			  printf("status = %i\n",status);
			  printf("Goddammit: %s\n",nc_strerror(status));
			  printf ("nc_def_var %s\n",varname[i2d]);
		}
	}
	status = nc_enddef(ncid);
	if (status != NC_NOERR)
	{
		  printf("status = %i\n",status);
		  printf("Goddammit: %s\n",nc_strerror(status));
		  printf ("nc_enddef\n");
	}

	status = nc_put_var_float (ncid, dx_id, &dx);
	status = nc_put_var_float (ncid, dy_id, &dy);
//	status = nc_put_var_float (ncid, angle_id, &angle);
//	status = nc_put_var_float (ncid, dl1_id, &dl1);
//	status = nc_put_var_float (ncid, dl2_id, &dl2);
//	status = nc_put_var_float (ncid, sl1_id, &sl1);
//	status = nc_put_var_float (ncid, sl2_id, &sl2);
//	status = nc_put_var_float (ncid, time_id, &time);

//	start[0] = 0; edges[0] = nz; status = nc_put_vara_float (ncid, zh_id, start, edges, zh);
//	start[0] = 0; edges[0] = nz+1; status = nc_put_vara_float (ncid, zf_id, start, edges, zf);
	start[0] = 0; edges[0] = nz; status = nc_put_vara_float (ncid, th0_id, start, edges, th0);
	start[0] = 0; edges[0] = nz; status = nc_put_vara_float (ncid, qv0_id, start, edges, qv0);
	start[0] = 0; edges[0] = nz; status = nc_put_vara_float (ncid, pi0_id, start, edges, pi0);
	start[0] = 0; edges[0] = nz; status = nc_put_vara_float (ncid, u0_id, start, edges, u0);
	start[0] = 0; edges[0] = nz; status = nc_put_vara_float (ncid, v0_id, start, edges, v0);
	start[0] = 0; edges[0] = ny; status = nc_put_vara_float (ncid, yh_id, start, edges, yh);
	start[0] = 0; edges[0] = nx; status = nc_put_vara_float (ncid, xh_id, start, edges, xh);
	start[0] = 0; edges[0] = ny+1; status = nc_put_vara_float (ncid, yf_id, start, edges, yf);
	start[0] = 0; edges[0] = nx+1; status = nc_put_vara_float (ncid, xf_id, start, edges, xf);

	start[0] = 0; start[1] = 0;
	edges[0] = ny; edges[1] = nx;

	// Loop over 2d vars, read & write
	for (i2d=0; i2d<n2d; i2d++)
	{
#ifdef OLD_FORMAT
		strcpy(fullvarname,"/2dfull/");
#else
		strcpy(fullvarname,"/2Dfull/");
#endif
		strcat(fullvarname,varname[i2d]);
		get2dfloat(file_id,fullvarname,twodvar,0,ny,0,nx); // ORF FIX BUG
		status = nc_put_vara_float (ncid, *(varid[i2d]),start,edges,twodvar);
		if (status != NC_NOERR)
		{
			  printf("status = %i\n",status);
			  printf("Goddammit: %s\n",nc_strerror(status));
			  printf ("nc_put_vara_float\n");
		}
	}

	H5Fclose(file_id);

	nc_close(ncid);
	fprintf(stderr,"\nDone.\n");
	exit (0);
}
