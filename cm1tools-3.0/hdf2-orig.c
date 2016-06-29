/* Leigh Orf cm1tools to convert buffered file-per-core cm1 hdf5 output
 * to other more useful formats
 *
 * Big rewrite June 2011. All top-level conversion code in one file, makes
 * maintenance easier. Symlinks to binary determine what is run.
 *
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "hdforf.h"
#include "errorf.h"

#define TRUE 1
#define FALSE 0


/* These vars are global, everyone uses them */
char topdir[512];
int dn;
char **timedir; 
char **nodedir;
int *times;
int ntimedirs;
int nx,ny,nz,nodex,nodey;
char firstfilename[512];

void grok_cm1hdf5_file_structure();
void hdf2nc(int argc, char *argv[], char *ncbase, int X0, int Y0, int X1, int Y1, int Z0, int Z1, int t0);
int hhmmss (int time);
void hdf2v5d();
void fix_fubar_v5d_coord (float *c, float *s,int snx, int sny, int snz);
void lflush (FILE * fp);

int main(int argc, char *argv[])
{
	char nstr[6];
	char ctime[6];
	char progname[512];
	int we_are_hdf2nc = FALSE;
	int we_are_hdf2v5d = FALSE;

	int time;
	int idum;

	strcpy(progname,argv[0]);
	if (strspn("hdf2nc",progname) == 6) we_are_hdf2nc = TRUE;
	else if (strspn("hdf2v5d",progname) == 7) we_are_hdf2v5d = TRUE;
	else
	{
		fprintf(stderr,"progname = %s\n",progname);
		ERROR_STOP("Must call program as either hdf2nc or hdf2v5d, should be symbolically linked");
	}


	if (we_are_hdf2nc) 
	{
		if (argc < 10)
		{
			fprintf(stderr,"Usage: %s [topdir] [ncbase] [X0] [Y0] [X1] [Y1] [Z0] [Z1] [time] [varname1 ... varnameX] \n",progname);
			ERROR_STOP("Insufficient arguments to hdf2nc");
		}
	}
	else if (we_are_hdf2v5d)
	{
		if (argc < 2)
		{
			fprintf(stderr,"Usage: %s [topdir] \n",progname);
			ERROR_STOP("Insufficient arguments to hdf2v5d");
		}
	}

	strcpy(topdir,argv[1]); // All programs must pass topdir as arg 1 

	grok_cm1hdf5_file_structure();
	get_first_hdf_file_name(topdir,timedir[0],nodedir[0],firstfilename); printf("First filename = %s\n",firstfilename);
	get_hdf_metadata(firstfilename,idum,idum,&nx,&ny,&nz,&nodex,&nodey); printf("ORF: nx = %i ny = %i nz = %i nodex = %i nodey = %i\n", nx,ny,nz,nodex,nodey);

	if (we_are_hdf2nc)
	{
		int X0,Y0,X1,Y1,Z0,Z1,t0;
		char ncbase[512];
		hdf2nc(argc,argv,ncbase,X0,Y0,X1,Y1,Z0,Z1,t0);
	}
	else if (we_are_hdf2v5d)
	{
		hdf2v5d();
	}
}

void grok_cm1hdf5_file_structure()
{
	int i;
	int nnodedirs;
	// ORF TEST 8/12/11
	int *alltimes;
	int ntottimes;

	ntimedirs = get_num_time_dirs(topdir); printf("ntimedirs: %i\n",ntimedirs);
	if (ntimedirs == 0) ERROR_STOP("No cm1 hdf5 files found");

	timedir = (char **)malloc(ntimedirs * sizeof(char *));
	for (i=0; i < ntimedirs; i++) timedir[i] = (char *)(malloc(256 * sizeof(char)));
	times = (int *)malloc(ntimedirs * sizeof(int));

	get_sorted_time_dirs(topdir,timedir,times,ntimedirs);

	for (i=0; i<ntimedirs; i++) printf("SORTED: %s\n",timedir[i]);

	nnodedirs =  get_num_node_dirs(topdir,timedir[0]);
	nodedir = (char **)malloc(nnodedirs * sizeof(char *));
	for (i=0; i < nnodedirs; i++) nodedir[i] = (char *)(malloc(7 * sizeof(char)));

	get_sorted_node_dirs(topdir,timedir[0],nodedir,&dn,nnodedirs);

	for (i=0; i<nnodedirs; i++) printf("SORTED: %s\n",nodedir[i]);
	printf("Spacing between node directories: %i\n",dn);

	printf("ORF: NEW: All times follow\n");
	alltimes = get_all_available_times(topdir,timedir,ntimedirs,nodedir,&ntottimes);
	for (i=0; i<ntottimes; i++)printf("TIMES (hdf2nc): %i\n",alltimes[i]);

}

void hdf2nc(int argc, char *argv[], char *ncbase, int X0, int Y0, int X1, int Y1, int Z0, int Z1, int t0)
{

#include <netcdf.h>
#define MAXVARS 100

	float *buffer,*ubuffer,*vbuffer;

	int i, nvar;
	char varname[MAXVARS][50];
	char ncfilename[512];

	int snx,sny,snz;

	int status;
	int ncid,nxid,nyid,nzid,x0id,y0id,z0id,x1id,y1id,z1id;
	int varnameid[MAXVARS];
	int dims[3];
	size_t start[3],edges[3];
	int ivar,bufsize;
	
	strcpy(ncbase,argv[2]);
	X0 = atoi(argv[3]);
	Y0 = atoi(argv[4]);
	X1 = atoi(argv[5]);
	Y1 = atoi(argv[6]);
	Z0 = atoi(argv[7]);
	Z1 = atoi(argv[8]);
	t0 = atoi(argv[9]);

	nvar = argc-10;

	for (i=0; i<nvar; i++)
	{
		strcpy(varname[i],argv[i+10]);
		fprintf(stderr,"%s \n",varname[i]);
	}
	sprintf(ncfilename,"%s.%05i.nc",ncbase,t0);
	
	snx = X1 - X0 + 1;
	sny = Y1 - Y0 + 1;
	snz = Z1 - Z0 + 1;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	edges[0] = snz;
	edges[1] = sny;
	edges[2] = snx;

	bufsize = snx * sny * snz * sizeof(float);
	ubuffer = (float *) malloc ((size_t) bufsize);
	vbuffer = (float *) malloc ((size_t) bufsize);
	buffer = (float *) malloc ((size_t) bufsize);
	//printf("snx = %i sny = %i snz = %i Bufsize = %i\n",snx,sny,snz,bufsize);
	if (buffer == NULL) ERROR_STOP("Cannot allocate buffer");

	status = nc_create (ncfilename, NC_CLOBBER, &ncid); if (status != NC_NOERR) ERROR_STOP ("nc_create failed");
	status = nc_def_dim (ncid, "nx", snx, &nxid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed"); 
	status = nc_def_dim (ncid, "ny", sny, &nyid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (ncid, "nz", snz, &nzid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_var (ncid, "X0", NC_INT, 0, dims, &x0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "Y0", NC_INT, 0, dims, &y0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "X1", NC_INT, 0, dims, &x1id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "Y1", NC_INT, 0, dims, &y1id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "Z0", NC_INT, 0, dims, &z0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "Z1", NC_INT, 0, dims, &z1id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");

	dims[0] = nzid;
	dims[1] = nyid;
	dims[2] = nxid;

	for (ivar = 0; ivar < nvar; ivar++)
	{
		status = nc_def_var (ncid, varname[ivar], NC_FLOAT, 3, dims, &(varnameid[ivar]));
		if (status != NC_NOERR) 
		{
			printf ("Cannot nc_def_var for var #%i %s\n", ivar, varname[ivar]);
			ERROR_STOP("nc_def_var failed");
		}
	}
	status = nc_enddef (ncid);
	if (status != NC_NOERR) ERROR_STOP("nc_enddef failed");

      status = nc_put_var_int (ncid,x0id,&X0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,y0id,&Y0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,x1id,&X1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,y1id,&Y1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,z0id,&Z0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,z1id,&Z1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");

	for (ivar = 0; ivar < nvar; ivar++)
	{
		fprintf(stderr,"Working on %s\n",varname[ivar]);

		if(!strcmp(varname[ivar],"hwin"))
		{
			read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,times,t0,"uinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,times,t0,"vinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			for(i=0; i<snx*sny*snz; i++)
			{
				buffer[i] = sqrt(ubuffer[i]*ubuffer[i]+vbuffer[i]*vbuffer[i]);
			}
			status = nc_put_vara_float (ncid, varnameid[ivar], start, edges, buffer);
		}
		else
		{
			read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,times,t0,varname[ivar],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			status = nc_put_vara_float (ncid, varnameid[ivar], start, edges, buffer);
		}
		if (status != NC_NOERR) 
		{
			fprintf(stderr,"Could not write variable %s at time %i to %s\n", varname[ivar],t0,ncfilename);
			ERROR_STOP("Write to netcdf file failed");
		}
	}
	status = nc_close(ncid); if (status != NC_NOERR) ERROR_STOP("Could not close netcdf file (nc_close)");
}

void hdf2v5d()
{
#define MAXVARS 100
	FILE *fp_data;
	char *datafile = "hdf2v5d.dat";
	char v5dfilename[512],cpcmd[512];
	int xoff,yoff,X0,Y0,X1,Y1,Z0,Z1,snx,sny,snz;
	int dt,t0,t1;
	int compression;
	int sbox;
	float umove,vmove;
	int nblocks;
	int *timestamp,*datestamp;
	int ivar,nv5dvars,ib,iz;
	float *buffer,*zh,*xyslice;
	float dx,dy;
	char v5dvarname[MAXVARS][10];
	char varname[MAXVARS][10];

	float dxkm,dykm;

	int projection = 1; // 1 = linear, rectangular, cylindrical equidistant (must represent dx, dy in degrees)
	int vertical = 2;   // 2 = Stretched vertical grid, explicitly indicate height AGL in array
	float proj_args[4]; // Projection arguments for our grid 
	int *nla;

	hid_t file_id,sd_id;

	if ((fp_data = fopen (datafile, "r")) == NULL) ERROR_STOP("Could not open hdf2v5d.dat");

	fscanf (fp_data, "%s", topdir); lflush (fp_data);
	fscanf (fp_data, "%s", v5dfilename); lflush (fp_data);
      fscanf (fp_data, "%d%d", &xoff, &yoff); lflush (fp_data);
	fscanf (fp_data, "%d%d%d%d%d%d", &X0, &Y0, &X1, &Y1, &Z0, &Z1); lflush (fp_data);
	fscanf (fp_data, "%i", &dt); lflush (fp_data);
	fscanf (fp_data, "%d%d", &t0, &t1); lflush (fp_data);
	fscanf (fp_data, "%d", &compression); lflush (fp_data);
	fscanf (fp_data, "%d %f %f", &sbox, &umove, &vmove); lflush (fp_data);
	ivar=0;
	while ((fscanf (fp_data, "%s", &(varname[ivar][0]))) != EOF)
	{
		if (!strcmp (varname[ivar], "uinterp"))
		{
			strcpy (v5dvarname[ivar], "U");
		}
		else if (!strcmp (varname[ivar], "vinterp"))
		{
			strcpy (v5dvarname[ivar], "V");
		}
		else if (!strcmp (varname[ivar], "winterp"))
		{
			strcpy (v5dvarname[ivar], "W");
		}
		else
		{
			strcpy (v5dvarname[ivar], varname[ivar]);
		}
		ivar++;
	}
	fclose (fp_data);
	nv5dvars = ivar;

      X0+=xoff; /* kludge - sort of */
      X1+=xoff;
      Y0+=yoff;
      Y1+=yoff;

	/* make backup dat file for future reference */
	sprintf (cpcmd, "cp hdf2v5d.dat hdf2v5d.dat.%s", v5dfilename); system (cpcmd);

	snx = X1 - X0 + 1;
	sny = Y1 - Y0 + 1;
	snz = Z1 - Z0 + 1;


	if ((zh = (float *) malloc (snz * sizeof (float))) == NULL) ERROR_STOP("Cannot allocate zh");
	if ((buffer = (float *) malloc (snx * sny * snz * sizeof (float))) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((xyslice = (float *) malloc (snx * sny * sizeof (float))) == NULL) ERROR_STOP("Cannot allocate 2D buffer slice");

	if ((file_id = H5Fopen (firstfilename, H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
	{
		fprintf(stderr,"Cannot open firstfilename which is %s, even though we have alredy opened it!\n",firstfilename);
		ERROR_STOP("Cannot open hdf file");
	}

	get1dfloat (file_id,"mesh/zh",zh,Z0,snz);
	get0dfloat (file_id, "mesh/dx", &dx);
	get0dfloat (file_id, "mesh/dy", &dy);

	if (H5Fclose(file_id) < 0) ERROR_STOP("Cannot close hdf file");

	dxkm = dx / 1000.0;
	dykm = dy / 1000.0;

	for (iz = 0; iz < snz; iz++)
	{
		zh[iz] *= 0.001;	/* km */
	}

	/* Note: proj_args 2 and 3 are increment between rows (colums) in
	 * degrees. We place our box on the equator where 1 degree lat = 1
	 * degree lon = 111 km. We could have chosen projection=0 but
	 * units are generic in that case and we cannot do trajectories */
	proj_args[0] = 0.0;
	proj_args[1] = 180.0;
	proj_args[2] = dykm / 111.0;
	proj_args[3] = dxkm / 111.0;

	/* number of levels for each variable can vary, ours never do*/
	if ((nla = (int *) malloc (nv5dvars * sizeof (int))) == NULL) ERROR_STOP("Cannot allocate nla array");
	for (ivar = 0; ivar < nv5dvars; ivar++) nla[ivar] = snz;

	nblocks = 1 + (t1 - t0) / dt;
	if ((timestamp = (int *) malloc (nblocks * sizeof (int))) == NULL) ERROR_STOP("Cannot allocate timestamp");
	if ((datestamp = (int *) malloc (nblocks * sizeof (int))) == NULL) ERROR_STOP("Cannot allocate datestamp")

	for (ib = 0; ib < nblocks; ib++)
	{
		timestamp[ib] = hhmmss (t0 + ib * dt);
		datestamp[ib] = 01001; //arbitrary
	}

	strcat(v5dfilename,".v5d");
	v5dCreate (v5dfilename, nblocks, nv5dvars, sny, snx, nla, v5dvarname, timestamp, datestamp, compression, projection, proj_args, vertical, zh);

	for (ib = 1; ib <= nblocks; ib++)
	{
		int itime;
		itime = (ib - 1) * dt + t0; printf("%05i: ",itime);
		for (ivar = 1; ivar <= nv5dvars; ivar++)
		{
			printf ("%s", varname[ivar - 1]); fflush (stdout);
			read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,times,itime,varname[ivar-1],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			fix_fubar_v5d_coord (buffer, xyslice,snx,sny,snz);
			v5dWrite (ib, ivar, buffer);
		}
		printf ("\n");
	}
	v5dClose ();
	printf("\nDone.\n");
}

int hhmmss (int time)
{
	int hh, mm, ss, val;
	ss = time % 60;
	mm = (time / 60) % 60;
	hh = time / 3600;
	val = ss + 100 * mm + 10000 * hh;
	return (val);
}

void fix_fubar_v5d_coord (float *c, float *s,int snx, int sny, int snz)
{

#define P2(x,y,mx) ((y)*(mx)+(x))

	int ix, iy, iz;
	float *cp1;
	float *s0;

	s0 = s;

	for (iz = 0; iz < snz; iz++)
	{
		cp1 = c;

		for (iy = 0; iy < sny; iy++)
			for (ix = 0; ix < snx; ix++)
				s[P2 (sny - 1 - iy, ix, sny)] = *c++;

		c = cp1;

		for (iy = 0; iy < sny; iy++)
			for (ix = 0; ix < snx; ix++)
				*c++ = *s++;

		s = s0;					/* restore slice pointer */
	}
}

void lflush (FILE * fp)
{
	while (fgetc (fp) != '\n');
}

