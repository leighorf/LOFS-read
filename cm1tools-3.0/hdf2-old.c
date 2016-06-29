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
#define MAXVARS 100

/* These vars are global, everyone uses them. Or I am lazy. */
char topdir[512];
int dn;
char **timedir; 
char **nodedir;
int *dirtimes;
int ntimedirs;
int nx,ny,nz,nodex,nodey;
char firstfilename[512];
char base[100];
int nnodedirs;
int *alltimes;
int ntottimes;

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
	int we_are_makevisit = FALSE;
	int we_are_linkfiles = FALSE;

	int time;
	int idum=0;

	strcpy(progname,argv[0]);
	if (strspn("hdf2nc",progname) == 6) we_are_hdf2nc = TRUE;
	else if (strspn("hdf2v5d",progname) == 7) we_are_hdf2v5d = TRUE;
	else if (strspn("makevisit",progname) == 9) we_are_makevisit = TRUE;
	else if (strspn("linkfiles",progname) == 9) we_are_linkfiles = TRUE;
	else
	{
		fprintf(stderr,"progname = %s\n",progname);
		ERROR_STOP("Must call program as either hdf2nc, hdf2v5d, makevisit, or linkfiles and should be symbolically linked");
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
	else if (we_are_makevisit)
	{
		if (argc < 2)
		{
			fprintf(stderr,"Usage: %s [topdir] \n",progname);
			ERROR_STOP("Insufficient arguments to makevisit");
		}
	}
	else if (we_are_linkfiles)
	{
		if (argc < 2)
		{
			fprintf(stderr,"Usage: %s [topdir] \n",progname);
			ERROR_STOP("Insufficient arguments to linkfiles");
		}
	}


	realpath(argv[1],topdir); // Need absolute path of 3D directory

	grok_cm1hdf5_file_structure();
	get_first_hdf_file_name(topdir,timedir[0],nodedir[0],firstfilename); printf("First filename = %s\n",firstfilename);
	get_hdf_metadata(firstfilename,idum,idum,&nx,&ny,&nz,&nodex,&nodey); //printf("ORF: nx = %i ny = %i nz = %i nodex = %i nodey = %i\n", nx,ny,nz,nodex,nodey);

	if (we_are_hdf2nc)
	{
		int X0,Y0,X1,Y1,Z0,Z1,t0;
		X0=Y0=X1=Y1=Z0=Z1=t0=0;//shut up compiler
		char ncbase[512];
		hdf2nc(argc,argv,ncbase,X0,Y0,X1,Y1,Z0,Z1,t0);
	}
	else if (we_are_hdf2v5d)
	{
		hdf2v5d(argc,argv);
	}
	else if (we_are_linkfiles)
	{
	  int snx0, snx1, sny0, sny1;
	  int idir,inode;
	  char dirname[512];
	  //adding subdomain
	  snx0 = sny0 = 0;
	  snx1 = nx-1; sny1 = ny-1;
	  if(argc == 4)
	  {
	      snx0 = atoi(argv[2]); snx1 = atoi(argv[3]);
	  }
	  else if(argc == 6)
	  {
	      snx0 = atoi(argv[2]); sny0 = atoi(argv[3]);
	      snx1 = atoi(argv[4]); sny1 = atoi(argv[5]);
	  }

		mkdir("sub3D",0755);
		for(idir=0;idir<ntimedirs;idir++)
		{
			sprintf(dirname,"sub3D/%s",timedir[idir]);
			mkdir (dirname,0755);
			for(inode=0;inode<nnodedirs;inode++)
			{
				sprintf(dirname,"sub3D/%s/%s",timedir[idir],nodedir[inode]);
				mkdir (dirname,0755);
				printf("Created %s\n",dirname);
			}
		}

		link_hdf_files (topdir, timedir, nodedir, ntimedirs, dn, dirtimes, alltimes, ntottimes,
			snx0, sny0, snx1, sny1, nx, ny, nz, nodex, nodey);

	}
	else if (we_are_makevisit)
	{
        int snx0, snx1, sny0, sny1, snz0, snz1;
        //adding subdomain
        snx0 = sny0 = snz0 = 0;
        snx1 = nx-1; sny1 = ny-1; snz1 = nz-1;
        if(argc == 4)
        {
            snx0 = atoi(argv[2]); snx1 = atoi(argv[3]);
        }
        else if(argc == 6)
        {
            snx0 = atoi(argv[2]); sny0 = atoi(argv[3]);
            snx1 = atoi(argv[4]); sny1 = atoi(argv[5]);
        }
        else if(argc == 8)
        {
            snx0 = atoi(argv[2]); sny0 = atoi(argv[3]);
            snx1 = atoi(argv[4]); sny1 = atoi(argv[5]);
            snz0 = atoi(argv[6]); snz1 = atoi(argv[7]);
        }

/* This code essentially reproduces what used to be in
 * avtcm1visitFileFormat, getting everything essential and stuffing it
 * into a file. We use hdf5 format because... we are a masochist(?).
 * Actually, because ascii text got cumbersome with large mesh arrays
 * etc.
 *
 * We need:
 * topdir
 * timedir
 * nodedir
 * ntimedirs
 * dn
 * dirtimes
 * alltimes array
 * varnames array
 * nx, ny, nz, nodex, nodey 
 *
 * Bunch of stuff already grokked once we are here
 *
 * */

		int i,rank,nvars;
		hid_t f_id,g_id,dset_id,dspace_id,strtype;
		H5G_info_t group_info;
		hsize_t dims[1];
		char dirbase[200];
		char visitfile[220];
		char groupname[220];
		char varname[MAXVARS][40]; //Yeah I know also hardcoded in avtcm1visitFileFormat.h

		float *xhfull,*yhfull,*zh;

		for (i=0; i < strlen(topdir)-2; i++) dirbase[i]=topdir[i];
		dirbase[i]='\0';

		sprintf(visitfile,"%s%s.cm1visit",dirbase,base);
		printf("visitfile = %s\n",visitfile);
		
		printf("topdir = %s\n",topdir);
		printf("ntimedirs = %i\n",ntimedirs);
		for (i=0; i<ntimedirs; i++)printf("%s ",timedir[i]);printf("\n");
		for (i=0; i<ntimedirs; i++)printf("%i ",dirtimes[i]);printf("\n");
		printf("nnodedirs = %i\n",nnodedirs);
		for (i=0; i<nnodedirs; i++)printf("%s ",nodedir[i]);printf("\n");
		printf("dn = %i\n",dn);
		printf("ntottimes = %i\n",ntottimes);
		for (i=0; i<ntottimes; i++)printf("%i ",alltimes[i]);printf("\n");
		printf("nx = %i ny = %i nz = %i nodex = %i nodey = %i\n",nx,ny,nz,nodex,nodey);

		// Done with file layout and basic metadata, now get mesh
		// and variable names

		if ((f_id = H5Fopen (firstfilename, H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
		{
			fprintf(stderr,"Cannot open firstfilename which is %s, even though we have alredy opened it!\n",firstfilename);
			ERROR_STOP("Cannot open hdf file");
		}

// Get full mesh data

		xhfull = (float *)malloc(nx * sizeof(float));
		yhfull = (float *)malloc(ny * sizeof(float));
		zh = (float *)malloc(nz * sizeof(float));
		get1dfloat (f_id,(char *)"mesh/xhfull",xhfull,0,nx);
		get1dfloat (f_id,(char *)"mesh/yhfull",yhfull,0,ny);
		get1dfloat (f_id,(char *)"mesh/zh",zh,0,nz);
		printf("Succesfully got mesh arrays\n");
// varnames
		sprintf(groupname,"%05i/3d",dirtimes[0]);
		g_id = H5Gopen(f_id,groupname,H5P_DEFAULT);
		H5Gget_info(g_id,&group_info);
		nvars = group_info.nlinks;
		for (i = 0; i < nvars; i++)
		{
		    H5Lget_name_by_idx(g_id,".",H5_INDEX_NAME,H5_ITER_INC,i,varname[i],40,H5P_DEFAULT); //40 characters per varname
		}
		H5Gclose(g_id);
		H5Fclose(f_id);
		printf("Succesfully got variable names\n");

		for (i = 0; i < nvars; i++) printf("%s ",varname[i]);printf("\n");

		// OK now we have everything, time to make hdf file

		if ((f_id = H5Fcreate (visitfile, H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT)) < 0)
		{
			fprintf(stderr,"Cannot create visit metadata hdf5 file %s\n",visitfile);
			ERROR_STOP("Cannot create hdf file");
		}

		// Try the lite interface

//		rank=1;dims[0]=1; strtype=H5Tcopy(H5T_C_S1); H5Tset_size(strtype,strlen(topdir));
//		H5LTmake_dataset     (f_id, "/topdir",    rank, dims, strtype, topdir);
		H5LTmake_dataset_string(f_id,"/topdir",topdir);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/ntimedirs", rank, dims, &ntimedirs);
		strtype=H5Tcopy(H5T_C_S1); H5Tset_size(strtype,H5T_VARIABLE); rank=1;dims[0]=ntimedirs;
		H5LTmake_dataset (f_id, "/timedir", rank, dims, strtype, timedir);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/nnodedirs", rank, dims, &nnodedirs);
		strtype=H5Tcopy(H5T_C_S1); H5Tset_size(strtype,H5T_VARIABLE); rank=1;dims[0]=nnodedirs;
		H5LTmake_dataset (f_id, "/nodedir", rank, dims, strtype, nodedir);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/dn", rank, dims, &dn);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/ntottimes", rank, dims, &ntottimes);
		rank=1;dims[0]=ntottimes; H5LTmake_dataset_int (f_id, "/alltimes", rank, dims, alltimes);
		rank=1;dims[0]=ntimedirs; H5LTmake_dataset_int (f_id, "/dirtimes", rank, dims, dirtimes);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/nx", rank, dims, &nx);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/ny", rank, dims, &ny);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/nz", rank, dims, &nz);
        
        //subdomains
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/snx0", rank, dims, &snx0);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/snx1", rank, dims, &snx1);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/sny0", rank, dims, &sny0);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/sny1", rank, dims, &sny1);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/snz0", rank, dims, &snz0);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/snz1", rank, dims, &snz1);

		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/nodex", rank, dims, &nodex);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/nodey", rank, dims, &nodey);
		rank=1;dims[0]=nx; H5LTmake_dataset_float (f_id, "/xhfull", rank, dims, xhfull);
		rank=1;dims[0]=ny; H5LTmake_dataset_float (f_id, "/yhfull", rank, dims, yhfull);
		rank=1;dims[0]=nz; H5LTmake_dataset_float (f_id, "/zh", rank, dims, zh);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/nvars", rank, dims, &nvars);
		//char varname[MAXVARS][40]; Making these fixed makes H5Lget_name_by_idx above easier
		strtype=H5Tcopy(H5T_C_S1); H5Tset_size(strtype,40); rank=1;dims[0]=nvars;
		H5LTmake_dataset (f_id, "/varname", rank, dims, strtype, varname[0]);
		if (H5Fclose(f_id) < 0) ERROR_STOP("Cannot close visit metadata hdf5 file");
		printf("Succesfully wrote %s which can now be read by VisIt using the cm1visit plugin.\n",visitfile);
	}

}

void grok_cm1hdf5_file_structure()
{
	int i;

	ntimedirs = get_num_time_dirs(topdir); printf("ntimedirs: %i\n",ntimedirs);
	if (ntimedirs == 0) ERROR_STOP("No cm1 hdf5 files found");

	timedir = (char **)malloc(ntimedirs * sizeof(char *));
	for (i=0; i < ntimedirs; i++) timedir[i] = (char *)(malloc(256 * sizeof(char)));
	dirtimes = (int *)malloc(ntimedirs * sizeof(int));

	get_sorted_time_dirs(topdir,timedir,dirtimes,ntimedirs,base);

//	for (i=0; i<ntimedirs; i++) printf("SORTED: %s\n",timedir[i]);

	nnodedirs =  get_num_node_dirs(topdir,timedir[0]);
	nodedir = (char **)malloc(nnodedirs * sizeof(char *));
	for (i=0; i < nnodedirs; i++) nodedir[i] = (char *)(malloc(7 * sizeof(char)));

	get_sorted_node_dirs(topdir,timedir[0],nodedir,&dn,nnodedirs);

//	for (i=0; i<nnodedirs; i++) printf("SORTED: %s\n",nodedir[i]);
	printf("Spacing between node directories: %i\n",dn);

	alltimes = get_all_available_times(topdir,timedir,ntimedirs,nodedir,&ntottimes);
//	printf("Times: ");
//	for (i=0; i<ntottimes; i++)printf("%i ",alltimes[i]);
//	printf("\n");

}

void hdf2nc(int argc, char *argv[], char *ncbase, int X0, int Y0, int X1, int Y1, int Z0, int Z1, int t0)
{

#include <netcdf.h>

	float *buffer,*ubuffer,*vbuffer,*wbuffer;
	float time_f;

	int i,ix,iy,iz,nvar;
	char varname[MAXVARS][50];
	char ncfilename[512];

	int snx,sny,snz;
	hid_t f_id;

	int status;
	//TRAJ2014
	int ncid;
	int nxh_dimid,nyh_dimid,nzh_dimid;
	int nxf_dimid,nyf_dimid,nzf_dimid,time_dimid,timeid;
	int x0id,y0id,z0id,x1id,y1id,z1id;
	int thsfcid;
	//TRAJ2014
	int xhid,yhid,zhid;
	int xfid,yfid,zfid;
	int x_dimid,y_dimid,z_dimid;
	int varnameid[MAXVARS];
	int dims[3];
	int d2[2];
	size_t start[3],edges[3];
	size_t s2[2],e2[2];
	int ivar,bufsize;
	float *xhfull,*yhfull,*zh,*xffull,*yffull,*zf;
	float *xhout,*yhout,*zhout,*xfout,*yfout,*zfout;
	//TRAJ2014
	
	strcpy(ncbase,argv[2]);
	X0 = atoi(argv[3]);
	Y0 = atoi(argv[4]);
	X1 = atoi(argv[5]);
	Y1 = atoi(argv[6]);
	Z0 = atoi(argv[7]);
	Z1 = atoi(argv[8]);
	t0 = atoi(argv[9]);

	time_f = (float) t0;

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

	s2[0]=0;
	s2[1]=0;
	e2[0]=sny;
	e2[1]=snx;

	bufsize = snx * sny * snz * sizeof(float);
	if ((buffer = (float *) malloc (bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((ubuffer = (float *) malloc (bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((vbuffer = (float *) malloc (bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((wbuffer = (float *) malloc (bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	//printf("snx = %i sny = %i snz = %i Bufsize = %i\n",snx,sny,snz,bufsize);
	if (buffer == NULL) ERROR_STOP("Cannot allocate buffer");
	if ((f_id = H5Fopen (firstfilename, H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
	{
		fprintf(stderr,"Cannot open firstfilename which is %s, even though we have alredy opened it!\n",firstfilename);
		ERROR_STOP("Cannot open hdf file");
	}
	xhfull = (float *)malloc(nx * sizeof(float));
	yhfull = (float *)malloc(ny * sizeof(float));
	zh = (float *)malloc(nz * sizeof(float));
	zf = (float *)malloc(nz * sizeof(float));
	get1dfloat (f_id,(char *)"mesh/xhfull",xhfull,0,nx);
	get1dfloat (f_id,(char *)"mesh/yhfull",yhfull,0,ny);
	get1dfloat (f_id,(char *)"mesh/zh",zh,0,nz);
	//TRAJ2014
	get1dfloat (f_id,(char *)"mesh/zf",zf,0,nz);

	printf("Succesfully got mesh arrays\n");
	xhout = (float *)malloc(snx * sizeof(float));
	yhout = (float *)malloc(sny * sizeof(float));
	zhout = (float *)malloc(snz * sizeof(float));
	//TRAJ2014
	xfout = (float *)malloc(snx * sizeof(float));
	yfout = (float *)malloc(sny * sizeof(float));
	zfout = (float *)malloc(snz * sizeof(float));

	for (iz=Z0; iz<=Z1; iz++) zhout[iz-Z0] = 0.001*zh[iz];
	for (iy=Y0; iy<=Y1; iy++) yhout[iy-Y0] = 0.001*yhfull[iy];
	for (ix=X0; ix<=X1; ix++) xhout[ix-X0] = 0.001*xhfull[ix];

	//ORF KLUDGE - WE NEED TO SAVE XFFULL YFFULL for u and v!!
	//TRAJ2014
	for (iz=Z0; iz<=Z1; iz++) zfout[iz-Z0] = 0.001*zf[iz];
	for (iy=Y0; iy<=Y1; iy++) yfout[iy-Y0] = 0.001*(yhfull[iy]-15.00);
	for (ix=X0; ix<=X1; ix++) xfout[ix-X0] = 0.001*(xhfull[ix]-15.00);


// ORF 1/25/13 put coordinate dimension in x y and z finally after all these years

	// netcdf v3
	status = nc_create (ncfilename, NC_CLOBBER|NC_64BIT_OFFSET, &ncid); if (status != NC_NOERR) ERROR_STOP ("nc_create failed");
	// netcdf v3
//	status = nc_create (ncfilename, NC_CLOBBER|NC_NETCDF4|NC_CLASSIC_MODEL, &ncid); if (status != NC_NOERR) ERROR_STOP ("nc_create failed");
	//TRAJ2014
	status = nc_def_dim (ncid, "xh", snx, &nxh_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed"); 
	status = nc_def_dim (ncid, "yh", sny, &nyh_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (ncid, "zh", snz, &nzh_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (ncid, "xf", snx, &nxf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed"); 
	status = nc_def_dim (ncid, "yf", sny, &nyf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (ncid, "zf", snz, &nzf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (ncid, "time", 1, &time_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
// Define Cartesian arrays
	//TRAJ2014
	status = nc_def_var (ncid, "xh", NC_FLOAT, 1, &nxh_dimid, &xhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "yh", NC_FLOAT, 1, &nyh_dimid, &yhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "zh", NC_FLOAT, 1, &nzh_dimid, &zhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");

	status = nc_def_var (ncid, "xf", NC_FLOAT, 1, &nxf_dimid, &xfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "yf", NC_FLOAT, 1, &nyf_dimid, &yfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "zf", NC_FLOAT, 1, &nzf_dimid, &zfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "time", NC_FLOAT, 1, &time_dimid, &timeid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");


	//TRAJ2014
	status = nc_put_att_text(ncid, xhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, yhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, zhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");

	status = nc_put_att_text(ncid, xfid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, yfid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, zfid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, timeid, "units", strlen("s"), "s");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");

	status = nc_def_var (ncid, "X0", NC_INT, 0, dims, &x0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "Y0", NC_INT, 0, dims, &y0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "X1", NC_INT, 0, dims, &x1id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "Y1", NC_INT, 0, dims, &y1id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "Z0", NC_INT, 0, dims, &z0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "Z1", NC_INT, 0, dims, &z1id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");

//	status = nc_def_var (ncid, "time", NC_FLOAT, 0, dims, &timeid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");



	for (ivar = 0; ivar < nvar; ivar++)
	{
		//TRAJ2014
		if(!strcmp(varname[ivar],"u"))
		{
			dims[0] = nzh_dimid;
			dims[1] = nyh_dimid;
			dims[2] = nxf_dimid;
		}
		else if (!strcmp(varname[ivar],"v"))
		{
			dims[0] = nzh_dimid;
			dims[1] = nyf_dimid;
			dims[2] = nxh_dimid;
		}
		else if (!strcmp(varname[ivar],"w"))
		{
			dims[0] = nzf_dimid;
			dims[1] = nyh_dimid;
			dims[2] = nxh_dimid;
		}
		else
		{
			dims[0] = nzh_dimid;
			dims[1] = nyh_dimid;
			dims[2] = nxh_dimid;
		}


		d2[0] = nyh_dimid;
		d2[1] = nxh_dimid;
		status = nc_def_var (ncid, "thpertsfc", NC_FLOAT, 2, d2, &thsfcid);
		status = nc_def_var (ncid, varname[ivar], NC_FLOAT, 3, dims, &(varnameid[ivar]));
		if (status != NC_NOERR) 
		{
			printf ("Cannot nc_def_var for var #%i %s\n", ivar, varname[ivar]);
			ERROR_STOP("nc_def_var failed");
		}
// time consuming and not really saving much space for interesting data
//                 status=nc_def_var_deflate(ncid, varnameid[ivar], 0, 1, 9);
//	            if (status != NC_NOERR) printf ("Cannot nc_def_var_deflate for var #%i %s\n", ivar, varname[ivar]);
	}
	status = nc_enddef (ncid);
	if (status != NC_NOERR) ERROR_STOP("nc_enddef failed");

	// For surface theta (Vapor wants 2D vars)
	read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"thpert",X0,Y0,X1,Y1,0,0,nx,ny,nz,nodex,nodey);
	status = nc_put_vara_float (ncid, thsfcid, s2, e2, buffer);

      status = nc_put_var_float (ncid,xhid,xhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_float (ncid,yhid,yhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_float (ncid,zhid,zhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	//TRAJ2014
      status = nc_put_var_float (ncid,xfid,xfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_float (ncid,yfid,yfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_float (ncid,zfid,zfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");

      status = nc_put_var_int (ncid,x0id,&X0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,y0id,&Y0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,x1id,&X1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,y1id,&Y1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,z0id,&Z0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,z1id,&Z1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_float (ncid,timeid,&time_f); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");

	for (ivar = 0; ivar < nvar; ivar++)
	{
		printf("Working on %s (",varname[ivar]);

		if(!strcmp(varname[ivar],"hwin"))
		{
			read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"uinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"vinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			for(i=0; i<snx*sny*snz; i++)
			{
				buffer[i] = sqrt(ubuffer[i]*ubuffer[i]+vbuffer[i]*vbuffer[i]);
			}
		}
		else if(!strcmp(varname[ivar],"hvort"))
		{
			read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"xvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"yvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			for(i=0; i<snx*sny*snz; i++)
			{
				buffer[i] = sqrt(ubuffer[i]*ubuffer[i]+vbuffer[i]*vbuffer[i]);
			}
		}
		else if(!strcmp(varname[ivar],"vortmag"))
		{
			read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"xvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"yvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(wbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"zvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			for(i=0; i<snx*sny*snz; i++)
			{
				buffer[i] = sqrt(ubuffer[i]*ubuffer[i]+vbuffer[i]*vbuffer[i] + wbuffer[i]*wbuffer[i]);
			}
		}
		else
		{
			read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,varname[ivar],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
		}
		status = nc_put_vara_float (ncid, varnameid[ivar], start, edges, buffer);
		if (status != NC_NOERR) 
		{
			fprintf(stderr,"Could not write variable %s at time %i to %s\n", varname[ivar],t0,ncfilename);
			ERROR_STOP("Write to netcdf file failed");
		}
		fprintf(stderr,")\n");
	}
	status = nc_close(ncid); if (status != NC_NOERR) ERROR_STOP("Could not close netcdf file (nc_close)");
}

void hdf2v5d(int argc, char *argv[])
{
	FILE *fp_data;
	char datafile[100];
	char v5dfilename[512],cpcmd[512];
	int xoff,yoff,X0,Y0,X1,Y1,Z0,Z1,snx,sny,snz;
	int nx,ny,ih;
	int dt,t0,t1;
	int compression;
	int sbox;
	float umove,vmove;
	int i, nblocks;
	int *timestamp,*datestamp;
	int ivar,nv5dvars,ib,iz;
	float *buffer,*ubuffer,*vbuffer,*wbuffer,*xh,*yh,*zh,*xyslice;
	float dx,dy;
	char v5dvarname[MAXVARS][10];
	char varname[MAXVARS][10];

	float dxkm,dykm;

	int projection = 1; // 1 = linear, rectangular, cylindrical equidistant (must represent dx, dy in degrees)
	int vertical = 2;   // 2 = Stretched vertical grid, explicitly indicate height AGL in array
	float proj_args[4]; // Projection arguments for our grid 
	int *nla;

	hid_t file_id,sd_id;


	if (argc == 3)	strcpy(datafile,argv[2]);
	else 			strcpy(datafile,"hdf2v5d.dat");


	if ((fp_data = fopen (datafile, "r")) == NULL)
	{
		char errorbuf[500];
		sprintf(errorbuf,"Could not open %s\n",datafile);
		ERROR_STOP(errorbuf);
	}

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
	sprintf (cpcmd, "cp %s hdf2v5d.dat.%s",datafile,v5dfilename); system (cpcmd);

	snx = X1 - X0 + 1;
	sny = Y1 - Y0 + 1;
	snz = Z1 - Z0 + 1;


	if ((zh = (float *) malloc (snz * sizeof (float))) == NULL) ERROR_STOP("Cannot allocate zh");
	if ((buffer = (float *) malloc (snx * sny * snz * sizeof (float))) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((ubuffer = (float *) malloc (snx * sny * snz * sizeof (float))) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((vbuffer = (float *) malloc (snx * sny * snz * sizeof (float))) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((wbuffer = (float *) malloc (snx * sny * snz * sizeof (float))) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((xyslice = (float *) malloc (snx * sny * sizeof (float))) == NULL) ERROR_STOP("Cannot allocate 2D buffer slice");

	if ((file_id = H5Fopen (firstfilename, H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
	{
		fprintf(stderr,"Cannot open firstfilename which is %s, even though we have alredy opened it!\n",firstfilename);
		ERROR_STOP("Cannot open hdf file");
	}

	get1dfloat (file_id,"mesh/zh",zh,Z0,snz);
//	get0dfloat (file_id, "mesh/dx", &dx);
//	get0dfloat (file_id, "mesh/dy", &dy);

// ARGH - if we use stretched horizontal grid this is WRONG-O!
// Should read in xh and yh mesh and take difference in center of domain
// This will ensure what is desired
//	get0dfloat (file_id, "mesh/dx", &dx);
//	get0dfloat (file_id, "mesh/dy", &dy);
	get0dint (file_id, "grid/nx", &nx);
	get0dint (file_id, "grid/ny", &ny);
	if ((xh = (float *) malloc (nx * sizeof (float))) == NULL) ERROR_STOP("Cannot allocate xh");
	if ((yh = (float *) malloc (ny * sizeof (float))) == NULL) ERROR_STOP("Cannot allocate yh");
	get1dfloat (file_id,"mesh/xhfull",xh,0,nx-1);
	get1dfloat (file_id,"mesh/yhfull",yh,0,ny-1);
	ih=nx/2; dx = xh[ih+1]-xh[ih];
	ih=ny/2; dy = yh[ih+1]-yh[ih];

	printf("Center of domain: dx = %f\t dy = %f\n",dx,dy);


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
			printf ("%s(", varname[ivar - 1]); fflush (stdout);
			if(!strcmp(varname[ivar-1],"hwin"))
			{
				read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,itime,"uinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
				read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,itime,"vinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
				for(i=0; i<snx*sny*snz; i++)
				{
					buffer[i] = sqrt(ubuffer[i]*ubuffer[i]+vbuffer[i]*vbuffer[i]);
				}
			}
			else if(!strcmp(varname[ivar-1],"vortmag"))
			{
				read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,itime,"xvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
				read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,itime,"yvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
				read_hdf_mult_md(wbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,itime,"zvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
				for(i=0; i<snx*sny*snz; i++)
				{
					buffer[i] = sqrt(ubuffer[i]*ubuffer[i]+vbuffer[i]*vbuffer[i]+wbuffer[i]*wbuffer[i]);
				}
			}
			else if(!strcmp(varname[ivar-1],"hvort"))
			{
				read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,itime,"xvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
				read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,itime,"yvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
				for(i=0; i<snx*sny*snz; i++)
				{
					buffer[i] = sqrt(ubuffer[i]*ubuffer[i]+vbuffer[i]*vbuffer[i]);
				}
			}
			else
			{
				read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,itime,varname[ivar-1],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			}
			printf (")"); fflush (stdout);
			fix_fubar_v5d_coord (buffer, xyslice,snx,sny,snz);
			v5dWrite (ib, ivar, buffer);
		}
		printf ("\n");
	}
	v5dClose ();
	free(buffer); free(ubuffer); free(vbuffer); free(wbuffer);
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
//fubar used to be another word, I'll let you guess what it was.
//Historical fact: Vis5d was written at SSEC initially to be compatible
//with satellite data - which scans top down, so it adoped a coordinate
//system where top left was (0,0), not bottom left like the flying
//spaghetti monster meant it to be.

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
