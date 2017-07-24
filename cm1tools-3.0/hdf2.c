/* Leigh Orf cm1tools to convert buffered file-per-core cm1 hdf5 output
 * to other more useful formats
 *
 * Big rewrite June 2011. All top-level conversion code in one file, makes
 * maintenance easier. Symlinks to binary determine what is run.
 *
 * Major cleanup 11/23/2016, cleared out all unused variables and got
 * gcc to mostly shut up. Getting this code stable before converting to
 * new cm1hdf5 (subsecond) format.
 *
 * Spring '17: Switched over to floating point format for times in LOFS
 *
 * Adding command line parsing 7/17
 */

#include <math.h>
#include "hdforf.h"
#include <netcdf.h>
#include "errorf.h"
#include <getopt.h>
void parse_cmdline_hdf2nc(int argc, char *argv[],
		char *histpath, int *got_histpath, char *ncbase, int *got_ncbase,
		double *time, int *got_time,
		int *X0, int *got_X0,
		int *Y0, int *got_Y0,
		int *X1, int *got_X1,
		int *Y1, int *got_Y1,
		int *Z0, int *got_Z0,
		int *Z1, int *got_Z1 );
#define TRUE 1
#define FALSE 0
#define MAXVARIABLES (100)
#define MAXSTR (512)

#define P3(x,y,z,mx,my) (((z)*(mx)*(my))+((y)*(mx))+(x))

/* These vars are global, everyone uses them. Or I am lazy. */
char topdir[PATH_MAX+1];
int dn;
char **timedir; 
char **nodedir;
double *dirtimes;
int ntimedirs;
int nx,ny,nz,nodex,nodey;
char firstfilename[MAXSTR];
char base[MAXSTR];
int nnodedirs;
double *alltimes;
int ntottimes;
int firsttimedirindex;

void grok_cm1hdf5_file_structure();
void hdf2nc(int argc, char *argv[], char *ncbase, int X0, int Y0, int X1, int Y1, int Z0, int Z1, double t0);
int hhmmss (int time);
void hdf2v5d();
void fix_fubar_v5d_coord (float *c, float *s,int snx, int sny, int snz);
void lflush (FILE * fp);

int main(int argc, char *argv[])
{
	char *cptr;
	char progname[MAXSTR];
	char histpath[MAXSTR];
	int we_are_hdf2nc = FALSE;
	int we_are_hdf2v5d = FALSE;
	int we_are_makevisit = FALSE;
	int we_are_linkfiles = FALSE;

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
		int X0,Y0,X1,Y1,Z0,Z1;
		double time;
		int got_histpath,got_ncbase,got_time;
		int got_X0,got_Y0,got_X1,got_Y1,got_Z0,got_Z1;

		got_X0=got_Y0=got_X1=got_Y1=got_Z0=got_Z1=0;
		X0=Y0=X1=Y1=Z0=Z1=0;
		time=0.0;
		char ncbase[MAXSTR];

		parse_cmdline_hdf2nc(argc, argv, histpath, &got_histpath, ncbase, &got_ncbase,
			&time, &got_time, &X0, &got_X0, &Y0, &got_Y0, &X1, &got_X1, &Y1, &got_Y1,
			&Z0, &got_Z0, &Z1, &got_Z1);

/* I want the full, absolute path to the top level directory, not the relative path. 
 * This is so I can find the files later on if I forget where I put
 * them, or at least to remind me where they were */

		if((cptr=realpath(histpath,topdir))==NULL)ERROR_STOP("realpath failed");
		grok_cm1hdf5_file_structure();
		get_hdf_metadata(firstfilename,&nx,&ny,&nz,&nodex,&nodey);
		//printf("ORF: nx = %i ny = %i nz = %i nodex = %i nodey = %i\n", nx,ny,nz,nodex,nodey);
		hdf2nc(argc,argv,ncbase,X0,Y0,X1,Y1,Z0,Z1,time);
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


	if (we_are_hdf2nc)
	{
		/*
		int X0,Y0,X1,Y1,Z0,Z1;
		double t0;
		X0=Y0=X1=Y1=Z0=Z1=0;
		t0=0.0;
		char ncbase[512];
		hdf2nc(argc,argv,ncbase,X0,Y0,X1,Y1,Z0,Z1,t0);
		*/
		exit(0); /* ORF we moved it all up */
	}
	else if (we_are_hdf2v5d)
	{
		hdf2v5d(argc,argv);
	}
	else if (we_are_linkfiles)
	{
	  int snx0, snx1, sny0, sny1;
	  int idir,inode;
	  char dirname[MAXSTR];
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
		hid_t f_id,g_id,strtype;
		H5G_info_t group_info;
		hsize_t dims[1];
		char dirbase[MAXSTR];
		char visitfile[MAXSTR];
		char groupname[MAXSTR];
		char varname[MAXVARIABLES][40]; //Yeah I know also hardcoded in avtcm1visitFileFormat.h

		float *xhfull,*yhfull,*zh;

		for (i=0; i < strlen(topdir)-2; i++) dirbase[i]=topdir[i];
		dirbase[i]='\0';

		sprintf(visitfile,"%s%s.cm1visit",dirbase,base);
		printf("visitfile = %s\n",visitfile);
		
		printf("topdir = %s\n",topdir);
		printf("ntimedirs = %i\n",ntimedirs);
		for (i=0; i<ntimedirs; i++)printf("%s ",timedir[i]);printf("\n");
		for (i=0; i<ntimedirs; i++)printf("%f ",dirtimes[i]);printf("\n");
		printf("nnodedirs = %i\n",nnodedirs);
		for (i=0; i<nnodedirs; i++)printf("%s ",nodedir[i]);printf("\n");
		printf("dn = %i\n",dn);
		printf("ntottimes = %i\n",ntottimes);
		for (i=0; i<ntottimes; i++)printf("%f ",alltimes[i]);printf("\n");
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
// varnames
//		sprintf(groupname,"%05i/3D",dirtimes[ntimedirs-1]); //is now last since we no longer have get_first_hdf_file_name routine;
// ORF: there will always be a group zero now that each file begins at
// zero
		sprintf(groupname,"%05i/3D",0); //is now last since we no longer have get_first_hdf_file_name routine;
		// first file name is the last "fist file name" from
		// sweeping through all times, so now "first file name" is
		// from last time direcotry 
		// This is getting convoluted now that I have to deal with
		// empty node dirs from saving subdomains.
		printf("MAKEVISIT: firstfilename = %s\t groupname = %s\n",firstfilename,groupname);
		g_id = H5Gopen(f_id,groupname,H5P_DEFAULT);
		H5Gget_info(g_id,&group_info);
		nvars = group_info.nlinks;
		for (i = 0; i < nvars; i++)
		{
		    H5Lget_name_by_idx(g_id,".",H5_INDEX_NAME,H5_ITER_INC,i,varname[i],40,H5P_DEFAULT); //40 characters per varname
		}
		H5Gclose(g_id);
		H5Fclose(f_id);

		for (i = 0; i < nvars; i++) printf("%s ",varname[i]);printf("\n");

		// OK now we have everything, time to make hdf file

		if ((f_id = H5Fcreate (visitfile, H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT)) < 0)
		{
			fprintf(stderr,"Cannot create visit metadata hdf5 file %s\n",visitfile);
			ERROR_STOP("Cannot create hdf file");
		}

		// Try the lite interface
		H5LTmake_dataset_string(f_id,"/topdir",topdir);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/ntimedirs", rank, dims, &ntimedirs);
		strtype=H5Tcopy(H5T_C_S1); H5Tset_size(strtype,H5T_VARIABLE); rank=1;dims[0]=ntimedirs;
		H5LTmake_dataset (f_id, "/timedir", rank, dims, strtype, timedir);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/nnodedirs", rank, dims, &nnodedirs);
		strtype=H5Tcopy(H5T_C_S1); H5Tset_size(strtype,H5T_VARIABLE); rank=1;dims[0]=nnodedirs;
		H5LTmake_dataset (f_id, "/nodedir", rank, dims, strtype, nodedir);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/dn", rank, dims, &dn);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/ntottimes", rank, dims, &ntottimes);
// ORF: 2017-01-25 old format was int, now we are float
//		rank=1;dims[0]=ntottimes; H5LTmake_dataset_int (f_id, "/alltimes", rank, dims, alltimes);
//		rank=1;dims[0]=ntimedirs; H5LTmake_dataset_int (f_id, "/dirtimes", rank, dims, dirtimes);
		rank=1;dims[0]=ntottimes; H5LTmake_dataset_double (f_id, "/alltimes", rank, dims, alltimes);
		rank=1;dims[0]=ntimedirs; H5LTmake_dataset_double (f_id, "/dirtimes", rank, dims, dirtimes);
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
		//char varname[MAXVARIABLES][40]; Making these fixed makes H5Lget_name_by_idx above easier
		strtype=H5Tcopy(H5T_C_S1); H5Tset_size(strtype,40); rank=1;dims[0]=nvars;
		H5LTmake_dataset (f_id, "/varname", rank, dims, strtype, varname[0]);
		if (H5Fclose(f_id) < 0) ERROR_STOP("Cannot close visit metadata hdf5 file");
		printf("Succesfully wrote %s which can now be read by VisIt using the cm1visit plugin.\n",visitfile);
	}

	exit(0);
}

void grok_cm1hdf5_file_structure()
{
	int i;
	ntimedirs = get_num_time_dirs(topdir); printf("ntimedirs: %i\n",ntimedirs);
	if (ntimedirs == 0) ERROR_STOP("No cm1 hdf5 files found");

	timedir = (char **)malloc(ntimedirs * sizeof(char *));
	for (i=0; i < ntimedirs; i++) timedir[i] = (char *)(malloc(MAXSTR * sizeof(char)));
	dirtimes = (double *)malloc(ntimedirs * sizeof(double)); // ORF NO LONGER INT, RIGHT???

	get_sorted_time_dirs(topdir,timedir,dirtimes,ntimedirs,base);

	nnodedirs =  get_num_node_dirs(topdir,timedir[0]);
	nodedir = (char **)malloc(nnodedirs * sizeof(char *));
	for (i=0; i < nnodedirs; i++) nodedir[i] = (char *)(malloc(8 * sizeof(char)));

	get_sorted_node_dirs(topdir,timedir[0],nodedir,&dn,nnodedirs);

	alltimes = get_all_available_times(topdir,timedir,ntimedirs,nodedir,nnodedirs,&ntottimes,firstfilename,&firsttimedirindex);
#ifdef DEBUG
	printf("Times: ");
	for (i=0; i<ntottimes; i++)printf("%lf ",alltimes[i]);
	printf("\n");
#endif
}

void hdf2nc(int argc, char *argv[], char *ncbase, int X0, int Y0, int X1, int Y1, int Z0, int Z1, double t0)
{
	float *buffer,*buf0,*ubuffer,*vbuffer,*wbuffer,*xvort,*yvort,*zvort;
	float *qvar1,*qvar2,*qvar3;
	float *thpert,*qvpert,*qtot;
	float reps;
	float *th0,*qv0;
	float *dum1,*dum2;
	double timearray[1];

	int i,ix,iy,iz,nvar;
	int j=0,k=0,ii;
	char varname[MAXVARIABLES][MAXSTR];
	char ncfilename[MAXSTR];

	int snx,sny,snz;
	hid_t f_id;

	int status;
	int ncid;
	int nxh_dimid,nyh_dimid,nzh_dimid;
	int nxf_dimid,nyf_dimid,nzf_dimid,time_dimid,timeid;
	int thsfcid,dbzsfcid;
	int x0id,y0id,z0id,x1id,y1id,z1id;
	int xhid,yhid,zhid;
	int xfid,yfid,zfid;
	int varnameid[MAXVARIABLES];
	int dims[4];
	int d2[3];
	size_t start[4],edges[4];
	size_t s2[3],e2[3];
	int ivar;
	long int bufsize;
	float *xhfull,*yhfull,*zh,*zf;
	float *xhout,*yhout,*zhout,*xfout,*yfout,*zfout;
	FILE *fp;
	char cmdfilename[512];
	//ORF for writing single time in unlimited time dimension/variable
	const size_t timestart = 0;
	const size_t timecount = 1;
	int n_cmdline_args=19;

	/* We used to set arguments based upon what order they were passed to
	 * hdf2nc. Now we use getopt with descriptive names that can be
	 * passed in any order, except the variable names all come last */
	
	/* This is kind of ugly but it works. All command line arguments
	 * to hdf2nc are mandatory, and we have 18 items before the list
	 * of variable names that brings up the rear of all the command
	 * line arguments (add one for the executable name!) */

	nvar = argc-n_cmdline_args;

	for (i=0; i<nvar; i++)
	{
		strcpy(varname[i],argv[i+n_cmdline_args]);
		fprintf(stdout,"%s \n",varname[i]);
	}
	sprintf(ncfilename,"%s.%012.6f.nc",ncbase,t0);
	
	snx = X1 - X0 + 1;
	sny = Y1 - Y0 + 1;
	snz = Z1 - Z0 + 1;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	start[3] = 0;
	edges[0] = 1;
	edges[1] = snz;
	edges[2] = sny;
	edges[3] = snx;
	//For 2D surface slices
	s2[0] = 0;
	s2[1] = 0;
	s2[2] = 0;
	e2[0] = 1;
	e2[1] = sny;
	e2[2] = snx;

	/* ORF LAZY allocate enough for any staggered combination, hence +1 for all three dimensions */
	bufsize = (long) (snx+1) * (long) (sny+1) * (long) (snz+1) * (long) sizeof(float);
	fprintf(stdout,"X0=%i Y0=%i X1=%i Y1=%i Z0=%i Z1=%i bufsize = %f GB\n",X0,Y0,X1,Y1,Z0,Z1,1.0e-9*bufsize);
	if ((buf0 = buffer = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((ubuffer = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((vbuffer = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((wbuffer = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((xvort = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((yvort = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((zvort = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((dum1 = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((dum2 = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	//printf("snx = %i sny = %i snz = %i Bufsize = %i\n",snx,sny,snz,bufsize);
	if (buffer == NULL) ERROR_STOP("Cannot allocate buffer");
	if ((f_id = H5Fopen (firstfilename, H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
	{
		fprintf(stderr,"Cannot open firstfilename which is %s, even though we have alredy opened it!\n",firstfilename);
		ERROR_STOP("Cannot open hdf file");
	}
	xhfull = (float *)malloc(nx * sizeof(float));
	yhfull = (float *)malloc(ny * sizeof(float));
	th0 = (float *)malloc(nz * sizeof(float));
	qv0 = (float *)malloc(nz * sizeof(float));
	zh = (float *)malloc(nz * sizeof(float));
	zf = (float *)malloc(nz * sizeof(float));
	get1dfloat (f_id,(char *)"mesh/xhfull",xhfull,0,nx);
	get1dfloat (f_id,(char *)"mesh/yhfull",yhfull,0,ny);
	/* ORF can't I just read xffull and yffull here?? */
	get1dfloat (f_id,(char *)"mesh/zh",zh,0,nz);
	get1dfloat (f_id,(char *)"mesh/zf",zf,0,nz);
	get1dfloat (f_id,(char *)"basestate/qv0",qv0,0,nz);
	get1dfloat (f_id,(char *)"basestate/th0",th0,0,nz);
	get1dfloat (f_id,(char *)"mesh/zf",zf,0,nz);

	xhout = (float *)malloc(snx * sizeof(float));
	yhout = (float *)malloc(sny * sizeof(float));
	zhout = (float *)malloc(snz * sizeof(float));
	xfout = (float *)malloc(snx * sizeof(float));
	yfout = (float *)malloc(sny * sizeof(float));
	zfout = (float *)malloc(snz * sizeof(float));

	for (iz=Z0; iz<=Z1; iz++) zhout[iz-Z0] = 0.001*zh[iz];
	for (iy=Y0; iy<=Y1; iy++) yhout[iy-Y0] = 0.001*yhfull[iy];
	for (ix=X0; ix<=X1; ix++) xhout[ix-X0] = 0.001*xhfull[ix];
	for (iz=Z0; iz<=Z1; iz++) zfout[iz-Z0] = 0.001*zf[iz];

	// Do we save these in the HDF5 files? We should just be able to read them directly.
	for (iy=Y0; iy<=Y1; iy++) yfout[iy-Y0] = 0.001*(yhfull[iy]-15.00); //ORF STUPID FIX BUG BUG BUG
	for (ix=X0; ix<=X1; ix++) xfout[ix-X0] = 0.001*(xhfull[ix]-15.00); //ORF STUPID FIX BUG BUG BUG

	H5Z_zfp_initialize();

	// ORF could make option: What kind of netCDF file to create?
	// If we add optional options to hdf2nc, we'll have to increment
	// n_cmd_line_args so we get all the varnames
	status = nc_create (ncfilename, NC_CLOBBER|NC_64BIT_OFFSET, &ncid); if (status != NC_NOERR) ERROR_STOP ("nc_create failed");

//	status = nc_create (ncfilename, NC_CLOBBER|NC_NETCDF4, &ncid); if (status != NC_NOERR)
//	{
//		printf("Error %i\n",status);
//		ERROR_STOP ("nc_create failed");
//	}
//	status = nc_create (ncfilename, NC_CLOBBER|NC_NETCDF4|NC_CLASSIC_MODEL, &ncid); if (status != NC_NOERR) ERROR_STOP ("nc_create failed");
//
//
	status = nc_def_dim (ncid, "xh", snx, &nxh_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed"); 
	status = nc_def_dim (ncid, "yh", sny, &nyh_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (ncid, "zh", snz, &nzh_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (ncid, "xf", snx, &nxf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed"); 
	status = nc_def_dim (ncid, "yf", sny, &nyf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (ncid, "zf", snz, &nzf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (ncid, "time", NC_UNLIMITED, &time_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_var (ncid, "xh", NC_FLOAT, 1, &nxh_dimid, &xhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "yh", NC_FLOAT, 1, &nyh_dimid, &yhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "zh", NC_FLOAT, 1, &nzh_dimid, &zhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "xf", NC_FLOAT, 1, &nxf_dimid, &xfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "yf", NC_FLOAT, 1, &nyf_dimid, &yfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "zf", NC_FLOAT, 1, &nzf_dimid, &zfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "time", NC_DOUBLE, 1, &time_dimid, &timeid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");

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

	d2[0] = time_dimid;
	d2[1] = nyh_dimid;
	d2[2] = nxh_dimid;
	status = nc_def_var (ncid, "thrhopert_sfc", NC_FLOAT, 3, d2, &thsfcid);
	status = nc_def_var (ncid, "dbz_sfc", NC_FLOAT, 3, d2, &dbzsfcid);

	for (ivar = 0; ivar < nvar; ivar++)
	{
		if(!strcmp(varname[ivar],"u"))
		{
			dims[0] = time_dimid;
			dims[1] = nzh_dimid;
			dims[2] = nyh_dimid;
			dims[3] = nxf_dimid;
		}
		else if (!strcmp(varname[ivar],"v"))
		{
			dims[0] = time_dimid;
			dims[1] = nzh_dimid;
			dims[2] = nyf_dimid;
			dims[3] = nxh_dimid;
		}
		else if (!strcmp(varname[ivar],"w"))
		{
			dims[0] = time_dimid;
			dims[1] = nzf_dimid;
			dims[2] = nyh_dimid;
			dims[3] = nxh_dimid;
		}
		else
		{
			dims[0] = time_dimid;
			dims[1] = nzh_dimid;
			dims[2] = nyh_dimid;
			dims[3] = nxh_dimid;
		}

		status = nc_def_var (ncid, varname[ivar], NC_FLOAT, 4, dims, &(varnameid[ivar]));
		if (status != NC_NOERR) 
		{
			printf ("Cannot nc_def_var for var #%i %s, status = %i, message = %s\n", ivar, varname[ivar],status,nc_strerror(status));
			ERROR_STOP("nc_def_var failed");
		}
// time consuming and not really saving much space for interesting data
//                 status=nc_def_var_deflate(ncid, varnameid[ivar], 0, 1, 9);
//	            if (status != NC_NOERR) printf ("Cannot nc_def_var_deflate for var #%i %s\n", ivar, varname[ivar]);
	}
	status = nc_enddef (ncid);

//	this is our 4d write:	status = nc_put_vara_float (ncid, varnameid[ivar], start, edges, buffer);

	if (status != NC_NOERR) ERROR_STOP("nc_enddef failed");

	// For surface theta (Vapor wants 2D vars)
	//
	// Save surface thrhopert and dbz for easy viewing as 2D vars in
	// Vapor
	// Should make this a command line option
	read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"thrhopert",X0,Y0,X1,Y1,0,0,nx,ny,nz,nodex,nodey);
	status = nc_put_vara_float (ncid, thsfcid, s2, e2, buffer);
	read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"dbz",X0,Y0,X1,Y1,0,0,nx,ny,nz,nodex,nodey);
	status = nc_put_vara_float (ncid, dbzsfcid, s2, e2, buffer);
	// Temporary: Read u, v, w, xvort,yvort,zvort up front
	// ORF 2016-03-07 ugh, we stored u v w for 20 meter runs
	/*
	*/
#ifdef READAHEAD

#ifndef UVW
	read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"uinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
	read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"vinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
	read_hdf_mult_md(wbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"winterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
#else
	read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"u",X0,Y0,X1+1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
	for(k=0; k<snz; k++)
	{
		for(j=0; j<sny; j++)
		{
			for(i=0; i<snx; i++)
			{
				ubuffer[P3(i,j,k,snx,sny)] = 0.5*(buffer[P3(i,j,k,snx+1,sny)] + buffer[P3(i+1,j,k,snx+1,sny)]);
			}
		}
	}
	read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"v",X0,Y0,X1,Y1+1,Z0,Z1,nx,ny,nz,nodex,nodey);
	for(k=0; k<snz; k++)
	{
		for(j=0; j<sny; j++)
		{
			for(i=0; i<snx; i++)
			{
				vbuffer[P3(i,j,k,snx,sny)] = 0.5*(buffer[P3(i,j,k,snx,sny+1)] + buffer[P3(i,j+1,k,snx,sny+1)]);
			}
		}
	}
	read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"w",X0,Y0,X1,Y1,Z0,Z1+1,nx,ny,nz,nodex,nodey);
	for(k=0; k<snz; k++)
	{
		for(j=0; j<sny; j++) 
		{
			for(i=0; i<snx; i++) 
			{
				wbuffer[P3(i,j,k,snx,sny)] = 0.5*(buffer[P3(i,j,k,snx,sny)] + buffer[P3(i,j,k+1,snx,sny)]);
			}
		}
	}
#endif
	read_hdf_mult_md(  xvort,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,  "xvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
	read_hdf_mult_md(  yvort,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,  "yvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
	read_hdf_mult_md(  zvort,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,  "zvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
#endif
      status = nc_put_var_float (ncid,xhid,xhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_float (ncid,yhid,yhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_float (ncid,zhid,zhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_float (ncid,xfid,xfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_float (ncid,yfid,yfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_float (ncid,zfid,zfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");

      status = nc_put_var_int (ncid,x0id,&X0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,y0id,&Y0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,x1id,&X1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,y1id,&Y1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,z0id,&Z0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,z1id,&Z1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	timearray[0] = t0;
      status = nc_put_vara_double (ncid,timeid,&timestart,&timecount,timearray); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");

	for (ivar = 0; ivar < nvar; ivar++)
	{
		printf("Working on %s (",varname[ivar]);
#ifdef READAHEAD
		if     (!strcmp(varname[ivar],"uinterp")) buffer = ubuffer;
		else if(!strcmp(varname[ivar],"vinterp")) buffer = vbuffer;
		else if(!strcmp(varname[ivar],"winterp")) buffer = wbuffer;
		else if(!strcmp(varname[ivar],"xvort")) buffer = xvort;
		else if(!strcmp(varname[ivar],"yvort")) buffer = yvort;
		else if(!strcmp(varname[ivar],"zvort")) buffer = zvort;
		else
#endif 
		if(!strcmp(varname[ivar],"hwin_sr"))
		{
			float usr,vsr;
#ifndef READAHEAD
			read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"uinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"vinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
#endif
			for(i=0; i<snx*sny*snz; i++)
			{
				usr = ubuffer[i];
				vsr = vbuffer[i];
				buffer[i] = sqrt(usr*usr+vsr*vsr);
			}
		}
		else if(!strcmp(varname[ivar],"hwin_gr"))
		{
			float umove = 15.1719; float vmove = 10.0781; /* ORF: KLUDGE: should save these in history files */
			float ugr, vgr;
#ifndef READAHEAD
			read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"uinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"vinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
#endif
			for(i=0; i<snx*sny*snz; i++)
			{
				ugr = ubuffer[i]+umove;
				vgr = vbuffer[i]+vmove;
				buffer[i] = sqrt(ugr*ugr+vgr*vgr);
			}
		}
		else if(!strcmp(varname[ivar],"hwin_sr_uv"))
		{
			float usr,vsr;
			read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"u",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"v",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			for(k=0; k<snz; k++)
			for(j=0; j<sny; j++)
			for(i=0; i<snx; i++)
			{
				usr = (i!=snx-1)?ubuffer[P3(i+1,j,k,snx,sny)]+ubuffer[P3(i,j,k,snx,sny)]:0.0; //for now just ignore last value
				vsr = (i!=sny-1)?vbuffer[P3(i,j+1,k,snx,sny)]+vbuffer[P3(i,j,k,snx,sny)]:0.0; //for now just ignore last value
				buffer[i] = sqrt(usr*usr+vsr*vsr);
			}
		}
		else if(!strcmp(varname[ivar],"hwin_gr_uv"))
		{
			float umove = 15.1719; float vmove = 10.0781; /* ORF: KLUDGE: should save these in history files */
			float ugr, vgr;
			read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"u",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"v",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			for(i=0; i<snx*sny*snz; i++)
			{
				ugr = (i!=snx-1)?ubuffer[P3(i+1,j,k,snx,sny)]+ubuffer[P3(i,j,k,snx,sny)]:0.0; //for now just ignore last value
				vgr = (i!=sny-1)?vbuffer[P3(i,j+1,k,snx,sny)]+vbuffer[P3(i,j,k,snx,sny)]:0.0; //for now just ignore last value
				ugr += umove;
				vgr += vmove;
				buffer[i] = sqrt(ugr*ugr+vgr*vgr);
			}
		}
		else if(!strcmp(varname[ivar],"hvort"))
		{
#ifndef READAHEAD
			read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"xvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"yvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
#endif
			for(i=0; i<snx*sny*snz; i++)
			{
#ifndef READAHEAD
				buffer[i] = sqrt(ubuffer[i]*ubuffer[i]+vbuffer[i]*vbuffer[i]);
#else
				buffer[i] = sqrt(xvort[i]*xvort[i]+yvort[i]*yvort[i]);
#endif
			}
		}
//#ifndef FUCKTHIS
// ORF because I'm saving a variable called vortmag!
		else if(!strcmp(varname[ivar],"vortmag"))
		{
#ifndef READAHEAD
			read_hdf_mult_md(xvort,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"xvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(yvort,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"yvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(zvort,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"zvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
#endif
			for(i=0; i<snx*sny*snz; i++)
			{
				buffer[i] = sqrt(xvort[i]*xvort[i]+yvort[i]*yvort[i]+zvort[i]*zvort[i]);
			}
		}
//#endif
		else if(!strcmp(varname[ivar],"streamvort")) /* streamwise vorticity */
		{
			for(i=0; i<snx*sny*snz; i++)
			{
				buffer[i] = (ubuffer[i]*xvort[i]+vbuffer[i]*yvort[i]+wbuffer[i]*zvort[i])/
					sqrt(ubuffer[i]*ubuffer[i]+vbuffer[i]*vbuffer[i]+wbuffer[i]*wbuffer[i]);
			}
		}
		else if(!strcmp(varname[ivar],"streamfrac")) /* streamwise vorticity */
		{
			for(i=0; i<snx*sny*snz; i++)
			{
				buffer[i] = (ubuffer[i]*xvort[i]+vbuffer[i]*yvort[i]+wbuffer[i]*zvort[i])/
					(sqrt(ubuffer[i]*ubuffer[i]+vbuffer[i]*vbuffer[i]+wbuffer[i]*wbuffer[i])*
					sqrt(xvort[i]*xvort[i]+yvort[i]*yvort[i]+zvort[i]*zvort[i]));
			}
		}
		else if(!strcmp(varname[ivar],"qcqi"))
		{
			read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qc",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qi",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			for(i=0; i<snx*sny*snz; i++)
			{
				buffer[i] = ubuffer[i]+vbuffer[i];
			}
		}
			//float *thpert,*qvpert,*qtot,*ql
		else if(!strcmp(varname[ivar],"thrhopert1"))
		{
			float *thpert_thrho,*qvpert_thrho; /* Local copies so as not to screw up other variables */
			if ((qvar1 = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate qvar1");
			if ((qvar2 = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate qvar2");
			if ((qvar3 = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate qvar2");
			reps = 461.5 / 287.04; //from constants.incl in cm1
			read_hdf_mult_md(qvar1,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qc",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(qvar2,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qg",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			for(i=0; i<snx*sny*snz; i++)
			{
				qvar3[i] = 0.001*(qvar1[i]+qvar2[i]);
			}
			read_hdf_mult_md(qvar2,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qr",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			for(i=0; i<snx*sny*snz; i++)
			{
				qvar3[i] = qvar3[i] + 0.001*qvar2[i];
			}
			qtot = qvar3;
			read_hdf_mult_md(qvar1,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"thpert",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			thpert_thrho = qvar1;
			read_hdf_mult_md(qvar2,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qvpert",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			/* NOTE ORF!!! qvpert now g/kg in some simulations!!!!! */
			qvpert_thrho = qvar2;
			ii = 0;
			for(k=0; k<snz; k++)
			{
				for(j=0; j<sny; j++)
				{
					for(i=0; i<snx; i++)
					{
						float theta,num1,denom1,thbar,num2,denom2;
						theta = th0[k]+thpert_thrho[P3(i,j,k,snx,sny)];
						num1 = 1.0+reps*(qv0[k]+0.001*qvpert_thrho[P3(i,j,k,snx,sny)]);
						denom1 = 1.0+qtot[P3(i,j,k,snx,sny)]+qv0[k]+0.001*qvpert_thrho[P3(i,j,k,snx,sny)];
						thbar = th0[k];
						num2 = 1.0+reps*qv0[k];
						denom2 = 1.0+qv0[k];

						buffer[ii] = theta*(num1/denom1) - thbar*(num2/denom2);

//						if(k==0)printf("theta = %f\tnum1 = %f\tdenom1 = %f\t thbar = %f\tnum2 = %f\tdenom2 = %f\tthrho_pert = %f\n",theta,num1,denom1,thbar,num2,denom2,buffer[ii]);
							
						ii++;
					}
				}
			}
			free(qvar1); free(qvar2); free(qvar3);
		}
		else if(!strcmp(varname[ivar],"thrhopert2"))
		{
			if ((qvar1 = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate qvar1");
			if ((qvar2 = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate qvar2");
			if ((qvar3 = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate qvar2");
			reps = 461.5 / 287.04; //from constants.incl in cm1
			read_hdf_mult_md(qvar1,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qc",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(qvar2,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qg",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			for(i=0; i<snx*sny*snz; i++)
			{
				qvar3[i] = qvar1[i]+qvar2[i];
			}
			read_hdf_mult_md(qvar2,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qr",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			for(i=0; i<snx*sny*snz; i++)
			{
				qvar3[i] = qvar3[i] + qvar2[i];
			}
			qtot = qvar3;
			read_hdf_mult_md(qvar1,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"thpert",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			thpert = qvar1;
			read_hdf_mult_md(qvar2,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qvpert",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			qvpert = qvar2;
			ii = 0;
			for(k=0; k<snz; k++)
			{
				for(j=0; j<sny; j++)
				{
					for(i=0; i<snx; i++)
					{
						float theta,qvactual,qcond,thbar,qvbar;
						theta = th0[k]+thpert[P3(i,j,k,snx,sny)];
						qvactual = qv0[k]+qvpert[P3(i,j,k,snx,sny)];
						qcond = qtot[P3(i,j,k,snx,sny)];
						thbar = th0[k];
						qvbar = qv0[k];

						buffer[ii] = theta*(1.0+0.61*qvactual-qcond) - thbar*(1.0+0.61*qvbar);
							
						ii++;
					}
				}
			}
			free(qvar1); free(qvar2); free(qvar3);
		}
		else if(!strcmp(varname[ivar],"thrhopert_diff_1e-1"))
		{
			read_hdf_mult_md(dum1,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"thrhopert_lossless",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(dum2,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"thrhopert_1e-1",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			ii=0;
			for(k=0; k<snz; k++)
			for(j=0; j<sny; j++)
			for(i=0; i<snx; i++)
			{
				ii=P3(i,j,k,snx,sny);
				buffer[ii] = dum1[ii]-dum2[ii];
			}
		}
		else if(!strcmp(varname[ivar],"thrhopert_diff_1e-2"))
		{
			read_hdf_mult_md(dum1,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"thrhopert_lossless",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(dum2,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"thrhopert_1e-2",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			ii=0;
			for(k=0; k<snz; k++)
			for(j=0; j<sny; j++)
			for(i=0; i<snx; i++)
			{
				ii=P3(i,j,k,snx,sny);
				buffer[ii] = dum1[ii]-dum2[ii];
			}
		}
		else if(!strcmp(varname[ivar],"thrhopert_diff_1e-3"))
		{
			read_hdf_mult_md(dum1,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"thrhopert_lossless",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(dum2,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"thrhopert_1e-3",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			ii=0;
			for(k=0; k<snz; k++)
			for(j=0; j<sny; j++)
			for(i=0; i<snx; i++)
			{
				ii=P3(i,j,k,snx,sny);
				buffer[ii] = dum1[ii]-dum2[ii];
			}
		}
		else if(!strcmp(varname[ivar],"qc_diff_1"))
		{
			read_hdf_mult_md(dum1,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qc_lossless",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(dum2,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qc_1",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			ii=0;
			for(k=0; k<snz; k++)
			for(j=0; j<sny; j++)
			for(i=0; i<snx; i++)
			{
				ii=P3(i,j,k,snx,sny);
				buffer[ii] = dum1[ii]-dum2[ii];
			}
		}
		else if(!strcmp(varname[ivar],"qc_diff_1e-1"))
		{
			read_hdf_mult_md(dum1,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qc_lossless",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(dum2,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qc_1e-1",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			ii=0;
			for(k=0; k<snz; k++)
			for(j=0; j<sny; j++)
			for(i=0; i<snx; i++)
			{
				ii=P3(i,j,k,snx,sny);
				buffer[ii] = dum1[ii]-dum2[ii];
			}
		}
		else if(!strcmp(varname[ivar],"qc_diff_1e-2"))
		{
			read_hdf_mult_md(dum1,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qc_lossless",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(dum2,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qc_1e-2",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			ii=0;
			for(k=0; k<snz; k++)
			for(j=0; j<sny; j++)
			for(i=0; i<snx; i++)
			{
				ii=P3(i,j,k,snx,sny);
				buffer[ii] = dum1[ii]-dum2[ii];
			}
		}
		else if(!strcmp(varname[ivar],"zeta_diff_1e-1"))
		{
			read_hdf_mult_md(dum1,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qc_lossless",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(dum2,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qc_1e-1",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			ii=0;
			for(k=0; k<snz; k++)
			for(j=0; j<sny; j++)
			for(i=0; i<snx; i++)
			{
				ii=P3(i,j,k,snx,sny);
				buffer[ii] = dum1[ii]-dum2[ii];
			}
		}
		else
		{
			read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,varname[ivar],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
		}
		status = nc_put_vara_float (ncid, varnameid[ivar], start, edges, buffer);
		buffer = buf0; /* I will never learn my lesson with pointers........ */
		if (status != NC_NOERR) 
		{
			fprintf(stdout,"Could not write variable %s at time %i to %s\n", varname[ivar],t0,ncfilename);
			ERROR_STOP("Write to netcdf file failed");
		}
		fprintf(stdout,")\n");
	}
	status = nc_close(ncid); if (status != NC_NOERR)
	{
//		printf("ncid = %i\n",ncid);
	    fprintf(stderr, "%s\n", nc_strerror(status));
		printf("status = %i\n",status);
		fprintf(stderr, "Warning: netcdf is throwing an hdf error but our file seems to be fine...\n");
//		ERROR_STOP("Could not close netcdf file (nc_close)");
	}

	// Write command line option into file for future reference

	sprintf(cmdfilename,"%s%s",ncfilename,".cmd");
	if ((fp = fopen(cmdfilename,"w")) != NULL)
	{
		for (i=0; i<argc; i++)
		{
			fprintf(fp,"%s ",argv[i]);
		}
		fprintf(fp,"\n");
		fclose(fp);
	}
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
	char v5dvarname[MAXVARIABLES][10];
	char varname[MAXVARIABLES][10];
	int iret;

	float dxkm,dykm;

	int projection = 1; // 1 = linear, rectangular, cylindrical equidistant (must represent dx, dy in degrees)
	int vertical = 2;   // 2 = Stretched vertical grid, explicitly indicate height AGL in array
	float proj_args[4]; // Projection arguments for our grid 
	int *nla;

	hid_t file_id;


	if (argc == 3)	strcpy(datafile,argv[2]);
	else 			strcpy(datafile,"hdf2v5d.dat");


	if ((fp_data = fopen (datafile, "r")) == NULL)
	{
		char errorbuf[500];
		sprintf(errorbuf,"Could not open %s\n",datafile);
		ERROR_STOP(errorbuf);
	}

	if ((iret=fscanf (fp_data, "%s", topdir))==EOF)ERROR_STOP("fscanf failed");; lflush (fp_data);
	if ((iret=fscanf (fp_data, "%s", v5dfilename))==EOF)ERROR_STOP("fscanf failed"); lflush (fp_data);
    if ((iret=fscanf (fp_data, "%d%d", &xoff, &yoff))==EOF)ERROR_STOP("fscanf failed"); lflush (fp_data);
	if ((iret=fscanf (fp_data, "%d%d%d%d%d%d", &X0, &Y0, &X1, &Y1, &Z0, &Z1))==EOF)ERROR_STOP("fscanf failed"); lflush (fp_data);
	if ((iret=fscanf (fp_data, "%i", &dt))==EOF)ERROR_STOP("fscanf failed"); lflush (fp_data);
	if ((iret=fscanf (fp_data, "%d%d", &t0, &t1))==EOF)ERROR_STOP("fscanf failed"); lflush (fp_data);
	if ((iret=fscanf (fp_data, "%d", &compression))==EOF)ERROR_STOP("fscanf failed"); lflush (fp_data);
	if ((iret=fscanf (fp_data, "%d %f %f", &sbox, &umove, &vmove))==EOF)ERROR_STOP("fscanf failed"); lflush (fp_data);
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
	if ((iret=sprintf (cpcmd, "cp %s hdf2v5d.dat.%s",datafile,v5dfilename))<0)ERROR_STOP("sprintf failed");
	if ((iret=system (cpcmd))==-1)ERROR_STOP("system copy command failed");

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
	v5dCreate (v5dfilename, nblocks, nv5dvars, sny, snx, nla, (const char (*)[10]) v5dvarname, timestamp, datestamp, compression, projection, proj_args, vertical, zh);

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

void	parse_cmdline_hdf2nc(int argc, char *argv[],
	char *histpath, int *got_histpath, char *ncbase, int *got_ncbase,
	double *time, int *got_time,
	int *X0, int *got_X0,
	int *Y0, int *got_Y0,
	int *X1, int *got_X1,
	int *Y1, int *got_Y1,
	int *Z0, int *got_Z0,
	int *Z1, int *got_Z1)
{
	enum { OPT_HISTPATH = 1000, OPT_NCBASE, OPT_TIME, OPT_X0, OPT_Y0, OPT_X1, OPT_Y1, OPT_Z0, OPT_Z1 };
	// see https://stackoverflow.com/questions/23758570/c-getopt-long-only-without-alias

	int bail = 0;

	if (argc == 1)
	{
		fprintf(stderr,
		"Usage: %s -histpath [histpath] -ncbase [ncbase] -x0 [X0] -y0 [Y0] -x1 [X1] -y1 [Y1] -z0 [Z0] -z1 [Z1] -time [time] [varname1 ... varnameX] \n",argv[0]);
		exit(0);
	}


	while (1)
	{
		static struct option long_options[] =
		{
			{"histpath", required_argument, 0, OPT_HISTPATH},
			{"ncbase",   required_argument, 0, OPT_NCBASE},
			{"time",     required_argument, 0, OPT_TIME},
			{"x0",       required_argument, 0, OPT_X0},
			{"y0",       required_argument, 0, OPT_Y0},
			{"x1",       required_argument, 0, OPT_X1},
			{"y1",       required_argument, 0, OPT_Y1},
			{"z0",       required_argument, 0, OPT_Z0},
			{"z1",       required_argument, 0, OPT_Z1}
		};

		int r;
		int option_index = 0;
		r = getopt_long_only (argc, argv,"",long_options,&option_index);
		if (r == -1) break;

		switch(r)
		{
			case OPT_HISTPATH:
				strcpy(histpath,optarg);
				*got_histpath=1;
				printf("histpath = %s\n",histpath);
				break;
			case OPT_NCBASE:
				strcpy(ncbase,optarg);
				*got_ncbase=1;
				printf("ncbase = %s\n",ncbase);
				break;
			case OPT_TIME:
				*time = atof(optarg);
				*got_time=1;
				printf("time = %f\n",*time);
				break;
			case OPT_X0:
				*X0 = atoi(optarg);
				*got_X0=1;
				printf("X0 = %i\n",*X0);
				break;
			case OPT_Y0:
				*Y0 = atoi(optarg);
				*got_Y0=1;
				printf("Y0 = %i\n",*Y0);
				break;
			case OPT_X1:
				*X1 = atoi(optarg);
				*got_X1=1;
				printf("X1 = %i\n",*X1);
				break;
			case OPT_Y1:
				*Y1 = atoi(optarg);
				*got_Y1=1;
				printf("Y1 = %i\n",*Y1);
				break;
			case OPT_Z0:
				*Z0 = atoi(optarg);
				*got_Z0=1;
				printf("Z0 = %i\n",*Z0);
				break;
			case OPT_Z1:
				*Z1 = atoi(optarg);
				*got_Z1=1;
				printf("Z1 = %i\n",*Z1);
				break;
		}
	}

		if (*got_histpath==0) { fprintf(stderr,"-histpath not specified\n"); bail = 1; }
		if (*got_ncbase==0)   { fprintf(stderr,"-ncbase not specified\n"); bail = 1; }
		if (*got_time==0)   { fprintf(stderr,"-time not specified\n"); bail = 1; }
		if (*got_X0==0)     { fprintf(stderr,"-x0 not specified\n"); bail = 1; }
		if (*got_Y0==0)     { fprintf(stderr,"-y0 not specified\n"); bail = 1; }
		if (*got_X1==0)     { fprintf(stderr,"-x1 not specified\n"); bail = 1; }
		if (*got_Y1==0)     { fprintf(stderr,"-y1 not specified\n"); bail = 1; }
		if (*got_Z0==0)     { fprintf(stderr,"-z0 not specified\n"); bail = 1; }
		if (*got_Z1==0)     { fprintf(stderr,"-z1 not specified\n"); bail = 1; }
		if (bail)           { fprintf(stderr,"Insufficient arguments to %s\n, exiting.\n",argv[0]); exit(-1); }
}
