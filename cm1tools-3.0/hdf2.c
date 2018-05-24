/* LOFS (Leigh Orf File System / Lack Of File System) tools to convert
 * to other perhaps more useful formats
 *
 * Big rewrite June 2011. All top-level conversion code in one file, makes
 * maintenance easier. Symlinks to binary determine what is run.
 *
 * Major cleanup 11/23/2016, cleared out all unused variables and got
 * gcc to mostly shut up. Getting this code stable before converting to
 * new cm1hdf5 (subsecond) LOFS format.
 *
 * 2/17: Switched over to floating point format for times in LOFS
 *
 * Added command line parsing 7/17
 *
 * Cleaned a bunch of stuff up 3/18, see git log
 *
 * 4/18: added OMP loop directives for our calculation loops
 */

#include "lofs-read.h"
#include <omp.h>
#include <time.h>

#define MAXVARIABLES (100)
#define MAXSTR (512)

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
int saved_X0,saved_Y0,saved_X1,saved_Y1;
const float MISSING=1.0E37;

int debug = 0;
int yes2d = 0;
int gzip = 0;
int filetype = NC_NETCDF4;
int saved_staggered_mesh_params = 0;
int nthreads = 1;


//Minimum number of required arguments to hdf2nc. Adding optional flags (to
//hdf2nc) will require incrementing in order to retrieve all the
//variable names. Make this a global just for simplicity.
int argc_hdf2nc_min=4;
int optcount=0;

void grok_cm1hdf5_file_structure();
void hdf2nc(int argc, char *argv[], char *base, int X0, int Y0, int X1, int Y1, int Z0, int Z1, double t0);
void makevisit(int argc, char *argv[], char *base,int X0,int Y0,int X1,int Y1,int Z0,int Z1);

void parse_cmdline_hdf2nc(int argc, char *argv[],
		char *histpath, char *base,
		double *time, int *X0, int *Y0, int *X1, int *Y1, int *Z0, int *Z1);

void parse_cmdline_makevisit(int argc, char *argv[],
		char *histpath, char *base,
		int *X0, int *Y0, int *X1, int *Y1, int *Z0, int *Z1);

extern char *optarg; /* This is handled by the getopt code */

int main(int argc, char *argv[])
{
	char *cptr;
	char progname[MAXSTR];
	char histpath[MAXSTR];
	int we_are_hdf2nc = FALSE;
	int we_are_makevisit = FALSE;
//	int we_are_hdf2v5d = FALSE;
//	int we_are_linkfiles = FALSE;
	int X0,Y0,X1,Y1,Z0,Z1;
	double time;

	strcpy(progname,argv[0]);
	if (strspn("hdf2nc",progname) == 6) we_are_hdf2nc = TRUE;
	else if (strspn("makevisit",progname) == 9) we_are_makevisit = TRUE;
	else
	{
		fprintf(stderr,"progname = %s\n",progname);
		ERROR_STOP("Must call program as either hdf2nc or makevisit");
	}

	X0=Y0=X1=Y1=Z0=Z1=-1;// set to bogus values
	time=0.0;
	if      (we_are_hdf2nc)    parse_cmdline_hdf2nc(argc, argv, histpath, base, &time, &X0, &Y0, &X1, &Y1, &Z0, &Z1 );
	else if (we_are_makevisit) parse_cmdline_makevisit(argc, argv, histpath, base, &X0, &Y0, &X1, &Y1, &Z0, &Z1);

	if((cptr=realpath(histpath,topdir))==NULL)ERROR_STOP("realpath failed");
	grok_cm1hdf5_file_structure();
	get_hdf_metadata(firstfilename,&nx,&ny,&nz,&nodex,&nodey);
	if(debug) printf("DEBUG: nx = %i ny = %i nz = %i nodex = %i nodey = %i\n", nx,ny,nz,nodex,nodey);
	/* If we didn't specify values at the command line, set them to values specifying all the saved data */
	if(X0<0)X0=saved_X0; if(Y0<0)Y0=saved_Y0; if(Z0<0)Z0=0;
	if(X1<0)X1=saved_X1; if(Y1<0)Y1=saved_Y1; if(Z1<0)Z1=nz-1;

	/* If our supplied indices are outside of the saved data,
	 * warn and adjust accordingly */

	/* First, look for idiocy */

	if (X0>X1||Y0>Y1||X1<saved_X0||Y1<saved_Y0||X0>saved_X1||Y0>saved_Y1)
	{
		printf(" *** X0=%i saved_X0=%i Y0=%i saved_Y0=%i X1=%i saved_X1=%i Y1=%i saved_Y1=%i\n",
				X0,saved_X0,Y0,saved_Y0,X1,saved_X1,Y1,saved_Y1);
		ERROR_STOP("Your requested indices are wack - check for weirdness at the command line\n");
	}
	if(X0<saved_X0)
	{
		printf("Oops: requested out of box: Adjusting X0 (%i) to saved_X0 (%i)\n",X0,saved_X0);
		X0=saved_X0;
	}
	if(Y0<saved_Y0)
	{
		printf("Oops: requested out of box: Adjusting Y0 (%i) to saved_Y0 (%i)\n",Y0,saved_Y0);
		Y0=saved_Y0;
	}
	if(X1>saved_X1)
	{
		printf("Oops: requested out of box: Adjusting X1 (%i) to saved_X1 (%i)\n",X1,saved_X1);
		X1=saved_X1;
	}
	if(Y1>saved_Y1)
	{
		printf("Oops: requested out of box: Adjusting Y1 (%i) to saved_Y1 (%i)\n",Y1,saved_Y1);
		X1=saved_X1;
	}

	if(we_are_hdf2nc)             hdf2nc(argc,argv,base,X0,Y0,X1,Y1,Z0,Z1,time);
	else if (we_are_makevisit) makevisit(argc,argv,base,X0,Y0,X1,Y1,Z0,Z1);
	exit(0);
}

void makevisit(int argc, char *argv[], char *base,int X0,int Y0,int X1,int Y1,int Z0,int Z1)
{
		int i,rank,nvars;
		hid_t f_id,g_id,strtype;
		H5G_info_t group_info;
		hsize_t dims[1];
		char visitfile[MAXSTR];
		char groupname[MAXSTR];
		char varname[MAXVARIABLES][40]; //Yeah I know also hardcoded in avtcm1visitFileFormat.h

		float *xhfull,*yhfull,*zh;

		sprintf(visitfile,"%s.cm1visit",base);
		printf("visitfile = %s\n",visitfile);
		printf("topdir = %s\n",topdir);
		printf("ntimedirs = %i\n",ntimedirs);
		if (debug) {for (i=0; i<ntimedirs; i++)printf("%s ",timedir[i]);printf("\n");}
		if (debug) {for (i=0; i<ntimedirs; i++)printf("%f ",dirtimes[i]);printf("\n");}
		printf("nnodedirs = %i\n",nnodedirs);
		if (debug) {for (i=0; i<nnodedirs; i++)printf("%s ",nodedir[i]);printf("\n");}
		printf("dn = %i\n",dn);
		printf("ntottimes = %i\n",ntottimes);
		if(debug) {for (i=0; i<ntottimes; i++)printf("%f ",alltimes[i]);printf("\n");}
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
		sprintf(groupname,"%05i/3D",0); //is now last since we no longer have get_first_hdf_file_name routine;

/* first file name is the last "fist file name" from sweeping through
all times, so now "first file name" is from last time direcotry. This is
a little convoluted now that I have to deal with empty node dirs from
saving subdomains. */

		printf("firstfilename = %s\n",firstfilename);
		printf("groupname = %s\n",groupname);
		g_id = H5Gopen(f_id,groupname,H5P_DEFAULT);
		H5Gget_info(g_id,&group_info);
		nvars = group_info.nlinks;
		for (i = 0; i < nvars; i++)
		{
		    H5Lget_name_by_idx(g_id,".",H5_INDEX_NAME,H5_ITER_INC,i,varname[i],40,H5P_DEFAULT); //40 characters per varname
		}
		H5Gclose(g_id);
		H5Fclose(f_id);

		printf("Variables retrieved: ");
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

//ORF: 2018-03-15 I have finally gotten rid of snx0 etc. and replaced
//with X0 etc. For reasons I no longer remember I used one variable name
//approach with hdf2nc and one with makevisit. However so I don't break
//things when reading old cm1visit files I will still write snx0 etc. so
//I don't have to deal with breakage. May never remove this, who knows :)
//Since I'm using the lite interface I can't do a link, so we just
//duplicate.

		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/X0",   rank, dims, &X0);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/snx0", rank, dims, &X0);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/X1",   rank, dims, &X1);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/snx1", rank, dims, &X1);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/Y0",   rank, dims, &Y0);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/sny0", rank, dims, &Y0);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/Y1",   rank, dims, &Y1);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/sny1", rank, dims, &Y1);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/Z0",   rank, dims, &Z0);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/snz0", rank, dims, &Z0);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/Z1",   rank, dims, &Z1);
		rank=1;dims[0]=1; H5LTmake_dataset_int (f_id, "/snz1", rank, dims, &Z1);

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

void grok_cm1hdf5_file_structure()
{
	int i;
	ntimedirs = get_num_time_dirs(topdir,debug); printf("ntimedirs: %i\n",ntimedirs);
	if (ntimedirs == 0) ERROR_STOP("No cm1 hdf5 files found");

	timedir = (char **)malloc(ntimedirs * sizeof(char *));
	for (i=0; i < ntimedirs; i++) timedir[i] = (char *)(malloc(MAXSTR * sizeof(char)));
	dirtimes = (double *)malloc(ntimedirs * sizeof(double)); // ORF NO LONGER INT, RIGHT???

	get_sorted_time_dirs(topdir,timedir,dirtimes,ntimedirs,base,debug);

	nnodedirs =  get_num_node_dirs(topdir,timedir[0],debug);
	nodedir = (char **)malloc(nnodedirs * sizeof(char *));
	// ORF 8 == 7 zero padded node number directory name plus 1 end of string char
	for (i=0; i < nnodedirs; i++) nodedir[i] = (char *)(malloc(8 * sizeof(char)));

	get_sorted_node_dirs(topdir,timedir[0],nodedir,&dn,nnodedirs,debug);

	alltimes = get_all_available_times(topdir,timedir,ntimedirs,nodedir,nnodedirs,&ntottimes,firstfilename,&firsttimedirindex,
			&saved_X0,&saved_Y0,&saved_X1,&saved_Y1,debug);
	if(debug)
	{
		printf("All available times: ");
		for (i=0; i<ntottimes; i++)printf("%lf ",alltimes[i]);
		printf("\n");
	}
}

void hdf2nc(int argc, char *argv[], char *base, int X0, int Y0, int X1, int Y1, int Z0, int Z1, double t0)
{
	float *buffer,*buf0,*ubuffer,*vbuffer,*wbuffer,*xvort,*yvort,*zvort;
	float *th0,*qv0;
	float *thrhopert;
	double timearray[1];

	int i,j,k,ix,iy,iz,nvar;
	char varname[MAXVARIABLES][MAXSTR];
	char ncfilename[MAXSTR];

	extern int H5Z_zfp_initialize(void);

	int NX,NY,NZ;
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
	float *xhfull,*yhfull,*xffull,*yffull,*zh,*zf;
	float *xhout,*yhout,*zhout,*xfout,*yfout,*zfout;
	FILE *fp;
	char cmdfilename[512];
	//ORF for writing single time in unlimited time dimension/variable
	const size_t timestart = 0;
	const size_t timecount = 1;
	/* readahead flags */
	int u_rh=0,v_rh=0,w_rh=0,xvort_rh=0,yvort_rh=0,zvort_rh=0,thrhopert_rh=0;

	float rv = 461.5;
	float rd = 287.04;
	float reps;

	reps = rv/rd;

	/* Since we tack on all the requested variables to the end of the
	 * command line, we have to find out where the 1st variable
	 * argument is. Since we have optional arguments, we keep track of
	 * things and then peel off the file name strings */

	nvar = argc - argc_hdf2nc_min - optcount;

	if (debug) printf("argc = %i, nvar = %i, optcount = %i\n",argc,nvar,optcount);

	printf("\nWe are requesting the following fields: ");
	for (i=0; i<nvar; i++)
	{
		strcpy(varname[i],argv[i+argc_hdf2nc_min+optcount]);
		printf("%s ",varname[i]);
	}
	printf("\n");
	sprintf(ncfilename,"%s.%012.6f.nc",base,t0);
	
	NX = X1 - X0 + 1;
	NY = Y1 - Y0 + 1;
	NZ = Z1 - Z0 + 1;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	start[3] = 0;
	edges[0] = 1;
	edges[1] = NZ;
	edges[2] = NY;
	edges[3] = NX;
	//For 2D surface slices
	s2[0] = 0;
	s2[1] = 0;
	s2[2] = 0;
	e2[0] = 1;
	e2[1] = NY;
	e2[2] = NX;

	/* ORF LAZY allocate enough for any staggered combination, hence +1 for all three dimensions */
	bufsize = (long) (NX+1) * (long) (NY+1) * (long) (NZ+1) * (long) sizeof(float);
	if(debug) fprintf(stdout,"X0=%i Y0=%i X1=%i Y1=%i Z0=%i Z1=%i bufsize = %f GB\n",X0,Y0,X1,Y1,Z0,Z1,1.0e-9*bufsize);
	if ((buf0 = buffer = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((ubuffer = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((vbuffer = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((wbuffer = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((xvort = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((yvort = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((zvort = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	if ((thrhopert = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate 3D buffer");
	//printf("NX = %i NY = %i NZ = %i Bufsize = %i\n",NX,NY,NZ,bufsize);
	if (buffer == NULL) ERROR_STOP("Cannot allocate buffer");
	if ((f_id = H5Fopen (firstfilename, H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
	{
		fprintf(stderr,"Cannot open firstfilename which is %s, even though we have alredy opened it!\n",firstfilename);
		ERROR_STOP("Cannot open hdf file");
	}
	xhfull = (float *)malloc(nx * sizeof(float));
	yhfull = (float *)malloc(ny * sizeof(float));
	xffull = (float *)malloc((nx+1) * sizeof(float));
	yffull = (float *)malloc((ny+1) * sizeof(float));
	th0 = (float *)malloc(nz * sizeof(float));
	qv0 = (float *)malloc(nz * sizeof(float));
	zh = (float *)malloc(nz * sizeof(float));
	zf = (float *)malloc(nz * sizeof(float));
	get1dfloat (f_id,(char *)"mesh/xhfull",xhfull,0,nx);
	get1dfloat (f_id,(char *)"mesh/yhfull",yhfull,0,ny);
/* THIS NEEDS TO BE TESTED IN CM1 FIRST... we haven't been saving these */
	if (saved_staggered_mesh_params)
	{
	 	get1dfloat (f_id,(char *)"mesh/xffull",xffull,0,nx+1);
		get1dfloat (f_id,(char *)"mesh/yffull",yffull,0,ny+1);
	}
	get1dfloat (f_id,(char *)"mesh/zh",zh,0,nz);
	get1dfloat (f_id,(char *)"mesh/zf",zf,0,nz);
	get1dfloat (f_id,(char *)"basestate/qv0",qv0,0,nz);
	get1dfloat (f_id,(char *)"basestate/th0",th0,0,nz);
	get1dfloat (f_id,(char *)"mesh/zf",zf,0,nz);

	xhout = (float *)malloc(NX * sizeof(float));
	yhout = (float *)malloc(NY * sizeof(float));
	zhout = (float *)malloc(NZ * sizeof(float));
	xfout = (float *)malloc(NX * sizeof(float));
	yfout = (float *)malloc(NY * sizeof(float)); // +1 when we do this right
	zfout = (float *)malloc(NZ * sizeof(float)); // +1 when we do this right

	for (iz=Z0; iz<=Z1; iz++) zhout[iz-Z0] = 0.001*zh[iz];
	for (iz=Z0; iz<=Z1; iz++) zfout[iz-Z0] = 0.001*zf[iz];     //NEED TO READ REAL ZF

	if (saved_staggered_mesh_params)
	{
		for (iy=Y0; iy<=Y1; iy++) yfout[iy-Y0] = 0.001*yffull[iy];
		for (ix=X0; ix<=X1; ix++) xfout[ix-X0] = 0.001*xffull[ix];
	}
	for (iy=Y0; iy<=Y1; iy++) yhout[iy-Y0] = 0.001*yhfull[iy];
	for (ix=X0; ix<=X1; ix++) xhout[ix-X0] = 0.001*xhfull[ix];

	H5Z_zfp_initialize();


/*

netCDF files come in several flavors... This is mostly due to the
passage of time - it's an old format! And files have gotten bigger!
And Leon has gotten larger!

The best is NC_NETCDF4, which is HDF5 but using the netCDF API...
With HDF5 you can make the hugest-assed files imaginable with nearly
unlimited variables/variable sizes. If netCDF4 doesn't work you probably
don't have it built into your netCDF distribution (it requires HDF5).

The second best is NC_64BIT_OFFSET (Can make files bigger than 2 GB)

See this hopefully not dead link for more:
http://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf/Large-File-Support.html#Large-File-Support

*/

	status = nc_create (ncfilename, NC_CLOBBER|filetype, &ncid); if (status != NC_NOERR) ERROR_STOP ("nc_create failed");


	status = nc_def_dim (ncid, "xh", NX, &nxh_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed"); 
	status = nc_def_dim (ncid, "yh", NY, &nyh_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (ncid, "zh", NZ, &nzh_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	if (saved_staggered_mesh_params)
	{
		status = nc_def_dim (ncid, "xf", NX, &nxf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed"); 
		status = nc_def_dim (ncid, "yf", NY, &nyf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	}
	status = nc_def_dim (ncid, "zf", NZ, &nzf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (ncid, "time", NC_UNLIMITED, &time_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_var (ncid, "xh", NC_FLOAT, 1, &nxh_dimid, &xhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "yh", NC_FLOAT, 1, &nyh_dimid, &yhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "zh", NC_FLOAT, 1, &nzh_dimid, &zhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	if (saved_staggered_mesh_params)
	{
		status = nc_def_var (ncid, "xf", NC_FLOAT, 1, &nxf_dimid, &xfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
		status = nc_def_var (ncid, "yf", NC_FLOAT, 1, &nyf_dimid, &yfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	}
	status = nc_def_var (ncid, "zf", NC_FLOAT, 1, &nzf_dimid, &zfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "time", NC_DOUBLE, 1, &time_dimid, &timeid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_put_att_text(ncid, xhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, yhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, zhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	if (saved_staggered_mesh_params)
	{
		status = nc_put_att_text(ncid, xfid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		status = nc_put_att_text(ncid, yfid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	}
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

	/* Save some surface 2d slices? */
	if (yes2d)
	{
		status = nc_def_var (ncid, "thrhopert_sfc", NC_FLOAT, 3, d2, &thsfcid);
		status = nc_def_var (ncid, "dbz_sfc", NC_FLOAT, 3, d2, &dbzsfcid);
	}

	for (ivar = 0; ivar < nvar; ivar++)
	{
		if(!strcmp(varname[ivar],"u"))
		{
			if(!saved_staggered_mesh_params) ERROR_STOP("We are asking for u but did not save xf mesh array, sorry!");
			dims[0] = time_dimid;
			dims[1] = nzh_dimid;
			dims[2] = nyh_dimid;
			dims[3] = nxf_dimid;
		}
		else if (!strcmp(varname[ivar],"v"))
		{
			if(!saved_staggered_mesh_params) ERROR_STOP("We are asking for v but did not save yf mesh array, sorry!");
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
		status = nc_put_att_float(ncid,varnameid[ivar],"missing_value",NC_FLOAT,1,&MISSING);
//		unfortunately this really slows things down
		if (gzip) status=nc_def_var_deflate(ncid, varnameid[ivar], 1, 1, 1);
	}
	status = nc_enddef (ncid);

	if (status != NC_NOERR) ERROR_STOP("nc_enddef failed");

	// For surface 2d fields (Vapor wants 2D vars)
	//
	// Save surface thrhopert and dbz for easy viewing as 2D vars in Vapor
	// Need also to get 2D data into cm1hdf5 files and then VisIt (TODO)
	// This is now only enabled by a command line option - default
	// is to not write any 2D fields.

	if (yes2d)
	{
		printf("\nWorking on surface 2D thrhopert and dbz (");
		read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"thrhopert",X0,Y0,X1,Y1,0,0,nx,ny,nz,nodex,nodey);
		status = nc_put_vara_float (ncid, thsfcid, s2, e2, buffer);
		read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"dbz",X0,Y0,X1,Y1,0,0,nx,ny,nz,nodex,nodey);
		status = nc_put_vara_float (ncid, dbzsfcid, s2, e2, buffer);
		printf(")\n");
	}

      status = nc_put_var_float (ncid,xhid,xhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_float (ncid,yhid,yhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_float (ncid,zhid,zhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	if (saved_staggered_mesh_params)
	{
		status = nc_put_var_float (ncid,xfid,xfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
		status = nc_put_var_float (ncid,yfid,yfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	}
	status = nc_put_var_float (ncid,zfid,zfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_int (ncid,x0id,&X0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,y0id,&Y0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,x1id,&X1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,y1id,&Y1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,z0id,&Z0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,z1id,&Z1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	timearray[0] = t0;
      status = nc_put_vara_double (ncid,timeid,&timestart,&timecount,timearray); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");

//	int u_rh=0,v_rh=0,w_rh=0,xvort_rh=0,yvort_rh=0,zvort_rh=0;

	/* So we don't read data more than once, flag stuff that we need
	 * ahead of time (read ahead) */
	for (ivar = 0; ivar < nvar; ivar++)
	{
		if(!strcmp(varname[ivar],"hwin_sr")) {u_rh=v_rh=1;}
		if(!strcmp(varname[ivar],"hwin_gr")) {u_rh=v_rh=1;}
		if(!strcmp(varname[ivar],"hdiv")) {u_rh=v_rh=1;}
		if(!strcmp(varname[ivar],"hvort")) {xvort_rh=yvort_rh=1;}
		if(!strcmp(varname[ivar],"vortmag")) {yvort_rh=yvort_rh=zvort_rh=1;}
		if(!strcmp(varname[ivar],"streamvort")) {u_rh=v_rh=w_rh=xvort_rh=yvort_rh=zvort_rh=1;}
		if(!strcmp(varname[ivar],"streamfrac")) {u_rh=v_rh=w_rh=xvort_rh=yvort_rh=zvort_rh=1;}
		if(!strcmp(varname[ivar],"xvort_baro")) {thrhopert_rh=1;}
		if(!strcmp(varname[ivar],"yvort_baro")) {thrhopert_rh=1;}
		if(!strcmp(varname[ivar],"yvort_stretch")) {u_rh=w_rh=yvort_rh=1;}
		if(!strcmp(varname[ivar],"yvort_tilt")) {v_rh=xvort_rh=yvort_rh=1;}
		if(!strcmp(varname[ivar],"xvort_stretch")) {v_rh=w_rh=xvort_rh=1;}
		if(!strcmp(varname[ivar],"xvort_tilt")) {u_rh=yvort_rh=zvort_rh=1;}
		if(!strcmp(varname[ivar],"zvort_stretch")) {u_rh=v_rh=zvort_rh=1;}
		if(!strcmp(varname[ivar],"zvort_tilt")) {w_rh=xvort_rh=zvort_rh=1;}
	}

	if (u_rh)
	{
		printf("Buffering uinterp:\n");
		read_hdf_mult_md(ubuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"uinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
	}
	if (v_rh)
	{
		printf("Buffering vinterp:\n");
		read_hdf_mult_md(vbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"vinterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
	}
	if (w_rh)
	{
		printf("Buffering winterp:\n");
		read_hdf_mult_md(wbuffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"winterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
	}
	if (xvort_rh)
	{
		printf("Buffering xvort:\n");
		read_hdf_mult_md(xvort,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"xvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
	}
	if (yvort_rh)
	{
		printf("Buffering yvort:\n");
		read_hdf_mult_md(yvort,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"yvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
	}
	if (zvort_rh)
	{
		printf("Buffering zvort:\n");
		read_hdf_mult_md(zvort,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"zvort",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
	}
	if (thrhopert_rh)
	{
		printf("Buffering thrhopert:\n");
		read_hdf_mult_md(thrhopert,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"thrhopert",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
	}


// Here is wher you can write your own code to calculate new fields based upon the fields you are reading in.
// Note, if uinterp, vinterp, and winterp were saved but not u, v, w, then you will lose accuracy if calculating
// derivatives of these variables since they have already been interpolated to the scalar grid.

	for (ivar = 0; ivar < nvar; ivar++)
	{
		printf("Working on %s (",varname[ivar]);
		fflush(stdout);

		if(!strcmp(varname[ivar],"hwin_sr")) //storm relative horizontal wind speed
		{
			float usr,vsr;
#pragma omp parallel for private(i)
			for(i=0; i<NX*NY*NZ; i++)
			{
				usr = ubuffer[i];
				vsr = vbuffer[i];
				buffer[i] = sqrt(usr*usr+vsr*vsr);
			}
		}
		else if(!strcmp(varname[ivar],"hwin_gr")) //ground relative horizontal wind speed
		{
			float usr,vsr;
			//TODO: put umove and vmove in the got dammed hdf5
			//files
			float umove=15.2;
			float vmove=10.5;
#pragma omp parallel for private(i)
			for(i=0; i<NX*NY*NZ; i++)
			{
				usr = ubuffer[i]+umove;
				vsr = vbuffer[i]+vmove;
				buffer[i] = sqrt(usr*usr+vsr*vsr);
			}
		}
		else if(!strcmp(varname[ivar],"hdiv")) // need to save this, more accurate with staggered vel. vars
		{
			float dxi,dyi,dzi;

			/* set edges to zero - should fill with missing (TODO)*/
			for(k=0; k<NZ; k++)
			for(j=0; j<NY; j++)
				buffer[P3(0,j,k,NX,NY)] = 
				buffer[P3(NX-1,j,k,NX,NY)] = MISSING;

			for(k=0; k<NZ; k++)
			for(i=0; i<NX; i++)
				buffer[P3(i,0,k,NX,NY)] = 
				buffer[P3(i,NY-1,k,NX,NY)] = MISSING;

#pragma omp parallel for private(i,j,k)
			for(k=0; k<NZ; k++)
			{
				for(j=1; j<NY-1; j++)
				{
					dyi=1.0/(yhfull[iy-Y0+1]-yhfull[iy-Y0-1]);
					for(i=1; i<NX-1; i++)
					{
						dxi=1.0/(xhfull[ix-X0+1]-xhfull[ix-X0-1]);
						buffer[P3(i,j,k,NX,NY)] =
						dxi * (ubuffer[P3(i+1,j,k,NX,NY)] - ubuffer[P3(i-1,j,k,NX,NY)]) +
						dyi * (vbuffer[P3(i,j+1,k,NX,NY)] - vbuffer[P3(i,j-1,k,NX,NY)]) ;
					}
				}
			}
		}
		else if(!strcmp(varname[ivar],"zvort_stretch")) 
		{
			float dxi,dyi,dzi;

			for(k=0; k<NZ; k++)
			for(j=0; j<NY; j++)
				buffer[P3(0,j,k,NX,NY)] = 
				buffer[P3(NX-1,j,k,NX,NY)] = MISSING;

			for(k=0; k<NZ; k++)
			for(i=0; i<NX; i++)
				buffer[P3(i,0,k,NX,NY)] = 
				buffer[P3(i,NY-1,k,NX,NY)] = MISSING;

#pragma omp parallel for private(i,j,k)
			for(k=1; k<NZ-1; k++)
			{
				for(j=1; j<NY-1; j++)
				{
					dyi=1.0/(yhfull[iy-Y0+1]-yhfull[iy-Y0-1]);
					for(i=0; i<NX; i++)
					{
						dxi=1.0/(xhfull[ix-X0+1]-xhfull[ix-X0-1]);
						buffer[P3(i,j,k,NX,NY)] =
							-1000.0*zvort[P3(i,j,k,NX,NY)] *
						(dxi * (ubuffer[P3(i+1,j,k,NX,NY)] - ubuffer[P3(i-1,j,k,NX,NY)]) +
						dyi * (vbuffer[P3(i,j+1,k,NX,NY)] - vbuffer[P3(i,j-1,k,NX,NY)])) ;
					}
				}
			}
		}
		else if(!strcmp(varname[ivar],"zvort_tilt")) 
		{
			float dxi,dyi,dzi;

			for(k=0; k<NZ; k++)
			for(j=0; j<NY; j++)
				buffer[P3(0,j,k,NX,NY)] = 
				buffer[P3(NX-1,j,k,NX,NY)] = MISSING;

			for(k=0; k<NZ; k++)
			for(i=0; i<NX; i++)
				buffer[P3(i,0,k,NX,NY)] = 
				buffer[P3(i,NY-1,k,NX,NY)] = MISSING;

#pragma omp parallel for private(i,j,k)
			for(k=1; k<NZ-1; k++)
			{
				for(j=1; j<NY-1; j++)
				{
					dyi=1.0/(yhfull[iy-Y0+1]-yhfull[iy-Y0-1]);
					for(i=0; i<NX; i++)
					{
						dxi=1.0/(xhfull[ix-X0+1]-xhfull[ix-X0-1]);
						buffer[P3(i,j,k,NX,NY)] =
							1000.0*(xvort[P3(i,j,k,NX,NY)] * dxi * (wbuffer[P3(i+1,j,k,NX,NY)] - wbuffer[P3(i-1,j,k,NX,NY)]) +
						        yvort[P3(i,j,k,NX,NY)] * dyi * (wbuffer[P3(i,j+1,k,NX,NY)] - wbuffer[P3(i,j-1,k,NX,NY)])) ;
					}
				}
			}

		}
		else if(!strcmp(varname[ivar],"xvort_baro")) // need to save this, more accurate with staggered vel. vars
		{
			float dxi,dyi,dzi;
//			float coeff = 1000.0 * 9.8/304.86; /* UGLY I know... this should be saved in hist files */
			float coeff;
			/* BUG:  Needs to be theta_bar[iz] */

			/* set edges to zero - should fill with missing (TODO)*/

			for(k=0; k<NZ; k++)
			for(i=0; i<NX; i++)
				buffer[P3(i,0,k,NX,NY)] = 
				buffer[P3(i,NY-1,k,NX,NY)] = MISSING;

#pragma omp parallel for private(i,j,k)
			for(k=0; k<NZ; k++)
			{
				coeff = 1000.0 * 9.81 / (th0[k]*(1.0+reps*qv0[iz]/(1.0+qv0[iz])));
				for(j=1; j<NY-1; j++)
				{
					for(i=0; i<NX; i++)
					{
						dyi=1.0/(yhfull[iy-Y0+1]-yhfull[iy-Y0-1]);
						buffer[P3(i,j,k,NX,NY)] =
							coeff * dyi * (thrhopert[P3(i,j+1,k,NX,NY)] - thrhopert[P3(i,j-1,k,NX,NY)]) ;
					}
				}
			}
		}
		else if(!strcmp(varname[ivar],"yvort_baro")) // need to save this, more accurate with staggered vel. vars
		{
			float dxi,dyi,dzi;
//			float coeff = 1000.0 * 9.8/304.86; /* UGLY I know... this should be saved in hist files */
			float coeff;
			/* BUG:  Needs to be theta_bar[iz] */

			/* set edges to zero - should fill with missing (TODO)*/

			for(k=0; k<NZ; k++)
			for(j=0; j<NY; j++)
				buffer[P3(0,j,k,NX,NY)] = 
				buffer[P3(NX-1,j,k,NX,NY)] = MISSING;

#pragma omp parallel for private(i,j,k)
			for(k=0; k<NZ; k++)
			{
				coeff = 1000.0 * 9.81 / (th0[k]*(1.0+reps*qv0[iz]/(1.0+qv0[iz])));
				for(j=0; j<NY; j++)
				{
					for(i=1; i<NX-1; i++)
					{
						dxi=1.0/(xhfull[ix-X0+1]-xhfull[ix-X0-1]);
						buffer[P3(i,j,k,NX,NY)] =
							-coeff * dxi * (thrhopert[P3(i+1,j,k,NX,NY)] - thrhopert[P3(i-1,j,k,NX,NY)]) ;
					}
				}
			}
		}
		else if(!strcmp(varname[ivar],"yvort_stretch")) // need to save this, more accurate with staggered vel. vars
		{
			float dxi,dyi,dzi;
			for(k=0; k<NZ; k++)
			for(j=0; j<NY; j++)
				buffer[P3(0,j,k,NX,NY)] = 
				buffer[P3(NX-1,j,k,NX,NY)] = MISSING;

			for(j=0; j<NY; j++)
			for(i=0; i<NX; i++)
				buffer[P3(i,j,0,NX,NY)] = 
				buffer[P3(i,j,NZ-1,NX,NY)] = MISSING;

#pragma omp parallel for private(i,j,k)
			for(k=1; k<NZ-1; k++)
			{
				dzi=1.0/(zh[iz-Z0+1]-zh[iz-Z0-1]);
				for(j=0; j<NY; j++)
				{
					for(i=1; i<NX-1; i++)
					{
						dxi=1.0/(xhfull[ix-X0+1]-xhfull[ix-X0-1]);
						buffer[P3(i,j,k,NX,NY)] = -1000.0 * yvort[P3(i,j,k,NX,NY)]*(dxi*(ubuffer[P3(i+1,j,k,NX,NY)]-ubuffer[P3(i-1,j,k,NX,NY)]) +
												dzi*(wbuffer[P3(i,j,k+1,NX,NY)]-wbuffer[P3(i,j,k-1,NX,NY)]));
					}
				}
			}
		}
		else if(!strcmp(varname[ivar],"yvort_tilt")) // need to save this, more accurate with staggered vel. vars
		{
			float dxi,dyi,dzi;
			for(k=0; k<NZ; k++)
			for(j=0; j<NY; j++)
				buffer[P3(0,j,k,NX,NY)] = 
				buffer[P3(NX-1,j,k,NX,NY)] = MISSING;

			for(j=0; j<NY; j++)
			for(i=0; i<NX; i++)
				buffer[P3(i,j,0,NX,NY)] = 
				buffer[P3(i,j,NZ-1,NX,NY)] = MISSING;

#pragma omp parallel for private(i,j,k)
			for(k=1; k<NZ; k++)
			{
				dzi=1.0/(zh[iz-Z0+1]-zh[iz-Z0-1]);
				for(j=0; j<NY; j++)
				{
					for(i=1; i<NX-1; i++)
					{
						dxi=1.0/(xhfull[ix-X0+1]-xhfull[ix-X0-1]);
						buffer[P3(i,j,k,NX,NY)] = 1000.0 * (xvort[P3(i,j,k,NX,NY)]*dxi*(vbuffer[P3(i+1,j,k,NX,NY)]-vbuffer[P3(i+1,j,k,NX,NY)]) +
										    zvort[P3(i,j,k,NX,NY)]*dzi*(vbuffer[P3(i,j,k+1,NX,NY)]-vbuffer[P3(i,j,k-1,NX,NY)]));
					}
				}
			}
		}
		else if(!strcmp(varname[ivar],"xvort_stretch")) // need to save this, more accurate with staggered vel. vars
		{
			float dxi,dyi,dzi;
			for(k=0; k<NZ; k++)
			for(i=0; i<NX; i++)
				buffer[P3(i,0,k,NX,NY)] = 
				buffer[P3(i,NY-1,k,NX,NY)] = MISSING;

			for(j=0; j<NY; j++)
			for(i=0; i<NX; i++)
				buffer[P3(i,j,0,NX,NY)] = 
				buffer[P3(i,j,NZ-1,NX,NY)] = MISSING;

#pragma omp parallel for private(i,j,k)
			for(k=1; k<NZ-1; k++)
			{
				dzi=1.0/(zh[iz-Z0+1]-zh[iz-Z0-1]);
				for(j=1; j<NY-1; j++)
				{
					dyi=1.0/(yhfull[iy-Y0+1]-yhfull[iy-Y0-1]);
					for(i=0; i<NX; i++)
					{
						buffer[P3(i,j,k,NX,NY)] = -1000.0 * xvort[P3(i,j,k,NX,NY)]*(dyi*(vbuffer[P3(i,j+1,k,NX,NY)]-vbuffer[P3(i,j-1,k,NX,NY)]) +
												dzi*(wbuffer[P3(i,j,k+1,NX,NY)]-wbuffer[P3(i,j,k-1,NX,NY)]));
					}
				}
			}
		}
		else if(!strcmp(varname[ivar],"xvort_tilt")) // need to save this, more accurate with staggered vel. vars
		{
			float dxi,dyi,dzi;
			for(k=0; k<NZ; k++)
			for(i=0; i<NX; i++)
				buffer[P3(i,0,k,NX,NY)] = 
				buffer[P3(i,NY-1,k,NX,NY)] = MISSING;

			for(j=0; j<NY; j++)
			for(i=0; i<NX; i++)
				buffer[P3(i,j,0,NX,NY)] = 
				buffer[P3(i,j,NZ-1,NX,NY)] = MISSING;

#pragma omp parallel for private(i,j,k)
			for(k=1; k<NZ-1; k++)
			{
				dzi=1.0/(zh[iz-Z0+1]-zh[iz-Z0-1]);
				for(j=1; j<NY-1; j++)
				{
					dyi=1.0/(yhfull[iy-Y0+1]-yhfull[iy-Y0-1]);
					for(i=0; i<NX; i++)
					{
						buffer[P3(i,j,k,NX,NY)] = 1000.0 * (yvort[P3(i,j,k,NX,NY)]*dyi*(ubuffer[P3(i,j+1,k,NX,NY)]-ubuffer[P3(i,j-1,k,NX,NY)]) +
										    zvort[P3(i,j,k,NX,NY)]*dzi*(ubuffer[P3(i,j,k+1,NX,NY)]-ubuffer[P3(i,j,k-1,NX,NY)]));
					}
				}
			}
		}
		else if(!strcmp(varname[ivar],"hvort")) //horizontal vorticity magnitude
		{
#pragma omp parallel for private(i)
			for(i=0; i<NX*NY*NZ; i++)
			{
				buffer[i] = sqrt(xvort[i]*xvort[i]+yvort[i]*yvort[i]);
			}
		}
		else if(!strcmp(varname[ivar],"vortmag")) //3D vorticity magnitude
		{
#pragma omp parallel for private(i)
			for(i=0; i<NX*NY*NZ; i++)
			{
				buffer[i] = sqrt(xvort[i]*xvort[i]+yvort[i]*yvort[i]+zvort[i]*zvort[i]);
			}
		}
		else if(!strcmp(varname[ivar],"streamvort")) // streamwise vorticity
		{
//			clock_t start = clock(), diff;
#pragma omp parallel for private(i)
			for(i=0; i<NX*NY*NZ; i++)
			{
				buffer[i] = (ubuffer[i]*xvort[i]+vbuffer[i]*yvort[i]+wbuffer[i]*zvort[i])/
					sqrt(ubuffer[i]*ubuffer[i]+vbuffer[i]*vbuffer[i]+wbuffer[i]*wbuffer[i]);
			}
//			diff = clock() - start;
//			int msec=diff*1000/CLOCKS_PER_SEC;
		}
		else if(!strcmp(varname[ivar],"streamfrac")) // streamwise vorticity fraction
		{
#pragma omp parallel for private(i)
			for(i=0; i<NX*NY*NZ; i++)
			{
				buffer[i] = (ubuffer[i]*xvort[i]+vbuffer[i]*yvort[i]+wbuffer[i]*zvort[i])/
					(sqrt(ubuffer[i]*ubuffer[i]+vbuffer[i]*vbuffer[i]+wbuffer[i]*wbuffer[i])*
					sqrt(xvort[i]*xvort[i]+yvort[i]*yvort[i]+zvort[i]*zvort[i]));
			}
		}
		else if(!strcmp(varname[ivar],"uinterp"))
		{
			if (!u_rh)read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,varname[ivar],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			else buffer=ubuffer;
		}
		else if(!strcmp(varname[ivar],"vinterp"))
		{
			if (!v_rh)read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,varname[ivar],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			else buffer=vbuffer;
		}
		else if(!strcmp(varname[ivar],"winterp"))
		{
			if (!w_rh)read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,varname[ivar],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			else buffer=wbuffer;
		}
		else if(!strcmp(varname[ivar],"xvort"))
		{
			if (!xvort_rh)read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,varname[ivar],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			else buffer=xvort;
		}
		else if(!strcmp(varname[ivar],"yvort"))
		{
			if (!yvort_rh)read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,varname[ivar],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			else buffer=yvort;
		}
		else if(!strcmp(varname[ivar],"zvort"))
		{
			if (!zvort_rh)read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,varname[ivar],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			else buffer=zvort;
		}
		else if(!strcmp(varname[ivar],"thrhopert"))
		{
			if (!thrhopert_rh)read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,varname[ivar],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			else buffer=thrhopert;
		}
		else // We have (hopefully) requested a variable that has been saved
		{
			read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,varname[ivar],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
		}
		status = nc_put_vara_float (ncid, varnameid[ivar], start, edges, buffer);
		buffer = buf0; /* I will never learn my lesson with pointers........ */
		if (status != NC_NOERR) 
		{
			printf("Could not write variable %s at time %f to %s\n", varname[ivar],t0,ncfilename);
			ERROR_STOP("Write to netcdf file failed");
		}
		printf(")\n");
	}
	status = nc_close(ncid); if (status != NC_NOERR)
	{
		fprintf(stderr, "%s\n", nc_strerror(status));
		printf("status = %i\n",status);
		fprintf(stderr, "Warning: netcdf is throwing an hdf error but our file seems to be fine...\n");
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

void parse_cmdline_makevisit(int argc, char *argv[],
		char *histpath, char *base,
		int *X0, int *Y0, int *X1, int *Y1, int *Z0, int *Z1)
{

	int got_X0,got_X1,got_Y0,got_Y1,got_Z0,got_Z1;
	int got_histpath,got_base;
	enum { OPT_HISTPATH = 1000, OPT_BASE, OPT_X0, OPT_Y0, OPT_X1, OPT_Y1, OPT_Z0, OPT_Z1, OPT_DEBUG };
	static struct option long_options[] =
	{
		{"histpath", required_argument, 0, OPT_HISTPATH},
		{"base",     required_argument, 0, OPT_BASE},
		{"x0",       optional_argument, 0, OPT_X0},
		{"y0",       optional_argument, 0, OPT_Y0},
		{"x1",       optional_argument, 0, OPT_X1},
		{"y1",       optional_argument, 0, OPT_Y1},
		{"z0",       optional_argument, 0, OPT_Z0},
		{"z1",       optional_argument, 0, OPT_Z1},
		{"debug",    optional_argument, 0, OPT_DEBUG},
		{0, 0, 0, 0}//sentinel, needed!
	};

	int bail = 0;
	got_histpath=got_base=got_X0=got_X1=got_Y0=got_Y1=got_Z0=got_Z1=0;

	if (argc == 1)
	{
		fprintf(stderr,
		"Usage: %s --histpath=[histpath] --base=[base] --x0=[X0] --y0=[Y0] --x1=[X1] --y1=[Y1] --z0=[Z0] --z1=[Z1]\n",argv[0]);
		exit(0);
	}

	while (1)
	{
		int r;
		int option_index = 0;
		r = getopt_long_only (argc, argv,"",long_options,&option_index);
//		printf("optarg = %s\n",optarg);
		if (r == -1) break;

		switch(r)
		{
			case OPT_HISTPATH:
				strcpy(histpath,optarg);
				got_histpath=1;
				printf("histpath = %s\n",histpath);
				break;
			case OPT_BASE:
				strcpy(base,optarg);
				got_base=1;
				printf("base = %s\n",base);
				break;
			case OPT_X0:
				*X0 = atoi(optarg);
				got_X0 = 1;
				optcount++;
				printf("X0 = %i\n",*X0);
				break;
			case OPT_Y0:
				*Y0 = atoi(optarg);
				got_Y0 = 1;
				optcount++;
				printf("Y0 = %i\n",*Y0);
				break;
			case OPT_X1:
				*X1 = atoi(optarg);
				got_X1 = 1;
				optcount++;
				printf("X1 = %i\n",*X1);
				break;
			case OPT_Y1:
				*Y1 = atoi(optarg);
				got_Y1 = 1;
				optcount++;
				printf("Y1 = %i\n",*Y1);
				break;
			case OPT_Z0:
				*Z0 = atoi(optarg);
				got_Z0 = 1;
				optcount++;
				printf("Z0 = %i\n",*Z0);
				break;
			case OPT_Z1:
				*Z1 = atoi(optarg);
				got_Z1 = 1;
				optcount++;
				printf("Z1 = %i\n",*Z1);
				break;
			case OPT_DEBUG:
				debug=1;
				optcount++;
				printf("debug = %i\n",debug);
				break;
			case '?':
				fprintf(stderr,"Exiting: unknown command line option.\n");
				exit(0);
				break;
		}
	}

		if (!got_histpath) { fprintf(stderr,"--histpath not specified\n"); bail = 1; }
		if (!got_base)   { fprintf(stderr,"--base not specified\n"); bail = 1; }

/* These are optional */
		if (!got_X0)      fprintf(stderr,"Setting x0 to default value of 0\n");
		if (!got_X1)      fprintf(stderr,"Setting x1 to default value of nx-1\n");
		if (!got_Y0)      fprintf(stderr,"Setting y0 to default value of 0\n");
		if (!got_Y1)      fprintf(stderr,"Setting y1 to default value of ny-1\n");
		if (!got_Z0)      fprintf(stderr,"Setting z0 to default value of 0\n");
		if (!got_Z1)      fprintf(stderr,"Setting z1 to default value of z1-1\n");
		if (bail)           { fprintf(stderr,"Insufficient arguments to %s, exiting.\n",argv[0]); exit(-1); }
}


void	parse_cmdline_hdf2nc(int argc, char *argv[],
	char *histpath, char *base, double *time,
	int *X0, int *Y0, int *X1, int *Y1, int *Z0, int *Z1 )
{
	int got_histpath,got_base,got_time,got_X0,got_X1,got_Y0,got_Y1,got_Z0,got_Z1;
	enum { OPT_HISTPATH = 1000, OPT_BASE, OPT_TIME, OPT_X0, OPT_Y0, OPT_X1, OPT_Y1, OPT_Z0, OPT_Z1,
		OPT_DEBUG, OPT_XYF, OPT_YES2D, OPT_NC3, OPT_COMPRESS, OPT_NTHREADS };
	// see https://stackoverflow.com/questions/23758570/c-getopt-long-only-without-alias
	static struct option long_options[] =
	{
		{"histpath", required_argument, 0, OPT_HISTPATH},
		{"base",     required_argument, 0, OPT_BASE},
		{"time",     required_argument, 0, OPT_TIME},
		{"x0",       optional_argument, 0, OPT_X0},
		{"y0",       optional_argument, 0, OPT_Y0},
		{"x1",       optional_argument, 0, OPT_X1},
		{"y1",       optional_argument, 0, OPT_Y1},
		{"z0",       optional_argument, 0, OPT_Z0},
		{"z1",       optional_argument, 0, OPT_Z1},
		{"debug",    optional_argument, 0, OPT_DEBUG},
		{"xyf",      optional_argument, 0, OPT_XYF},
		{"yes2d",    optional_argument, 0, OPT_YES2D},
		{"nc3",      optional_argument, 0, OPT_NC3},
		{"compress", optional_argument, 0, OPT_COMPRESS},
		{"nthreads", optional_argument, 0, OPT_NTHREADS},
		{0, 0, 0, 0}//sentinel, needed!
	};

	got_histpath=got_base=got_time=got_X0=got_X1=got_Y0=got_Y1=got_Z0=got_Z1=0;

	int bail = 0;

	if (argc == 1)
	{
		fprintf(stderr,
		"Usage: %s --histpath=[histpath] --base=[base] --x0=[X0] --y0=[Y0] --x1=[X1] --y1=[Y1] --z0=[Z0] --z1=[Z1] --time=[time] [varname1 ... varnameN] \n",argv[0]);
		exit(0);
	}

	while (1)
	{
		int r;
		int option_index = 0;
		r = getopt_long_only (argc, argv,"",long_options,&option_index);
		if (r == -1) break;

		switch(r)
		{
			case OPT_HISTPATH:
				strcpy(histpath,optarg);
				got_histpath=1;
				printf("histpath = %s\n",histpath);
				break;
			case OPT_BASE:
				strcpy(base,optarg);
				got_base=1;
				printf("base = %s\n",base);
				break;
			case OPT_TIME:
				*time = atof(optarg);
				got_time=1;
				printf("time = %f\n",*time);
				break;
			case OPT_X0:
				*X0 = atoi(optarg);
				got_X0=1;
				optcount++;
				printf("X0 = %i\n",*X0);
				break;
			case OPT_Y0:
				*Y0 = atoi(optarg);
				got_Y0=1;
				optcount++;
				printf("Y0 = %i\n",*Y0);
				break;
			case OPT_X1:
				*X1 = atoi(optarg);
				got_X1=1;
				optcount++;
				printf("X1 = %i\n",*X1);
				break;
			case OPT_Y1:
				*Y1 = atoi(optarg);
				got_Y1=1;
				optcount++;
				printf("Y1 = %i\n",*Y1);
				break;
			case OPT_Z0:
				*Z0 = atoi(optarg);
				got_Z0=1;
				optcount++;
				printf("Z0 = %i\n",*Z0);
				break;
			case OPT_Z1:
				*Z1 = atoi(optarg);
				got_Z1=1;
				optcount++;
				printf("Z1 = %i\n",*Z1);
				break;
			case OPT_DEBUG:
				debug=1;
				optcount++;
				break;
			case OPT_XYF:
				saved_staggered_mesh_params=1;
				optcount++;
				break;
			case OPT_YES2D:
				yes2d=1;
				optcount++;
				break;
			case OPT_COMPRESS:
				gzip=1;
				optcount++;
				break;
			case OPT_NC3:
				filetype=NC_64BIT_OFFSET;
				optcount++;
				break;
			case OPT_NTHREADS:
				nthreads=atoi(optarg);
				omp_set_num_threads(nthreads);
				optcount++;
				break;
			case '?':
				fprintf(stderr,"Exiting: unknown command line option.\n");
				exit(0);
				break;
		}
	}

		if (got_histpath==0) { fprintf(stderr,"--histpath not specified\n"); bail = 1; }
		if (got_base==0)   { fprintf(stderr,"--base not specified\n"); bail = 1; }
		if (got_time==0)   { fprintf(stderr,"--time not specified\n"); bail = 1; }

/* These are now optional */
		if (!got_X0)      fprintf(stderr,"Will set X0 to saved_X0\n");
		if (!got_Y0)      fprintf(stderr,"Will set Y0 to saved_Y0\n");
		if (!got_X1)      fprintf(stderr,"Will set X1 to saved_X1\n");
		if (!got_Y1)      fprintf(stderr,"Will set Y1 to saved_Y1\n");
		if (!got_Z0)      fprintf(stderr,"Setting Z0 to default value of 0\n");
		if (!got_Z1)      fprintf(stderr,"Setting Z1 to default value of nz-1\n");

		if (bail)           { fprintf(stderr,"Insufficient arguments to %s, exiting.\n",argv[0]); exit(-1); }
}
