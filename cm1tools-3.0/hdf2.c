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
 *
 * 4/30/19: Added --swaths flag to convert all the swaths (now saved in
 * 3D files by CM1/LOFS)
 *
 * 5/1/10: Added --allvars flag to convert all the 3D fields
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
int nnodedirs;
double *alltimes;
int ntottimes;
int firsttimedirindex;
int saved_X0,saved_Y0,saved_X1,saved_Y1;
float umove = 0.0, vmove = 0.0; /* Need to save these in the history files dammit! */
//const float MISSING=1.0E37;
const float MISSING=0.0; //Ugh deal with these later

int debug = 0;
int do_swaths = 0;
int do_allvars = 0;
int gzip = 0;
int use_interp = 0;
int use_box_offset = 0;
int filetype = NC_NETCDF4;
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
		double *time, int *X0, int *Y0, int *X1, int *Y1, int *Z0, int *Z1, int *regenerate_cache);

void parse_cmdline_makevisit(int argc, char *argv[],
		char *histpath, char *base,
		int *X0, int *Y0, int *X1, int *Y1, int *Z0, int *Z1);

extern char *optarg; /* This is handled by the getopt code */

int n2d,i2d;
const char **twodvarname;
int *twodvarid;
int ncid;
int d2[3];

// For the swath stuff, we iterate through all appropriate groups rather
// than select them one by one, these two routines enable us to do that
// along with the H5Giterate function

herr_t twod_first_pass_hdf2nc(hid_t loc_id, const char *name, void *opdata)
{
    H5G_stat_t statbuf;
    H5Gget_objinfo(loc_id, name, FALSE, &statbuf);
    n2d++;
    return 0;
}

herr_t twod_second_pass_hdf2nc(hid_t loc_id, const char *name, void *opdata)
{
    H5G_stat_t statbuf;
    H5Gget_objinfo(loc_id, name, FALSE, &statbuf);
    strcpy((char *)twodvarname[n2d],name);
    nc_def_var (ncid, twodvarname[n2d], NC_FLOAT, 3, d2, &(twodvarid[n2d]));
    n2d++;
    return 0;
}

// These are defined in parsedir.c
extern int cmpstringp (const void *p1, const void *p2);
extern void sortchararray (char **strarray, int nel);

int main(int argc, char *argv[])
{
	char base[MAXSTR];
	char *cptr;
	char progname[MAXSTR];
	char histpath[MAXSTR];
	int we_are_hdf2nc = FALSE;
	int we_are_makevisit = FALSE;
//	int we_are_hdf2v5d = FALSE;
//	int we_are_linkfiles = FALSE;
	int X0,Y0,X1,Y1,Z0,Z1;
	int i2d;
	double time;
	int regenerate_cache = 0; //We pass this to our parsedir routines in the case that we want to force cache file refreshing for whatever reason

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
	if      (we_are_hdf2nc)    parse_cmdline_hdf2nc(argc, argv, histpath, base, &time, &X0, &Y0, &X1, &Y1, &Z0, &Z1, &regenerate_cache );
	else if (we_are_makevisit) parse_cmdline_makevisit(argc, argv, histpath, base, &X0, &Y0, &X1, &Y1, &Z0, &Z1);

	if((cptr=realpath(histpath,topdir))==NULL)ERROR_STOP("realpath failed");
	grok_cm1hdf5_file_structure(base,regenerate_cache);
	get_hdf_metadata(firstfilename,&nx,&ny,&nz,&nodex,&nodey); //NOTE: saved_X0 etc. are set now
	if(debug) printf("DEBUG: nx = %i ny = %i nz = %i nodex = %i nodey = %i\n", nx,ny,nz,nodex,nodey);


//Often I want to subset from already subsetted LOFS data to make a
//netCDF file. Easiest way is to first make a full subsetted netcdf file
//then use ncview to get the i,j indices for the subset (rather than
//having to do the math by hand to find the new indices). This allows
//that with the --offset option

	if(use_box_offset)
	{
		if(X0<0||Y0<0||X1<0||Y1<0)
			ERROR_STOP ("X0,Y0,X1,Y1 must be specified at command line with --offset option\n");
		X0+=saved_X0;
		X1+=saved_X0;
		Y0+=saved_Y0;
		Y1+=saved_Y0;
	}

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
		ERROR_STOP("Your requested indices are wack, or you have missing cm1hdf5 files, goodbye!\n");
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

void grok_cm1hdf5_file_structure(char *base,int regenerate_cache)
{
	int i;
	ntimedirs = get_num_time_dirs(topdir,debug,regenerate_cache); printf("ntimedirs: %i\n",ntimedirs);
	if (ntimedirs == 0) ERROR_STOP("No cm1 hdf5 files found");

	timedir = (char **)malloc(ntimedirs * sizeof(char *));
	for (i=0; i < ntimedirs; i++) timedir[i] = (char *)(malloc(MAXSTR * sizeof(char)));
	dirtimes = (double *)malloc(ntimedirs * sizeof(double));//times are float not int

	get_sorted_time_dirs(topdir,timedir,dirtimes,ntimedirs,debug,regenerate_cache);

	nnodedirs =  get_num_node_dirs(topdir,timedir[0],debug,regenerate_cache);
	nodedir = (char **)malloc(nnodedirs * sizeof(char *));
	// ORF 8 == 7 zero padded node number directory name plus 1 end of string char
	for (i=0; i < nnodedirs; i++) nodedir[i] = (char *)(malloc(8 * sizeof(char)));

	get_sorted_node_dirs(topdir,timedir[0],nodedir,&dn,nnodedirs,debug,regenerate_cache);

	alltimes = get_all_available_times(topdir,timedir,ntimedirs,nodedir,nnodedirs,&ntottimes,firstfilename,&firsttimedirindex,
			&saved_X0,&saved_Y0,&saved_X1,&saved_Y1,debug,regenerate_cache);
	if(debug)
	{
		printf("All available times: ");
		for (i=0; i<ntottimes; i++)printf("%lf ",alltimes[i]);
		printf("\n");
	}
}

// OK being a bit clever here ... fun with macros. This will make the
// code a lot easier to compare to native CM1 Fortran90 code that we are
// copying anyway. I adopt TEM for his tem array, UA for ua etc.

#define BUF(x,y,z) buf0[P3(x,y,z,NX,NY)]
#define TEM(x,y,z) dum0[P3(x,y,z,NX+1,NY+1)]
#define TEM1(x,y,z) dum1[P3(x,y,z,NX+1,NY+1)]
#define UA(x,y,z) ustag[P3(x+1,y+1,z,NX+2,NY+2)]
#define VA(x,y,z) vstag[P3(x+1,y+1,z,NX+2,NY+2)]
#define WA(x,y,z) wstag[P3(x+1,y+1,z,NX+2,NY+2)]

//hdf2ncdammit
void hdf2nc(int argc, char *argv[], char *base, int X0, int Y0, int X1, int Y1, int Z0, int Z1, double t0)
{
	float *buffer,*buf0,*ustag,*vstag,*wstag,*xvort,*yvort,*zvort;
	float *dum0,*dum00,*dum1,*dum10;
	float *twodbuffer,*twodbuf0;
	float *writeptr;
	float *qc,*qi,*qs;
	float *u0,*v0,*pres0,*pi0,*th0,*qv0;
	float *twodfield;
	int u0id,v0id,pres0id,pi0id,th0id,qv0id;
	float *thrhopert;
	double timearray[1];

	int i,j,k,ix,iy,iz,nvar;
	char varname[MAXVARIABLES][MAXSTR];
	char varname_available[MAXVARIABLES][MAXSTR];
	char varname_cmdline[MAXVARIABLES][MAXSTR];
	char **varname_tmp;
	char ncfilename[MAXSTR];
	char groupname[MAXSTR];

	extern int H5Z_zfp_initialize(void);

	int NX,NY,NZ;
	int ni,nj,nk; //For cm1-like code
// 2019-05-01 new, for getting all 3d variables
	hid_t f_id,g_id;
	H5G_info_t group_info;
	int nvar_available;
	int nvar_cmdline;

	int status;
	int nxh_dimid,nyh_dimid,nzh_dimid;
	int nxf_dimid,nyf_dimid,nzf_dimid,time_dimid,timeid;
	int thsfcid,dbzsfcid;
	int x0id,y0id,z0id,x1id,y1id,z1id;
	int xhid,yhid,zhid;
	int xfid,yfid,zfid;
	int varnameid[MAXVARIABLES];
	int dims[4];
	size_t start[4],edges[4];
	size_t s2[3],e2[3];
	int ivar;
	long int bufsize;
	float bufsize_gb;
	float *xhfull,*yhfull,*xffull,*yffull,*zh,*zf;
	float *uh,*uf,*vh,*vf,*mh,*mf;
	float *xhout,*yhout,*zhout,*xfout,*yfout,*zfout;
	FILE *fp;
	char cmdfilename[512];
	//ORF for writing single time in unlimited time dimension/variable
	const size_t timestart = 0;
	const size_t timecount = 1;
	/* readahead flags */
	int u_rh=0,v_rh=0,w_rh=0,xvort_rh=0,yvort_rh=0,zvort_rh=0,thrhopert_rh=0;
	int qc_rh=0,qi_rh=0,qr_rh=0,qs_rh=0,qg_rh=0;
	int twodslice=0;

	float dx,dy,dz,rdx,rdy,rdz; // reproduce CM1 approach
	float dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz;

	float rv = 461.5;
	float rd = 287.04;
	float reps;



	reps = rv/rd;

	/* Since we tack on all the requested variables to the end of the
	 * command line, we have to find out where the 1st variable
	 * argument is. Since we have optional arguments, we keep track of
	 * things and then peel off the file name strings */

	nvar_cmdline = argc - argc_hdf2nc_min - optcount;

	printf("argc = %i, nvar_cmdline = %i, optcount = %i\n",argc,nvar_cmdline,optcount);

	if ((f_id = H5Fopen (firstfilename, H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
	{
		fprintf(stderr,"Cannot open firstfilename which is %s, even though we have alredy opened it!\n",firstfilename);
		ERROR_STOP("Cannot open hdf file");
	}

	// Get a list of all saved variables - and make it easy to just
	// convert all available variables with another command line
	// option. TODO: cache these?

	sprintf(groupname,"%05i/3D",0);//All vars in 3D group 00000 are always available
	g_id = H5Gopen(f_id,groupname,H5P_DEFAULT);
	H5Gget_info(g_id,&group_info);
	nvar_available = group_info.nlinks;
	for (i = 0; i < nvar_available; i++)
	{
	    H5Lget_name_by_idx(g_id,".",H5_INDEX_NAME,H5_ITER_INC,i,varname_available[i],40,H5P_DEFAULT); //40 characters per varname
	}
	H5Gclose(g_id);

	printf("Variables available: ");
	for (i = 0; i < nvar_available; i++) printf("%s ",varname_available[i]);printf("\n");

	printf("\nWe are requesting the following fields: ");
	for (i=0; i<nvar_cmdline; i++)
	{
		strcpy(varname_cmdline[i],argv[i+argc_hdf2nc_min+optcount]);//HERE IS WHERE WE POPULATE VARNAME
		printf("%s ",varname_cmdline[i]);
	}
	printf("\n");
	sprintf(ncfilename,"%s.%012.6f.nc",base,t0);
	
	NX = X1 - X0 + 1;
	NY = Y1 - Y0 + 1;
	NZ = Z1 - Z0 + 1;

	nk=NZ;nj=NY;ni=NX; //For cm1-like code


/* These are standard for on the scalar mesh and requesting 3D data */
	/* Set below in loop not here anymore */

	/*
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	start[3] = 0;
	edges[0] = 1;
	edges[1] = NZ;
	edges[2] = NY;
	edges[3] = NX;
	*/

//For 2D surface slices and swaths/static slices
//See d2 definition and how it is used in iteration function
	s2[0] = 0;
	s2[1] = 0;
	s2[2] = 0;
	e2[0] = 1;
	e2[1] = NY;
	e2[2] = NX;

	xhfull = (float *)malloc(nx * sizeof(float));
	yhfull = (float *)malloc(ny * sizeof(float));
	xffull = (float *)malloc((nx+1) * sizeof(float));
	yffull = (float *)malloc((ny+1) * sizeof(float));
	th0 = (float *)malloc(nz * sizeof(float));
	qv0 = (float *)malloc(nz * sizeof(float));
	u0 = (float *)malloc(nz * sizeof(float));
	v0 = (float *)malloc(nz * sizeof(float));
	pres0 = (float *)malloc(nz * sizeof(float));
	pi0 = (float *)malloc(nz * sizeof(float));
	zh = (float *)malloc((nz+2) * sizeof(float));
	zf = (float *)malloc((nz+2) * sizeof(float));
	get0dfloat (f_id,(char *)"mesh/dx",&dx); rdx=1.0/dx;
	get0dfloat (f_id,(char *)"mesh/dy",&dy); rdy=1.0/dy;
//	well hell we need to store this
//	get0dfloat (f_id,(char *)"mesh/dz",&dz); rdz=1.0/dz;
	dz=dx;rdz=rdx;
	get1dfloat (f_id,(char *)"mesh/xhfull",xhfull,0,nx);
	get1dfloat (f_id,(char *)"mesh/yhfull",yhfull,0,ny);
	get1dfloat (f_id,(char *)"mesh/xffull",xffull,0,nx+1);
	get1dfloat (f_id,(char *)"mesh/yffull",yffull,0,ny+1);
	get1dfloat (f_id,(char *)"mesh/zh",zh,0,nz);
	get1dfloat (f_id,(char *)"mesh/zf",zf,0,nz);
	get1dfloat (f_id,(char *)"basestate/qv0",qv0,0,nz);
	get1dfloat (f_id,(char *)"basestate/th0",th0,0,nz);
	get1dfloat (f_id,(char *)"basestate/u0",u0,0,nz);
	get1dfloat (f_id,(char *)"basestate/v0",v0,0,nz);
	get1dfloat (f_id,(char *)"basestate/pres0",pres0,0,nz);
	get1dfloat (f_id,(char *)"basestate/pi0",pi0,0,nz);

	xhout = (float *)malloc(NX * sizeof(float));
	yhout = (float *)malloc(NY * sizeof(float));
	zhout = (float *)malloc(NZ * sizeof(float));

	xfout = (float *)malloc((NX+1) * sizeof(float));
	yfout = (float *)malloc((NY+1) * sizeof(float));
	zfout = (float *)malloc((NZ+1) * sizeof(float));

// We recreate George's mesh/derivative calculation paradigm even though
// we are usually isotropic. We need to have our code here match what
// CM1 does internally for stretched and isotropic meshes.
//
// Becuase C cannot do have negative array indices (i.e., uh[-1]) like
// F90 can, we have to offset everything to keep the same CM1-like code
// We malloc enough space for the "ghost zones" and then make sure we
// offset by the correct amount on each side. The macros take care of
// the offsetting.

	uh = (float *)malloc((NX+2) * sizeof(float));
	uf = (float *)malloc((NX+2) * sizeof(float));
	vh = (float *)malloc((NY+2) * sizeof(float));
	vf = (float *)malloc((NY+2) * sizeof(float));
	mh = (float *)malloc((NZ+2) * sizeof(float));
	mf = (float *)malloc((NZ+2) * sizeof(float));
#define UH(ix) uh[ix+1]
#define UF(ix) uf[ix+1]
#define VH(iy) vh[iy+1]
#define VF(iy) vf[iy+1]
#define MH(iz) mh[iz+1]
#define MF(iz) mf[iz+1]
	for (ix=X0-1; ix<X1+1; ix++) UH(ix-X0) = dx/(xffull[ix+1]-xffull[ix]);
	for (ix=X0-1; ix<X1+1; ix++) UF(ix-X0) = dx/(xhfull[ix]-xhfull[ix-1]);
	for (iy=Y0-1; iy<Y1+1; iy++) VH(iy-Y0) = dy/(yffull[iy+1]-yffull[iy]);
	for (iy=Y0-1; iy<Y1+1; iy++) VF(iy-Y0) = dy/(yhfull[iy]-yhfull[iy-1]);
	zf[0] = -zf[2]; //param.F... WE NEED TO LOOK VERY CLOSELY AT NEAR SFC CALCS/DERIVATIVES
	for (iz=Z0-1; iz<Z1+1; iz++) MH(iz-Z0) = dz/(zf[iz+1]-zf[iz]);
	for (iz=Z0-1; iz<Z1+1; iz++) MF(iz-Z0) = dz/(zh[iz]-zf[iz-1]);

	for (iz=Z0; iz<=Z1+1; iz++) zfout[iz-Z0] = 0.001*zf[iz]; 
	for (iy=Y0; iy<=Y1+1; iy++) yfout[iy-Y0] = 0.001*yffull[iy];
	for (ix=X0; ix<=X1+1; ix++) xfout[ix-X0] = 0.001*xffull[ix];

	for (iz=Z0; iz<=Z1; iz++)   zhout[iz-Z0] = 0.001*zh[iz];
	for (iy=Y0; iy<=Y1; iy++)   yhout[iy-Y0] = 0.001*yhfull[iy];
	for (ix=X0; ix<=X1; ix++)   xhout[ix-X0] = 0.001*xhfull[ix];

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
	status = nc_def_dim (ncid, "xf", NX+1, &nxf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed"); 
	status = nc_def_dim (ncid, "yf", NY+1, &nyf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (ncid, "zf", NZ+1, &nzf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (ncid, "time", NC_UNLIMITED, &time_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_var (ncid, "xh", NC_FLOAT, 1, &nxh_dimid, &xhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "yh", NC_FLOAT, 1, &nyh_dimid, &yhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "zh", NC_FLOAT, 1, &nzh_dimid, &zhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "xf", NC_FLOAT, 1, &nxf_dimid, &xfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "yf", NC_FLOAT, 1, &nyf_dimid, &yfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "zf", NC_FLOAT, 1, &nzf_dimid, &zfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "time", NC_DOUBLE, 1, &time_dimid, &timeid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_put_att_text(ncid, xhid, "long_name", strlen("x-coordinate in Cartesian system"), "x-coordinate in Cartesian system");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, xhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, xhid, "axis", strlen("X"), "X");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, yhid, "long_name", strlen("y-coordinate in Cartesian system"), "y-coordinate in Cartesian system");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, yhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, yhid, "axis", strlen("Y"), "Y");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, zhid, "long_name", strlen("z-coordinate in Cartesian system"), "z-coordinate in Cartesian system");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, zhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, zhid, "axis", strlen("Z"), "Z");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, xfid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, yfid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, zfid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, timeid, "units", strlen("s"), "s");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, timeid, "axis", strlen("T"), "T");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, timeid, "long_name", strlen("time"), "time");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_def_var (ncid, "X0", NC_INT, 0, dims, &x0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "Y0", NC_INT, 0, dims, &y0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "X1", NC_INT, 0, dims, &x1id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "Y1", NC_INT, 0, dims, &y1id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "Z0", NC_INT, 0, dims, &z0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "Z1", NC_INT, 0, dims, &z1id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_put_att_text(ncid, x0id, "long_name", strlen("westmost grid index from LOFS data"), "westmost grid index from LOFS data");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, x1id, "long_name", strlen("eastmost grid index from LOFS data"), "eastmost grid index from LOFS data");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, y0id, "long_name", strlen("southmost grid index from LOFS data"), "southmost grid index from LOFS data");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, y1id, "long_name", strlen("northmost grid index from LOFS data"), "northmost grid index from LOFS data");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, z0id, "long_name", strlen("bottom grid index from LOFS data"), "bottom grid index from LOFS data");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(ncid, z1id, "long_name", strlen("top grid index from LOFS data"), "top grid index from LOFS data");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");

	d2[0] = time_dimid;
	d2[1] = nyh_dimid;
	d2[2] = nxh_dimid;

	if (do_swaths)
	{
		bufsize = (long) (NX) * (long) (NY) * (long) sizeof(float);
		if ((twodfield = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our 3D variable buffer array");

		n2d = 0;
		H5Giterate(f_id, "/00000/2D/static",NULL,twod_first_pass_hdf2nc,NULL);
		H5Giterate(f_id, "/00000/2D/swath",NULL,twod_first_pass_hdf2nc,NULL);

		twodvarname = (const char **)malloc(n2d*sizeof(char *));
		twodvarid =   (int *)        malloc(n2d*sizeof(int));

		printf("There are %i 2D static/swath fields (the former domain of the 2D files).\n",n2d);
		bufsize = (long) (NX) * (long) (NY) * (long) (n2d) * (long) sizeof(float);
		if ((twodbuf0 = twodbuffer = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our 3D variable buffer array");

		for (i2d=0; i2d<n2d; i2d++)
		{
			twodvarname[i2d] = (char *)malloc(50*sizeof(char)); // 50 characters per variable
		}

		n2d = 0;
		H5Giterate(f_id, "/00000/2D/static",NULL,twod_second_pass_hdf2nc,NULL);
		H5Giterate(f_id, "/00000/2D/swath",NULL,twod_second_pass_hdf2nc,NULL);

		/* And, like magic, we have populated our netcdf id arrays for all the swath slices */
	}


//OK we are going for simple here. We create a single varname
//array that contains all available variables, plus any
//additional requested variables. We then check for duplicates
//and remove them - this alphabetizes the order of the variables,
//however, which could be annoying and/or useful

	if (!do_allvars)
	{
		nvar = nvar_cmdline;
		nvar_available=0; /* liar! */
	}
//With faking nvar_available to zero this loop now sorts/uniqs
//all possible variables lists (just command line, just allvars,
//or a combination of the two)
	if (nvar !=0||do_allvars) 
	{
		int ndupes = 0;
		varname_tmp = (char **)malloc(MAXVARIABLES * sizeof(char *));
		for (i=0; i < MAXVARIABLES; i++) varname_tmp[i] = (char *)(malloc(MAXSTR * sizeof(char)));

		for (i=0; i<nvar_available; i++)
		{
			strcpy(varname_tmp[i],varname_available[i]);
		}
		for (i=0;i<nvar_cmdline; i++)
		{
			strcpy(varname_tmp[i+nvar_available],varname_cmdline[i]);
		}

		sortchararray (varname_tmp,nvar_available+nvar_cmdline);

		strcpy(varname[0],varname_tmp[0]);
		j=1;
		/* Get rid of accidental duplicates */
		for (i=0; i<nvar_available+nvar_cmdline-1; i++)
		{
			if(strcmp(varname_tmp[i],varname_tmp[i+1]))
			{
				strcpy(varname[j],varname_tmp[i+1]);
				j++;
			}
		}
		nvar = j;
		ndupes = nvar_available+nvar_cmdline-nvar;
		if(ndupes!=0)printf("We got rid of %i duplicate requested variables\n",ndupes);
		for (i=0; i < MAXVARIABLES; i++) free(varname_tmp[i]);
		free(varname_tmp);
	}

// This is our main "loop over all requested variable names" loop that
// sets all the metadata shit.

	for (ivar = 0; ivar < nvar; ivar++)
	{
		/* u v and w live on their own mesh (Arakawa C grid)*/

		/* Recommend, however, for making netcdf files, to just
		 * request uinterp vinterp winterp which NOW are cacluated
		 * HERE rather than saved in LOFS */

		/* We are going to preserve u v w with the extra point for
		 * saving u v and w easily while facilitating averaging
		 * which requires 1 fewer point */

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
		else /* scalar grid */
		{
			dims[0] = time_dimid;
			dims[1] = nzh_dimid;
			dims[2] = nyh_dimid;
			dims[3] = nxh_dimid;
		}

// I'm now going to create truly 2D files when we ask for them rather
//than putting them into a 3D container as I have done in the past

		if(X0==X1)
		{
			twodslice = TRUE;
			dims[0] = time_dimid;
			dims[1] = nzh_dimid;
			dims[2] = nyh_dimid;
			start[0] = 0;
			start[1] = 0;
			start[2] = 0;
			edges[0] = 1;
			edges[1] = NZ;
			edges[2] = NY;
			status = nc_def_var (ncid, varname[ivar], NC_FLOAT, 3, dims, &(varnameid[ivar]));
		}
		else if(Y0==Y1)
		{
			twodslice = TRUE;
			dims[0] = time_dimid;
			dims[1] = nzh_dimid;
			dims[2] = nxh_dimid;
			start[0] = 0;
			start[1] = 0;
			start[2] = 0;
			edges[0] = 1;
			edges[1] = NZ;
			edges[2] = NX;
			status = nc_def_var (ncid, varname[ivar], NC_FLOAT, 3, dims, &(varnameid[ivar]));
		}
		else if(Z0==Z1)
		{
			twodslice = TRUE;
			dims[0] = time_dimid;
			dims[1] = nyh_dimid;
			dims[2] = nxh_dimid;
			start[0] = 0;
			start[1] = 0;
			start[2] = 0;
			edges[0] = 1;
			edges[1] = NY;
			edges[2] = NX;
			status = nc_def_var (ncid, varname[ivar], NC_FLOAT, 3, dims, &(varnameid[ivar]));
		}
		else status = nc_def_var (ncid, varname[ivar], NC_FLOAT, 4, dims, &(varnameid[ivar]));

//		printf("dims:  %5i %5i %5i %5i %s\n",dims[0],dims[1],dims[2],dims[3],varname[ivar]);
//		printf("start: %5i %5i %5i %5i %s\n",start[0],start[1],start[2],start[3],varname[ivar]);
//		printf("edges: %5i %5i %5i %5i %s\n",edges[0],edges[1],edges[2],edges[3],varname[ivar]);
//		printf("\n");

		if (status != NC_NOERR) 
		{
			printf ("Cannot nc_def_var for var #%i %s, status = %i, message = %s\n", ivar, varname[ivar],status,nc_strerror(status));
			ERROR_STOP("nc_def_var failed");
		}

// Here is where we just go through the different CM1 variable names
// that George (and I) use, adding variable attributes merrily. C really
// needs a 'case' statement with character data... instead we will have
// to settle with the if-then-else ladder.
//
// We are still before the nc_enddef call, in case you are lost

		if(!strcmp(varname[ivar],"uinterp"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name", strlen("eastward_wind"), "eastward_wind");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("m/s"), "m/s");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"vinterp"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name", strlen("northward_wind"), "northward_wind");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("m/s"), "m/s");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"winterp"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name", strlen("upward_air_velocity"), "upward_air_velocity");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("m/s"), "m/s");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"prespert"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "long_name", strlen("perturbation_pressure_from_base_state"), "perturbation_pressure_from_base_state");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("hPa"), "hPa");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
				
		else if(!strcmp(varname[ivar],"thpert"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "long_name", strlen("perturbation_potential_temperature_from_base_state"),
					"perturbation_potential_temperature_from_base_state");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("K"), "K");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"thrhopert"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "long_name", strlen("perturbation_density_potential_temperature_from_base_state"),
					"perturbation_density_potential_temperature_from_base_state");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("K"), "K");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"xvort"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name", strlen("x_vorticity"), "x_vorticity");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("s^-1"), "s^-1");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"yvort"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name", strlen("y_vorticity"), "y_vorticity");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("s^-1"), "s^-1");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"zvort"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name", strlen("z_vorticity"), "z_vorticity");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("s^-1"), "s^-1");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"vortmag"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name", strlen("vorticity_magnitude"), "vorticity_magnitude");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("s^-1"), "s^-1");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"dbz"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name", strlen("simulated_radar_reflectivity"), "simulated_radar_reflectivity");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("dBZ"), "dBZ");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"rhopert"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name", strlen("perturbation_density_from_base_state"), "perturbation_density_from_base_state");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("kg/m^3"), "kg/m^3");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"khh"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("horizontal_subgrid_eddy_scalar_diffusivity"), "horizontal_subgrid_eddy_scalar_diffusivity");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("m^2/s"), "m^2/s");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"khv"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("vertical_subgrid_eddy_scalar_diffusivity"), "vertical_subgrid_eddy_scalar_diffusivity");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("m^2/s"), "m^2/s");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"kmh"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("horizontal_subgrid_eddy_momentum_viscosity"), "horizontal_subgrid_eddy_momentum_viscosity");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("m^2/s"), "m^2/s");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"kmv"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("vertical_subgrid_eddy_momentum_viscosity"), "vertical_subgrid_eddy_momentum_viscosity");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("m^2/s"), "m^2/s");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}

		else if(!strcmp(varname[ivar],"nci"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("number_concentration_of_cloud_ice_crystals"), "number_concentration_of_cloud_ice_crystals");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("cm^-3"), "cm^-3");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"ncr"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("number_concentration_of_rain_droplets"), "number_concentration_of_rain_droplets");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("m^-3"), "m^-3");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"ncs"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("number_concentration_of_snow_crystals"), "number_concentration_of_snow_crystals");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("m^-3"), "m^-3");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"ncg"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("number_concentration_of_graupel_particles"), "number_concentration_of_graupel_particles");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("m^-3"), "m^-3");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"qvpert"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("perturbation_mixing_ratio_from_base_state"), "perturbation_mixing_ratio_from_base_state");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("g/kg"), "g/kg");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"qc"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("cloud_water_mixing_ratio"), "cloud_water_mixing_ratio");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("g/kg"), "g/kg");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"qi"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("cloud_ice_mixing_ratio"), "cloud_ice_mixing_ratio");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("g/kg"), "g/kg");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"qr"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("rain_mixing_ratio"), "rain_mixing_ratio");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("g/kg"), "g/kg");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"qs"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("snow_mixing_ratio"), "snow_mixing_ratio");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("g/kg"), "g/kg");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"qg"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("graupel_mixing_ratio"), "graupel_mixing_ratio");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("g/kg"), "g/kg");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
// some diagnostics
		else if(!strcmp(varname[ivar],"qcloud"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("sum_of_qc_and_qi"), "sum_of_qc_and_qi");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("g/kg"), "g/kg");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"qprecip"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("sum_of_qr_qs_and_qg"), "sum_of_qr_qs_and_qg");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("g/kg"), "g/kg");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"streamvort"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("storm_relative_streamwise_vorticity_magnitude"), "storm_relative_streamwise_vorticity_magnitude");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("s^-1"), "s^-1");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"streamfrac"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("fraction_of_vorticity_that_is_streamwise"), "fraction_of_vorticity_that_is_streamwise");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("#"), "#");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"hwin_sr"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("storm_relative_horizontal_wind"), "storm_relative_horizontal_wind");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("m/s"), "m/s");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"windmag_sr"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("storm_relative_3D_wind"), "storm_relative_3D_wind");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("m/s"), "m/s");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}
		else if(!strcmp(varname[ivar],"hwin_gr"))
		{
			status = nc_put_att_text(ncid, varnameid[ivar], "standard_name",
					strlen("ground_relative_horizontal_wind"), "ground_relative_horizontal_wind");
					if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
			status = nc_put_att_text(ncid, varnameid[ivar], "units", strlen("m/s"), "m/s");
			if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
		}

//  UGH FINALLY OVER

//  ONLY DO THIS if I actually have missing values, it really degrades
//  performance of Vapor (and pehraps other software)
//
//  status =  nc_put_att_float(ncid,varnameid[ivar],"missing_value",NC_FLOAT,1,&MISSING);
//
//  Unfortunately this really slows things down. WE NEED
//  ZFP HERE DAMMIT (although that would mean uncompressing and
//  recompressing ZFP data)


//  gzip, lowest level, is good for getting something like 2:1 compression with normal use
//  For large files, you will feel the cost of doing the compression
		if (gzip) status=nc_def_var_deflate(ncid, varnameid[ivar], 1, 1, 1);
	}

	//end of looping over variables before nc_enddef

// But wait! There's more! Let's actually do the world a favor and save
// the sounding data into the netCDF files mmmmkay???

	status = nc_def_var (ncid, "u0", NC_FLOAT, 1, &nzh_dimid, &u0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "v0", NC_FLOAT, 1, &nzh_dimid, &v0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "th0", NC_FLOAT, 1, &nzh_dimid, &th0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "pres0", NC_FLOAT, 1, &nzh_dimid, &pres0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "pi0", NC_FLOAT, 1, &nzh_dimid, &pi0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (ncid, "qv0", NC_FLOAT, 1, &nzh_dimid, &qv0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");

	status = nc_enddef (ncid); if (status != NC_NOERR) ERROR_STOP("nc_enddef failed");

// WE HAVE ENDED OUR FREAKING DEFINITIONS
// NOW LET'S DO A BUNCH OF ACTUAL CRAP

      status = nc_put_var_float (ncid,xhid,xhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_float (ncid,yhid,yhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_float (ncid,zhid,zhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (ncid,xfid,xfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (ncid,yfid,yfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (ncid,zfid,zfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (ncid,u0id,u0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (ncid,v0id,v0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (ncid,th0id,th0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (ncid,pres0id,pres0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (ncid,pi0id,pi0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (ncid,qv0id,qv0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
      status = nc_put_var_int (ncid,x0id,&X0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,y0id,&Y0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,x1id,&X1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,y1id,&Y1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,z0id,&Z0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
      status = nc_put_var_int (ncid,z1id,&Z1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	timearray[0] = t0;
      status = nc_put_vara_double (ncid,timeid,&timestart,&timecount,timearray); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");

	/* So we don't read data more than once, flag stuff that we need
	 * ahead of time (read ahead) */
	for (ivar = 0; ivar < nvar; ivar++)
	{
		if(!strcmp(varname[ivar],"uinterp")&&!use_interp) {u_rh=1;}
		if(!strcmp(varname[ivar],"vinterp")&&!use_interp) {v_rh=1;}
		if(!strcmp(varname[ivar],"winterp")&&!use_interp) {w_rh=1;}
		if(!strcmp(varname[ivar],"hwin_sr")) {u_rh=v_rh=1;}
		if(!strcmp(varname[ivar],"hdiv")) {u_rh=v_rh=1;}
		if(!strcmp(varname[ivar],"windmag_sr")) {u_rh=v_rh=w_rh=1;}
		if(!strcmp(varname[ivar],"zvort")) {u_rh=v_rh=1;}
		if(!strcmp(varname[ivar],"xvort")) {v_rh=w_rh=1;}
		if(!strcmp(varname[ivar],"yvort")) {u_rh=w_rh=1;}
		if(!strcmp(varname[ivar],"vortmag")) {u_rh=v_rh=w_rh=1;}
//		if(!strcmp(varname[ivar],"thrhopert")) {q_liq_solid_rh=1;thpert_rh=1;}
// THESE DO NOT WORK, pop up as they do
		if(!strcmp(varname[ivar],"hwin_gr")) {u_rh=v_rh=1;}
		if(!strcmp(varname[ivar],"rotvortx")) {v_rh=w_rh=xvort_rh=1;}
		if(!strcmp(varname[ivar],"rotvorty")) {u_rh=w_rh=yvort_rh=1;}
		if(!strcmp(varname[ivar],"rotvortz")) {u_rh=v_rh=zvort_rh=1;}
		if(!strcmp(varname[ivar],"rotvortmag")) {u_rh=v_rh=w_rh=xvort_rh=yvort_rh=zvort_rh=1;}
		if(!strcmp(varname[ivar],"roheldens")) {u_rh=v_rh=w_rh=xvort_rh=yvort_rh=zvort_rh=1;}
		if(!strcmp(varname[ivar],"nroheldens")) {u_rh=v_rh=w_rh=xvort_rh=yvort_rh=zvort_rh=1;}
		if(!strcmp(varname[ivar],"vorheldens")) {u_rh=v_rh=w_rh=xvort_rh=yvort_rh=zvort_rh=1;}
		if(!strcmp(varname[ivar],"nvorheldens")) {u_rh=v_rh=w_rh=xvort_rh=yvort_rh=zvort_rh=1;}
		if(!strcmp(varname[ivar],"isquared")) {u_rh=v_rh=w_rh=1;}
		if(!strcmp(varname[ivar],"hvort")) {xvort_rh=yvort_rh=1;}
		if(!strcmp(varname[ivar],"vortmag")) {xvort_rh=yvort_rh=zvort_rh=1;}
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
//Might use these, but I usually only read cloud stuff once, so this
//isn't used... yet
		if(!strcmp(varname[ivar],"qcloud")) {qc_rh=qi_rh=1;}
		if(!strcmp(varname[ivar],"qprecip")) {qr_rh=qs_rh=qg_rh=1;}
	}

// HERE IS WHERE WE DO OUR FREAKING MALLOCS

// ORF 2018-12-04: We only malloc beyond our single 3D buffer when
// requesting on of our diagnostics that require calculations that are a
// function of other saved 3D arrays
//
// Extra point in each dimension is for staggered variables and/or
// diagnostics that require differencing and we therefore assign
// boundary values to missing. Since we are (almost) always far enough
// from the true model lateral boundaries, we would be smarter to just
// read in one extra value in each dimension for ALL variables (except
// for the z=-1!) to serve as a "ghost zone" for the "true" calculation
// array; plus, not everything supports missing value attributes and my
// missing value is freaking huge

// We only malloc here if we are requesting at least one 3D field! */

	if (nvar>0)
	{
		bufsize = (long) (NX+2) * (long) (NY+2) * (long) (NZ+1) * (long) sizeof(float);
		bufsize_gb = 1.0e-9*bufsize;
		printf("\nAllocating %8.5f GB of memory for our main 3D variable array\n",bufsize_gb);
		if ((buf0 = buffer = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our 3D variable buffer array");
		printf("\nAllocating %8.5f GB of memory for our main 3D temporary array\n",bufsize_gb);
		// TODO We should eventually test for only variables that need
		// these extra arrays and only allocate if necessary
		if ((dum00 = dum0 = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our first 3D temp calculation array");
		if ((dum10 = dum1 = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our second 3D temp calculation array");
	}

	//Grab stack o' swaths and blast them home
	if (do_swaths)
	{
		printf("\nWorking on 2D static fields and swaths ("); 

		//Note: nd2 is interpreted as the number of swaths when varname=="swaths"

		// This is supposed to read all of them with n2d 2D slices. Holy crap it actually works. BEEEEERRRRR
		read_hdf_mult_md(twodbuf0,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"swaths",X0,Y0,X1,Y1,0,n2d,nx,ny,nz,nodex,nodey);

		for (i2d=0;i2d<n2d;i2d++)
		{
			for (iy=0; iy<NY; iy++)
				for (ix=0; ix<NX; ix++)
					twodfield[P2(ix,iy,NX)] = twodbuffer[P3(ix,iy,i2d,NX,NY)];
			writeptr = twodfield;
			status = nc_put_vara_float (ncid, twodvarid[i2d], s2, e2, writeptr);
		}
		free(twodbuf0);
		printf(")\n");
	}//So we are all done with swaths.

	if (u_rh)
	{
		printf("\nAllocating %8.5f GB of memory and buffering ustag:\n",bufsize_gb);
		if ((ustag = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate ustag");
		read_hdf_mult_md(ustag,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"u",X0-1,Y0-1,X1+1,Y1+1,Z0,Z1,nx,ny,nz,nodex,nodey);
	}
	if (v_rh)
	{
		printf("\nAllocating %8.5f GB of memory and buffering vstag:\n",bufsize_gb);
		if ((vstag = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate vstag");
// OK we need to be careful, assume we are far enough from lateral
// boundaries here. Cannot combine --allvars with diagnostics if not
// careful.
		read_hdf_mult_md(vstag,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"v",X0-1,Y0-1,X1+1,Y1+1,Z0,Z1,nx,ny,nz,nodex,nodey);
	}
	if (w_rh)
	{
		printf("\nAllocating %8.5f GB of memory and buffering wstag:\n",bufsize_gb);
		if ((wstag = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate wstag");
//		read_hdf_mult_md(wstag,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"winterp",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
		read_hdf_mult_md(wstag,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"w",X0-1,Y0-1,X1+1,Y1+1,Z0,Z1+1,nx,ny,nz,nodex,nodey);
	}

// You could write your own code to calculate new fields based upon the
// fields you are reading in. A bunch of examples follow.

	for (ivar = 0; ivar < nvar; ivar++)
	{
		// Default; only change edges if u, v, w. Note, dims is
		// linked via the variable id, set in the first-pass
		// metadata/before nc_enddef loop, but start and edges are
		// passed to nc_writef, and must be adjustd within this loop
		// for any staggered variables (u, v, w which I call ustag,
		// vstag, wstag in this code, not just u, v, w in cm1)
		//
		// 2019-05-29 Added twodslice boolean for dealing with XY XZ
		// YZ slices. Only works with non-derived variables. If you
		// pass --interp and have uinterp vinterp winterp then you
		// get the desired result. TODO: Make 2D slices work with
		// ustag vstag wstag

		if (!twodslice){
			start[0]=0;start[1]=0;start[2]=0;start[3]=0;
			edges[0]=1;edges[1]=NZ;edges[2]=NY;edges[3]=NX;
		}

		// 2019-05-06: All diagnostic calculations that used
		// [uvw] interp variables with centered differences for
		// derivates have been nuked below.
		//
		// TODO: recalculate these things properly using the
		// staggered velocity variables, with same approach as
		// native CM1

		printf("Working on %s (",varname[ivar]);
		fflush(stdout);

/************************** BEGINNING OF ROT VORT STUFF ************************/
//		if(!strcmp(varname[ivar],"rotvortmag")||!strcmp(varname[ivar],"roheldens")||!strcmp(varname[ivar],"nroheldens"))
//		else if(!strcmp(varname[ivar],"rotvortx")) 
//		else if(!strcmp(varname[ivar],"rotvorty")) 
//		else if(!strcmp(varname[ivar],"rotvortz"))
//		else if(!strcmp(varname[ivar],"isquared")||!strcmp(varname[ivar],"vorheldens")||!strcmp(varname[ivar],"nvorheldens"))
/************************** END OF ROT VORT STUFF ************************/

/*
#define BUF(x,y,z) buf0[P3(x,y,z,NX,NY)]
#define TEM(x,y,z) dum0[P3(x,y,z,NX+1,NY+1)]
#define UA(x,y,z) ustag[P3(x+1,y+1,z,NX+2,NY+2)]
#define VA(x,y,z) vstag[P3(x+1,y+1,z,NX+2,NY+2)]
#define WA(x,y,z) wstag[P3(x+1,y+1,z,NX+2,NY+2)]
*/
		if(!strcmp(varname[ivar],"uinterp")&&!use_interp) //We now calculate this, do not save it any more
		{
#define UINTERP BUF
#pragma omp parallel for private(ix,iy,iz)
			for(iz=0; iz<NZ; iz++)
			for(iy=0; iy<NY; iy++)
			for(ix=0; ix<NX; ix++)
			{
				UINTERP(ix,iy,iz) = 0.5*(UA(ix,iy,iz)+UA(ix+1,iy,iz));
			}
			writeptr = buffer;
		}
		else if(!strcmp(varname[ivar],"vinterp")&&!use_interp) //We now calculate this, do not save it any more
		{
#define VINTERP BUF
#pragma omp parallel for private(ix,iy,iz)
			for(iz=0; iz<NZ; iz++)
			for(iy=0; iy<NY; iy++)
			for(ix=0; ix<NX; ix++)
			{
//				buffer[P3(ix,iy,iz,NX,NY)] = 0.5*(vstag[P3(ix,iy,iz,NX,NY+1)]+vstag[P3(ix,iy+1,iz,NX,NY+1)]);
				VINTERP(ix,iy,iz) = 0.5*(VA(ix,iy,iz)+VA(ix,iy+1,iz));
			}
			writeptr = buffer;
		}
		else if(!strcmp(varname[ivar],"winterp")&&!use_interp) //We now calculate this, do not save it any more
		{
#define WINTERP BUF
#pragma omp parallel for private(ix,iy,iz)
			for(iz=0; iz<NZ; iz++)
			for(iy=0; iy<NY; iy++)
			for(ix=0; ix<NX; ix++)
			{
//				buffer[P3(ix,iy,iz,NX,NY)] = 0.5*(wstag[P3(ix,iy,iz,NX,NY)]+wstag[P3(ix,iy,iz+1,NX,NY)]);
//		bufsize = (long) (NX+2) * (long) (NY+2) * (long) (NZ+1) * (long) sizeof(float);
// #define WA(x,y,z) wstag[P3(x+1,y+1,z,NX+2,NY+2)]
// #define P3(x,y,z,mx,my) (((z)*(mx)*(my))+((y)*(mx))+(x))

				WINTERP(ix,iy,iz) = 0.5*(WA(ix,iy,iz)+WA(ix,iy,iz+1));
			}
			writeptr = buffer;
		}
		else if(!strcmp(varname[ivar],"hwin_sr")) //storm relative horizontal wind speed
		{
			float usr,vsr;
#define HWIN_SR BUF
#pragma omp parallel for private(ix,iy,iz,usr,vsr)
			for(iz=0; iz<NZ; iz++)
			for(iy=0; iy<NY; iy++)
			for(ix=0; ix<NX; ix++)
			{
				usr = 0.5*(UA(ix,iy,iz)+UA(ix+1,iy,iz));
				vsr = 0.5*(VA(ix,iy,iz)+VA(ix,iy+1,iz));
				HWIN_SR(ix,iy,iz) = sqrt(usr*usr+vsr*vsr);
			}
			writeptr = buffer;
		}
		else if(!strcmp(varname[ivar],"windmag_sr")) //storm relative 3D wind speed
		{
			float usr,vsr,wsr; //wsr is dumb but whatevah
#define WINDMAG_SR BUF
#pragma omp parallel for private(ix,iy,iz,usr,vsr,wsr)
			for(iz=0; iz<NZ; iz++)
			for(iy=0; iy<NY; iy++)
			for(ix=0; ix<NX; ix++)
			{
				usr = 0.5*(UA(ix,iy,iz)+UA(ix+1,iy,iz));
				vsr = 0.5*(VA(ix,iy,iz)+VA(ix,iy+1,iz));
				wsr = 0.5*(WA(ix,iy,iz)+WA(ix,iy,iz+1));
				WINDMAG_SR(ix,iy,iz) = sqrt(usr*usr+vsr*vsr+wsr*wsr);
			}
			writeptr = buffer;
		}
//		else if(!strcmp(varname[ivar],"hwin_gr")) //ground relative horizontal wind speed
//		{
//			float usr,vsr;
//			//TODO: put umove and vmove in the got dammed hdf5
//			//files
//			float umove=15.2;
//			float vmove=10.5;
//			These are now command line options until we store
//			these in the cm1hdf5 files
//#pragma omp parallel for private(i,usr,vsr)
//			for(i=0; i<NX*NY*NZ; i++)
//			{
//				usr = ustag[i]+umove;
//				vsr = vstag[i]+vmove;
//				buffer[i] = sqrt(usr*usr+vsr*vsr);
//			}
//			writeptr = buffer;
//		}
		else if(!strcmp(varname[ivar],"hdiv")) // uses staggered velocity variables!
		{
#define HDIV BUF
#pragma omp parallel for private(i,j,k,dudx,dvdy)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
			{
				dudx = (UA(i+1,j,k)-UA(i,j,k))*rdx*UH(i);
				dvdy = (VA(i,j+1,k)-VA(i,j,k))*rdy*VH(j);
				HDIV(i,j,k) = dudx + dvdy;
			}
			writeptr = buffer;
		}
		else if(!strcmp(varname[ivar],"vortmag"))
		{
#pragma omp parallel for private(i,j,k,dwdy,dvdz)
			for(k=1; k<nk; k++)
			for(j=0; j<nj+1; j++)
			for(i=0; i<ni; i++)
			{
				dwdy = (WA(i,j,k)-WA(i,j-1,k))*rdy*VF(i);
				dvdz = (VA(i,j,k)-VA(i,j,k-1))*rdz*MF(k);
				TEM(i,j,k) = dwdy - dvdz;
			}
//This is dependent upon our current free slip bc, see CM1 for other
//decisions
			for(j=0; j<nj+1; j++)
			for(i=0; i<ni; i++)
			{
				TEM(i,j,0)=TEM(i,j,1);
				TEM(i,j,nk)=TEM(i,j,nk-1);
			}

#define XVORT BUF
#pragma omp parallel for private(i,j,k)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				XVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i,j+1,k)+TEM(i,j,k+1)+TEM(i,j+1,k+1));

#pragma omp parallel for private(i,j,k)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				TEM1(i,j,k) = XVORT(i,j,k)*XVORT(i,j,k);


//		}
//		{
#pragma omp parallel for private(i,j,k,dudz,dwdx)
			for(k=1; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni+1; i++)
			{
				dudz = (UA(i,j,k)-UA(i,j,k-1))*rdz*MF(k);
				dwdx = (WA(i,j,k)-WA(i-1,j,k))*rdx*UF(i);
				TEM(i,j,k) = dudz - dwdx;
			}
//This is dependent upon our current free slip bc, see CM1 for other
//decisions
			for(j=0; j<nj; j++)
			for(i=0; i<ni+1; i++)
			{
				TEM(i,j,0)=TEM(i,j,1);
				TEM(i,j,nk)=TEM(i,j,nk-1);
			}
#define YVORT BUF
#pragma omp parallel for private(i,j,k)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				YVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i+1,j,k)+TEM(i,j,k+1)+TEM(i+1,j,k+1));

			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				TEM1(i,j,k) += YVORT(i,j,k)*YVORT(i,j,k);

//		}
//		{
#pragma omp parallel for private(i,j,k,dvdx,dudy)
			for(k=0; k<nk; k++)
			for(j=0; j<nj+1; j++)
			for(i=0; i<ni+1; i++)
			{
				dvdx = (VA(i,j,k)-VA(i-1,j,k))*rdx*UF(i);
				dudy = (UA(i,j,k)-UA(i,j-1,k))*rdy*VF(j);
				TEM(i,j,k) = dvdx - dudy;
			}
#define ZVORT BUF
#pragma omp parallel for private(i,j,k)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				ZVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i+1,j,k)+TEM(i,j+1,k)+TEM(i+1,j+1,k));

#pragma omp parallel for private(i,j,k)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				TEM1(i,j,k) += ZVORT(i,j,k)*ZVORT(i,j,k);

#define VORTMAG BUF
#pragma omp parallel for private(i,j,k)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				VORTMAG(i,j,k) = sqrt(TEM1(i,j,k));

			writeptr=buffer;

		}
		else if(!strcmp(varname[ivar],"xvort"))
		{
#pragma omp parallel for private(i,j,k,dwdy,dvdz)
			for(k=1; k<nk; k++)
			for(j=0; j<nj+1; j++)
			for(i=0; i<ni; i++)
			{
				dwdy = (WA(i,j,k)-WA(i,j-1,k))*rdy*VF(i);
				dvdz = (VA(i,j,k)-VA(i,j,k-1))*rdz*MF(k);
				TEM(i,j,k) = dwdy - dvdz;
			}
//This is dependent upon our current free slip bc, see CM1 for other
//decisions
			for(j=0; j<nj+1; j++)
			for(i=0; i<ni; i++)
			{
				TEM(i,j,0)=TEM(i,j,1);
				TEM(i,j,nk)=TEM(i,j,nk-1);
			}

#define XVORT BUF
#pragma omp parallel for private(i,j,k)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				XVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i,j+1,k)+TEM(i,j,k+1)+TEM(i,j+1,k+1));

			writeptr=buffer;

		}
		else if(!strcmp(varname[ivar],"yvort"))
		{
#pragma omp parallel for private(i,j,k,dudz,dwdx)
			for(k=1; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni+1; i++)
			{
				dudz = (UA(i,j,k)-UA(i,j,k-1))*rdz*MF(k);
				dwdx = (WA(i,j,k)-WA(i-1,j,k))*rdx*UF(i);
				TEM(i,j,k) = dudz - dwdx;
			}
//This is dependent upon our current free slip bc, see CM1 for other
//decisions
			for(j=0; j<nj; j++)
			for(i=0; i<ni+1; i++)
			{
				TEM(i,j,0)=TEM(i,j,1);
				TEM(i,j,nk)=TEM(i,j,nk-1);
			}
#define YVORT BUF
#pragma omp parallel for private(i,j,k)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				YVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i+1,j,k)+TEM(i,j,k+1)+TEM(i+1,j,k+1));

			writeptr=buffer;

		}
		else if(!strcmp(varname[ivar],"zvort")) // uses staggered velocity variables!
		{
#pragma omp parallel for private(i,j,k,dvdx,dudy)
			for(k=0; k<nk; k++)
			for(j=0; j<nj+1; j++)
			for(i=0; i<ni+1; i++)
			{
				dvdx = (VA(i,j,k)-VA(i-1,j,k))*rdx*UF(i);
				dudy = (UA(i,j,k)-UA(i,j-1,k))*rdy*VF(j);
				TEM(i,j,k) = dvdx - dudy;
			}
#define ZVORT BUF
#pragma omp parallel for private(i,j,k)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				ZVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i+1,j,k)+TEM(i,j+1,k)+TEM(i+1,j+1,k));

			writeptr=buffer;

		}
		else if(!strcmp(varname[ivar],"thrhopert")) //We now calcluate here
		{
			read_hdf_mult_md(dum0,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qtot",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			read_hdf_mult_md(dum0,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"thpert",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
#pragma omp parallel for private(i,j,k,dvdx,dudy)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
			{
//do vort first :)
			}

		}
//		else if(!strcmp(varname[ivar],"zvort_stretch")) 
//		else if(!strcmp(varname[ivar],"zvort_tilt")) 
//		else if(!strcmp(varname[ivar],"xvort_baro")) 
//		else if(!strcmp(varname[ivar],"yvort_baro")) 
//		else if(!strcmp(varname[ivar],"yvort_stretch"))
//		else if(!strcmp(varname[ivar],"yvort_tilt"))
//		else if(!strcmp(varname[ivar],"xvort_stretch"))
//		else if(!strcmp(varname[ivar],"xvort_tilt"))
//		else if(!strcmp(varname[ivar],"hvort"))
//		else if(!strcmp(varname[ivar],"vortmag"))
//		else if(!strcmp(varname[ivar],"streamvort"))
//		else if(!strcmp(varname[ivar],"streamfrac"))
		else if(!strcmp(varname[ivar],"hvort"))
		{
#pragma omp parallel for private(i,j,k,dwdy,dvdz)
			for(k=1; k<nk; k++)
			for(j=0; j<nj+1; j++)
			for(i=0; i<ni; i++)
			{
				dwdy = (WA(i,j,k)-WA(i,j-1,k))*rdy*VF(i);
				dvdz = (VA(i,j,k)-VA(i,j,k-1))*rdz*MF(k);
				TEM(i,j,k) = dwdy - dvdz;
			}
//This is dependent upon our current free slip bc, see CM1 for other
//decisions
			for(j=0; j<nj+1; j++)
			for(i=0; i<ni; i++)
			{
				TEM(i,j,0)=TEM(i,j,1);
				TEM(i,j,nk)=TEM(i,j,nk-1);
			}

#define XVORT BUF
#pragma omp parallel for private(i,j,k)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				XVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i,j+1,k)+TEM(i,j,k+1)+TEM(i,j+1,k+1));

#pragma omp parallel for private(i,j,k)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				TEM1(i,j,k) = XVORT(i,j,k)*XVORT(i,j,k);


//		}
//		{
#pragma omp parallel for private(i,j,k,dudz,dwdx)
			for(k=1; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni+1; i++)
			{
				dudz = (UA(i,j,k)-UA(i,j,k-1))*rdz*MF(k);
				dwdx = (WA(i,j,k)-WA(i-1,j,k))*rdx*UF(i);
				TEM(i,j,k) = dudz - dwdx;
			}
//This is dependent upon our current free slip bc, see CM1 for other
//decisions
			for(j=0; j<nj; j++)
			for(i=0; i<ni+1; i++)
			{
				TEM(i,j,0)=TEM(i,j,1);
				TEM(i,j,nk)=TEM(i,j,nk-1);
			}
#define YVORT BUF
#pragma omp parallel for private(i,j,k)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				YVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i+1,j,k)+TEM(i,j,k+1)+TEM(i+1,j,k+1));

			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				TEM1(i,j,k) += YVORT(i,j,k)*YVORT(i,j,k);

#define HVORT BUF
#pragma omp parallel for private(i,j,k)
			for(k=0; k<nk; k++)
			for(j=0; j<nj; j++)
			for(i=0; i<ni; i++)
				HVORT(i,j,k) = sqrt(TEM1(i,j,k));

			writeptr=buffer;

		}

		else if(!strcmp(varname[ivar],"u"))
		{
			if (!u_rh){read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,
					varname[ivar],X0,Y0,X1+1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);writeptr=buffer;}
			else writeptr=ustag;
		}
		else if(!strcmp(varname[ivar],"v"))
		{
			if (!v_rh){read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,
					varname[ivar],X0,Y0,X1,Y1+1,Z0,Z1,nx,ny,nz,nodex,nodey);writeptr=buffer;}
			else writeptr=vstag;
		}
		else if(!strcmp(varname[ivar],"w"))
		{
			if (!w_rh){read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,
					varname[ivar],X0,Y0,X1,Y1,Z0,Z1+1,nx,ny,nz,nodex,nodey);writeptr=buffer;}
			else writeptr=wstag;
		}
//		else if(!strcmp(varname[ivar],"xvort"))
//		else if(!strcmp(varname[ivar],"yvort"))
//		else if(!strcmp(varname[ivar],"zvort"))
//		else if(!strcmp(varname[ivar],"qcloud"))
//		{
//			read_hdf_mult_md(dumarray,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qc",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
//			read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qi",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
//			for(i=0; i<NX*NY*NZ; i++) *buffer++ += *dumarray++;
//			buffer=buf0; dumarray=dum0;
//			writeptr = buffer;
//		}
//		else if(!strcmp(varname[ivar],"qprecip"))
//		{
//			read_hdf_mult_md(dumarray,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qr",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
//			read_hdf_mult_md(buffer,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qs",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
//			for(i=0; i<NX*NY*NZ; i++) *buffer++ += *dumarray++;
//			buffer=buf0; dumarray=dum0;
//			read_hdf_mult_md(dumarray,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"qg",X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
//			for(i=0; i<NX*NY*NZ; i++) *buffer++ += *dumarray++;
//			buffer=buf0; dumarray=dum0;
//			writeptr = buffer;
//		}
//		else if(!strcmp(varname[ivar],"qvstupid"))
		else // We have (hopefully) requested a variable that has been saved
		{
			read_hdf_mult_md(buf0,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,varname[ivar],X0,Y0,X1,Y1,Z0,Z1,nx,ny,nz,nodex,nodey);
			writeptr = buf0;
		}
		if(!strcmp(varname[ivar],"u")) edges[3] = NX+1;
		if(!strcmp(varname[ivar],"v")) edges[2] = NY+1;
		if(!strcmp(varname[ivar],"w")) edges[1] = NZ+1;
dumpit:	status = nc_put_vara_float (ncid, varnameid[ivar], start, edges, writeptr);
		if (status != NC_NOERR) 
		{
			printf("dims:  %5i %5i %5i %5i %s\n",dims[0],dims[1],dims[2],dims[3],varname[ivar]);
			printf("start: %5i %5i %5i %5i %s\n",start[0],start[1],start[2],start[3],varname[ivar]);
			printf("edges: %5i %5i %5i %5i %s\n",edges[0],edges[1],edges[2],edges[3],varname[ivar]);
			printf("\n");

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
	int *X0, int *Y0, int *X1, int *Y1, int *Z0, int *Z1, int *regenerate_cache )
{
	int got_histpath,got_base,got_time,got_X0,got_X1,got_Y0,got_Y1,got_Z0,got_Z1;
	enum { OPT_HISTPATH = 1000, OPT_BASE, OPT_TIME, OPT_X0, OPT_Y0, OPT_X1, OPT_Y1, OPT_Z0, OPT_Z1,
		OPT_DEBUG, OPT_REGENERATECACHE, OPT_ALLVARS, OPT_SWATHS, OPT_NC3, OPT_COMPRESS, OPT_NTHREADS, OPT_UMOVE, OPT_VMOVE, OPT_OFFSET, OPT_INTERP };
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
		{"recache",  optional_argument, 0, OPT_REGENERATECACHE},
		{"allvars",  optional_argument, 0, OPT_ALLVARS},
		{"swaths",   optional_argument, 0, OPT_SWATHS},
		{"nc3",      optional_argument, 0, OPT_NC3},
		{"compress", optional_argument, 0, OPT_COMPRESS},
		{"nthreads", optional_argument, 0, OPT_NTHREADS},
		{"umove", optional_argument, 0, OPT_UMOVE},
		{"vmove", optional_argument, 0, OPT_VMOVE},
		{"offset", optional_argument, 0, OPT_OFFSET},
		{"interp", optional_argument, 0, OPT_INTERP},
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
			case OPT_REGENERATECACHE:
				*regenerate_cache=1;
				optcount++;
				break;
			case OPT_SWATHS:
				do_swaths=1;
				optcount++;
				break;
			case OPT_ALLVARS:
				do_allvars=1;
				optcount++;
				break;
			case OPT_COMPRESS:
				gzip=1;
				optcount++;
				break;
			case OPT_INTERP:
				use_interp=1;
				optcount++;
				break;
			case OPT_OFFSET:
				use_box_offset=1;
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
			case OPT_UMOVE:
				umove=atof(optarg);
				optcount++;
				break;
			case OPT_VMOVE:
				vmove=atof(optarg);
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
//Some comments that I had higher up describing our methodology with calculating things like zvort
//
// 2019-05-06 THIS LOOKS GOOD but should be tested against CM1 zvort!
//
// I have copied George's approach here for doing vorticity. See
// misclibs/calcvort in cm1r16. The "trick" is to read in "ghost
// zones" of ustag and vstag since to calculate vortz on the scalar
// mesh you must first calcluate zeta (using "upwind" nearest-neighbor
// differencing) but those calculations go a bit to the left and a
// bit to the right of where we would prefer. So we have one extra
// point to the west(south) and one extra point to the east(north)
// for ustag(vstag). Because C does not have the "handy" way to index
// negative values like Fortran, we just start with zero as usual such
// that the smallest index is always 0.
//
// Here is the CM1 Fortran version.
//
// Keep in mind the ghost zones are actual ghost zones for the model;
// for us they are "lateral boundaries" for the full cube, and we are
// staying away from the model's true boundaries. So the idea is,
// horizntally, we should always be requesting data with enough padding
// to the left and right so as to not fail. We will eventually probably
// want to do this right someday.

/* CM1 fortran:

      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua

	ngxy=3 (number of ghost zones in x and y)
	ib=1-ngxy
	jb=1-ngxy
	kb=1-ngz

	ie=ni+ngxy
	je=nj+ngxy
	ke=nk+ngz

So, CM1 can safely reference ua(-2:ni+3,-2:nj+2,0:nz+1) whereas I am
only reading in ustag(0:nx+2,0:ny+2,0:0:nz+1) [in C-ese]. So the array
indexing will be different in the C code than with George's. I am
keeping the same loop bounds and just changing the indexing.

        do k=1,nk
          do j=1,nj+1
          do i=1,ni+1
            tem(i,j,k) = (va(i,j,k)-va(i-1,j,k))*rdx*uf(i)   &
                        -(ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)
          enddo
          enddo
          do j=1,nj
          do i=1,ni
            zvort(i,j,k) = 0.25*(tem(i,j,k)+tem(i+1,j,k)+tem(i,j+1,k)+tem(i+1,j+1,k))
          enddo
          enddo
        enddo
*/

//#define P3(x,y,z,mx,my) (((z)*(mx)*(my))+((y)*(mx))+(x))
