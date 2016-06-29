/* Test reader code for the new .cm1visit file. This will be tested and
 * then stuck in avtcm1visitFileFormat.C */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "hdforf.h"
#include "errorf.h"

#define TRUE 1
#define FALSE 0
#define MAXVARS 100

int main(int argc, char *argv[])
{
	int i,idum=0,nnodedirs;
	hid_t file_id, group_id;
	H5G_info_t group_info;
	char groupname[512];
	char base[512];
	hid_t f_id,g_id,dset_id,dspace_id,strtype;
	hsize_t dims[1];
	char *topdir;
	char **timedir,**nodedir;
	int *dirtimes;
	int ntimedirs;
	int ntottimes;
	int dn;
	int *alltimes;
	int nx,ny,nz,nodex,nodey;
	float *xhfull,*yhfull,*zh;
	int nvars;
	char varname[100][40]; // 100 vars 40 chars long

	if ((file_id = H5Fopen (argv[1], H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
	{
		fprintf(stderr,"Cannot open %s, something is seriously amiss!\n",argv[1]);
		ERROR_STOP("Cannot open hdf file");
	}

	topdir = (char *)malloc(512*sizeof(char));


	H5LTread_dataset_string(file_id,"/topdir",topdir); printf("topdir = %s\n",topdir);
	H5LTread_dataset_int(file_id,"/ntimedirs",&ntimedirs); printf("ntimedirs = %i\n",ntimedirs);
	dirtimes = (int *)malloc(ntimedirs * sizeof(int));
	timedir = (char **)malloc(ntimedirs * sizeof(char *)); for (i=0; i < ntimedirs; i++) timedir[i] = (char *)(malloc(256 * sizeof(char)));
	strtype=H5Tcopy(H5T_C_S1); H5Tset_size(strtype,H5T_VARIABLE);
	H5LTread_dataset(file_id,"/timedir",strtype,timedir); for (i=0; i < ntimedirs; i++) printf("%s ",timedir[i]);printf("\n");
	H5LTread_dataset_int(file_id,"/nnodedirs",&nnodedirs); printf("nnodedirs = %i\n",nnodedirs);
	nodedir = (char **)malloc(nnodedirs * sizeof(char *)); for (i=0; i < nnodedirs; i++) nodedir[i] = (char *)(malloc(256 * sizeof(char)));
	strtype=H5Tcopy(H5T_C_S1); H5Tset_size(strtype,H5T_VARIABLE);
	H5LTread_dataset(file_id,"/nodedir",strtype,nodedir); for (i=0; i < nnodedirs; i++) printf("%s ",nodedir[i]);printf("\n");
	H5LTread_dataset_int(file_id,"/dn",&dn); printf("dn = %i\n",dn);
	H5LTread_dataset_int(file_id,"/ntottimes",&ntottimes); printf("ntottimes = %i\n",ntottimes);
	alltimes = (int *)malloc(ntottimes * sizeof(int));
	H5LTread_dataset_int(file_id,"/alltimes",alltimes); for (i=0; i < ntottimes; i++) printf("%i ",alltimes[i]);printf("\n");
	H5LTread_dataset_int(file_id,"/nx",&nx); printf("nx = %i\n",nx);
	H5LTread_dataset_int(file_id,"/ny",&ny); printf("ny = %i\n",ny);
	H5LTread_dataset_int(file_id,"/nz",&nz); printf("nz = %i\n",nz);
	H5LTread_dataset_int(file_id,"/nodex",&nodex); printf("nodex = %i\n",nodex);
	H5LTread_dataset_int(file_id,"/nodey",&nodey); printf("nodey = %i\n",nodey);
	xhfull = (float *)malloc(nx * sizeof(float)); yhfull = (float *)malloc(ny * sizeof(float)); zh = (float *)malloc(nz * sizeof(float));
	H5LTread_dataset_float(file_id,"/xhfull",xhfull); for (i=0; i < nx; i++) printf("%f ",xhfull[i]);printf("\n");
	H5LTread_dataset_float(file_id,"/yhfull",yhfull); for (i=0; i < ny; i++) printf("%f ",yhfull[i]);printf("\n");
	H5LTread_dataset_float(file_id,"/zh",zh); for (i=0; i < nz; i++) printf("%f ",zh[i]);printf("\n");
	H5LTread_dataset_int(file_id,"/nvars",&nvars); printf("nvars = %i\n",nvars);
	strtype=H5Tcopy(H5T_C_S1); H5Tset_size(strtype,40);
	H5LTread_dataset (file_id, "/varname", strtype, *varname);
	for(i=0;i<nvars;i++)printf("%s ",varname[i]);printf("\n");



}
