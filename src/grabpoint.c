#include "../include/lofs-read.h"
#include "../include/lofs-dirstruct.h"
#include "../include/lofs-hdf2nc.h"
#include "../include/lofs-limits.h"

//ORF dealing with floating point names in our netcdf files is a pain
//smalleps needs to be set carefully, it's only for making netcdf file names
#define smalleps ((double)1.0e-3)

int main(int argc, char *argv[])
{
	int i,nzfpacc_LOFS,status;

	dir_meta dm;
	hdf_meta hm;
	grid gd;
	cmdline cmd;
	ncstruct nc;
	mesh msh;
	sounding snd;
	readahead rh;
	buffers b;
	zfpacc zfpacc;

	hid_t hdf_file_id;

	/* begin */

	cmd.argc_hdf2nc_min = 3; /* minimum arguments to this routine */

	init_structs(&cmd,&dm,&gd,&nc,&rh,&hm,&zfpacc);

	parse_cmdline_grabpoint(argc,argv,&cmd,&dm,&gd,&zfpacc);
	
	cmd.nvar_cmdline = argc - cmd.argc_hdf2nc_min - cmd.optcount;
		
	if((realpath(cmd.histpath,dm.topdir))==NULL)
	{
		fprintf(stderr,"%s: No such directory\n",cmd.histpath);
		ERROR_STOP("realpath failed");
	}

	get_num_time_dirs(&dm,cmd); //Sets dm.ntimedirs

	/* Malloc our time directory arrays */
	dm.timedir = (char **)malloc(dm.ntimedirs * sizeof(char *));
	for (i=0; i < dm.ntimedirs; i++)
		dm.timedir[i] = (char *)(malloc(MAXSTR * sizeof(char)));
	dm.dirtimes = (double *)malloc(dm.ntimedirs * sizeof(double));

	get_sorted_time_dirs(&dm,cmd); //Sets dm.timedir char array
	get_num_node_dirs(&dm,cmd);    //Sets dm.nnodedirs

	dm.nodedir = (char **)malloc(dm.nnodedirs * sizeof(char *));
	for (i=0; i < dm.nnodedirs; i++)
		dm.nodedir[i] = (char *)(malloc(8 * sizeof(char)));
	// ORF 8 is 7 zero padded node number directory name plus 1 end of string char
	// TODO: make these constants/macros

	get_sorted_node_dirs(&dm,cmd); //Sets dm.nodedir char array

	get_saved_base(dm.timedir[0],dm.saved_base);

//	printf("Simulation identifier (saved_base) = %s\n",dm.saved_base);

	get_all_available_times(&dm,&gd,cmd); //Gets all times in one double precision array

	if(cmd.debug)
	{
		printf("All available times: ");
		for (i=0; i<dm.ntottimes; i++)
		{
			printf("%lf ",dm.alltimes[i]);
			printf("\n");
		}
	}

	if ((hdf_file_id = H5Fopen (dm.firstfilename, H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
	{
		fprintf(stderr,"Unable to open %s, bailing!\n", dm.firstfilename);
		ERROR_STOP("Can't open firstfilename! Weird...");
	} // Keep open as we need metadata, 1d sounding data, etc.

	//ORF 2021-03-26 we now also collect ZFP accuracy attributes for all 3D LOFS vars.
	//These are written as global attributes to the netCDF file below.

	get_hdf_metadata(dm,&hm,&cmd,&nc,argv,&hdf_file_id,&zfpacc);

	if(cmd.verbose)
	for (i = 0; i < cmd.nvar_cmdline; i++)
	{
		printf("Index %i var3dvarname=%s zfp_accuracy_LOFS = %f\n",i , nc.var3d[i].varname, nc.var3d[i].zfpacc_LOFS);
	}

//	printf("3D variables available: ");
//	for (i = 0; i < hm.nvar_available; i++) printf("%s ",hm.varname_available[i]);
//	printf("\n");

	if(cmd.verbose&&cmd.nvar_cmdline > 0)
	{
		printf("We are requesting the following variables: ");
		for (i=0; i<cmd.nvar_cmdline; i++) printf("%s ",cmd.varname_cmdline[i]);
		printf("\n");
	}

	if (cmd.debug)
		printf("nx = %i ny = %i nz = %i rankx = %i rankx = %i\n",hm.nx,hm.ny,hm.nz,hm.rankx,hm.ranky);


	/* Allocate memory for 1d mesh and sounding arrays */

	gd.NX=2;gd.NY=2;gd.NZ=2;
	allocate_1d_arrays(hm, gd, &msh, &snd);

	get1dfloat (hdf_file_id,(char *)"mesh/xhfull",msh.xhfull,0,hm.nx);
	get1dfloat (hdf_file_id,(char *)"mesh/yhfull",msh.yhfull,0,hm.ny);
	get1dfloat (hdf_file_id,(char *)"mesh/xffull",msh.xffull,0,hm.nx+1);
	get1dfloat (hdf_file_id,(char *)"mesh/yffull",msh.yffull,0,hm.ny+1);
	get1dfloat (hdf_file_id,(char *)"mesh/zh",msh.zh,0,hm.nz);
	get1dfloat (hdf_file_id,(char *)"mesh/zf",msh.zf,0,hm.nz);
	get0dfloat (hdf_file_id,(char *)"mesh/umove",&msh.umove);
	get0dfloat (hdf_file_id,(char *)"mesh/vmove",&msh.vmove);

//All righty. This ugly code makes beautiful output. Only the 1st
//call to grabpoint should have the --header argument so the
//resultant file is pure 'csv' that can be read by pandas. All of the
//requested LOFS variables will be interpolated and put in columns.
	{
		float *interpval;
		char header[MAXSTR];
		char values[MAXSTR];
		char tmpstr[MAXSTR];
		interpval = (float *)malloc(cmd.nvar_cmdline*sizeof(float));
		for (i=0; i<cmd.nvar_cmdline; i++)
		{
			interpval[i] = grabpoint(&gd,hm,dm,cmd,msh,cmd.varname_cmdline[i]);
		}
		if(cmd.header)
		{
			sprintf(header,"%20s%20s%20s%20s","time","xpos","ypos","zpos");
			for (i=0; i<cmd.nvar_cmdline; i++)
			{
				sprintf(tmpstr,"%20s",cmd.varname_cmdline[i]);
				strcat(header,tmpstr);
			}
			printf("%s\n",header);
		}
		sprintf(values,"%20.7f%20.7f%20.7f%20.7f",cmd.time,gd.XC,gd.YC,gd.ZC);
		for (i=0; i<cmd.nvar_cmdline; i++)
		{
			sprintf(tmpstr,"%20.7f",interpval[i]);
			strcat(values,tmpstr);
		}
		printf("%s\n",values);
	}
	exit(0);
}
