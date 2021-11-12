#include "../include/lofs-read.h"
#include "../include/lofs-dirstruct.h"
#include "../include/lofs-hdf2nc.h"
#include "../include/lofs-limits.h"

//ORF dealing with floating point names in our netcdf files is a pain
//smalleps needs to be set carefully, it's only for making netcdf file names
#define smalleps ((double)1.0e-3)

int main(int argc, char *argv[])
{
	int i,status;

	dir_meta dm;
	hdf_meta hm;
	grid gd;
	cmdline cmd;
	ncstruct nc;
	mesh msh;
	sounding snd;
	readahead rh;
	buffers b;

	hid_t hdf_file_id;

	/* begin */

	cmd.argc_hdf2nc_min = 3; /* minimum arguments to this routine */

	init_structs(&cmd,&dm,&gd,&nc,&rh);

	parse_cmdline_hdf2nc(argc,argv,&cmd,&dm,&gd);
	
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

	printf("Simulation identifier (saved_base) = %s\n",dm.saved_base);

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

	//ORF 2021-03-26 we now also collect have ZFP accuracy attributes for all 3D LOFS vars
	get_hdf_metadata(dm,&hm,&cmd,&nc,argv,&hdf_file_id);

	for (i = 0; i < cmd.nvar_cmdline; i++)
	{
		printf("Index %i var3dvarname=%s zfp_accuracy_LOFS = %f\n",i , nc.var3d[i].varname, nc.var3d[i].zfpacc_LOFS);
	}

	printf("3D variables available: ");
	for (i = 0; i < hm.nvar_available; i++) printf("%s ",hm.varname_available[i]);
	printf("\n");

	if(cmd.verbose&&cmd.nvar_cmdline > 0)
	{
		printf("We are requesting the following variables: ");
		for (i=0; i<cmd.nvar_cmdline; i++) printf("%s ",cmd.varname_cmdline[i]);
		printf("\n");
	}

	if (cmd.debug)
		printf("nx = %i ny = %i nz = %i rankx = %i rankx = %i\n",hm.nx,hm.ny,hm.nz,hm.rankx,hm.ranky);

	/* Check for idiocy and tweak the span (X0-X1/Y0-Y1/Z0-Z1) as necessary */

	set_span(&gd,hm,cmd);

	/* Here is where hdf2nc began in old code */

	/* Set base if not at command line and create netcdf file name */

	if (!cmd.got_base) strcpy(cmd.base,dm.saved_base);

 /* ORF 2021-03-02 We can now specify a directory for the netCDF files.
 By default it is the directory from which you execute hdf2nc. If the
 ncdir directory does not exist it will be created. This is optional but
 will be useful for things like ensembles where we are chewing through
 many simulations */

//ORF 2021-11-10
//Should handle this better.
//WE don't want file names with stuff like file.1234.1999998.nc we want file.1234.02000034 or something
//This comes from our subsecond time steps which can be 0.200000000000002134231, 0.333333333333304975, etc...
//(you get the picture)


//  TODO NOW: Change netcdf file name to be in centiseconds (hundredths of a second)

/* We now pass --centiseconds flag to hdf2nc
 * This should be the data save time step in centiseconds (integer value from 0 to 99)
 * Will help us contsruct less weird file names, in
 * cases where we have sequential subsecond saves.
 *
 * if --centiseconds is not passed to the command line it will be set to
 * zero and the code will assume integer second time steps.
 *
 * Regardless the netcdf file times are in centiseconds now and forever amen.
 */

	{
		int itime,ifrac;
		int cs;
		itime = (int)(cmd.time+smalleps);
		ifrac = (int)(100.0*(cmd.time+smalleps-itime));
		cs=cmd.centiseconds;
		printf("cmd.time = %f,ifrac = %i cs = %i\n",cmd.time,ifrac,cs);
//		exit(0);

		if (cmd.got_ncdir)
		{
			DIR* dir; int stat;
			dir = opendir(cmd.ncdir);
			if(dir) // Already exists, do nothing
			{
				closedir(dir);
			}
			else
			{
				//Below fails for mkdir -p type request. Using mkdir_p
				//routine originating from stack overflow
				//stat = mkdir(cmd.ncdir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
				stat = mkdir_p(cmd.ncdir);
				if(stat==-1)
				{
					fprintf(stderr,"%s: Cannot create directory\n",cmd.ncdir);
					ERROR_STOP("mkdir failed");
				}
			}
//			sprintf(nc.ncfilename,"%s/%s.%012.2f.nc",cmd.ncdir,cmd.base,cmd.time);
			sprintf(nc.ncfilename,"%s/%s.%06i%02i.nc",cmd.ncdir,cmd.base,itime,ifrac);
		}
		else
		{
//			sprintf(nc.ncfilename,"%s.%012.2f.nc",cmd.base,cmd.time);
			sprintf(nc.ncfilename,"%s.%06i%02i.nc",cmd.base,itime,ifrac);
		}
		printf("nc.ncfilename = %s\n",nc.ncfilename);//exit(0);
	}

// ORF 2021-11-12 smalleps above keeps our floating point file names from
// having lots of 9s (being a tad 'too small')

	//ORF 2021-07-16
	//ZFP needs chunk dimensions evenly divisible by four
	//As well as our horizontal dimensions divisible by four
	printf("Original: gd.X0=%5i gd.X1=%5i gd.NX=%5i\n",gd.X0,gd.X1,gd.X1-gd.X0+1);
	printf("Original: gd.Y0=%5i gd.Y1=%5i gd.NY=%5i\n",gd.Y0,gd.Y1,gd.Y1-gd.Y0+1);
	printf("Original: gd.Z0=%5i gd.Z1=%5i gd.NZ=%5i\n",gd.Z0,gd.Z1,gd.Z1-gd.Z0+1);
	if(cmd.zfp)
	{
		int x1a,y1a,z1a;
		x1a=gd.X1;
		y1a=gd.Y1;
		z1a=gd.Z1;
		while((gd.X1-gd.X0+1)%4!=0) gd.X1--; 
		while((gd.Y1-gd.Y0+1)%4!=0) gd.Y1--; 
		while((gd.Z1-gd.Z0+1)%4!=0) gd.Z1--; 
		if(x1a-gd.X1 !=0) printf("Adjusted for ZFP writes: gd.X0=%5i gd.X1=%5i gd.NX=%5i\n",gd.X0,gd.X1,gd.X1-gd.X0+1);
		if(y1a-gd.Y1 !=0) printf("Adjusted for ZFP writes: gd.Y0=%5i gd.Y1=%5i gd.NY=%5i\n",gd.Y0,gd.Y1,gd.Y1-gd.Y0+1);
		if(z1a-gd.Z1 !=0) printf("Adjusted for ZFP writes: gd.Z0=%5i gd.Z1=%5i gd.NZ=%5i\n",gd.Z0,gd.Z1,gd.Z1-gd.Z0+1);
	}

	gd.NX = gd.X1 - gd.X0 + 1;
	gd.NY = gd.Y1 - gd.Y0 + 1;
	gd.NZ = gd.Z1 - gd.Z0 + 1;

	/* Allocate memory for 1d mesh and sounding arrays */
	allocate_1d_arrays(hm, gd, &msh, &snd);

	set_1d_arrays(hm,gd,&msh,&snd,&hdf_file_id);

	status = nc_create (nc.ncfilename, NC_CLOBBER|cmd.filetype, &nc.ncid);
	if (status != NC_NOERR) ERROR_STOP ("nc_create failed");

	set_netcdf_attributes(&nc,gd,&cmd,&b,&hm,&hdf_file_id);

	status = nc_enddef (nc.ncid);
	if (status != NC_NOERR) ERROR_STOP("nc_enddef failed");

	nc_write_1d_data(nc,gd,msh,snd,cmd);

	set_readahead(&rh,nc,cmd);

	malloc_3D_arrays(&b,gd,rh,cmd);

	H5Z_zfp_initialize();

	if (cmd.do_swaths) do_the_swaths(hm,nc,dm,gd,cmd);

	do_readahead(&b,gd,rh,dm,hm,cmd);

	do_requested_variables(&b,nc,gd,msh,&snd,rh,dm,hm,cmd);

	status = nc_close(nc.ncid);  if (status != NC_NOERR)
	{
		printf("******nc.ncid = %i\n",nc.ncid);
		fprintf(stderr, "%s\n", nc_strerror(status));
		printf("status = %i\n",status);
		fprintf(stderr, "Warning: netcdf is throwing an error when we close...\n");
	}

	H5Z_zfp_finalize();

	if(cmd.write_cmd_file)write_hdf2nc_command_txtfile(argc,argv,nc);

	//ORF FREE ALL 3D ARRAYS HERE

	free_3D_arrays(&b,gd,rh,cmd);

	if(cmd.gzip) compress_with_nccopy(nc,cmd);

	exit(0);

}
