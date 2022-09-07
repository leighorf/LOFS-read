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

	parse_cmdline_hdf2nc(argc,argv,&cmd,&dm,&gd,&zfpacc);
	
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
	get_hdf_metadata(dm,&hm,&cmd,&nc,argv,&hdf_file_id,&zfpacc);

	//ORF 2022-08-23 Here is where we can compare saved LOFS zfp accuracy parameter to the
	//ones we chose to write for the netcdf file. Because we uncompress and recompress we
	//must make sure we do not recompress with more accuracy than what was saved if we are
	//just passing LOFS variables to the netCDF file

	if(cmd.verbose)
	for (i = 0; i < cmd.nvar_cmdline; i++)
	{
		printf("Index %i var3dvarname=%s zfp_accuracy_LOFS = %f\n",i , nc.var3d[i].varname, nc.var3d[i].zfpacc_LOFS);
	}

	printf("3D variables available: ");
	for (i = 0; i < hm.nvar_available; i++) printf("%s ",hm.varname_available[i]);
	printf("\n");

//STOPPED HERE 2022-08-24
//Had to regress, I sprintf'ed the ZFP strings so I could write them as
//global metadata to the netCDF files, but some bad shit started
//happening. See saved diffs for what I was doing.
	nzfpacc_LOFS=list_LOFS_zfpacc(&hm,&hdf_file_id); //List all of the ZFP accuracy values for each available LOFS var
	for (i = 0; i < nzfpacc_LOFS; i++)
	{
		printf(hm.zfpacc_LOFS_all[i]);
	}


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

//  Change netcdf file name to be in centiseconds (hundredths of a second)

/* We now can pass --centiseconds flag to hdf2nc
 * This should be the data save time step in centiseconds (integer value from 0 to 99)
 * Will help us contsruct less weird file names, in
 * cases where we have sequential subsecond saves.
 *
 * if --centiseconds is not passed to the command line and you have subsecond data it may
 * or may not work.
*/

/*
 * ORF 2022-08-19 actually the above is not true. cs is not used in this
 * code, see below. However with the 0.2 second ER10 run, everything
 * just works. If you specify the time at the command line with enough
 * accuracy/precision it works. Will revisit this later if we get to
 * weird subsecond time steps (like 1/3, 1/7 etc.)
 *
 * Regardless the netcdf file times are in centiseconds now and forever
 * amen.
 */

	{
		int itime,ifrac;
		int cs;
		itime = (int)(cmd.time+smalleps);
		ifrac = (int)(100.0*(cmd.time+smalleps-itime));
		cs=cmd.centiseconds;
		printf("cmd.time = %f,ifrac = %i cs = %i\n",cmd.time,ifrac,cs);

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
	if(cmd.verbose)
	{
		printf("Original: gd.X0=%5i gd.X1=%5i gd.NX=%5i\n",gd.X0,gd.X1,gd.X1-gd.X0+1);
		printf("Original: gd.Y0=%5i gd.Y1=%5i gd.NY=%5i\n",gd.Y0,gd.Y1,gd.Y1-gd.Y0+1);
		printf("Original: gd.Z0=%5i gd.Z1=%5i gd.NZ=%5i\n",gd.Z0,gd.Z1,gd.Z1-gd.Z0+1);
	}
	if(cmd.zfp)
	{
		int x1a,y1a,z1a;
		x1a=gd.X1;
		y1a=gd.Y1;
		z1a=gd.Z1;
		while((gd.X1-gd.X0+1)%4!=0) gd.X1--; 
		while((gd.Y1-gd.Y0+1)%4!=0) gd.Y1--; 
		while((gd.Z1-gd.Z0+1)%4!=0) gd.Z1--; 
		if(x1a-gd.X1 !=0) printf("***Adjusted for ZFP writes (nx%%4=0): gd.X0=%5i gd.X1=%5i gd.NX=%5i\n",gd.X0,gd.X1,gd.X1-gd.X0+1);
		if(y1a-gd.Y1 !=0) printf("***Adjusted for ZFP writes (ny%%4=0): gd.Y0=%5i gd.Y1=%5i gd.NY=%5i\n",gd.Y0,gd.Y1,gd.Y1-gd.Y0+1);
		if(z1a-gd.Z1 !=0) printf("***Adjusted for ZFP writes (nz%%4=0): gd.Z0=%5i gd.Z1=%5i gd.NZ=%5i\n",gd.Z0,gd.Z1,gd.Z1-gd.Z0+1);
	}

	//ORF 2022-08-23-TODO: If you do not pass any x0,y0,x1,y1,z0,z1 location data
	//at the command line, for calculated variables, there may be zeroes
	//along the border (try calculating tempC for instance). Some check
	//isn't being made I think. See: rc vs gd in the verbose output...
	//It's OK to read in any data, but writing it out we still need to
	//be in multiples of 4 in x,y,z. If you make sure --x0, --x1 etc.
	//line up (divisible by 4) then no weird border issue.

	gd.NX = gd.X1 - gd.X0 + 1;
	gd.NY = gd.Y1 - gd.Y0 + 1;
	gd.NZ = gd.Z1 - gd.Z0 + 1;

	/* Allocate memory for 1d mesh and sounding arrays */
	allocate_1d_arrays(hm, gd, &msh, &snd);

	set_1d_arrays(hm,gd,&msh,&snd,&hdf_file_id);

	status = nc_create (nc.ncfilename, NC_CLOBBER|cmd.filetype, &nc.ncid);
	if (status != NC_NOERR) ERROR_STOP ("nc_create failed");

	set_netcdf_attributes(&nc,gd,&cmd,&b,&hm,&hdf_file_id,&zfpacc);

	status = nc_enddef (nc.ncid);
	if (status != NC_NOERR) ERROR_STOP("nc_enddef failed");

	nc_write_1d_data(nc,gd,msh,snd,cmd);

	set_readahead(&rh,nc,cmd);

	malloc_3D_arrays(&b,gd,rh,cmd);

	H5Z_zfp_initialize();

	if (cmd.do_swaths) do_the_swaths(hm,nc,dm,gd,cmd);

	do_readahead(&b,gd,rh,dm,hm,cmd);

	do_requested_variables(&b,nc,gd,msh,&snd,rh,dm,hm,cmd);

	status = nc_close(nc.ncid);

	if (status != NC_NOERR)
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
