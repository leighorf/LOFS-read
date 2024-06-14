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

	parse_cmdline_nukefiles(argc,argv,&cmd,&dm,&gd,&zfpacc);
	
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
	nuke_lofs_files(dm,hm,rc,cmd);

	exit(0);
}
