#include "../include/lofs-read.h"
#include "../include/lofs-dirstruct.h"
#include "../include/lofs-hdf2nc.h"
#include "../include/lofs-limits.h"

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

	get_hdf_metadata(dm,&hm,&cmd,argv,&hdf_file_id);

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
		printf("nx = %i ny = %i nz = %i nodex = %i nodey = %i\n",hm.nx,hm.ny,hm.nz,hm.nodex,hm.nodey);

	/* Check for idiocy and tweak the span (X0-X1/Y0-Y1/Z0-Z1) as necessary */

	set_span(&gd,hm,cmd);

	/* Here is where hdf2nc began in old code */

	/* Set base if not at command line and create netcdf file name */

	if (!cmd.got_base) strcpy(cmd.base,dm.saved_base);
	sprintf(nc.ncfilename,"%s.%012.6f.nc",cmd.base,cmd.time);

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

	do_requested_variables(&b,nc,gd,msh,rh,dm,hm,cmd);

	status = nc_close(nc.ncid);  if (status != NC_NOERR)
	{
		printf("******nc.ncid = %i\n",nc.ncid);
		fprintf(stderr, "%s\n", nc_strerror(status));
		printf("status = %i\n",status);
		fprintf(stderr, "Warning: netcdf is throwing an error when we close...\n");
	}

	H5Z_zfp_finalize();

	write_hdf2nc_command_txtfile(argc,argv,nc);

	//ORF FREE ALL 3D ARRAYS HERE

	free_3D_arrays(&b,gd,rh,cmd);

	if(cmd.gzip) compress_with_nccopy(nc,cmd);


	exit(0);

}
