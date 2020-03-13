#include <omp.h>
#include "include/dirstruct.h"
#include "include/limits.h"
#include "include/lofs-read.h"
#include "include/hdf2nc.h"

int main(int argc, char *argv[])
{
	int i;

	dir_meta dm;
	hdf_meta hm;
	grid gd;
	cmdline cmd;

	/* begin */

	cmd.argc_hdf2nc_min = 3; /* minimum arguments to this routine */

	init_structs(&cmd,&dm,&gd);

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
	for (i=0; i < dm.ntimedirs; i++) dm.timedir[i] = (char *)(malloc(MAXSTR * sizeof(char)));
	dm.dirtimes = (double *)malloc(dm.ntimedirs * sizeof(double));

	get_sorted_time_dirs(&dm,cmd); //Sets dm.timedir char array
	get_num_node_dirs(&dm,cmd);    //Sets dm.nnodedirs

	dm.nodedir = (char **)malloc(dm.nnodedirs * sizeof(char *));
	for (i=0; i < dm.nnodedirs; i++) dm.nodedir[i] = (char *)(malloc(8 * sizeof(char)));
	// ORF 8 is 7 zero padded node number directory name plus 1 end of string char
	// TODO: make these constants/macros

	get_sorted_node_dirs(&dm,cmd); //Sets dm.nodedir char array

	get_saved_base(dm.timedir[0],dm.saved_base);

	printf("Simulation identifier (saved_base) = %s\n",dm.saved_base);

	get_all_available_times(&dm,&gd,cmd); //Gets all times in one double precision array

	if(cmd.debug) { printf("All available times: "); for (i=0; i<dm.ntottimes; i++)printf("%lf ",dm.alltimes[i]); printf("\n"); }

	get_hdf_metadata(dm,&hm,&cmd,argv);

	printf("Variables available: "); for (i = 0; i < hm.nvar_available; i++) printf("%s ",hm.varname_available[i]);printf("\n");
	if(cmd.verbose&&cmd.nvar_cmdline > 0)
	{
		printf("We are requesting the following variables: ");
		for (i=0; i<cmd.nvar_cmdline; i++) printf("%s ",cmd.varname_cmdline[i]);
	} printf("\n");

	if (cmd.debug) printf("nx = %i ny = %i nz = %i nodex = %i nodey = %i\n",hm.nx,hm.ny,hm.nz,hm.nodex,hm.nodey);

	/* We have all our metadata (have filled all our structures) */

	/* Check for idiocy and tweak the span (X0-X1/Y0-Y1/Z0-Z1) as necessary */

	set_span(&gd,hm,cmd);

	/* Here is where hdf2nc began in old code */

}
