#include <omp.h>
#include "include/dirstruct.h"
#include "include/limits.h"
#include "include/hdf2nc.h"
#include "include/lofs-read.h"

void init_struct1(cmdline *cmd,dir_meta *dm, grid *gd)
{
	cmd->histpath        = (char *) malloc(MAXSTR*sizeof(char));
	cmd->base            = (char *) malloc(MAXSTR*sizeof(char));
	dm-> firstfilename   = (char *) malloc(MAXSTR*sizeof(char));
	dm-> saved_base      = (char *) malloc(MAXSTR*sizeof(char));
	dm-> topdir          = (char *) malloc(MAXSTR*sizeof(char));

	dm-> regenerate_cache = 0;

	gd->X0=gd->Y0=gd->X1=gd->Y1=gd->Z0=gd->Z1=-1;
	gd->saved_X0=gd->saved_X1=0;
	gd->saved_Y0=gd->saved_Y1=0;
	gd->saved_Z0=gd->saved_Z1=0;
	cmd->time=0.0; cmd->got_base=0; cmd->optcount=0;
	cmd->debug=0;
}

void get_saved_base(char *timedir, char *saved_base)
{
	// Just grab the basename so we can have it set automatically
	// for the netcdf files if so desired. Apparenlty the easiest
	// way to do this is to just remove the last 14 characters of
	// one of the timedir directories

	int tdlen;
//ORF TODO: make 14 a constant somewhere
//We may wish to allow flexibility in the floating point
//left and right of the decimal.
	tdlen=strlen(timedir)-14;
	strncpy(saved_base,timedir,tdlen);
	saved_base[tdlen]='\0';
}


int main(int argc, char *argv[])
{
	int i;

	dir_meta dm;
	hdf_meta hm;
	grid gd;
	cmdline cmd;

	init_struct1(&cmd,&dm,&gd);

	parse_cmdline_lofs2nc(argc,argv,&cmd,&dm,&gd);

	if((realpath(cmd.histpath,dm.topdir))==NULL)
	{
		fprintf(stderr,"%s: No such directory\n",cmd.histpath);
		ERROR_STOP("realpath failed");
	}

	get_num_time_dirs(&dm,cmd); //Sets dm.ntimedirs

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

	printf("ORF DEBUG timedir[0] = %s saved_base = %s\n",dm.timedir[0],dm.saved_base);

	get_saved_base(dm.timedir[0],dm.saved_base);

	printf("Simulation identifier (saved_base) = %s\n",dm.saved_base);

	get_all_available_times(&dm,&gd,cmd);

	if(1)
	{
		printf("All available times: ");
		for (i=0; i<dm.ntottimes; i++)printf("%lf ",dm.alltimes[i]);
		printf("\n");
	}

}
