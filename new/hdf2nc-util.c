#include "include/dirstruct.h"
#include "include/limits.h"
#include "include/hdf2nc.h"
#include "include/lofs-read.h"

void init_structs(cmdline *cmd,dir_meta *dm, grid *gd)
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
	cmd->verbose=0;
	cmd->use_box_offset=0;
	cmd->use_interp=0;
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

void set_span(grid *gd,hdf_meta hm,cmdline cmd)
{
	/* We do this because of stencils */

	gd->saved_X0+=1; gd->saved_X1-=1;
	gd->saved_Y0+=1; gd->saved_Y1-=1;

	/* X and Y values are with respect to (0,0) of the actual saved data, not with respect
	 * to actual (0,0) (done because I often choose new netCDF subdomains from ncview).
	 * --offset allows this
	 */

	if(cmd.use_box_offset)
	{
		if(gd->X0<0||gd->Y0<0||gd->X1<0||gd->Y1<0) ERROR_STOP ("X0,Y0,X1,Y1 must be specified at command line with --offset option\n");
		gd->X0+=gd->saved_X0;
		gd->X1+=gd->saved_X0;
		gd->Y0+=gd->saved_Y0;
		gd->Y1+=gd->saved_Y0;
	}

	/* If we didn't specify values at the command line, set them to values specifying all the saved data */
	if(gd->X0<0)gd->X0=gd->saved_X0; if(gd->Y0<0)gd->Y0=gd->saved_Y0; if(gd->Z0<0)gd->Z0=0;
	if(gd->X1<0)gd->X1=gd->saved_X1; if(gd->Y1<0)gd->Y1=gd->saved_Y1; if(gd->Z1<0)gd->Z1=hm.nz-2; //ORF -2 not -1 now


	if (gd->X0>gd->X1||gd->Y0>gd->Y1||gd->X1<gd->saved_X0||gd->Y1<gd->saved_Y0||gd->X0>gd->saved_X1||gd->Y0>gd->saved_Y1)
	{
		printf(" *** X0=%i saved_X0=%i Y0=%i saved_Y0=%i X1=%i saved_X1=%i Y1=%i saved_Y1=%i\n",
				gd->X0,gd->saved_X0,gd->Y0,gd->saved_Y0,gd->X1,gd->saved_X1,gd->Y1,gd->saved_Y1);
		ERROR_STOP("Your requested indices are wack, or you have missing cm1hdf5 files, goodbye!\n");
	}
	if(gd->X0<gd->saved_X0)
	{
		printf("Oops: requested out of box: Adjusting X0 (%i) to saved_X0 (%i)\n",gd->X0,gd->saved_X0);
		gd->X0=gd->saved_X0;
	}
	if(gd->Y0<gd->saved_Y0)
	{
		printf("Oops: requested out of box: Adjusting Y0 (%i) to saved_Y0 (%i)\n",gd->Y0,gd->saved_Y0);
		gd->Y0=gd->saved_Y0;
	}
	if(gd->X1>gd->saved_X1)
	{
		printf("Oops: requested out of box: Adjusting X1 (%i) to saved_X1 (%i)\n",gd->X1,gd->saved_X1);
		gd->X1=gd->saved_X1;
	}
	if(gd->Y1>gd->saved_Y1)
	{
		printf("Oops: requested out of box: Adjusting Y1 (%i) to saved_Y1 (%i)\n",gd->Y1,gd->saved_Y1);
		gd->X1=gd->saved_X1;
	}
}

