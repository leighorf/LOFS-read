//What is this, the fifth iteration of this code?
#include <omp.h>

#include "include/dirstruct.h"
#include "include/limits.h"
#include "include/hdf2nc.h"
#include "include/lofs-read.h"

//Use pointers only when modifying otherwiwe pass copies? Or just use pointers everwhere?


void init_struct1(cmdline *cmd,dir_meta *dm, grid *gd)
{
	cmd->histpath        = (char *) malloc(MAXSTR*sizeof(char));
	cmd->base            = (char *) malloc(MAXSTR*sizeof(char));
	dm-> firstfilename   = (char *) malloc(MAXSTR*sizeof(char));
	dm-> saved_base      = (char *) malloc(MAXSTR*sizeof(char));
	dm-> topdir          = (char *) malloc(MAXSTR*sizeof(char));

	gd->X0=gd->Y0=gd->X1=gd->Y1=gd->Z0=gd->Z1=-1;
	cmd->time=0.0; cmd->got_base=0; cmd->optcount=0;
}


int main(int argc, char *argv[])
{
	dir_meta dm;
	hdf_meta hm;
	grid gd;
	cmdline cmd;

	init_struct1(&cmd,&dm,&gd);

	parse_cmdline_lofs2nc(argc,argv,&cmd,&dm,&gd);

	if((cptr=realpath(cmd->histpath,dm->topdir))==NULL)ERROR_STOP("realpath failed");


	dm->timedir = (char **)malloc(ntimedirs * sizeof(char *));
	for (i=0; i < ntimedirs; i++) timedir[i] = (char *)(malloc(MAXSTR * sizeof(char)));
	dirtimes = (double *)malloc(ntimedirs * sizeof(double));//times are float not int
}
