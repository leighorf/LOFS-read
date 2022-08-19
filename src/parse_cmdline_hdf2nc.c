#include <omp.h>
#include "../include/lofs-dirstruct.h"
#include "../include/lofs-limits.h"
#include "../include/lofs-hdf2nc.h"
#include "../include/lofs-read.h"

void parse_cmdline_hdf2nc(int argc, char *argv[], cmdline *cmd, dir_meta *dm, grid *gd)
{
	int got_histpath,got_time,got_X0,got_X1,got_Y0,got_Y1,got_Z0,got_Z1;
	enum { OPT_HISTPATH = 1000, OPT_NCDIR, OPT_BASE, OPT_TIME, OPT_X0, OPT_Y0, OPT_X1, OPT_Y1, OPT_Z0, OPT_Z1,
		OPT_DEBUG, OPT_VERBOSE, OPT_REGENERATECACHE, OPT_ALLVARS, OPT_SWATHS, OPT_NC3, OPT_COMPRESS_GZIP,
		OPT_COMPRESS_ZFP, OPT_COMPRESS_ZFP_LOSSLESS, OPT_NTHREADS, OPT_OFFSET, OPT_NOCMD, OPT_INTERP, OPT_CENTISECONDS, OPT_TWODWRITE };

	static struct option long_options[] =
	{
		{"histpath", required_argument, 0, OPT_HISTPATH},
		{"ncdir", optional_argument, 0, OPT_NCDIR},
		{"base",     optional_argument, 0, OPT_BASE},
		{"time",     required_argument, 0, OPT_TIME},
		{"x0",       optional_argument, 0, OPT_X0},
		{"y0",       optional_argument, 0, OPT_Y0},
		{"x1",       optional_argument, 0, OPT_X1},
		{"y1",       optional_argument, 0, OPT_Y1},
		{"z0",       optional_argument, 0, OPT_Z0},
		{"z1",       optional_argument, 0, OPT_Z1},
		{"centiseconds", optional_argument, 0, OPT_CENTISECONDS},
		{"debug",    optional_argument, 0, OPT_DEBUG},
		{"verbose",  optional_argument, 0, OPT_VERBOSE},
		{"recache",  optional_argument, 0, OPT_REGENERATECACHE},
		{"allvars",  optional_argument, 0, OPT_ALLVARS},
		{"swaths",   optional_argument, 0, OPT_SWATHS},
		{"nc3",      optional_argument, 0, OPT_NC3},
		{"gzip",     optional_argument, 0, OPT_COMPRESS_GZIP},
		{"zfp",      optional_argument, 0, OPT_COMPRESS_ZFP},
		{"zfplossless",      optional_argument, 0, OPT_COMPRESS_ZFP_LOSSLESS},
		{"nthreads", optional_argument, 0, OPT_NTHREADS},
		{"twodwrite",optional_argument, 0, OPT_TWODWRITE},
		{"offset",   optional_argument, 0, OPT_OFFSET},
		{"nocmd",   optional_argument, 0, OPT_NOCMD},
		{"interp", optional_argument, 0, OPT_INTERP},
		{0, 0, 0, 0}//sentinel, needed!
	};

	got_histpath=got_time=got_X0=got_X1=got_Y0=got_Y1=got_Z0=got_Z1=0;

	int bail = 0;
	if (argc == 1)
	{
		fprintf(stderr,
		"Usage: %s --histpath=[histpath] --base=[base] --x0=[X0] --y0=[Y0] --x1=[X1] --y1=[Y1] --z0=[Z0] --z1=[Z1] --time=[time] [varname1 ... varnameN] \n",argv[0]);
		exit(0);
	}

	while (1)
	{
		int r;
		int option_index = 0;
		r = getopt_long_only (argc, argv,"",long_options,&option_index);
		if (r == -1) break;

		switch(r)
		{
			case OPT_HISTPATH:
				strcpy(cmd->histpath,optarg);
				got_histpath=1;
				break;
			case OPT_NCDIR:
				strcpy(cmd->ncdir,optarg);
				cmd->got_ncdir=1;
				cmd->optcount++;
				break;
			case OPT_BASE:
				strcpy(cmd->base,optarg);
				cmd->got_base=1;
				cmd->optcount++;
				break;
			case OPT_TIME:
				cmd->time = atof(optarg);
//				printf("cmd->time = %s %12.6f\n",optarg,cmd->time);exit(0);
				got_time=1;
				break;
			case OPT_X0:
				gd->X0 = atoi(optarg);
				got_X0=1;
				cmd->optcount++;
				break;
			case OPT_Y0:
				gd->Y0 = atoi(optarg);
				got_Y0=1;
				cmd->optcount++;
				break;
			case OPT_X1:
				gd->X1 = atoi(optarg);
				got_X1=1;
				cmd->optcount++;
				break;
			case OPT_Y1:
				gd->Y1 = atoi(optarg);
				got_Y1=1;
				cmd->optcount++;
				break;
			case OPT_Z0:
				gd->Z0 = atoi(optarg);
				got_Z0=1;
				cmd->optcount++;
				break;
			case OPT_Z1:
				gd->Z1 = atoi(optarg);
				got_Z1=1;
				cmd->optcount++;
				break;
			case OPT_CENTISECONDS:
				cmd->centiseconds = atoi(optarg);
				cmd->optcount++;
				break;
			case OPT_DEBUG:
				cmd->debug=1;
				cmd->optcount++;
				break;
			case OPT_VERBOSE:
				cmd->verbose=1;
				cmd->optcount++;
				break;
			case OPT_REGENERATECACHE:
				dm->regenerate_cache=1;
				cmd->optcount++;
				break;
			case OPT_SWATHS:
				cmd->do_swaths=1;
				cmd->optcount++;
				break;
			case OPT_ALLVARS:
				cmd->do_allvars=1;
				cmd->optcount++;
				break;
			case OPT_COMPRESS_GZIP:
				cmd->gzip = 1;
				cmd->optcount++;
				break;
			case OPT_COMPRESS_ZFP:
				cmd->zfp = 1;
				cmd->optcount++;
				printf("*** ZFP accuracy factors set in hdf2nc-util.c ***\n");
				break;
			case OPT_COMPRESS_ZFP_LOSSLESS:
				cmd->zfplossless = 1;
				cmd->optcount++;
				printf("*** ZFP LOSSLESS data chosen ***\n");
				break;
			case OPT_INTERP:
				cmd->use_interp=1;
				cmd->optcount++;
				break;
			case OPT_OFFSET:
				cmd->use_box_offset=1;
				cmd->optcount++;
				break;
			case OPT_NOCMD:
				cmd->write_cmd_file=0;
				cmd->optcount++;
				break;
			case OPT_TWODWRITE:
				cmd->twodwrite=1;
				cmd->optcount++;
				break;
			case OPT_NC3:
				cmd->filetype=NC_64BIT_OFFSET;
				cmd->optcount++;
				break;
			case OPT_NTHREADS:
				cmd->nthreads=atoi(optarg);
				omp_set_num_threads(cmd->nthreads);
				cmd->optcount++;
				break;
			case '?':
				fprintf(stderr,"Exiting: unknown command line option.\n");
				exit(0);
				break;
		}
	}
	if (cmd->debug==1)cmd->verbose=1; //show everything

	if (!got_histpath) { fprintf(stderr,"--histpath not specified\n"); bail = 1; }
	if (!got_time)     { fprintf(stderr,"--time not specified\n");     bail = 1; }

	if (!got_X0&&cmd->debug)      fprintf(stderr,"Will set X0 to saved_X0\n");
	if (!got_Y0&&cmd->debug)      fprintf(stderr,"Will set Y0 to saved_Y0\n");
	if (!got_X1&&cmd->debug)      fprintf(stderr,"Will set X1 to saved_X1\n");
	if (!got_Y1&&cmd->debug)      fprintf(stderr,"Will set Y1 to saved_Y1\n");
	if (!got_Z0&&cmd->debug)      fprintf(stderr,"Setting Z0 to default value of 0\n");
	if (!got_Z1&&cmd->debug)      fprintf(stderr,"Setting Z1 to default value of nz-2\n");


	if (bail)   { fprintf(stderr,"Insufficient arguments to %s, exiting.\n",argv[0]); exit(-1); }
}
