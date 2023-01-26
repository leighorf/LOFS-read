#include <omp.h>
#include "../include/lofs-dirstruct.h"
#include "../include/lofs-limits.h"
#include "../include/lofs-hdf2nc.h"
#include "../include/lofs-read.h"

void parse_cmdline_hdf2nc(int argc, char *argv[], cmdline *cmd, dir_meta *dm, grid *gd, zfpacc *zfpacc)
{
	int got_histpath,got_time,got_X0,got_X1,got_Y0,got_Y1,got_Z0,got_Z1;
	enum { OPT_HISTPATH = 1000, OPT_NCDIR, OPT_BASE, OPT_TIME, OPT_X0, OPT_Y0, OPT_X1, OPT_Y1, OPT_Z0, OPT_Z1,
		OPT_DEBUG, OPT_VERBOSE, OPT_REGENERATECACHE, OPT_SWATHS, OPT_NC3, OPT_COMPRESS_GZIP,
		OPT_COMPRESS_ZFP, OPT_COMPRESS_ZFP_LOSSLESS, OPT_NTHREADS, OPT_OFFSET, OPT_NOCMD, OPT_INTERP, OPT_CENTISECONDS, OPT_TWODWRITE,
		OPT_DEV_SHM_CACHE, OPT_CHECK_CMD_FILE,
        OPT_UINTERP_ZFPACC, OPT_VINTERP_ZFPACC, OPT_WINTERP_ZFPACC, OPT_U_ZFPACC, OPT_V_ZFPACC, OPT_W_ZFPACC,
		OPT_HWIN_SR_ACC, OPT_HWIN_GR_ACC, OPT_WINDMAG_SR_ACC, OPT_XVORT_ZFPACC, OPT_YVORT_ZFPACC, OPT_ZVORT_ZFPACC, OPT_VORTMAG_ZFPACC,
        OPT_QC_ZFPACC, OPT_QI_ZFPACC, OPT_QS_ZFPACC, OPT_QR_ZFPACC, OPT_QG_ZFPACC, OPT_QV_ZFPACC, OPT_QVPERT_ZFPACC,
        OPT_DBZ_ZFPACC, OPT_NCI_ZFPACC, OPT_NCG_ZFPACC, OPT_NCR_ZFPACC, OPT_NCS_ZFPACC, OPT_PRESPERT_ZFPACC,
		OPT_QHL_ZFPACC,OPT_CCN_ZFPACC,OPT_CCW_ZFPACC,OPT_CRW_ZFPACC,OPT_CCI_ZFPACC,OPT_CSW_ZFPACC,OPT_CHW_ZFPACC,OPT_CHL_ZFPACC,OPT_VHW_ZFPACC,OPT_VHL_ZFPACC,
        OPT_THRHOPERT_ZFPACC, OPT_RHO_ZFPACC, OPT_RHOPERT_ZFPACC, OPT_THPERT_ZFPACC, OPT_TH_ZFPACC,
        OPT_PI_ZFPACC, OPT_PRS_ZFPACC, OPT_PIPERT_ZFPACC, OPT_TKE_SG_ZFPACC, OPT_KH_ZFPACC, OPT_KM_ZFPACC,
		OPT_WB_BUOY_ACC, OPT_UB_PGRAD_ACC, OPT_VB_PGRAD_ACC, OPT_WB_PGRAD_ACC, OPT_XVORT_STRETCH_ACC,
		OPT_YVORT_STRETCH_ACC, OPT_ZVORT_STRETCH_ACC, OPT_XVORT_BARO_ACC, OPT_YVORT_BARO_ACC, OPT_XVORT_SOLENOID_ACC,
		OPT_YVORT_SOLENOID_ACC, OPT_ZVORT_SOLENOID_ACC, OPT_HVORT_ACC, OPT_STREAMVORT_ACC, OPT_QIQVPERT_ACC,
		OPT_QTOT_ACC, OPT_TEMPC_ACC, OPT_HDIV_ACC
	};

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
		{"devshmcache", optional_argument, 0, OPT_DEV_SHM_CACHE},
		{"checkcmd", optional_argument, 0, OPT_CHECK_CMD_FILE},
		/* All the ZFP stuff for netcdf output now */
		{"u_acc",optional_argument, 0,    OPT_U_ZFPACC},
		{"v_acc",optional_argument, 0,    OPT_V_ZFPACC},
		{"w_acc",optional_argument, 0,    OPT_W_ZFPACC},
		{"uinterp_acc",optional_argument, 0,        OPT_UINTERP_ZFPACC},
		{"vinterp_acc",optional_argument, 0,        OPT_VINTERP_ZFPACC},
		{"winterp_acc",optional_argument, 0,        OPT_WINTERP_ZFPACC},
		{"hwin_sr_acc",optional_argument, 0,           OPT_HWIN_SR_ACC},
		{"hwin_gr_acc",optional_argument, 0,           OPT_HWIN_GR_ACC},
		{"windmag_sr_acc",optional_argument, 0,        OPT_WINDMAG_SR_ACC},
		{"xvort_acc",optional_argument, 0,          OPT_XVORT_ZFPACC},
		{"yvort_acc",optional_argument, 0,          OPT_YVORT_ZFPACC},
		{"zvort_acc",optional_argument, 0,          OPT_ZVORT_ZFPACC},
		{"vortmag_acc",optional_argument, 0,        OPT_VORTMAG_ZFPACC},
		{"qc_acc",optional_argument, 0,             OPT_QC_ZFPACC},
		{"qi_acc",optional_argument, 0,             OPT_QI_ZFPACC},
		{"qs_acc",optional_argument, 0,             OPT_QS_ZFPACC},
		{"qr_acc",optional_argument, 0,             OPT_QR_ZFPACC},
		{"qg_acc",optional_argument, 0,             OPT_QG_ZFPACC},
		{"qv_acc",optional_argument, 0,             OPT_QV_ZFPACC},
		{"qvpert_acc",optional_argument, 0,         OPT_QVPERT_ZFPACC},
		{"dbz_acc",optional_argument, 0,            OPT_DBZ_ZFPACC},
		{"nci_acc",optional_argument, 0,            OPT_NCI_ZFPACC},
		{"ncg_acc",optional_argument, 0,            OPT_NCG_ZFPACC},
		{"ncr_acc",optional_argument, 0,            OPT_NCR_ZFPACC},
		{"ncs_acc",optional_argument, 0,            OPT_NCS_ZFPACC},
		{"qhl_acc",optional_argument, 0,            OPT_QHL_ZFPACC},
		{"ccn_acc",optional_argument, 0,            OPT_CCN_ZFPACC},
		{"ccw_acc",optional_argument, 0,            OPT_CCW_ZFPACC},
		{"crw_acc",optional_argument, 0,            OPT_CRW_ZFPACC},
		{"cci_acc",optional_argument, 0,            OPT_CCI_ZFPACC},
		{"chw_acc",optional_argument, 0,            OPT_CHW_ZFPACC},
		{"chl_acc",optional_argument, 0,            OPT_CHL_ZFPACC},
		{"chl_vhw",optional_argument, 0,            OPT_VHW_ZFPACC},
		{"chl_vhl",optional_argument, 0,            OPT_VHL_ZFPACC},
		{"prespert_acc",optional_argument, 0,       OPT_PRESPERT_ZFPACC},
		{"thrhopert_acc",optional_argument, 0,      OPT_THRHOPERT_ZFPACC},
		{"rho_acc",optional_argument, 0,               OPT_RHO_ZFPACC},
		{"rhopert_acc",optional_argument, 0,           OPT_RHOPERT_ZFPACC},
		{"thpert_acc",optional_argument, 0,         OPT_THPERT_ZFPACC},
		{"th_acc",optional_argument, 0,             OPT_TH_ZFPACC},
		{"pi_acc",optional_argument, 0,             OPT_PI_ZFPACC},
		{"prs_acc",optional_argument, 0,            OPT_PRS_ZFPACC},
		{"pipert_acc",optional_argument, 0,         OPT_PIPERT_ZFPACC},
		{"tke_sg_acc",optional_argument, 0,         OPT_TKE_SG_ZFPACC},
		{"kh_acc",optional_argument, 0,             OPT_KH_ZFPACC},
		{"km_acc",optional_argument, 0,             OPT_KM_ZFPACC},
		{"wb_buoy_acc",optional_argument, 0,        OPT_WB_BUOY_ACC},
		{"ub_pgrad_acc",optional_argument, 0,       OPT_UB_PGRAD_ACC},
		{"vb_pgrad_acc",optional_argument, 0,       OPT_VB_PGRAD_ACC},
		{"wb_pgrad_acc",optional_argument, 0,       OPT_WB_PGRAD_ACC},
		{"xvort_stretch_acc",optional_argument, 0,  OPT_XVORT_STRETCH_ACC},
		{"yvort_stretch_acc",optional_argument, 0,  OPT_YVORT_STRETCH_ACC},
		{"zvort_stretch_acc",optional_argument, 0,  OPT_ZVORT_STRETCH_ACC},
		{"xvort_baro_acc",optional_argument, 0,        OPT_XVORT_BARO_ACC},
		{"yvort_baro_acc",optional_argument, 0,        OPT_YVORT_BARO_ACC},
		{"xvort_solenoid_acc",optional_argument, 0,    OPT_XVORT_SOLENOID_ACC},
		{"yvort_solenoid_acc",optional_argument, 0,    OPT_YVORT_SOLENOID_ACC},
		{"zvort_solenoid_acc",optional_argument, 0,    OPT_ZVORT_SOLENOID_ACC},
		{"hvort_acc",optional_argument, 0,             OPT_HVORT_ACC},
		{"streamvort_acc",optional_argument, 0,        OPT_STREAMVORT_ACC},
		{"qiqvpert_acc",optional_argument, 0,          OPT_QIQVPERT_ACC},
		{"qtot_acc",optional_argument, 0,              OPT_QTOT_ACC},
		{"tempC_acc",optional_argument, 0,             OPT_TEMPC_ACC},
		{"hdiv_acc",optional_argument, 0,              OPT_HDIV_ACC},

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
//			case OPT_ALLVARS:
//				cmd->do_allvars=1;
//				cmd->optcount++;
//				break;
			case OPT_COMPRESS_GZIP:
				cmd->gzip = 1;
				cmd->optcount++;
				break;
			case OPT_COMPRESS_ZFP:
				cmd->zfp = 1;
				cmd->optcount++;
				printf("*** Default ZFP accuracy factors set in hdf2nc-util.c ***\n");
				break;
			case OPT_COMPRESS_ZFP_LOSSLESS:
				cmd->zfplossless = 1;
				cmd->optcount++;
				printf("*** ZFP lossless data chosen ***\n");
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
			case OPT_DEV_SHM_CACHE:
				cmd->devshmcache=1;
				cmd->optcount++;
				break;
			case OPT_CHECK_CMD_FILE:
				cmd->checkcmd=1;
				cmd->optcount++;
				break;
            case OPT_U_ZFPACC:
				zfpacc->netcdf->u = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_V_ZFPACC:
				zfpacc->netcdf->v = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_W_ZFPACC:
				zfpacc->netcdf->w = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_UINTERP_ZFPACC:
				zfpacc->netcdf->uinterp = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_VINTERP_ZFPACC:
				zfpacc->netcdf->vinterp = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_WINTERP_ZFPACC:
				zfpacc->netcdf->winterp = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_HWIN_SR_ACC:
				zfpacc->netcdf->hwin_sr = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_HWIN_GR_ACC:
				zfpacc->netcdf->hwin_gr = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_WINDMAG_SR_ACC:
				zfpacc->netcdf->windmag_sr = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_XVORT_ZFPACC:
				zfpacc->netcdf->xvort = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_YVORT_ZFPACC:
				zfpacc->netcdf->yvort = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_ZVORT_ZFPACC:
				zfpacc->netcdf->zvort = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_VORTMAG_ZFPACC:
				zfpacc->netcdf->vortmag = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_QC_ZFPACC:
				zfpacc->netcdf->qc = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_QI_ZFPACC:
				zfpacc->netcdf->qi = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_QS_ZFPACC:
				zfpacc->netcdf->qs = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_QR_ZFPACC:
				zfpacc->netcdf->qr = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_QG_ZFPACC:
				zfpacc->netcdf->qg = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_QV_ZFPACC:
				zfpacc->netcdf->qv = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_QVPERT_ZFPACC:
				zfpacc->netcdf->qvpert = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_DBZ_ZFPACC:
				zfpacc->netcdf->dbz = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_NCI_ZFPACC:
				zfpacc->netcdf->nci = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_NCG_ZFPACC:
				zfpacc->netcdf->ncg = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_NCR_ZFPACC:
				zfpacc->netcdf->ncr = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_NCS_ZFPACC:
				zfpacc->netcdf->ncs = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_PRESPERT_ZFPACC:
				zfpacc->netcdf->prespert = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_THRHOPERT_ZFPACC:
				zfpacc->netcdf->thrhopert = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_RHO_ZFPACC:
				zfpacc->netcdf->rho = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_RHOPERT_ZFPACC:
				zfpacc->netcdf->rhopert = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_THPERT_ZFPACC:
				zfpacc->netcdf->thpert = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_TH_ZFPACC:
				zfpacc->netcdf->th = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_PI_ZFPACC:
				zfpacc->netcdf->pi = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_PRS_ZFPACC:
				zfpacc->netcdf->prs = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_PIPERT_ZFPACC:
				zfpacc->netcdf->pipert = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_TKE_SG_ZFPACC:
				zfpacc->netcdf->tke_sg = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_KH_ZFPACC:
				zfpacc->netcdf->kh = atof(optarg);
				cmd->optcount++;
				break;
            case OPT_KM_ZFPACC:
				zfpacc->netcdf->km = atof(optarg);
				cmd->optcount++;
				break;

			case OPT_WB_BUOY_ACC:
				zfpacc->netcdf->wb_buoy = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_UB_PGRAD_ACC:
				zfpacc->netcdf->ub_pgrad = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_VB_PGRAD_ACC:
				zfpacc->netcdf->vb_pgrad = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_WB_PGRAD_ACC:
				zfpacc->netcdf->wb_pgrad = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_XVORT_STRETCH_ACC:
				zfpacc->netcdf->xvort_stretch = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_YVORT_STRETCH_ACC:
				zfpacc->netcdf->yvort_stretch = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_ZVORT_STRETCH_ACC:
				zfpacc->netcdf->zvort_stretch = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_XVORT_BARO_ACC:
				zfpacc->netcdf->xvort_baro = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_YVORT_BARO_ACC:
				zfpacc->netcdf->yvort_baro = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_XVORT_SOLENOID_ACC:
				zfpacc->netcdf->xvort_solenoid = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_YVORT_SOLENOID_ACC:
				zfpacc->netcdf->yvort_solenoid = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_ZVORT_SOLENOID_ACC:
				zfpacc->netcdf->zvort_solenoid = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_HVORT_ACC:
				zfpacc->netcdf->hvort = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_STREAMVORT_ACC:
				zfpacc->netcdf->streamvort = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_QIQVPERT_ACC:
				zfpacc->netcdf->qiqvpert = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_QTOT_ACC:
				zfpacc->netcdf->qtot = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_TEMPC_ACC:
				zfpacc->netcdf->tempC = atof(optarg);
				cmd->optcount++;
				break;
			case OPT_HDIV_ACC:
				zfpacc->netcdf->hdiv = atof(optarg);
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

void parse_cmdline_grabpoint(int argc, char *argv[], cmdline *cmd, dir_meta *dm, grid *gd, zfpacc *zfpacc)
{
	int got_histpath,got_time,got_XC,got_YC,got_ZC;
	enum { OPT_HISTPATH = 1000, OPT_TIME, OPT_XC, OPT_YC, OPT_ZC, OPT_HEADER };

	static struct option long_options[] =
	{
		{"histpath", required_argument, 0, OPT_HISTPATH},
		{"time",     required_argument, 0, OPT_TIME},
		{"xc",       required_argument, 0, OPT_XC},
		{"yc",       required_argument, 0, OPT_YC},
		{"zc",       required_argument, 0, OPT_ZC},
		{"header",   optional_argument, 0, OPT_HEADER},
		{0, 0, 0, 0}//sentinel, needed!
	};

	got_histpath=got_time=got_XC=got_YC=got_ZC=0;

	int bail = 0;
	if (argc == 1)
	{
		fprintf(stderr,
		"Usage: %s --histpath=[histpath] --xc=[XC] --yc=[YC] --zc=[ZC] --time=[time] varname\n",argv[0]);
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
			case OPT_TIME:
				cmd->time = atof(optarg);
//				printf("cmd->time = %s %12.6f\n",optarg,cmd->time);exit(0);
				got_time=1;
				break;
			case OPT_XC:
				gd->XC = atof(optarg);
				got_XC=1;
				cmd->optcount++;
				break;
			case OPT_YC:
				gd->YC = atof(optarg);
				got_YC=1;
				cmd->optcount++;
				break;
			case OPT_ZC:
				gd->ZC = atof(optarg);
				got_ZC=1;
				cmd->optcount++;
				break;
			case OPT_HEADER:
				cmd->header=1;
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
	if (!got_XC)     { fprintf(stderr,"--xc not specified\n");     bail = 1; }
	if (!got_YC)     { fprintf(stderr,"--yc not specified\n");     bail = 1; }
	if (!got_ZC)     { fprintf(stderr,"--zc not specified\n");     bail = 1; }

	if (bail)   { fprintf(stderr,"Insufficient arguments to %s, exiting.\n",argv[0]); exit(-1); }
}
