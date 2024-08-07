#include "../include/lofs-dirstruct.h"
#include "../include/lofs-limits.h"
#include "../include/lofs-hdf2nc.h"
#include "../include/lofs-read.h"
#include "../include/lofs-macros.h"

void init_structs(cmdline *cmd,dir_meta *dm, grid *gd,ncstruct *nc, readahead *rh, hdf_meta *hm, zfpacc *zfpacc)
{
	int i;
//	char syscmd[MAXSTR];

	cmd->histpath        = (char *) malloc(MAXSTR*sizeof(char));
	cmd->base            = (char *) malloc(MAXSTR*sizeof(char));
	dm->firstfilename    = (char *) malloc(MAXSTR*sizeof(char));
	dm->saved_base       = (char *) malloc(MAXSTR*sizeof(char));
	dm->topdir           = (char *) malloc(MAXSTR*sizeof(char));
	dm->devshmdir        = (char *) malloc(MAXSTR*sizeof(char));
	dm->cachedir         = (char *) malloc(MAXSTR*sizeof(char));
	nc->ncfilename       = (char *) malloc(MAXSTR*sizeof(char));
	zfpacc->lofs         = (lofs *) malloc(sizeof(lofs));
	nc->var3d            = (var3dstruct *) malloc(MAXVARIABLES*sizeof(var3dstruct));
	zfpacc->netcdf       = (netcdf *) malloc(sizeof(netcdf));

	cmd->ncdir = (char *)(malloc(MAXSTR * sizeof(char)));

	strcpy(dm->devshmdir,"/dev/shm/orf");

	dm->regenerate_cache = 0;

	gd->X0=-1;
	gd->Y0=-1;
	gd->X1=-1;
	gd->Y1=-1;
	gd->Z0=-1;
	gd->Z1=-1;
	gd->saved_X0=0;
	gd->saved_X1=0;
	gd->saved_Y0=0;
	gd->saved_Y1=0;
	gd->saved_Z0=0;
	gd->saved_Z1=0;
	cmd->er10=0;
	cmd->tusc30=0;
	cmd->time=0.0;
	cmd->got_base=0;
	cmd->got_ncdir=0;
	cmd->optcount=0;
	cmd->nthreads=8; //Changed 2/15/24
	cmd->twodwrite=0;
	cmd->inprogress=0;
	cmd->write_cmd_file=0; //CHANGED 2023-08-18. Defaults to no file.
	cmd->debug=0;
	cmd->gzip=0;
	cmd->zfp=0;
	cmd->zfplossless=0;
	cmd->bitgroom1=0;
	cmd->bitgroom2=0;
	cmd->bitgroom3=0;
	cmd->bitgroom_nsd=8; // Eight significant digits/bits
	cmd->devshmcache=0;
	cmd->checkcmd=0;
	cmd->centiseconds=0;
	cmd->verbose=0;
	cmd->header=0;
//	cmd->do_allvars=0;
	cmd->use_box_offset=0;
	cmd->do_swaths=0;
	cmd->filetype=NC_NETCDF4;
	nc->twodslice=0;
	for (i=0; i < MAXVARIABLES; i++)nc->var3d[i].zfpacc_LOFS=(double) -9999.0;

	rh->u=0; rh->v=0; rh->w=0;
	rh->ppert=0; rh->thrhopert=0;

//Below not really readahead, used for mallocs
	rh->qiqvpert=0;
	rh->qcqi=0;
	rh->qgqhqr=0;
	rh->qtot=0;
	rh->tempC=0;
	rh->budgets=0;
	rh->vortmag=0;
	rh->hvort=0;
	rh->streamvort=0;
	rh->interp=0;

//Block of text that will be saved as metadata to netCDF file - contains
//all LOFS variable zfp accuracy parameters
	hm->zfpacc_LOFS_all        = (char **) malloc(MAXVARIABLES*sizeof(char *));
	for (i=0; i < MAXVARIABLES; i++) hm->zfpacc_LOFS_all[i] = (char *)(malloc(MAXSTR * sizeof(char)));

/*
 2022-08-23 Here are our default zfp accuracy parameters for saved
 netcdf variables - these are reasonably sane choices but watch out for
 post-processing issues! Either twiddle these to your heart's content,
 or use command line options to save them (see parse_command_line), or
 create a routine that overwrites them. Could 'bulk set' them with say
 --zfp-set1 --zfp-set2 etc. new options, if you want to have different
 "tiers" of zfp or different research projects, etc...
*/

/* Leigh Orf's default ZFP accuracy parameters */

/* Each can be overridden on the command line i.e.: hdf2nc --thrhopert_acc=0.001 --prespert_acc=0.01 ... etc */

	zfpacc->netcdf->u       =        1.0e-2;
	zfpacc->netcdf->v       =        1.0e-2;
	zfpacc->netcdf->w       =        1.0e-2;
	zfpacc->netcdf->uinterp =        1.0e-2;
	zfpacc->netcdf->vinterp =        1.0e-2;
	zfpacc->netcdf->winterp =        1.0e-2;
	zfpacc->netcdf->windmag_sr =     1.0e-2;
	zfpacc->netcdf->hwin_sr =        1.0e-2;
	zfpacc->netcdf->hwin_gr =        1.0e-2;
	zfpacc->netcdf->u_gr =           1.0e-2;
	zfpacc->netcdf->v_gr =           1.0e-2;
	zfpacc->netcdf->thrhopert =      1.0e-2;
	zfpacc->netcdf->prespert =       1.0e-1;
	zfpacc->netcdf->rhopert =        1.0e-5;
	zfpacc->netcdf->xvort =          1.0e-2;
	zfpacc->netcdf->yvort =          1.0e-2;
	zfpacc->netcdf->zvort =          1.0e-2;
	zfpacc->netcdf->vortmag =        1.0e-2;
	zfpacc->netcdf->streamvort =     1.0e-2;
/* Keep in mind all mixing ratios are g/kg here*/
/* Accuracy parameters less than 0.0 result in LOSSLESS ZFP */
	zfpacc->netcdf->qc =             1.0e-3;
	zfpacc->netcdf->qi =             1.0e-3;
	zfpacc->netcdf->qvpert =         1.0e-3;
	zfpacc->netcdf->ncr =            -1.0;
	zfpacc->netcdf->nci =            -1.0;
	zfpacc->netcdf->ncs =            -1.0;
	zfpacc->netcdf->ncg =            -1.0;
	zfpacc->netcdf->dbz =            5.0;

	zfpacc->netcdf->qr =             -1.0;
	zfpacc->netcdf->qs =             -1.0;
	zfpacc->netcdf->qg =             -1.0;
/* NSSL microphysics only */
	zfpacc->netcdf->qhl =            -1.0;
	zfpacc->netcdf->crw =            -1.0;
	zfpacc->netcdf->csw =            -1.0;
	zfpacc->netcdf->chl =            -1.0;
	zfpacc->netcdf->chw =            -1.0;
/* End Rachael lossless vars */

/* Don't save these NSSL vars unless you really need them*/
	zfpacc->netcdf->ccn =            -1.0;
	zfpacc->netcdf->ccw =            -1.0;
	zfpacc->netcdf->cci =            -1.0;
	zfpacc->netcdf->vhw =            -1.0;
	zfpacc->netcdf->vhl =            -1.0;
	zfpacc->netcdf->zhl =            -1.0;
	zfpacc->netcdf->zhw =            -1.0;
	zfpacc->netcdf->zrw =            -1.0;

	zfpacc->netcdf->tke_sg =         1.0e-2;
	zfpacc->netcdf->khh =            1.0e-2;
	zfpacc->netcdf->khh_interp =     1.0e-2;
	zfpacc->netcdf->kmh =            1.0e-2;
	zfpacc->netcdf->kmh_interp =     1.0e-2;
	zfpacc->netcdf->khv =            1.0e-2;
	zfpacc->netcdf->khv_interp =     1.0e-2;
	zfpacc->netcdf->kmv =            1.0e-2;
	zfpacc->netcdf->kmv_interp =     1.0e-2;

	zfpacc->netcdf->wb_buoy =        1.0e-3;
	zfpacc->netcdf->wb_buoy_interp = 1.0e-3;
	zfpacc->netcdf->ub_pgrad =       1.0e-2;
	zfpacc->netcdf->vb_pgrad =       1.0e-2;
	zfpacc->netcdf->wb_pgrad =       1.0e-2;
	zfpacc->netcdf->ub_pgrad_interp =1.0e-2;
	zfpacc->netcdf->vb_pgrad_interp =1.0e-2;
	zfpacc->netcdf->wb_pgrad_interp =1.0e-2;
	zfpacc->netcdf->xvort_stretch =  1.0e-4;
	zfpacc->netcdf->yvort_stretch =  1.0e-4;
	zfpacc->netcdf->zvort_stretch =  1.0e-4;
	zfpacc->netcdf->xvort_tilt =     1.0e-4;
	zfpacc->netcdf->yvort_tilt =     1.0e-4;
	zfpacc->netcdf->zvort_tilt =     1.0e-4;
	zfpacc->netcdf->xvort_baro =     1.0e-4;
	zfpacc->netcdf->yvort_baro =     1.0e-4;
	zfpacc->netcdf->xvort_solenoid = 1.0e-5; //Solenoidal terms get a bit more accuracy - they are small!
	zfpacc->netcdf->yvort_solenoid = 1.0e-5;
	zfpacc->netcdf->zvort_solenoid = 1.0e-5;
	zfpacc->netcdf->hvort =          1.0e-2;
	zfpacc->netcdf->qiqvpert =       1.0e-4;
	zfpacc->netcdf->qcqi =           1.0e-4;
	zfpacc->netcdf->qgqhqr =         1.0e-4;
	zfpacc->netcdf->qtot =           1.0e-4;
	zfpacc->netcdf->tempC =          1.0e-1;
	zfpacc->netcdf->hdiv =           1.0e-3;
	zfpacc->netcdf->liutexmag =      1.0e-3;
	zfpacc->netcdf->liutex_x =      1.0e-3;
	zfpacc->netcdf->liutex_y =      1.0e-3;
	zfpacc->netcdf->liutex_z =      1.0e-3;

	/* These are read only in the sense that we read them from LOFS attributes. These
	 * are the zfp accuracy values that were written with CM1-LOFS. Initialize
	 * negative for existence testing */

	zfpacc->lofs->u = -1.0;
	zfpacc->lofs->v = -1.0;
	zfpacc->lofs->w = -1.0;
	zfpacc->lofs->uinterp = -1.0;
	zfpacc->lofs->vinterp = -1.0;
	zfpacc->lofs->winterp = -1.0;
	zfpacc->lofs->prespert = -1.0;
	zfpacc->lofs->thrhopert = -1.0;
	zfpacc->lofs->dbz = -1.0;
	zfpacc->lofs->qc = -1.0;
	zfpacc->lofs->qr = -1.0;
	zfpacc->lofs->qi = -1.0;
	zfpacc->lofs->qs = -1.0;
	zfpacc->lofs->qg = -1.0;
	zfpacc->lofs->nci = -1.0;
	zfpacc->lofs->ncg = -1.0;
	zfpacc->lofs->ncr = -1.0;
	zfpacc->lofs->ncs = -1.0;
	zfpacc->lofs->qvpert = -1.0;
	zfpacc->lofs->thpert = -1.0;
	zfpacc->lofs->th = -1.0;
	zfpacc->lofs->prs = -1.0;
	zfpacc->lofs->pi = -1.0;
	zfpacc->lofs->pipert = -1.0;
	zfpacc->lofs->rho = -1.0;
	zfpacc->lofs->rhopert = -1.0;
	zfpacc->lofs->tke_sg = -1.0;
	zfpacc->lofs->km = -1.0;
	zfpacc->lofs->kh = -1.0;
	zfpacc->lofs->qv = -1.0;
}

void dealloc_structs(cmdline *cmd,dir_meta *dm, grid *gd,ncstruct *nc, readahead *rh) {
}


void copy_grid_to_requested_cube (requested_cube *rc, grid gd)
{
	rc->X0=gd.X0; rc->X1=gd.X1;
	rc->Y0=gd.Y0; rc->Y1=gd.Y1;
	rc->Z0=gd.Z0; rc->Z1=gd.Z1;
	rc->NX=gd.NX;
	rc->NY=gd.NY;
	rc->NZ=gd.NZ;
}

void write_hdf2nc_command_txtfile(int argc, char *argv[],ncstruct nc)
{
	int i;
	FILE *fp;
	char *cmdfilename;
	cmdfilename = (char *)malloc(MAXSTR*sizeof(char));
	sprintf(cmdfilename,"%s%s",nc.ncfilename,".cmd");
	if ((fp = fopen(cmdfilename,"w")) != NULL)
	{
		for (i=0; i<argc; i++)
		{
			fprintf(fp,"%s ",argv[i]);
		}
	fprintf(fp,"\n");
	fclose(fp);
    }
	free(cmdfilename);
}


void get_saved_base(char *timedir, char *saved_base)
{
	// Just grab the basename so we can have it set automatically
	// for the netcdf files if so desired. Apparenlty the easiest
	// way to do this is to just remove the last 14 characters of
	// one of the timedir directories

	int i,tdlen;
//ORF TODO: make 14 a constant somewhere
//We may wish to allow flexibility in the floating point
//left and right of the decimal.
	tdlen=strlen(timedir)-14;
	for(i=0;i<tdlen;i++)saved_base[i]=timedir[i];
	//strncpy(saved_base,timedir,tdlen);
	saved_base[tdlen]='\0';
}

void set_span(grid *gd,hdf_meta hm,cmdline cmd)
{
	/* We do this because of stencils. Note, saved_Z0 is always 0 in LOFS.
	 * If we want to save only elevated data, you could set lower data to float(0.0) and compress*/

	gd->saved_X0+=1; gd->saved_X1-=1;
	gd->saved_Y0+=1; gd->saved_Y1-=1;
	gd->saved_Z1-=2;

	/* X and Y values are with respect to (0,0) of the actual saved data, not with respect
	 * to actual (0,0) (done because I often choose new netCDF subdomains from ncview).
	 * --offset allows this
	 */

	/* If we didn't specify values at the command line, value is -1
	 * Set them to values specifying all the saved data */

	if(cmd.use_box_offset)
	{
		if (gd->X0<0) gd->X0 = gd->saved_X0; else gd->X0 = gd->X0 + gd->saved_X0;
		if (gd->X1<0) gd->X1 = gd->saved_X1; else gd->X1 = gd->X1 + gd->saved_X0;
		if (gd->Y0<0) gd->Y0 = gd->saved_Y0; else gd->Y0 = gd->Y0 + gd->saved_Y0;
		if (gd->Y1<0) gd->Y1 = gd->saved_Y1; else gd->Y1 = gd->Y1 + gd->saved_Y0;
	}
	else
	{
		if (gd->X0<0) gd->X0 = gd->saved_X0;
		if (gd->X1<0) gd->X1 = gd->saved_X1;
		if (gd->Y0<0) gd->Y0 = gd->saved_Y0;
		if (gd->Y1<0) gd->Y1 = gd->saved_Y1;
	}
	if (gd->Z0<0) gd->Z0 = gd->saved_Z0;
	if (gd->Z1<0) gd->Z1 = gd->saved_Z1; //ORF -2 not -1 now

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
		gd->Y1=gd->saved_Y1;
	}
}

float grabpoint(grid *gd,hdf_meta hm,dir_meta dm,cmdline cmd, mesh msh, char *varname)
{
	int ix,iy,iz,i,j,k;
	gd->saved_X0+=1; gd->saved_X1-=1;
	gd->saved_Y0+=1; gd->saved_Y1-=1;
	gd->saved_Z1-=2;
	float xc,yc,zc; //these are requested and do not need to lie on the grid
	int ix0,iy0,iz0;
	int nx,ny,nz;
	float *xh,*yh,*zh,*xf,*yf,*zf;
	float *buf,*b0;
	float interpval;
	requested_cube rc;
	float w0,w1,w2,p0,p1,p2,p3,p4,p5,p6;
	float dx,dy,dz;
	float avg=0;
	point p[2][2][2];

	//Here we figure out the grid indices we need
	//Need to determine if we are u,v,w, or scalar to choose correct
	//mesh

	buf = (float *)malloc(8*sizeof(float));
	
	xc=gd->XC;
	yc=gd->YC;
	zc=gd->ZC;
	
	nx = hm.nx;
	ny = hm.ny;
	nz = hm.nz;

	xh = msh.xhfull;
	yh = msh.yhfull;
	zh = msh.zh;

	xf = msh.xffull;
	yf = msh.yffull;
	zf = msh.zf;

	//More annoying code, this to fill our 2x2x2 cube of structures so
	//we can then do the final interpolation. Lots of checks for what
	//mesh we are on, which depends on whether you are asking for
	//u,v,w,or a scalar

	if(same(varname,"u")||same(varname,"u_gr"))
	{
		for(ix=0;ix<nx;ix++)
		{
			if(xf[ix] <= xc && xf[ix+1] > xc)
			{
				ix0=ix;
				for (k=0;k<2;k++)
				{
					for(j=0;j<2;j++)
					{
						p[k][j][0].XC = xf[ix0];
						p[k][j][1].XC = xf[ix0+1];
					}
				}
				break;
			}
		}
		for(iy=0;iy<ny;iy++)
		{
			if(yh[iy] <= yc && yh[iy+1] > yc)
			{
				iy0=iy;
				for (k=0;k<2;k++)
				{
					for(i=0;i<2;i++)
					{
						p[k][0][i].YC = yh[iy0];
						p[k][1][i].YC = yh[iy0+1];
					}
				}
				break;
			}
		}
		for(iz=0;iz<nz;iz++)
		{
			if(zh[iz] <= zc && zh[iz+1] > zc)
			{
				iz0=iz;
				for (j=0;j<2;j++)
				{
					for(i=0;i<2;i++)
					{
						p[0][j][i].ZC = zh[iz0];
						p[1][j][i].ZC = zh[iz0+1];
					}
				}
				break;
			}
		}
	}
	else if(same(varname,"v")||same(varname,"v_gr"))
	{
		for(ix=0;ix<nx;ix++)
		{
			if(xh[ix] <= xc && xh[ix+1] > xc)
			{
				ix0=ix;
				for (k=0;k<2;k++)
				{
					for(j=0;j<2;j++)
					{
						p[k][j][0].XC = xh[ix0];
						p[k][j][1].XC = xh[ix0+1];
					}
				}
				break;
			}
		}
		for(iy=0;iy<ny;iy++)
		{
			if(yf[iy] <= yc && yf[iy+1] > yc)
			{
				iy0=iy;
				for (k=0;k<2;k++)
				{
					for(i=0;i<2;i++)
					{
						p[k][0][i].YC = yf[iy0];
						p[k][1][i].YC = yf[iy0+1];
					}
				}
				break;
			}
		}
		for(iz=0;iz<nz;iz++)
		{
			if(zh[iz] <= zc && zh[iz+1] > zc)
			{
				iz0=iz;
				for (j=0;j<2;j++)
				{
					for(i=0;i<2;i++)
					{
						p[0][j][i].ZC = zh[iz0];
						p[1][j][i].ZC = zh[iz0+1];
					}
				}
				break;
			}
		}
	}
	else if(same(varname,"w"))
	{
		for(ix=0;ix<nx;ix++)
		{
			if(xh[ix] <= xc && xh[ix+1] > xc)
			{
				ix0=ix;
				for (k=0;k<2;k++)
				{
					for(j=0;j<2;j++)
					{
						p[k][j][0].XC = xh[ix0];
						p[k][j][1].XC = xh[ix0+1];
					}
				}
				break;
			}
		}
		for(iy=0;iy<ny;iy++)
		{
			if(yh[iy] <= yc && yh[iy+1] > yc)
			{
				iy0=iy;
				for (k=0;k<2;k++)
				{
					for(i=0;i<2;i++)
					{
						p[k][0][i].YC = yh[iy0];
						p[k][1][i].YC = yh[iy0+1];
					}
				}
				break;
			}
		}
		for(iz=0;iz<nz;iz++)
		{
			if(zf[iz] <= zc && zf[iz+1] > zc)
			{
				iz0=iz;
				for (j=0;j<2;j++)
				{
					for(i=0;i<2;i++)
					{
						p[0][j][i].ZC = zf[iz0];
						p[1][j][i].ZC = zf[iz0+1];
					}
				}
				break;
			}
		}
	}
	else
	{
		for(ix=0;ix<nx;ix++)
		{
			if(xh[ix] <= xc && xh[ix+1] > xc)
			{
				ix0=ix;
				for (k=0;k<2;k++)
				{
					for(j=0;j<2;j++)
					{
						p[k][j][0].XC = xh[ix0];
						p[k][j][1].XC = xh[ix0+1];
					}
				}
				break;
			}
		}
		for(iy=0;iy<ny;iy++)
		{
			if(yh[iy] <= yc && yh[iy+1] > yc)
			{
				iy0=iy;
				for (k=0;k<2;k++)
				{
					for(i=0;i<2;i++)
					{
						p[k][0][i].YC = yh[iy0];
						p[k][1][i].YC = yh[iy0+1];
					}
				}
				break;
			}
		}
		for(iz=0;iz<nz;iz++)
		{
			if(zh[iz] <= zc && zh[iz+1] > zc)
			{
				iz0=iz;
				for (j=0;j<2;j++)
				{
					for(i=0;i<2;i++)
					{
						p[0][j][i].ZC = zh[iz0];
						p[1][j][i].ZC = zh[iz0+1];
					}
				}
				break;
			}
		}
	}

	rc.X0=gd->X0=ix0; rc.X1=gd->X1=ix0+1;
	rc.Y0=gd->Y0=iy0; rc.Y1=gd->Y1=iy0+1;
	rc.Z0=gd->Z0=iz0; rc.Z1=gd->Z1=iz0+1;
	rc.NX=2;
	rc.NY=2;
	rc.NZ=2;

	if (gd->X0>gd->X1||gd->Y0>gd->Y1||gd->X1<gd->saved_X0||gd->Y1<gd->saved_Y0||gd->X0>gd->saved_X1||gd->Y0>gd->saved_Y1)
	{
		fprintf(stderr,"X: %i %i %i %i\n",gd->saved_X0,gd->X0,gd->X1,gd->saved_X1);
		fprintf(stderr,"Y: %i %i %i %i\n",gd->saved_X0,gd->Y0,gd->Y1,gd->saved_Y1);
		fprintf(stderr,"xc: %f\n",xc);
		fprintf(stderr,"yc: %f\n",yc);
		fprintf(stderr,"zc: %f\n",zc);
		ERROR_STOP("Above numbers should be monotonically increasing. Your requested data does not fit within the saved data!\n");
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
		gd->Y1=gd->saved_Y1;
	}

	b0=buf;
	/* 2023-03-02 ORF TODO: If we want to track derived stuff we're gonna have to
	 * do all that extra checking shit.. for now we assume all
	 * interpolated things are in the LOFS files... except for u_gr and
	 * v_gr which are handled below */

	if(same(varname,"u_gr"))
	{
		read_lofs_buffer(b0,"u",dm,hm,rc,cmd);
	}
	else if(same(varname,"v_gr"))
	{
		read_lofs_buffer(b0,"v",dm,hm,rc,cmd);
	}
	else
	{
		read_lofs_buffer(b0,varname,dm,hm,rc,cmd);
	}

	for(k=0;k<2;k++)
	for(j=0;j<2;j++)
	for(i=0;i<2;i++)
	{
		p[k][j][i].val=*buf++;
		avg+=p[k][j][i].val;
//		printf("%s %04i %04i %04i %14.8f %14.8f %14.8f %14.8f\n",varname,rc.X0+i,rc.Y0+j,rc.Z0+k,p[k][j][i].XC,p[k][j][i].YC,p[k][j][i].ZC,p[k][j][i].val);
	}
	avg/=8.0;

	dx = p[0][0][1].XC-p[0][0][0].XC;
	dy = p[0][1][0].YC-p[0][0][0].YC;
	dz = p[1][0][0].ZC-p[0][0][0].ZC;

	w0 = (p[0][0][1].XC-xc)/dx;
	w1 = (p[0][1][0].YC-yc)/dy;
	w2 = (p[1][0][0].ZC-zc)/dz;

	p0 = p[0][0][0].val*w0+p[0][0][1].val*(1.0-w0);
	p1 = p[0][1][0].val*w0+p[0][1][1].val*(1.0-w0);
	p2 = p0*w1+p1*(1.0-w1);
	p3 = p[1][0][0].val*w0+p[1][0][1].val*(1.0-w0);
	p4 = p[1][1][0].val*w0+p[1][1][1].val*(1.0-w0);
	p5 = p3*w1+p4*(1.0-w1);
	p6 = p2*w2+p5*(1.0-w2);
	interpval = p6;

	free(b0);

	if(same(varname,"u_gr")) interpval += msh.umove;
	if(same(varname,"v_gr")) interpval += msh.vmove;

	return interpval;
}
void allocate_1d_arrays(hdf_meta hm, grid gd, mesh *msh, sounding *snd) {

	msh->xhfull = (float *)malloc(hm.nx * sizeof(float));
	msh->yhfull = (float *)malloc(hm.ny * sizeof(float));
	msh->xffull = (float *)malloc((hm.nx+1) * sizeof(float));
	msh->yffull = (float *)malloc((hm.ny+1) * sizeof(float));
	msh->zh = (float *)malloc(hm.nz * sizeof(float));
	msh->zf = (float *)malloc(hm.nz * sizeof(float));

	msh->xhout = (float *)malloc(gd.NX * sizeof(float));
	msh->yhout = (float *)malloc(gd.NY * sizeof(float));
	msh->zhout = (float *)malloc(gd.NZ * sizeof(float));

	msh->xfout = (float *)malloc(gd.NX * sizeof(float));
	msh->yfout = (float *)malloc(gd.NY * sizeof(float));
	msh->zfout = (float *)malloc(gd.NZ * sizeof(float));

	msh->uh = (float *)malloc((gd.NX) * sizeof(float));
	msh->uf = (float *)malloc((gd.NX+1) * sizeof(float));
	msh->vh = (float *)malloc((gd.NY) * sizeof(float));
	msh->vf = (float *)malloc((gd.NY+1) * sizeof(float));
	msh->mh = (float *)malloc((gd.NZ) * sizeof(float));
	msh->mf = (float *)malloc((gd.NZ+1) * sizeof(float));


	snd->th0 = (float *)malloc((gd.NZ+1) * sizeof(float));
	snd->thv0 = (float *)malloc((gd.NZ+1) * sizeof(float));
	snd->qv0 = (float *)malloc((gd.NZ+1) * sizeof(float));
	snd->u0 = (float *)malloc((gd.NZ+1) * sizeof(float));
	snd->v0 = (float *)malloc((gd.NZ+1) * sizeof(float));
	snd->pres0 = (float *)malloc((gd.NZ+1) * sizeof(float));
	snd->pi0 = (float *)malloc((gd.NZ+1) * sizeof(float));
	snd->rho0 = (float *)malloc((gd.NZ+1) * sizeof(float));
}

void set_1d_arrays(hdf_meta hm, grid gd, mesh *msh, sounding *snd, cmdline cmd, hid_t *f_id)
{
#include "../include/lofs-constants.h"
	int ix,iy,iz,k;

	printf("gd.X0 = %i gd.X1 = %i gd.Y0 = %i gd.Y1 = %i gd.Z0 = %i gd.Z1 = %i\n",gd.X0,gd.X1,gd.Y0,gd.Y1,gd.Z0,gd.Z1);

	get0dfloat (*f_id,(char *)"mesh/dx",&msh->dx); msh->rdx=1.0/msh->dx;
	get0dfloat (*f_id,(char *)"mesh/dy",&msh->dy); msh->rdy=1.0/msh->dy;
	if (H5Lexists(*f_id,"mesh/dz",H5P_DEFAULT)>0)
	{
		get0dfloat (*f_id,(char *)"mesh/dz",&msh->dz); msh->rdz=1.0/msh->dz;
	}
	else
	{
		printf("********* dz not found in LOFS metadata; assuming isotropic, setting dz=dx\n");
		msh->dz = msh->dx; msh->rdz=1.0/msh->dz;
	}
//	printf("set_1d_arrays: cmd.tusc30 = %i\n", cmd.tusc30);
	/* Before we stored umove and vmove in the cm1hdf5 files */
	if (cmd.tusc30)
	{
		msh->umove = 18.75;
		msh->vmove = 9.0;
	}
	/* There is a bafflingly werid thing with the ER10 dataset... when I
	 * run on concurrent files I get either:
			umove = 15.200000 vmove = 10.500000 (correct)
			umove = 18.750000 vmove = 9.000000 (incorrect)
     I suspect it's someting with firstfilename but I am really quite baffled.
	 I did not change box speed for the 10 meter run! Certainly not efery 0.2 seconds LOL */
	else if (cmd.er10)
	{
		msh->umove = 15.2;
		msh->vmove = 10.5;
	}
	else //The normal case!
	{
		get0dfloat (*f_id,(char *)"mesh/umove",&msh->umove);
		get0dfloat (*f_id,(char *)"mesh/vmove",&msh->vmove);
	}

	printf("\numove = %f vmove = %f\n",msh->umove,msh->vmove);
	get1dfloat (*f_id,(char *)"mesh/xhfull",msh->xhfull,0,hm.nx);
	get1dfloat (*f_id,(char *)"mesh/yhfull",msh->yhfull,0,hm.ny);
	get1dfloat (*f_id,(char *)"mesh/xffull",msh->xffull,0,hm.nx+1);
	get1dfloat (*f_id,(char *)"mesh/yffull",msh->yffull,0,hm.ny+1);
	get1dfloat (*f_id,(char *)"mesh/zh",msh->zh,0,hm.nz);
	get1dfloat (*f_id,(char *)"mesh/zf",msh->zf,0,hm.nz);
	get1dfloat (*f_id,(char *)"basestate/qv0",snd->qv0,gd.Z0,gd.NZ+1);
	//ASSUMES Z0=0!!
//	for (k=gd.Z0; k<gd.NZ; k++) snd->qv0[k-gd.Z0] *= 1000.0; // g/kg now
	get1dfloat (*f_id,(char *)"basestate/th0",snd->th0,gd.Z0,gd.NZ+1);
	get1dfloat (*f_id,(char *)"basestate/u0",snd->u0,gd.Z0,gd.NZ+1);
	get1dfloat (*f_id,(char *)"basestate/v0",snd->v0,gd.Z0,gd.NZ+1);
	get1dfloat (*f_id,(char *)"basestate/pres0",snd->pres0,gd.Z0,gd.NZ+1);
	get1dfloat (*f_id,(char *)"basestate/rh0",snd->rho0,gd.Z0,gd.NZ+1);//NOTE: rh0 is kind of a typo, should be rho0; rh means rel. hum. in CM1
	get1dfloat (*f_id,(char *)"basestate/pi0",snd->pi0,gd.Z0,gd.NZ+1);
	//Temporary - we now save this in LOFS
	for (iz=0; iz < gd.NZ+1; iz++) snd->thv0[iz] = snd->th0[iz] * (1.0+reps*snd->qv0[iz])/(1.0+snd->qv0[iz]);

// We recreate George's mesh/derivative calculation paradigm even though
// we are usually isotropic. We need to have our code here match what
// CM1 does internally for stretched and isotropic meshes.
//
// Becuase C cannot do have negative array indices (i.e., uh[-1]) like
// F90 can, we have to offset everything to keep the same CM1-like code
// We malloc enough space for the "ghost zones" and then make sure we
// offset by the correct amount on each side. The macros take care of
// the offsetting.

	// Carefully consider the UH,UF etc. macros and make sure they are appropriate for where you are in
	// the code (should be in pointer land)
	//

	for (ix=gd.X0; ix<=gd.X1; ix++) UHp(ix-gd.X0) = msh->dx/(msh->xffull[ix+1]-msh->xffull[ix]);
	for (ix=gd.X0; ix<=gd.X1+1; ix++) UFp(ix-gd.X0) = msh->dx/(msh->xhfull[ix]-msh->xhfull[ix-1]);
	for (iy=gd.Y0; iy<=gd.Y1; iy++) VHp(iy-gd.Y0) = msh->dy/(msh->yffull[iy+1]-msh->yffull[iy]);
	for (iy=gd.Y0; iy<=gd.Y1+1; iy++) VFp(iy-gd.Y0) = msh->dy/(msh->yhfull[iy]-msh->yhfull[iy-1]);
	for (iz=gd.Z0; iz<=gd.Z1; iz++) MHp(iz-gd.Z0) = msh->dz/(msh->zf[iz+1]-msh->zf[iz]);
	// the iz-1 will index to -1 if we don't start from iz=1
	// this is how things are set in CM1 Param.F line 6874
	// ORF this only makes sense when Z0=0 for real
	if(gd.Z0==0)
	{
		for (iz=gd.Z0+1; iz<=gd.Z1+1; iz++) MFp(iz-gd.Z0) = msh->dz/(msh->zh[iz]-msh->zf[iz-1]);
		MFp(0) = MFp(1);
	}
	else
		for (iz=gd.Z0; iz<=gd.Z1+1; iz++) MFp(iz-gd.Z0) = msh->dz/(msh->zh[iz]-msh->zf[iz-1]);
	
	for (iz=gd.Z0; iz<=gd.Z1; iz++) msh->zfout[iz-gd.Z0] = msh->zf[iz]; 
	for (iz=gd.Z0; iz<=gd.Z1; iz++) msh->zhout[iz-gd.Z0] = msh->zh[iz];
	for (iy=gd.Y0; iy<=gd.Y1; iy++) msh->yfout[iy-gd.Y0] = msh->yffull[iy];
	for (ix=gd.X0; ix<=gd.X1; ix++) msh->xfout[ix-gd.X0] = msh->xffull[ix];
	for (iy=gd.Y0; iy<=gd.Y1; iy++) msh->yhout[iy-gd.Y0] = msh->yhfull[iy];
	for (ix=gd.X0; ix<=gd.X1; ix++) msh->xhout[ix-gd.X0] = msh->xhfull[ix];
}


void set_nc_meta_global_string(int ncid, char *name, char *value)
{
	int status;
	size_t len;
	const void *v;

	v = value;

	len=strlen(value);
	status = nc_put_att_text(ncid,NC_GLOBAL,name,len,v);
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_txt failed");
}

void set_nc_meta_global_integer(int ncid, char *name, int *value)
{
	int status;
	size_t len;
	const void *v;

	v = value;
	len=1;
	status = nc_put_att_int(ncid,NC_GLOBAL,name,NC_INT,len,value);
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_int failed");
}


//set_nc_meta(ncid_g, twodvarid[n2d_hdf2nc],"long_name",description_string_filtered,units_string);
void set_nc_meta_name_units(int ncid, int varnameid, char *lnstring, char *long_name, char *units)
{
	int len,status,i;

	len=strlen(long_name);
	status = nc_put_att_text(ncid, varnameid,lnstring,len,long_name);
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");

	len=strlen(units);
	status = nc_put_att_text(ncid,varnameid, "units",len,units);
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
}

/* basically lifted from HDF5-ZFP code, keep it simple */

#define ZFP_ID 32013

#define H5Z_ZFP_MODE_RATE      1
#define H5Z_ZFP_MODE_PRECISION 2
#define H5Z_ZFP_MODE_ACCURACY  3
#define H5Z_ZFP_MODE_EXPERT    4
#define H5Z_ZFP_MODE_REVERSIBLE 5

#define set_zfp_accuracy_cdata(A, CD) \
{ double *p = (double *) &CD[2];      \
CD[0]=CD[1]=CD[2]=CD[3]=0;            \
CD[0]=H5Z_ZFP_MODE_ACCURACY; *p=A;}

#define set_zfp_lossless(CD) \
{ CD[0]=H5Z_ZFP_MODE_REVERSIBLE; }



//		We are looping over total requested variables (index ivar)
//		nid = nc->ncid;
//		v3did = &(nc->var3d[ivar]);
//else if(same(var,"thrhopert"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->thrhopert,       cmd->zfp,nid,hm,v3did,"long_name",var,"K");
//void set_nc_meta_name_units_compression(double zfpacc_netcdf,int do_zfp,int do_bg1, int do_bg2, int do_bg3, int ncid, hdf_meta *hm, var3dstruct *v3d, char *lnstring, char *long_name, char *units)
void set_nc_meta_name_units_compression(double zfpacc_netcdf,cmdline cmd, int ncid, hdf_meta *hm, var3dstruct *v3d, char *lnstring, char *long_name, char *units)
{
	int len,status;
	int varnameid;
	int do_zfp,do_bg1,do_bg2,do_bg3;
	int i,do_zfp_lossless=0,flag_adjust=0;
	unsigned int cdata[4]; /* for the ZFP stuff */
	char attstr[MAXSTR];
	float zfpacc_netcdf_old;

	do_zfp=cmd.zfp;
	do_bg1=cmd.bitgroom1;
	do_bg2=cmd.bitgroom2;
	do_bg3=cmd.bitgroom3;

	i=0;
	if(do_zfp)i++;
	if(do_bg1)i++;
	if(do_bg2)i++;
	if(do_bg3)i++;
	if(i>1) ERROR_STOP("Check your compression options, you have chosen more than one!");

	if(zfpacc_netcdf < 0.0) do_zfp_lossless=1;

	len=strlen(long_name);
	status = nc_put_att_text(ncid,v3d->varnameid,lnstring,len,long_name);
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");

	len=strlen(units);
	status = nc_put_att_text(ncid,v3d->varnameid, "units",len,units);
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");

/* 2022-08-23 TODO ORF: Simply interpret negative accuracy parameters as requesting
 * ZFP_REVERSIBLE (lossless) like I do with CM1-LOFS. */


	if (v3d->is_LOFS_var)
	{
		if (do_zfp && zfpacc_netcdf>0.0 && (zfpacc_netcdf < v3d->zfpacc_LOFS))
		{
			zfpacc_netcdf_old = zfpacc_netcdf;
			zfpacc_netcdf = v3d->zfpacc_LOFS;
			flag_adjust = 1;
		}
	}


	if (do_zfp_lossless)
	{
		if (!H5Zfilter_avail(ZFP_ID))
		{
			char *hdf5_plugin_path;
			printf("ZFP filter not available!\n");
			hdf5_plugin_path = getenv("HDF5_PLUGIN_PATH");
			printf("Check your HDF5_PLUGIN_PATH; it is currently %s\n",hdf5_plugin_path);
			ERROR_STOP("ZFP filter not available");
		}
		sprintf(attstr,"zfpacc_netcdf_%s",long_name);
		printf("%30s = %14.7f",attstr,zfpacc_netcdf);
		if(flag_adjust) printf(" **** ADJUSTED UPWARDS TO LOFS VALUE (from %f to %f)\n",zfpacc_netcdf_old,zfpacc_netcdf); else printf("\n");

		status = nc_put_att_double(ncid,v3d->varnameid, "zfp_accuracy_netcdf",NC_DOUBLE,1,&zfpacc_netcdf);
		if (status != NC_NOERR) ERROR_STOP("nc_put_att_double failed");

		set_zfp_lossless(cdata);
		status = nc_def_var_filter(ncid,v3d->varnameid,ZFP_ID,4,cdata);
		if(status != NC_NOERR)
		{
			ERROR_STOP("nc_def_var_filter failed");
		}
	}
	else if (do_zfp)
	{
		v3d->zfpacc_netcdf = zfpacc_netcdf; // This is the value we set in this code, not what was saved in LOFS

		sprintf(attstr,"zfpacc_netcdf_%s",long_name);
		printf("%30s = %14.7f",attstr,zfpacc_netcdf);
		if(flag_adjust) printf(" **** ADJUSTED UPWARDS TO LOFS VALUE (from %f to %f)\n",zfpacc_netcdf_old,zfpacc_netcdf); else printf("\n");

		status = nc_put_att_double(ncid,v3d->varnameid, "zfp_accuracy_netcdf",NC_DOUBLE,1,&zfpacc_netcdf);
		if (status != NC_NOERR) ERROR_STOP("nc_put_att_double failed");

		set_zfp_accuracy_cdata(zfpacc_netcdf,cdata);
		if (!H5Zfilter_avail(ZFP_ID))
		{
			char *hdf5_plugin_path;
			printf("ZFP filter not available!\n");
			hdf5_plugin_path = getenv("HDF5_PLUGIN_PATH");
			printf("Check your HDF5_PLUGIN_PATH; it is currently %s\n",hdf5_plugin_path);
			ERROR_STOP("ZFP filter not available");
		}

		status = nc_def_var_filter(ncid,v3d->varnameid,ZFP_ID,4,cdata);
		if(status != NC_NOERR)
		{
			ERROR_STOP("nc_def_var_filter failed");
		}
	}

/* BitGroom, Granular BitRound, and Bit Round
   Requires NetCDF4 4.9.0 or greater
	 
	Selectively zeroes out least significant bits in floating point
	data, and these zeroes can then be gzip compressed efficiently,
	giving us the desired lossy compression ratios that exceed lossless.
	
	
	Documented as of 2023 here: https://docs.unidata.ucar.edu/netcdf-c/current/md_quantize.html
	Important definitions for vars here: https://docs.unidata.ucar.edu/netcdf-c/current/netcdf_8h.html
	
	We have, from the second link:
	
	#define 	NC_QUANTIZE_BITGROOM   1
	#define 	NC_QUANTIZE_GRANULARBR 2
	#define 	NC_QUANTIZE_BITROUND   3
	
	Using the quantize example from the first link above:
	 
		if (nc_def_var_quantize(ncid, varid1, NC_QUANTIZE_BITGROOM, NSD_3)) ERR;

	The last argument is number of significant digits, that is the only knob we can turn.
	Then, gzip or you don't see any compression benefits (the whole point basicaly!!).

	By default we stick with 3 significant digits and gzip level 1. Just change it here rather than muck with 
	more command line options for now.

*/

/* Last option is number of significant digits, but it's interpeted
 * differently by each of the 3 algorithms */

	else if (do_bg1)
	{
//		printf("Choosing BitGroom Option 1\n");
		status = nc_def_var_quantize(ncid,v3d->varnameid,NC_QUANTIZE_BITGROOM, cmd.bitgroom_nsd);
		if(status != NC_NOERR)
		{
			ERROR_STOP("nc_def_var_quantize bg1 failed");
		}
		status = nc_def_var_deflate(ncid, v3d->varnameid, 0, 1, 1);
		if(status != NC_NOERR)
		{
			ERROR_STOP("nc_def_var_deflate bg1 failed");
		}
	}
	else if (do_bg2)
	{
//		printf("Choosing BitGroom Option 2\n");
		status = nc_def_var_quantize(ncid,v3d->varnameid,NC_QUANTIZE_GRANULARBR, cmd.bitgroom_nsd);
		if(status != NC_NOERR)
		{
			ERROR_STOP("nc_def_var_quantize bg2 failed");
		}
		status = nc_def_var_deflate(ncid, v3d->varnameid, 0, 1, 1);
		if(status != NC_NOERR)
		{
			ERROR_STOP("nc_def_var_deflate bg2 failed");
		}
	}
	else if (do_bg3)
	{
//		printf("Choosing BitGroom Option 3\n");
//
//		ORF 2023-02-28: This is the option that is easiest I think to envision in
//		terms of ACTUAL significant digits (after which they are all
//		zeroed out). But you need to crank NSD up for this one. The
//		other ones you can get away with smaller NSD but you get
//		crappier compression ratios! Kind of an apples to oranges
//		comparison, really.
		status = nc_def_var_quantize(ncid,v3d->varnameid,NC_QUANTIZE_BITROUND, cmd.bitgroom_nsd);
		if(status != NC_NOERR)
		{
			ERROR_STOP("nc_def_var_quantize bg3 failed");
		}
		status = nc_def_var_deflate(ncid, v3d->varnameid, 0, 1, 1);
		if(status != NC_NOERR)
		{
			ERROR_STOP("nc_def_var_deflate bg3 failed");
		}
	}
}

int n2d_hdf2nc;
const char **twodvarname_hdf2nc;
int *twodvarid;
int ncid_g;
int d2[3];

herr_t twod_first_pass_hdf2nc(hid_t loc_id, const char *name, void *opdata)
{
    H5G_stat_t statbuf;
    H5Gget_objinfo(loc_id, name, FALSE, &statbuf);
    n2d_hdf2nc++;
    return 0;
}

herr_t twod_second_pass_hdf2nc(hid_t loc_id, const char *name, void *opdata)
{
//loc_id is just f_id, the hdf5 file id 
    char attname[50];//ORF FIX 50
    char **units_string,**description_string,*description_string_filtered,*units_for_nc;
    H5G_stat_t statbuf;
    H5A_info_t ainfo;
    hid_t attr,filetype,sdim,space,memtype,lapl_id,dset;
    hsize_t dims[1];
    herr_t status;
    int i,ndims,stringlen;

    dims[0]=1;
    H5Gget_objinfo(loc_id, name, FALSE, &statbuf);
    strcpy((char *)twodvarname_hdf2nc[n2d_hdf2nc],name);

// The following h5_wank grabs the units and description strings
// We then process the description a bit, as the netcdf long_name
// attribute does not contain spaces or commas

    strcpy(attname,"units");
    dset = H5Dopen(loc_id,name,H5P_DEFAULT);
    attr = H5Aopen(dset,attname,H5P_DEFAULT);
    filetype = H5Aget_type(attr);
    sdim = H5Tget_size(filetype);sdim++;
    space = H5Aget_space(attr);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    units_string = (char **) malloc (dims[0] * sizeof (char *));
    units_string[0] = (char *) malloc (dims[0] * sdim * sizeof(char));
    for (i=1; i<dims[0]; i++) units_string[i] = units_string[0] + i * sdim;
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype,sdim);
    status = H5Aread(attr,memtype,units_string[0]);
    status = H5Aclose(attr);
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(filetype);
    status = H5Tclose(memtype);

    strcpy(attname,"description");
    dset = H5Dopen(loc_id,name,H5P_DEFAULT);
    attr = H5Aopen(dset,attname,H5P_DEFAULT);
    filetype = H5Aget_type(attr);
    sdim = H5Tget_size(filetype);sdim++;
    space = H5Aget_space(attr);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    description_string = (char **) malloc (dims[0] * sizeof (char *));
    description_string[0] = (char *) malloc (dims[0] * sdim * sizeof(char));
    for (i=1; i<dims[0]; i++) description_string[i] = description_string[0] + i * sdim;
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype,sdim);
    status = H5Aread(attr,memtype,description_string[0]);
    status = H5Aclose(attr);
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(filetype);
    status = H5Tclose(memtype);

//end h5wank. begin ncwank.
//ORF: What is with the stringlen exit?? Buffer overflow...

    stringlen=strlen(description_string[0]); if(stringlen > 500) exit(0);
    description_string_filtered = (char *) malloc ((stringlen+1) * sizeof (char));
    for(i=0; i<stringlen; i++) description_string_filtered[i] = (description_string[0][i] == ' ' || description_string[0][i] == ',') ? '_' : description_string[0][i];
    description_string_filtered[stringlen]='\0';
    stringlen=strlen(units_string[0]); if(stringlen > 500) exit(0);
    units_for_nc = (char *) malloc ((stringlen+1) * sizeof (char));
    strcpy(units_for_nc,units_string[0]);
//ORF NEW gzip compress swaths
    nc_def_var (ncid_g, twodvarname_hdf2nc[n2d_hdf2nc], NC_FLOAT, 3, d2, &(twodvarid[n2d_hdf2nc]));
	nc_def_var_deflate(ncid_g,twodvarid[n2d_hdf2nc], 1, 1, 1);
    set_nc_meta_name_units(ncid_g, twodvarid[n2d_hdf2nc],"long_name",description_string_filtered,units_for_nc);

    free(units_string[0]); free(units_string);
    free(description_string[0]); free(description_string);
    free(description_string_filtered);
    free(units_for_nc);

    n2d_hdf2nc++;
    return 0;
}

void set_netcdf_attributes(ncstruct *nc, grid gd, cmdline *cmd, buffers *b, hdf_meta *hm, hid_t *f_id, zfpacc *zfpacc, int argc, char *argv[] )
{
	int status;
	int nv;
	int ivar;
	int i,isLOFS;
	char var[MAXSTR];
	int nid;
	var3dstruct *v3did;

	status = nc_def_dim (nc->ncid, "xh", gd.NX, &nc->nxh_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed"); 
	status = nc_def_dim (nc->ncid, "yh", gd.NY, &nc->nyh_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (nc->ncid, "zh", gd.NZ, &nc->nzh_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	//ORF shave off one point (was NX+1 etc.)
	status = nc_def_dim (nc->ncid, "xf", gd.NX, &nc->nxf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed"); 
	status = nc_def_dim (nc->ncid, "yf", gd.NY, &nc->nyf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (nc->ncid, "zf", gd.NZ, &nc->nzf_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_dim (nc->ncid, "time", NC_UNLIMITED, &nc->time_dimid); if (status != NC_NOERR) ERROR_STOP("nc_def_dim failed");
	status = nc_def_var (nc->ncid, "xh", NC_FLOAT, 1, &nc->nxh_dimid, &nc->xhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (nc->ncid, "yh", NC_FLOAT, 1, &nc->nyh_dimid, &nc->yhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (nc->ncid, "zh", NC_FLOAT, 1, &nc->nzh_dimid, &nc->zhid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (nc->ncid, "xf", NC_FLOAT, 1, &nc->nxf_dimid, &nc->xfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (nc->ncid, "yf", NC_FLOAT, 1, &nc->nyf_dimid, &nc->yfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (nc->ncid, "zf", NC_FLOAT, 1, &nc->nzf_dimid, &nc->zfid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (nc->ncid, "time", NC_DOUBLE, 1, &nc->time_dimid, &nc->timeid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_put_att_text(nc->ncid, nc->xhid, "long_name", strlen("x-coordinate in Cartesian system"), "x-coordinate in Cartesian system");
	status = nc_put_att_text(nc->ncid, nc->xfid, "long_name", strlen("x-coordinate for u in Cartesian system"), "x-coordinate for u in Cartesian system"); if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->xfid, "axis", strlen("X"), "X");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->xhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->xhid, "axis", strlen("X"), "X");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->yhid, "long_name", strlen("y-coordinate in Cartesian system"), "y-coordinate in Cartesian system");
	status = nc_put_att_text(nc->ncid, nc->yfid, "long_name", strlen("y-coordinate for v in Cartesian system"), "y-coordinate for v in Cartesian system"); if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->yfid, "axis", strlen("Y"), "Y");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->yhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->yhid, "axis", strlen("Y"), "Y");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->zhid, "long_name", strlen("z-coordinate in Cartesian system"), "z-coordinate in Cartesian system");
	status = nc_put_att_text(nc->ncid, nc->zfid, "long_name", strlen("z-coordinate for w in Cartesian system"), "z-coordinate for w in Cartesian system"); if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->zfid, "axis", strlen("Z"), "Z");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->zhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->zhid, "axis", strlen("Z"), "Z");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->xfid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->yfid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->zfid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->timeid, "units", strlen("s"), "s");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->timeid, "axis", strlen("T"), "T");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->timeid, "long_name", strlen("time"), "time");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_def_var (nc->ncid, "X0", NC_INT, 0, nc->dims, &nc->x0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (nc->ncid, "Y0", NC_INT, 0, nc->dims, &nc->y0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (nc->ncid, "X1", NC_INT, 0, nc->dims, &nc->x1id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (nc->ncid, "Y1", NC_INT, 0, nc->dims, &nc->y1id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (nc->ncid, "Z0", NC_INT, 0, nc->dims, &nc->z0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_def_var (nc->ncid, "Z1", NC_INT, 0, nc->dims, &nc->z1id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	status = nc_put_att_text(nc->ncid, nc->x0id, "long_name", strlen("westmost grid index from LOFS data"), "westmost grid index from LOFS data"); if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->x1id, "long_name", strlen("eastmost grid index from LOFS data"), "eastmost grid index from LOFS data"); if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->y0id, "long_name", strlen("southmost grid index from LOFS data"), "southmost grid index from LOFS data"); if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->y1id, "long_name", strlen("northmost grid index from LOFS data"), "northmost grid index from LOFS data"); if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->z0id, "long_name", strlen("bottom grid index from LOFS data"), "bottom grid index from LOFS data"); if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->z1id, "long_name", strlen("top grid index from LOFS data"), "top grid index from LOFS data"); if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");

	nc->d2[0] = nc->time_dimid;
	nc->d2[1] = nc->nyh_dimid;
	nc->d2[2] = nc->nxh_dimid;

	//d2 global (to this file) routine for H5Giter
	d2[0]=nc->d2[0];
	d2[1]=nc->d2[1];
	d2[2]=nc->d2[2];

	if (cmd->do_swaths)
	{
		int i2d,bufsize;

		n2d_hdf2nc = 0;
		H5Giterate(*f_id, "/00000/2D/static",NULL,twod_first_pass_hdf2nc,NULL);
		H5Giterate(*f_id, "/00000/2D/swath",NULL,twod_first_pass_hdf2nc,NULL);

		twodvarname_hdf2nc = (const char **)malloc(n2d_hdf2nc*sizeof(char *));
		twodvarid =   (int *)        malloc(n2d_hdf2nc*sizeof(int));

		hm->n2dswaths = n2d_hdf2nc;
		printf("There are %i 2D static/swath fields.\n",n2d_hdf2nc);

		for (i2d=0; i2d<n2d_hdf2nc; i2d++)
		{
			twodvarname_hdf2nc[i2d] = (char *)malloc(50*sizeof(char)); // 50 characters per variable
		}

		n2d_hdf2nc = 0;
		ncid_g=nc->ncid; //Needed for iteration func.
		H5Giterate(*f_id, "/00000/2D/static",NULL,twod_second_pass_hdf2nc,NULL);
		H5Giterate(*f_id, "/00000/2D/swath",NULL,twod_second_pass_hdf2nc,NULL);

		for (i2d=0; i2d<n2d_hdf2nc; i2d++) free((void *)twodvarname_hdf2nc[i2d]); //shaddap compiler
		free(twodvarname_hdf2nc);

		/* And, like magic, we have populated our netcdf id arrays for all the swath slices */
	}


/* I had the best intentions here, but we can not do this. By
 * sorting the variable names alphabetically, but not being able to keep
 * the other items linked, everything gets scrambled when trying to
 * reference things by variable name, which is how this damned code
 * works. I will still check for dupes, but will exit if I find them.
 * This is where I really need a higher level language (C++, Python) so
 * I can 'sort by keys' and still keep everything linked.
 *
 * This means the variables in your netCDF file will be in the order you
 * have them at the command line!
 *
 * Also, I am removing the --allvars option, and checking for duplicates,
 * because it relies on my sort/merge approach that breaks things.
 */

// TODO: Create code that checks for dupes and exits if it finds them.

	cmd->nvar = cmd->nvar_cmdline;


/* Begin block of code that wrecked my brain */

//	if (!cmd->do_allvars)
//	{
//		cmd->nvar = cmd->nvar_cmdline;
//		nv=0; /* liar! */
//	}
//	else nv=hm->nvar_available;
//With faking nvar_available to zero this loop now sorts/uniqs
//all possible variables lists (just command line, just allvars,
//or a combination of the two)
//	if (cmd->nvar !=0||cmd->do_allvars) 
//	{
//		char **varname_tmp;
//		int ndupes = 0;
//		int i,j;

//		varname_tmp = (char **)malloc(MAXVARIABLES * sizeof(char *));
//		for (i=0; i < MAXVARIABLES; i++) varname_tmp[i] = (char *)(malloc(MAXSTR * sizeof(char)));

//		for (i=0; i<nv; i++)
//		{
//			strcpy(varname_tmp[i],hm->varname_available[i]);
//		}
//		for (i=0;i<cmd->nvar_cmdline; i++)
//		{
//			strcpy(varname_tmp[i+nv],cmd->varname_cmdline[i]);
//		}

//		cmd->nvar = cmd->nvar_cmdline;
//		strcpy(nc->var3d[j].varname,varname_tmp[i+1]); //THIS IS WHERE VARNAME IS SET!

//		sortchararray (varname_tmp,nv+cmd->nvar_cmdline); //This right here breaks parity with stuff like zfp_acc and is_LOFS_var

//		strcpy(nc->var3d[0].varname,varname_tmp[0]);
//		j=1;
//		 Get rid of accidental duplicates 
//		for (i=0; i<nv+cmd->nvar_cmdline-1; i++)
//		{
//			if(strcmp(varname_tmp[i],varname_tmp[i+1]))
//			{
//				strcpy(nc->varname[j],varname_tmp[i+1]); //THIS IS WHERE VARNAME IS SET!
//				strcpy(nc->var3d[j].varname,varname_tmp[i+1]); //THIS IS WHERE VARNAME IS SET!
//				j++;
//			}
//		}
//		cmd->nvar = j;
//		ndupes = nv+cmd->nvar_cmdline-cmd->nvar;
//		if(ndupes!=0)printf("We got rid of %i duplicate requested variables\n",ndupes);
//		for (i=0; i < MAXVARIABLES; i++) free(varname_tmp[i]);
//		free(varname_tmp);
//	}


/* End block of code that wrecked my brain */

/* Below is handled within set_nc_meta_name_units_compression now and actually works */

/*
	for (ivar = 0; ivar < cmd->nvar; ivar++)
	{
		double lofsacc,netcdfacc;
		lofsacc = nc->var3d[ivar].zfpacc_LOFS;
		netcdfacc = nc->var3d[ivar].zfpacc_netcdf;
		isLOFS = nc->var3d[ivar].is_LOFS_var;
		strcpy(var,nc->var3d[ivar].varname); //var set here
		if (isLOFS && (netcdfacc < lofsacc))
		{
			printf("WARNING WARNING WARNING: %s is an LOFS variable with zfp_acc = %f, but the netCDF zfp_acc is %f, which is smaller!\n",var,lofsacc,netcdfacc);
			printf("Setting %s zfpacc_netdcf to zfpacc_LOFS (Was: %f. Is now: %f)\n",var,netcdfacc,lofsacc);
			nc->var3d[ivar].zfpacc_netcdf = nc->var3d[ivar].zfpacc_LOFS; //does not work
		}
	}
*/


// This is our main "loop over all requested variable names" loop that
// sets all the metadata shit.
//
// Set some global metadata
	 	set_nc_meta_global_string(nc->ncid,"cm1_lofs_version", "1.0");
		i=1; set_nc_meta_global_integer(nc->ncid,"uniform_mesh",&i);

		{
			int j,k;
			k=0;
			char *cmdstring;
			cmdstring = (char *)(malloc(MAXSTR * sizeof(char)));

			for (i=0; i<argc; i++)
			{
				for (j=0; j<strlen(argv[i]); j++)
					{
						cmdstring[k]=argv[i][j];
						k++;
					}
				cmdstring[k]=' ';k++;
			}
			cmdstring[k-1]='\0';
			if(cmd->debug==1) printf("cmdstring = %s\n",cmdstring);
			status = nc_put_att_text(nc->ncid,NC_GLOBAL,"commandline",strlen(cmdstring),cmdstring);
		}


	for (ivar = 0; ivar < cmd->nvar; ivar++)
	{
		int mesh_is_u = 0, mesh_is_v = 0, mesh_is_w = 0;

		strcpy(var,nc->var3d[ivar].varname); //var set here

		if (same(var,"u")) mesh_is_u=1;
		if (same(var,"v")) mesh_is_v=1;
		if (same(var,"w")) mesh_is_w=1;
		/* ORF these guys live on the w grid actually, time to
		 * get this right. TODO: Have interp option available. */
		if (same(var,"khh")) mesh_is_w=1;
		if (same(var,"khv")) mesh_is_w=1;
		if (same(var,"kmh")) mesh_is_w=1;
		if (same(var,"kmv")) mesh_is_w=1;

		/* ORF 2020-06-25
		 * All of the budget stuff now matched to the right
		 * mesh. */

		if (same(var,"ub_cor"))    mesh_is_u=1;
		if (same(var,"ub_hadv"))   mesh_is_u=1;
		if (same(var,"ub_hediff")) mesh_is_u=1;
		if (same(var,"ub_hidiff")) mesh_is_u=1;
		if (same(var,"ub_hturb"))  mesh_is_u=1;
		if (same(var,"ub_pbl"))    mesh_is_u=1;
		if (same(var,"ub_pgrad"))  mesh_is_u=1;
		if (same(var,"ub_rdamp"))  mesh_is_u=1;
		if (same(var,"ub_subs"))   mesh_is_u=1;
		if (same(var,"ub_vadv"))   mesh_is_u=1;
		if (same(var,"ub_vediff")) mesh_is_u=1;
		if (same(var,"ub_vidiff")) mesh_is_u=1;
		if (same(var,"ub_vturb"))  mesh_is_u=1;

		if (same(var,"vb_cor"))    mesh_is_v=1;
		if (same(var,"vb_hadv"))   mesh_is_v=1;
		if (same(var,"vb_hediff")) mesh_is_v=1;
		if (same(var,"vb_hidiff")) mesh_is_v=1;
		if (same(var,"vb_hturb"))  mesh_is_v=1;
		if (same(var,"vb_pbl"))    mesh_is_v=1;
		if (same(var,"vb_pgrad"))  mesh_is_v=1;
		if (same(var,"vb_rdamp"))  mesh_is_v=1;
		if (same(var,"vb_subs"))   mesh_is_v=1;
		if (same(var,"vb_vadv"))   mesh_is_v=1;
		if (same(var,"vb_vediff")) mesh_is_v=1;
		if (same(var,"vb_vidiff")) mesh_is_v=1;
		if (same(var,"vb_vturb"))  mesh_is_v=1;

		if (same(var,"wb_buoy"))   mesh_is_w=1;
		if (same(var,"wb_hadv"))   mesh_is_w=1;
		if (same(var,"wb_hediff")) mesh_is_w=1;
		if (same(var,"wb_hidiff")) mesh_is_w=1;
		if (same(var,"wb_hturb"))  mesh_is_w=1;
		if (same(var,"wb_pbl"))    mesh_is_w=1;
		if (same(var,"wb_pgrad"))  mesh_is_w=1;
		if (same(var,"wb_rdamp"))  mesh_is_w=1;
		if (same(var,"wb_subs"))   mesh_is_w=1;
		if (same(var,"wb_vadv"))   mesh_is_w=1;
		if (same(var,"wb_vediff")) mesh_is_w=1;
		if (same(var,"wb_vidiff")) mesh_is_w=1;
		if (same(var,"wb_vturb"))  mesh_is_w=1;

		/* u v and w live on their own mesh (Arakawa C grid)*/

		/* Recommend, however, for making netcdf files, to just
		 * request uinterp vinterp winterp which NOW are cacluated
		 * HERE rather than saved in LOFS */

		/* We are going to preserve u v w with the extra point for
		 * saving u v and w easily while facilitating averaging
		 * which requires 1 fewer point */


		/* ORF 2020-04-16
		 * NOTE start/edges that are set here are ignored. I've changed over to writing
		 * netCDF files in Z slices, and set the appropirate start and edges values in
		 * do_requested_variables. However the dimid stuff must still be done here and is
		 * used */

		/* ORF 2020-06
		 * I don't do Z slices any more, performance sucked
		 */

		if(mesh_is_u)
		{
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = nc->nzh_dimid;
			nc->dims[2] = nc->nyh_dimid;
			nc->dims[3] = nc->nxf_dimid;

			nc->start[0] = 0;
			nc->start[1] = 0;
			nc->start[2] = 0;
			nc->start[3] = 0;

			nc->edges[0] = 1;
			nc->edges[1] = gd.NZ;
			nc->edges[2] = gd.NY;
			nc->edges[3] = gd.NX;
		}
		else if (mesh_is_v)
		{
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = nc->nzh_dimid;
			nc->dims[2] = nc->nyf_dimid;
			nc->dims[3] = nc->nxh_dimid;

			nc->start[0] = 0;
			nc->start[1] = 0;
			nc->start[2] = 0;
			nc->start[3] = 0;

			nc->edges[0] = 1;
			nc->edges[1] = gd.NZ;
			nc->edges[2] = gd.NY;
			nc->edges[3] = gd.NX;
		}
		else if (mesh_is_w)
		{
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = nc->nzf_dimid;
			nc->dims[2] = nc->nyh_dimid;
			nc->dims[3] = nc->nxh_dimid;

			nc->start[0] = 0;
			nc->start[1] = 0;
			nc->start[2] = 0;
			nc->start[3] = 0;

			nc->edges[0] = 1;
			nc->edges[1] = gd.NZ;
			nc->edges[2] = gd.NY;
			nc->edges[3] = gd.NX;
		}
		else /* scalar grid */
		{
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = nc->nzh_dimid;
			nc->dims[2] = nc->nyh_dimid;
			nc->dims[3] = nc->nxh_dimid;

			nc->start[0] = 0;
			nc->start[1] = 0;
			nc->start[2] = 0;
			nc->start[3] = 0;

			nc->edges[0] = 1;
			nc->edges[1] = gd.NZ;
			nc->edges[2] = gd.NY;
			nc->edges[3] = gd.NX;
		}

// If we request a 2D slice set parameters accordingly

		if(gd.X0==gd.X1)
		{
			nc->twodslice = TRUE;
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = (mesh_is_w)?nc->nzf_dimid:nc->nzh_dimid;
			nc->dims[2] = (mesh_is_v)?nc->nyf_dimid:nc->nyh_dimid;
			nc->start[0] = 0;
			nc->start[1] = 0;
			nc->start[2] = 0;
			nc->edges[0] = 1;
			nc->edges[1] = gd.NZ;
			nc->edges[2] = gd.NY;
			status = nc_def_var (nc->ncid, nc->var3d[ivar].varname, NC_FLOAT, 3, nc->dims, &(nc->var3d[ivar].varnameid));
		}
		else if(gd.Y0==gd.Y1)
		{
			nc->twodslice = TRUE;
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = (mesh_is_w)?nc->nzf_dimid:nc->nzh_dimid;
			nc->dims[2] = (mesh_is_u)?nc->nxf_dimid:nc->nxh_dimid;
			nc->start[0] = 0;
			nc->start[1] = 0;
			nc->start[2] = 0;
			nc->edges[0] = 1;
			nc->edges[1] = gd.NZ;
			nc->edges[2] = gd.NX;
			status = nc_def_var (nc->ncid, nc->var3d[ivar].varname, NC_FLOAT, 3, nc->dims, &(nc->var3d[ivar].varnameid));
		}
		else if(gd.Z0==gd.Z1)
		{
			nc->twodslice = TRUE;
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = (mesh_is_v)?nc->nyf_dimid:nc->nyh_dimid;
			nc->dims[2] = (mesh_is_u)?nc->nxf_dimid:nc->nxh_dimid;
			nc->start[0] = 0;
			nc->start[1] = 0;
			nc->start[2] = 0;
			nc->edges[0] = 1;
			nc->edges[1] = gd.NY;
			nc->edges[2] = gd.NX;
			status =  nc_def_var (nc->ncid, nc->var3d[ivar].varname, NC_FLOAT, 3, nc->dims, &(nc->var3d[ivar].varnameid));
		}
		else
		{
			//If using ZFP, must make sure chunk dimensions work correctly for the
			//requested data, meaning NX%chunkx==0 etc. If ZFP is chosen
			//even if NX%4==0 you will get bad values at boundaries.
			//TODO: We will give the option of choosing chunking parameters at the command line.
			//We will by default choose chunking parameters as follows.
			//This will always work because we have already forced
			//NX,NY,NZ % 4 to be 0
			//However, you may wish to choose your own frigging HDF chunking parameters bud!
			int xchunk,ychunk,zchunk;

			xchunk=gd.NX/4;
			ychunk=gd.NY/4;
			zchunk=gd.NZ/2;

			size_t chunkdims[4] = {1,zchunk,ychunk,xchunk}; 

			status = nc_def_var (nc->ncid, nc->var3d[ivar].varname, NC_FLOAT, 4, nc->dims, &(nc->var3d[ivar].varnameid));
			if (cmd->zfp==1) status = nc_def_var_chunking(nc->ncid,nc->var3d[ivar].varnameid,NC_CHUNKED,chunkdims);
		}

		if (status != NC_NOERR) 
		{
			printf ("Cannot nc_def_var for var #%i %s, status = %i, message = %s\n", ivar, nc->var3d[ivar].varname,status,nc_strerror(status));
			ERROR_STOP("nc_def_var failed");
		}


// Set all the ZFP LOFS metadata, finally!
//		for (i=0; i<hm->nzfplofs;i++)
//		{
//			char delimiter[]=" ";
//			char *a,*b,*c;
//			status = nc_put_att_string(nc->ncid,NC_GLOBAL,"LOFS_CM1_ZFP_SAVED",hm->nzfplofs,hm->zfpacc_LOFS_all);
//		}

// We are still before the nc_enddef call, in case you are lost

		nid = nc->ncid;
		v3did = &(nc->var3d[ivar]);

// I really hate this section but at some fundamental level, you need to just go step by
// step through each variable and set your metadata - and compression - accordingly

		if(same(var,"u"))				    set_nc_meta_name_units_compression(zfpacc->netcdf->u,               *cmd,nid,hm,v3did,"long_name",var,"m/s");
		else if(same(var,"v"))			    set_nc_meta_name_units_compression(zfpacc->netcdf->v,               *cmd,nid,hm,v3did,"long_name",var,"m/s");
		else if(same(var,"w"))			    set_nc_meta_name_units_compression(zfpacc->netcdf->w,               *cmd,nid,hm,v3did,"long_name",var,"m/s");
		else if(same(var,"hwin_sr"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->hwin_sr,         *cmd,nid,hm,v3did,"long_name",var,"m/s");
		else if(same(var,"hwin_gr"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->hwin_gr,         *cmd,nid,hm,v3did,"long_name",var,"m/s");
		else if(same(var,"u_gr"))	        set_nc_meta_name_units_compression(zfpacc->netcdf->u_gr,            *cmd,nid,hm,v3did,"long_name",var,"m/s");
		else if(same(var,"v_gr"))	        set_nc_meta_name_units_compression(zfpacc->netcdf->v_gr,            *cmd,nid,hm,v3did,"long_name",var,"m/s");
		else if(same(var,"windmag_sr"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->windmag_sr,      *cmd,nid,hm,v3did,"long_name",var,"m/s");
		else if(same(var,"xvort"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->xvort,           *cmd,nid,hm,v3did,"long_name",var,"s^-1");
		else if(same(var,"yvort"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->yvort,           *cmd,nid,hm,v3did,"long_name",var,"s^-1");
		else if(same(var,"zvort"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->zvort,           *cmd,nid,hm,v3did,"long_name",var,"s^-1");
		else if(same(var,"vortmag"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->vortmag,         *cmd,nid,hm,v3did,"long_name",var,"s^-1");
		else if(same(var,"qvpert"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->qvpert,          *cmd,nid,hm,v3did,"long_name",var,"g/kg");
		else if(same(var,"dbz"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->dbz,             *cmd,nid,hm,v3did,"long_name",var,"dBZ");
		else if(same(var,"uinterp"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->uinterp,         *cmd,nid,hm,v3did,"long_name",var,"m/s");
		else if(same(var,"vinterp"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->vinterp,         *cmd,nid,hm,v3did,"long_name",var,"m/s");
		else if(same(var,"winterp"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->winterp,         *cmd,nid,hm,v3did,"long_name",var,"m/s");
		else if(same(var,"thrhopert"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->thrhopert,       *cmd,nid,hm,v3did,"long_name",var,"K");
		else if(same(var,"prespert"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->prespert,        *cmd,nid,hm,v3did,"long_name",var,"hPa");
		else if(same(var,"tke_sg"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->tke_sg,          *cmd,nid,hm,v3did,"long_name",var,"m^2/s^2");
		else if(same(var,"qc"))			    set_nc_meta_name_units_compression(zfpacc->netcdf->qc,              *cmd,nid,hm,v3did,"long_name",var,"g/kg");
		else if(same(var,"qr"))			    set_nc_meta_name_units_compression(zfpacc->netcdf->qr,              *cmd,nid,hm,v3did,"long_name",var,"g/kg");
		else if(same(var,"ncr"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->ncr,             *cmd,nid,hm,v3did,"long_name",var,"cm^-3");
		else if(same(var,"qg"))			    set_nc_meta_name_units_compression(zfpacc->netcdf->qg,              *cmd,nid,hm,v3did,"long_name",var,"g/kg");
		else if(same(var,"qhl"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->qhl,             *cmd,nid,hm,v3did,"long_name",var,"cm^-3");
		else if(same(var,"ncg"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->ncg,             *cmd,nid,hm,v3did,"long_name",var,"cm^-3");
		else if(same(var,"qi"))			    set_nc_meta_name_units_compression(zfpacc->netcdf->qi,              *cmd,nid,hm,v3did,"long_name",var,"g/kg");
		else if(same(var,"nci"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->nci,             *cmd,nid,hm,v3did,"long_name",var,"cm^-3");
		else if(same(var,"qs"))			    set_nc_meta_name_units_compression(zfpacc->netcdf->qs,              *cmd,nid,hm,v3did,"long_name",var,"g/kg");
		else if(same(var,"ncs"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->ncs,             *cmd,nid,hm,v3did,"long_name",var,"cm^-3");
		else if(same(var,"cci"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->cci,             *cmd,nid,hm,v3did,"long_name",var,"cm^-3");
		else if(same(var,"ccn"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->ccn,             *cmd,nid,hm,v3did,"long_name",var,"#");
		else if(same(var,"ccw"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->ccw,             *cmd,nid,hm,v3did,"long_name",var,"cm^-3");
		else if(same(var,"chl"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->chl,             *cmd,nid,hm,v3did,"long_name",var,"cm^-3");
		else if(same(var,"vhl"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->vhl,             *cmd,nid,hm,v3did,"long_name",var,"cm^-3");
		else if(same(var,"vhw"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->vhw,             *cmd,nid,hm,v3did,"long_name",var,"cm^-3");
		else if(same(var,"chw"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->chw,             *cmd,nid,hm,v3did,"long_name",var,"cm^-3");
		else if(same(var,"crw"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->crw,             *cmd,nid,hm,v3did,"long_name",var,"cm^-3");
		else if(same(var,"csw"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->csw,             *cmd,nid,hm,v3did,"long_name",var,"cm^-3");
		else if(same(var,"zhl"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->zhl,             *cmd,nid,hm,v3did,"long_name",var,"Z");
		else if(same(var,"zhw"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->zhw,             *cmd,nid,hm,v3did,"long_name",var,"Z");
		else if(same(var,"zrw"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->zrw,             *cmd,nid,hm,v3did,"long_name",var,"Z");
		else if(same(var,"rho"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->rho,             *cmd,nid,hm,v3did,"long_name",var,"kg/m^3");
		else if(same(var,"qv"))		        set_nc_meta_name_units_compression(zfpacc->netcdf->qv,              *cmd,nid,hm,v3did,"long_name",var,"g/kg");
		else if(same(var,"wb_buoy"))        set_nc_meta_name_units_compression(zfpacc->netcdf->wb_buoy,         *cmd,nid,hm,v3did,"long_name",var,"m/s^2");
		else if(same(var,"wb_buoy_interp")) set_nc_meta_name_units_compression(zfpacc->netcdf->wb_buoy_interp,         *cmd,nid,hm,v3did,"long_name",var,"m/s^2");
		else if(same(var,"ub_pgrad"))       set_nc_meta_name_units_compression(zfpacc->netcdf->ub_pgrad,        *cmd,nid,hm,v3did,"long_name",var,"m/s^2");
		else if(same(var,"vb_pgrad"))       set_nc_meta_name_units_compression(zfpacc->netcdf->vb_pgrad,        *cmd,nid,hm,v3did,"long_name",var,"m/s^2");
		else if(same(var,"wb_pgrad"))       set_nc_meta_name_units_compression(zfpacc->netcdf->wb_pgrad,        *cmd,nid,hm,v3did,"long_name",var,"m/s^2");
		else if(same(var,"ub_pgrad_interp"))set_nc_meta_name_units_compression(zfpacc->netcdf->ub_pgrad_interp,        *cmd,nid,hm,v3did,"long_name",var,"m/s^2");
		else if(same(var,"vb_pgrad_interp"))set_nc_meta_name_units_compression(zfpacc->netcdf->vb_pgrad_interp,        *cmd,nid,hm,v3did,"long_name",var,"m/s^2");
		else if(same(var,"wb_pgrad_interp"))set_nc_meta_name_units_compression(zfpacc->netcdf->wb_pgrad_interp,        *cmd,nid,hm,v3did,"long_name",var,"m/s^2");
		else if(same(var,"pipert"))	        set_nc_meta_name_units_compression(zfpacc->netcdf->pipert,          *cmd,nid,hm,v3did,"long_name",var,"None");
		else if(same(var,"thpert"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->thpert,          *cmd,nid,hm,v3did,"long_name",var,"K");
		else if(same(var,"rhopert"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->rhopert,         *cmd,nid,hm,v3did,"long_name",var,"kg/m^3");
		else if(same(var,"khh"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->khh,              *cmd,nid,hm,v3did,"long_name",var,"m^2/s");
		else if(same(var,"khv"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->khv,              *cmd,nid,hm,v3did,"long_name",var,"m^2/s");
		else if(same(var,"kmh"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->kmh,              *cmd,nid,hm,v3did,"long_name",var,"m^2/s");
		else if(same(var,"kmv"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->kmv,              *cmd,nid,hm,v3did,"long_name",var,"m^2/s");
		else if(same(var,"khh_interp"))		set_nc_meta_name_units_compression(zfpacc->netcdf->khh_interp,              *cmd,nid,hm,v3did,"long_name",var,"m^2/s");
		else if(same(var,"khv_interp"))		set_nc_meta_name_units_compression(zfpacc->netcdf->khv_interp,              *cmd,nid,hm,v3did,"long_name",var,"m^2/s");
		else if(same(var,"kmh_interp"))		set_nc_meta_name_units_compression(zfpacc->netcdf->kmh_interp,              *cmd,nid,hm,v3did,"long_name",var,"m^2/s");
		else if(same(var,"kmv_interp"))		set_nc_meta_name_units_compression(zfpacc->netcdf->kmv_interp,              *cmd,nid,hm,v3did,"long_name",var,"m^2/s");
		else if(same(var,"xvort_stretch"))  set_nc_meta_name_units_compression(zfpacc->netcdf->xvort_stretch,   *cmd,nid,hm,v3did,"long_name",var,"s^-2");
		else if(same(var,"yvort_stretch"))  set_nc_meta_name_units_compression(zfpacc->netcdf->yvort_stretch,   *cmd,nid,hm,v3did,"long_name",var,"s^-2");
		else if(same(var,"zvort_stretch"))  set_nc_meta_name_units_compression(zfpacc->netcdf->zvort_stretch,   *cmd,nid,hm,v3did,"long_name",var,"s^-2");
		else if(same(var,"xvort_tilt"))  set_nc_meta_name_units_compression(zfpacc->netcdf->xvort_tilt,   *cmd,nid,hm,v3did,"long_name",var,"s^-2");
		else if(same(var,"yvort_tilt"))  set_nc_meta_name_units_compression(zfpacc->netcdf->yvort_tilt,   *cmd,nid,hm,v3did,"long_name",var,"s^-2");
		else if(same(var,"zvort_tilt"))  set_nc_meta_name_units_compression(zfpacc->netcdf->zvort_tilt,   *cmd,nid,hm,v3did,"long_name",var,"s^-2");
		else if(same(var,"xvort_baro"))     set_nc_meta_name_units_compression(zfpacc->netcdf->xvort_baro,      *cmd,nid,hm,v3did,"long_name",var,"s^-2");
		else if(same(var,"yvort_baro"))     set_nc_meta_name_units_compression(zfpacc->netcdf->yvort_baro,      *cmd,nid,hm,v3did,"long_name",var,"s^-2");
		else if(same(var,"xvort_solenoid")) set_nc_meta_name_units_compression(zfpacc->netcdf->xvort_solenoid,  *cmd,nid,hm,v3did,"long_name",var,"s^-2");
		else if(same(var,"yvort_solenoid")) set_nc_meta_name_units_compression(zfpacc->netcdf->yvort_solenoid,  *cmd,nid,hm,v3did,"long_name",var,"s^-2");
		else if(same(var,"zvort_solenoid")) set_nc_meta_name_units_compression(zfpacc->netcdf->zvort_solenoid,  *cmd,nid,hm,v3did,"long_name",var,"s^-2");
		else if(same(var,"hvort"))		    set_nc_meta_name_units_compression(zfpacc->netcdf->hvort,           *cmd,nid,hm,v3did,"long_name",var,"s^-1");
		else if(same(var,"streamvort"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->streamvort,      *cmd,nid,hm,v3did,"long_name",var,"s^-1");
		else if(same(var,"qiqvpert"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->qiqvpert,        *cmd,nid,hm,v3did,"long_name",var,"g/kg");
		else if(same(var,"qcqi"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->qcqi,        *cmd,nid,hm,v3did,"long_name",var,"g/kg");
		else if(same(var,"qgqhqr"))	    set_nc_meta_name_units_compression(zfpacc->netcdf->qgqhqr,        *cmd,nid,hm,v3did,"long_name",var,"g/kg");
		else if(same(var,"qtot"))	    	set_nc_meta_name_units_compression(zfpacc->netcdf->qtot,            *cmd,nid,hm,v3did,"long_name",var,"g/kg");
		else if(same(var,"tempC"))	    	set_nc_meta_name_units_compression(zfpacc->netcdf->tempC,           *cmd,nid,hm,v3did,"long_name",var,"degC");
		else if(same(var,"hdiv"))	    	set_nc_meta_name_units_compression(zfpacc->netcdf->hdiv,            *cmd,nid,hm,v3did,"long_name",var,"s^-1");
		else if(same(var,"liutexmag"))	   	set_nc_meta_name_units_compression(zfpacc->netcdf->liutexmag,       *cmd,nid,hm,v3did,"long_name",var,"s^-1");
		else if(same(var,"liutex_x"))	   	set_nc_meta_name_units_compression(zfpacc->netcdf->liutex_x,       *cmd,nid,hm,v3did,"long_name",var,"s^-1");
		else if(same(var,"liutex_y"))	   	set_nc_meta_name_units_compression(zfpacc->netcdf->liutex_y,       *cmd,nid,hm,v3did,"long_name",var,"s^-1");
		else if(same(var,"liutex_z"))	   	set_nc_meta_name_units_compression(zfpacc->netcdf->liutex_z,       *cmd,nid,hm,v3did,"long_name",var,"s^-1");

	} // End of big ivar loop


//	exit(0);

/* Write sounding data attributes, also umove and vmove */

	status = nc_def_var (nc->ncid, "u0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->u0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta_name_units(nc->ncid,nc->u0id,"long_name","base_state_u","m/s");
	status = nc_def_var (nc->ncid, "v0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->v0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta_name_units(nc->ncid,nc->v0id,"long_name","base_state_v","m/s");
	status = nc_def_var (nc->ncid, "th0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->th0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta_name_units(nc->ncid,nc->th0id,"long_name","base_state_potential_temperature","K");
	status = nc_def_var (nc->ncid, "thv0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->thv0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta_name_units(nc->ncid,nc->thv0id,"long_name","base_state_virtual_potential_temperature","K");
	status = nc_def_var (nc->ncid, "pres0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->pres0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta_name_units(nc->ncid,nc->pres0id,"long_name","base_state_pressure","hPa");
	status = nc_def_var (nc->ncid, "pi0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->pi0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta_name_units(nc->ncid,nc->pi0id,"long_name","base_state_nondimensional_pressure","#");
	status = nc_def_var (nc->ncid, "qv0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->qv0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta_name_units(nc->ncid,nc->qv0id,"long_name","base_state_water_vapor_mixing_ratio","g/kg");
	status = nc_def_var (nc->ncid, "rho0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->rho0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta_name_units(nc->ncid,nc->rho0id,"long_name","base_state_density","kg/m^3");

	status = nc_def_var (nc->ncid, "umove", NC_FLOAT, 0, nc->dims, &nc->umoveid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta_name_units(nc->ncid,nc->umoveid,"long_name","x_box_translation_component","m/s");
	status = nc_def_var (nc->ncid, "vmove", NC_FLOAT, 0, nc->dims, &nc->vmoveid); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta_name_units(nc->ncid,nc->vmoveid,"long_name","y_box_translation_component","m/s");
}

void nc_write_1d_data (ncstruct nc, grid gd, mesh msh, sounding snd, cmdline cmd)
{
	double timearray[1];
	float *tmparray;
	//ORF for writing single time in unlimited time dimension/variable
	const size_t timestart = 0;
	const size_t timecount = 1;
	int status;
	int i;

	/* I do not want all MKS units in my netcdf files. But I set all
	 * the arrays to MKS units for consistency. Here I make my
	 * adjustments - mesh in km, mixing ratios in g/kg */

	/* TODO: Make command line option */

	tmparray = (float *)malloc (gd.NX*sizeof(float));
	for (i=gd.X0; i<=gd.X1; i++) tmparray[i-gd.X0] = 0.001*msh.xhout[i-gd.X0];
	status = nc_put_var_float (nc.ncid,nc.xhid,tmparray); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	for (i=gd.X0; i<=gd.X1; i++) tmparray[i-gd.X0] = 0.001*msh.xfout[i-gd.X0];
	status = nc_put_var_float (nc.ncid,nc.xfid,tmparray); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	free(tmparray);
	tmparray = (float *)malloc (gd.NY*sizeof(float));
	for (i=gd.Y0; i<=gd.Y1; i++) tmparray[i-gd.Y0] = 0.001*msh.yhout[i-gd.Y0];
	status = nc_put_var_float (nc.ncid,nc.yhid,tmparray); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	for (i=gd.Y0; i<=gd.Y1; i++) tmparray[i-gd.Y0] = 0.001*msh.yfout[i-gd.Y0];
	status = nc_put_var_float (nc.ncid,nc.yfid,tmparray); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	free(tmparray);
	tmparray = (float *)malloc (gd.NZ*sizeof(float));
	for (i=gd.Z0; i<=gd.Z1; i++) tmparray[i-gd.Z0] = 0.001*msh.zhout[i-gd.Z0];
	status = nc_put_var_float (nc.ncid,nc.zhid,tmparray); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	for (i=gd.Z0; i<=gd.Z1; i++) tmparray[i-gd.Z0] = 0.001*msh.zfout[i-gd.Z0];
	status = nc_put_var_float (nc.ncid,nc.zfid,tmparray); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	for (i=gd.Z0; i<=gd.Z1; i++) tmparray[i-gd.Z0] = 0.01*snd.pres0[i-gd.Z0];
	status = nc_put_var_float (nc.ncid,nc.pres0id,tmparray); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	for (i=gd.Z0; i<=gd.Z1; i++) tmparray[i-gd.Z0] = 1000.0*snd.qv0[i-gd.Z0];
	status = nc_put_var_float (nc.ncid,nc.qv0id,tmparray); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	free(tmparray);

//	status = nc_put_var_float (nc.ncid,nc.xhid,msh.xhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
//	status = nc_put_var_float (nc.ncid,nc.yhid,msh.yhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
//	status = nc_put_var_float (nc.ncid,nc.zhid,msh.zhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
//	status = nc_put_var_float (nc.ncid,nc.xfid,msh.xfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
//	status = nc_put_var_float (nc.ncid,nc.yfid,msh.yfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
//	status = nc_put_var_float (nc.ncid,nc.zfid,msh.zfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
//	status = nc_put_var_float (nc.ncid,nc.pres0id,snd.pres0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
//	status = nc_put_var_float (nc.ncid,nc.qv0id,snd.qv0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.u0id,snd.u0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.v0id,snd.v0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.th0id,snd.th0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.thv0id,snd.thv0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.pi0id,snd.pi0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.rho0id,snd.rho0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_int (nc.ncid,nc.x0id,&gd.X0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	status = nc_put_var_int (nc.ncid,nc.y0id,&gd.Y0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	status = nc_put_var_int (nc.ncid,nc.x1id,&gd.X1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	status = nc_put_var_int (nc.ncid,nc.y1id,&gd.Y1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	status = nc_put_var_int (nc.ncid,nc.z0id,&gd.Z0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	status = nc_put_var_int (nc.ncid,nc.z1id,&gd.Z1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	status = nc_put_var_float (nc.ncid,nc.umoveid,&msh.umove); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.vmoveid,&msh.vmove); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	timearray[0] = cmd.time;
	status = nc_put_vara_double (nc.ncid,nc.timeid,&timestart,&timecount,timearray);
	if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
}

void set_readahead(readahead *rh,ncstruct nc, cmdline cmd)
{
	int ivar;
	char *var;

/*
 *
 *
 * In this routine we set flags for "readahead" variables (actual LOFS
 * variables that we read before doing calculations, and read into
 * buffers named after the actual LOFS variable). This is for doing
 * memory management. We only want to malloc what is the minimum for
 * doing our calculations because a lot of these files are BIG.
 *
 * So for u, v, w, thrhopert, and prespert, they are true "readahead"
 * variables.
 *
 * FOR THE REST: They are not LOFS variables but we need generic temporary
 * arrays for caluclations. Whereas u,v,w,thrhopert and prespert will
 * NEVER be re-used, the other ones CAN be re-used for each variable.
 * Some calculations require more temporary arrays than others.
 *
 * I should probably have a separate structure for "dovar" as opposed to
 * "readahead" but the readahead structure is pretty entangled and it's
 * just not worth it. So long as the memory management works (and boy
 * you'll know when it doesn't) and doesn't allocate tons of unused
 * memory, we are good.
 *
 */ 

	var = (char *) malloc (MAXSTR * sizeof (char));

	for (ivar = 0; ivar < cmd.nvar; ivar++)
	{
		var=nc.var3d[ivar].varname;

		if(same(var,"pipert"))  {rh->ppert=1;}
		if(same(var,"wb_buoy")) {rh->thrhopert=1; rh->budgets=1;}
		if(same(var,"wb_buoy_interp")) {rh->thrhopert=1; rh->budgets=1;}
		if(same(var,"ub_pgrad")) {rh->ppert=1; rh->thrhopert=1; rh->budgets=1;}
		if(same(var,"vb_pgrad")) {rh->ppert=1; rh->thrhopert=1; rh->budgets=1;}
		if(same(var,"wb_pgrad")) {rh->ppert=1; rh->thrhopert=1; rh->budgets=1;}
		if(same(var,"ub_pgrad_interp")) {rh->ppert=1; rh->thrhopert=1; rh->budgets=1;}
		if(same(var,"vb_pgrad_interp")) {rh->ppert=1; rh->thrhopert=1; rh->budgets=1;}
		if(same(var,"wb_pgrad_interp")) {rh->ppert=1; rh->thrhopert=1; rh->budgets=1;}
//		if(same(var,"uinterp")) {rh->interp=1;}
//		if(same(var,"vinterp")) {rh->interp=1;}
//		if(same(var,"winterp")) {rh->interp=1;}
		if(same(var,"hwin_sr")) {rh->u=1;rh->v=1;}
		if(same(var,"hwin_gr")) {rh->u=1;rh->v=1;}
		if(same(var,"u_gr")) {rh->u=1;}
		if(same(var,"v_gr")) {rh->v=1;}
		if(same(var,"windmag_sr")) {rh->u=1;rh->v=1;rh->w=1;}
		if(same(var,"hdiv")) {rh->u=1;rh->v=1;}
		if(same(var,"xvort")) {rh->v=1;rh->w=1;}
		if(same(var,"yvort")) {rh->u=1;rh->w=1;}
		if(same(var,"zvort")) {rh->u=1;rh->v=1;}
		if(same(var,"xvort_stretch")) {rh->v=1;rh->w=1;rh->budgets=1;}
		if(same(var,"xvort_tilt")) {rh->u=1;rh->v=1;rh->w=1;rh->budgets=1;}
		if(same(var,"yvort_stretch")) {rh->u=1;rh->w=1;rh->budgets=1;}
		if(same(var,"yvort_tilt")) {rh->u=1;rh->v=1;rh->w=1;rh->budgets=1;}
		if(same(var,"zvort_stretch")) {rh->u=1;rh->v=1;rh->budgets=1;}
		if(same(var,"zvort_tilt")) {rh->u=1;rh->v=1;rh->w=1;rh->budgets=1;}
		if(same(var,"xvort_baro")) {rh->thrhopert=1;rh->budgets=1;}
		if(same(var,"yvort_baro")) {rh->thrhopert=1;rh->budgets=1;}
		if(same(var,"xvort_solenoid")) {rh->ppert=1;rh->thrhopert=1;rh->budgets=1;}
		if(same(var,"yvort_solenoid")) {rh->ppert=1;rh->thrhopert=1;rh->budgets=1;}
		if(same(var,"zvort_solenoid")) {rh->ppert=1;rh->thrhopert=1;rh->budgets=1;}
		if(same(var,"hvort")) {rh->u=1;rh->v=1;rh->w=1;rh->hvort=1;}
		if(same(var,"vortmag")) {rh->u=1;rh->v=1;rh->w=1;rh->vortmag=1;}
		if(same(var,"streamvort")) {rh->u=1;rh->v=1;rh->w=1;rh->streamvort=1;}
		if(same(var,"qiqvpert")) {rh->qiqvpert=1;}
		if(same(var,"qcqi")) {rh->qcqi=1;}
		if(same(var,"qgqhqr")) {rh->qgqhqr=1;}
		if(same(var,"tempC")) {rh->tempC=1;}
		if(same(var,"liutexmag")) {rh->u=1;rh->v=1;rh->w=1;}
		if(same(var,"liutex_x")) {rh->u=1;rh->v=1;rh->w=1;}
		if(same(var,"liutex_y")) {rh->u=1;rh->v=1;rh->w=1;}
		if(same(var,"liutex_z")) {rh->u=1;rh->v=1;rh->w=1;}
	}
	//free(var);
}

void malloc_3D_arrays (buffers *b, grid gd, readahead rh,cmdline cmd)
{
	if (cmd.nvar>0)
	{
		long bufsize,bswrite,totbufsize;
		int ibuf=0;

		printf("***********malloc bufsize X = %i bufsize Y = %i bufsize Z = %i tot=%i\n",gd.NX+2,gd.NY+2, gd.NZ+1,((gd.NX+2)*(gd.NY+2)*(gd.NZ+1)));
		bufsize = (long) (gd.NX+2) * (long) (gd.NY+2) * (long) (gd.NZ+1) * (long) sizeof(float);
		bswrite = (long) (gd.NX) * (long) (gd.NY) * (long) (gd.NZ) * (long) sizeof(float);
		totbufsize = bufsize;

		/* We always allocate 1 array (buf0) with bufsize bytes and 1
		 * array (threedbuf) with bswrite bytes. bswrite is the size we
		 * want. threedbuf is used ONLY ONCE being filled at the last
		 * moment before being written to disk. The other array is what
		 * our LOFS variable, at the very least, is initially read into.
		 * The arrays are larger than bswrite because of the
		 * calculations involving derivatives that we must often do */

		if(cmd.verbose)printf("b->buf0: Attempting to allocate %6.2f GB of memory...\n",1.0e-9*bufsize);
		if ((b->buf0 = b->buf = (float *) malloc ((size_t)bufsize)) == NULL)
			ERROR_STOP("Cannot allocate our 3D variable buffer array");
		ibuf++;

		if(cmd.verbose)printf("b->threedbuf: Attempting to allocate %6.2f GB of memory...\n",1.0e-9*bufsize);
		if ((b->threedbuf = (float *) malloc ((size_t)bswrite)) == NULL)
			ERROR_STOP("Cannot allocate our 3D variable write array");
		totbufsize+=bswrite;
		ibuf++;
		if (rh.ppert)
		{
			if(cmd.verbose)printf("b->ppert: Attempting to allocate %6.2f GB of memory...\n",1.0e-9*bufsize);
			if((b->ppert = (float *) malloc ((size_t)bufsize)) == NULL)
				ERROR_STOP("Cannot allocate our prespert/pipert buffer array");
			totbufsize+=bufsize;
			ibuf++;
		}
		if (rh.thrhopert)
		{
			if(cmd.verbose)printf("b->thrhopert: Attempting to allocate %6.2f GB of memory...\n",1.0e-9*bufsize);
			if((b->thrhopert = (float *) malloc ((size_t)bufsize)) == NULL)
				ERROR_STOP("Cannot allocate our thrhopert buffer array");
			totbufsize+=bufsize;
			ibuf++;
		}
		if (rh.u)
		{
			if(cmd.verbose)printf("b->ustag: Attempting to allocate %6.2f GB of memory...\n",1.0e-9*bufsize);
			if ((b->ustag = (float *) malloc ((size_t)bufsize)) == NULL)
				ERROR_STOP("Cannot allocate our ustag buffer array");
			totbufsize+=bufsize;
			ibuf++;
		}
		if (rh.v)
		{
			if(cmd.verbose)printf("b->vstag: Attempting to allocate %6.2f GB of memory...\n",1.0e-9*bufsize);
			if ((b->vstag = (float *) malloc ((size_t)bufsize)) == NULL)
				ERROR_STOP("Cannot allocate our vstag buffer array");
			totbufsize+=bufsize;
			ibuf++;
		}
		if (rh.w)
		{
			if(cmd.verbose)printf("b->wstag: Attempting to allocate %6.2f GB of memory...\n",1.0e-9*bufsize);
			if ((b->wstag = (float *) malloc ((size_t)bufsize)) == NULL)
				ERROR_STOP("Cannot allocate our wstag buffer array");
			totbufsize+=bufsize;
			ibuf++;
		}
		if (rh.vortmag||rh.hvort||rh.streamvort||rh.budgets||rh.qiqvpert||rh.qtot||rh.qcqi||rh.qgqhqr||rh.tempC)//Not really readahead, but if we calculated these we need another array
		{
			if(cmd.verbose)printf("b->dum0: Attempting to allocate %6.2f GB of memory...\n",1.0e-9*bufsize);
			if ((b->dum0 = (float *) malloc ((size_t)bufsize)) == NULL)
				ERROR_STOP("Cannot allocate our first 3D temp calculation array");
			totbufsize+=bufsize;
			ibuf++;
		}
		if (rh.vortmag||rh.hvort||rh.streamvort||rh.budgets||rh.tempC)//Not really readahead, but if we calculated these we need another array
		{
			if(cmd.verbose)printf("b->dum1: Attempting to allocate %6.2f GB of memory...\n",1.0e-9*bufsize);
			if ((b->dum1 = (float *) malloc ((size_t)bufsize)) == NULL)
				ERROR_STOP("Cannot allocate our second 3D temp calculation array");
			totbufsize+=bufsize;
			ibuf++;
		}
		printf("We allocated a total of %6.2f GB of memory for %i floating point buffers\n",1.0e-9*totbufsize,ibuf);
	}
}

void free_3D_arrays (buffers *b, grid gd, readahead rh,cmdline cmd)
{
	if (cmd.nvar>0)
	{
		free (b->buf);
		if(!cmd.twodwrite) free (b->threedbuf);
		if (rh.ppert) free (b->ppert);
		if (rh.thrhopert) free (b->thrhopert);
		if (rh.u) free (b->ustag);
		if (rh.v) free (b->vstag);
		if (rh.w) free (b->wstag);
		if (rh.vortmag||rh.hvort||rh.streamvort||rh.budgets||rh.qiqvpert||rh.qtot||rh.qcqi||rh.qgqhqr||rh.tempC) free(b->dum0);
		if (rh.vortmag||rh.hvort||rh.streamvort||rh.budgets||rh.tempC) free(b->dum1);
		//TODO more checks required here.  We want to be absolutely
		//sure to free all memory before doing external compression,
		//should we choose the --gzip option
	}
}

void do_the_swaths(hdf_meta hm, ncstruct nc, dir_meta dm, grid gd, cmdline cmd)
{
	int i2d,ix,iy,status;
	float *twodfield,*swathbuf,*writeptr;
	float bufsize;
	requested_cube rc;
	size_t s2[3] = {0,0,0};
	size_t e2[3] = {1,gd.NY,gd.NX};

	copy_grid_to_requested_cube(&rc,gd);

	printf("swaths: reading...");FL; 

	bufsize = (long) (rc.NX) * (long) (rc.NY) * (long) sizeof(float);
	if ((twodfield = (float *) malloc ((size_t)bufsize)) == NULL)
		ERROR_STOP("Cannot allocate our 2D swath buffer array");
	bufsize = (long) (rc.NX) * (long) (rc.NY) * (long) (hm.n2dswaths) * (long) sizeof(float);
	if ((swathbuf = (float *) malloc ((size_t)bufsize)) == NULL)
		ERROR_STOP("Cannot allocate our 3D swaths buffer array");

	read_lofs_buffer(swathbuf,"swaths",dm,hm,rc,cmd);

	printf("writing...");FL;
	for (i2d=0;i2d<hm.n2dswaths;i2d++)
	{
		for (iy=0; iy<rc.NY; iy++)
			for (ix=0; ix<rc.NX; ix++)
				twodfield[P2(ix,iy,rc.NX)] = swathbuf[P3(ix,iy,i2d,rc.NX,rc.NY)];
		writeptr = twodfield;
		status = nc_put_vara_float (nc.ncid, twodvarid[i2d], s2, e2, writeptr);
	}
	free(swathbuf);
	free(twodfield);
	BL;
}

void do_readahead(buffers *b,grid gd,readahead rh,dir_meta dm,hdf_meta hm,cmdline cmd)
{
	requested_cube rc;

	/* We select extra data for doing spatial averaging */
	/* By shrinking in saved_x0,saved_x1 etc by 1 point on either side (see set_span), this will not fail
	 * if we do not specify X0, X1 etc.*/

	if (rh.ppert)
	{
		rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
		rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1+1;
		rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;
		printf("readahead: prespert...");
		read_lofs_buffer(b->ppert,"prespert",dm,hm,rc,cmd);
		BL;
	}
	if (rh.thrhopert)
	{
		rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
		rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1+1;
		rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;
		printf("readahead: thrhopert...");
		read_lofs_buffer(b->thrhopert,"thrhopert",dm,hm,rc,cmd);
//ORF TEMPORARY FOR MICROBURST
//      printf("MICROBURST: SUBSTITUTE THPERT FOR THRHOPERT\n");
//		read_lofs_buffer(b->thrhopert,"thpert",dm,hm,rc,cmd);
		BL;
	}
	if (rh.u)
	{
		rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
		rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1+1;
		rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;
		printf("readahead: u...");
		read_lofs_buffer(b->ustag,"u",dm,hm,rc,cmd);
		BL;
	}
	if (rh.v)
	{
		rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
		rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1+1;
		rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;
		printf("readahead: v...");
		read_lofs_buffer(b->vstag,"v",dm,hm,rc,cmd);
		BL;
	}
	if (rh.w)
	{
		rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
		rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1+1;
		rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;
		printf("readahead: w...");
		read_lofs_buffer(b->wstag,"w",dm,hm,rc,cmd);
		BL;
	}

}

// It is much, much faster to do this externally and also has the benefit of compressing
// all arrays, not just the 3D ones. So for the best performance, just make sure nccopy is
// in your path!

void compress_with_nccopy(ncstruct nc,cmdline cmd)
{
	char strbuf[MAXSTR];
	off_t unc_fsize,comp_fsize;
	float ratio;
	struct stat st;
	int retval;

	printf("compressing...\n");fflush(stdout);
	sprintf(strbuf,"mv %s %s.uncompressed",nc.ncfilename,nc.ncfilename);
	retval=system(strbuf);
	if(retval!=0)
	{
		fprintf(stderr,"Command: %s ...failed!\n",strbuf);
		return;
	}
	sprintf(strbuf,"nccopy -d9 -s %s.uncompressed %s",nc.ncfilename,nc.ncfilename);
	if(cmd.verbose)
	{
		printf("Calling external compression program: executing:\n %s ...\n",strbuf);
		fflush(stdout);
	}
	retval=system(strbuf);
	if(retval!=0)
	{
		fprintf(stderr,"Command: %s ...failed!\n",strbuf);
		//nccopy failed so move file name back to .nc, warn and bail
		sprintf(strbuf,"mv %s.uncompressed %s",nc.ncfilename,nc.ncfilename);
		retval=system(strbuf);
		fprintf(stderr,"%s will not be compressed.\n",nc.ncfilename);
		return;
	}

//We have compressed the file. Get file savings info.

	sprintf(strbuf,"%s.uncompressed",nc.ncfilename);
	stat(strbuf,&st);
	unc_fsize=st.st_size;

	sprintf(strbuf,"%s",nc.ncfilename);
	stat(strbuf,&st);
	comp_fsize=st.st_size;

	ratio = (float)unc_fsize/(float)comp_fsize;

	sprintf(strbuf,"\n%12li bytes %s.uncompressed\n%12li bytes %s\nFile compression ratio of %7.2f:1\n",
			unc_fsize,nc.ncfilename,comp_fsize,nc.ncfilename,ratio);
	printf("%s",strbuf);

	sprintf(strbuf,"rm %s.uncompressed",nc.ncfilename);
	retval=system(strbuf);
	if(retval!=0)
	{
		fprintf(stderr,"Command: %s ...failed!\n",strbuf);
		return;
	}
}


void add_CM1_LOFS_zfp_metadata_to_netcdf_file (hdf_meta *hm, hid_t *f_id, ncstruct nc)
{
	hid_t g_id,d_id,attr_id,attr_memtype;
	herr_t status;
	H5O_info_t dset_info;
	H5G_info_t group_info;
	int i,j,k,nattr;
	double zval;
	char groupname[MAXSTR];
	char attrname[MAXATTR][MAXSTR]; //ORF FIX TODO
	char attstr[MAXSTR];
	htri_t existence;
	hid_t lapl_id;

	sprintf(groupname,"%05i/3D",0);//All vars in 3D group 00000 are always available
	g_id = H5Gopen(*f_id,groupname,H5P_DEFAULT);
	H5Gget_info(g_id,&group_info);
	hm->nvar_available = group_info.nlinks;
//	printf("nvar avail = %i\n",hm->nvar_available);
//	printf("g_id = %i varname = %s\n",g_id,hm->varname_available[0]);

	k=0;
	for (i = 0; i < hm->nvar_available; i++)
	{
		if ((d_id = H5Dopen (g_id,hm->varname_available[i],H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
		if ((status = H5Oget_info(d_id,&dset_info,H5O_INFO_NUM_ATTRS)) < 0) ERROR_STOP("Could not H5Oget_info");
		nattr=dset_info.num_attrs;
//		printf("nattr = %i\n",nattr);
		for (j=0; j<nattr; j++)
		{
			if ((status=H5Aget_name_by_idx(d_id,".",H5_INDEX_NAME,H5_ITER_INC,j,attrname[j],40,H5P_DEFAULT))<0) ERROR_STOP("Could not H5Aget_name_by_idex");
			if ((attr_id=H5Aopen(d_id,attrname[j],H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Aopen");
//			printf("attrname[%i] = %s\n",j,attrname[j]);
			if (same(attrname[j],"zfp_accuracy"))
			{
				if ((attr_memtype=H5Aget_type(attr_id)) < 0) ERROR_STOP("Could not H5Aget_type");
				if ((status=H5Aread(attr_id,attr_memtype,&zval)) < 0) ERROR_STOP("Could not H5Aread"); //grab ZFP accuracy from LOFS 3dvar
//				printf("LOFS variable %10s: zfpacc_LOFS = %14.7f\n",hm->varname_available[i],zval);
				sprintf(attstr,"zfpacc_LOFS_%s",hm->varname_available[i]);
				sprintf(hm->zfpacc_LOFS_all[k],"%30s = %14.7f\n",attstr,zval);
				status = nc_put_att_double(nc.ncid,NC_GLOBAL,attstr,NC_DOUBLE,1,&zval);
				if (status != NC_NOERR) ERROR_STOP("nc_put_att_double failed");
				k++;
				break;
			}
		}
		H5Dclose(d_id);
	}
	hm->nzfplofs=k;
}
