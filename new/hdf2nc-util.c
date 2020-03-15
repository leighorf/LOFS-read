#include "include/dirstruct.h"
#include "include/limits.h"
#include "include/hdf2nc.h"
#include "include/lofs-read.h"
#include "include/macros.h"

void init_structs(cmdline *cmd,dir_meta *dm, grid *gd,ncstruct *nc, readahead *rh)
{
	int i;

	cmd->histpath        = (char *) malloc(MAXSTR*sizeof(char));
	cmd->base            = (char *) malloc(MAXSTR*sizeof(char));
	dm-> firstfilename   = (char *) malloc(MAXSTR*sizeof(char));
	dm-> saved_base      = (char *) malloc(MAXSTR*sizeof(char));
	dm-> topdir          = (char *) malloc(MAXSTR*sizeof(char));
	nc->ncfilename       = (char *) malloc(MAXSTR*sizeof(char));
	nc->varnameid        = (int  *) malloc(MAXVARIABLES*sizeof(int));

	dm-> regenerate_cache = 0;

	gd->X0=gd->Y0=gd->X1=gd->Y1=gd->Z0=gd->Z1=-1;
	gd->saved_X0=gd->saved_X1=0;
	gd->saved_Y0=gd->saved_Y1=0;
	gd->saved_Z0=gd->saved_Z1=0;
	cmd->time=0.0; cmd->got_base=0; cmd->optcount=0;
	cmd->debug=0;
	cmd->verbose=0;
	cmd->do_allvars=0;
	cmd->use_box_offset=0;
	cmd->use_interp=0;
	cmd->do_swaths=0;
	nc->twodslice=0;

	nc->varname = (char **)malloc(MAXVARIABLES * sizeof(char *));
	for (i=0; i < MAXVARIABLES; i++) nc->varname[i] = (char *)(malloc(MAXSTR * sizeof(char)));

	rh->u=0; rh->v=0; rh->w=0;
	rh->vortmag=0; rh->hvort=0; rh->streamvort=0;//Not really readahead, used for mallocs
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

void set_1d_arrays(hdf_meta hm, grid gd, mesh *msh, sounding *snd, hid_t *f_id)
{
	int ix,iy,iz,k;

	msh->xhfull = (float *)malloc(hm.nx * sizeof(float));
	msh->yhfull = (float *)malloc(hm.ny * sizeof(float));
	msh->xffull = (float *)malloc((hm.nx+1) * sizeof(float));
	msh->yffull = (float *)malloc((hm.ny+1) * sizeof(float));
	msh->zh = (float *)malloc((hm.nz) * sizeof(float));
	msh->zf = (float *)malloc((hm.nz) * sizeof(float));
	snd->th0 = (float *)malloc(gd.NZ * sizeof(float));
	snd->qv0 = (float *)malloc(gd.NZ * sizeof(float));
	snd->u0 = (float *)malloc(gd.NZ * sizeof(float));
	snd->v0 = (float *)malloc(gd.NZ * sizeof(float));
	snd->pres0 = (float *)malloc(gd.NZ * sizeof(float));
	snd->pi0 = (float *)malloc(gd.NZ * sizeof(float));

	get0dfloat (*f_id,(char *)"mesh/dx",&msh->dx); //rdx=1.0/msh->dx;
	get0dfloat (*f_id,(char *)"mesh/dy",&msh->dy); //rdy=1.0/msh->dy;
	msh->dz=msh->dx;//rdz=rdx; //BAD BAD BAD
	get1dfloat (*f_id,(char *)"mesh/xhfull",msh->xhfull,0,hm.nx);
	get1dfloat (*f_id,(char *)"mesh/yhfull",msh->yhfull,0,hm.ny);
	get1dfloat (*f_id,(char *)"mesh/xffull",msh->xffull,0,hm.nx+1);
	get1dfloat (*f_id,(char *)"mesh/yffull",msh->yffull,0,hm.ny+1);
	get1dfloat (*f_id,(char *)"mesh/zh",msh->zh,0,hm.nz);
	get1dfloat (*f_id,(char *)"mesh/zf",msh->zf,0,hm.nz);
	get1dfloat (*f_id,(char *)"basestate/qv0",snd->qv0,gd.Z0,gd.NZ);
	for (k=0; k<gd.NZ; k++) snd->qv0[k] *= 1000.0; // g/kg now
	get1dfloat (*f_id,(char *)"basestate/th0",snd->th0,gd.Z0,gd.NZ);
	get1dfloat (*f_id,(char *)"basestate/u0",snd->u0,gd.Z0,gd.NZ);
	get1dfloat (*f_id,(char *)"basestate/v0",snd->v0,gd.Z0,gd.NZ);
	get1dfloat (*f_id,(char *)"basestate/pres0",snd->pres0,gd.Z0,gd.NZ);
	get1dfloat (*f_id,(char *)"basestate/pi0",snd->pi0,gd.Z0,gd.NZ);

	msh->xhout = (float *)malloc(gd.NX * sizeof(float));
	msh->yhout = (float *)malloc(gd.NY * sizeof(float));
	msh->zhout = (float *)malloc(gd.NZ * sizeof(float));

	msh->xfout = (float *)malloc(gd.NX * sizeof(float));
	msh->yfout = (float *)malloc(gd.NY * sizeof(float));
	msh->zfout = (float *)malloc(gd.NZ * sizeof(float));

// We recreate George's mesh/derivative calculation paradigm even though
// we are usually isotropic. We need to have our code here match what
// CM1 does internally for stretched and isotropic meshes.
//
// Becuase C cannot do have negative array indices (i.e., uh[-1]) like
// F90 can, we have to offset everything to keep the same CM1-like code
// We malloc enough space for the "ghost zones" and then make sure we
// offset by the correct amount on each side. The macros take care of
// the offsetting.

	msh->uh = (float *)malloc((gd.NX+2) * sizeof(float));
	msh->uf = (float *)malloc((gd.NX+2) * sizeof(float));
	msh->vh = (float *)malloc((gd.NY+2) * sizeof(float));
	msh->vf = (float *)malloc((gd.NY+2) * sizeof(float));
	msh->mh = (float *)malloc((gd.NZ+2) * sizeof(float));
	msh->mf = (float *)malloc((gd.NZ+2) * sizeof(float));

	// Carefully consider the UH,UF etc. macros and make sure they are appropriate for where you are in
	// the code (should be in pointer land)

	for (ix=gd.X0-1; ix<gd.X1+1; ix++) UH(ix-gd.X0) = msh->dx/(msh->xffull[ix+1]-msh->xffull[ix]);
	for (ix=gd.X0-1; ix<gd.X1+1; ix++) UF(ix-gd.X0) = msh->dx/(msh->xhfull[ix]-msh->xhfull[ix-1]);
	for (iy=gd.Y0-1; iy<gd.Y1+1; iy++) VH(iy-gd.Y0) = msh->dy/(msh->yffull[iy+1]-msh->yffull[iy]);
	for (iy=gd.Y0-1; iy<gd.Y1+1; iy++) VF(iy-gd.Y0) = msh->dy/(msh->yhfull[iy]-msh->yhfull[iy-1]);
	for (iz=gd.Z0-1; iz<gd.Z1+1; iz++) MH(iz-gd.Z0) = msh->dz/(msh->zf[iz+1]-msh->zf[iz]);
	for (iz=gd.Z0-1; iz<gd.Z1+1; iz++) MF(iz-gd.Z0) = msh->dz/(msh->zh[iz]-msh->zf[iz-1]);
	for (iz=gd.Z0; iz<=gd.Z1; iz++) msh->zfout[iz-gd.Z0] = 0.001*msh->zf[iz]; 
	for (iy=gd.Y0; iy<=gd.Y1; iy++) msh->yfout[iy-gd.Y0] = 0.001*msh->yffull[iy];
	for (ix=gd.X0; ix<=gd.X1; ix++) msh->xfout[ix-gd.X0] = 0.001*msh->xffull[ix];
	for (iz=gd.Z0; iz<=gd.Z1; iz++) msh->zhout[iz-gd.Z0] = 0.001*msh->zh[iz];
	for (iy=gd.Y0; iy<=gd.Y1; iy++) msh->yhout[iy-gd.Y0] = 0.001*msh->yhfull[iy];
	for (ix=gd.X0; ix<=gd.X1; ix++) msh->xhout[ix-gd.X0] = 0.001*msh->xhfull[ix];
}
/*
void set_nc_meta(ncstruct nc, int ivar, char *namestandard, char *name, char *units)
{
	int len;

	len=strlen(name);
	status = nc_put_att_text(nc->ncid, nc->varnameid[ivar],namestandard,len,name); if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	len=strlen(units); status = nc_put_att_text(nc->ncid, nc->varnameid[ivar], "units",len,units);
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
}
*/

void set_nc_meta(int ncid, int varnameid, char *namestandard, char *name, char *units)
{
	int len,status;

	len=strlen(name);
	status = nc_put_att_text(ncid, varnameid,namestandard,len,name);
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");

	len=strlen(units);
	status = nc_put_att_text(ncid,varnameid, "units",len,units);
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
}

// We need these global variables only for the H5Giter function 
int n2d;
const char **twodvarname;
int *twodvarid;
int ncid;
int d2[3];

herr_t twod_first_pass_hdf2nc(hid_t loc_id, const char *name, void *opdata)
{
    H5G_stat_t statbuf;
    H5Gget_objinfo(loc_id, name, FALSE, &statbuf);
    n2d++;
    return 0;
}

herr_t twod_second_pass_hdf2nc(hid_t loc_id, const char *name, void *opdata)
{
    H5G_stat_t statbuf;
    H5Gget_objinfo(loc_id, name, FALSE, &statbuf);
    strcpy((char *)twodvarname[n2d],name);
    nc_def_var (ncid, twodvarname[n2d], NC_FLOAT, 3, d2, &(twodvarid[n2d]));
    n2d++;
    return 0;
}

void set_netcdf_attributes(ncstruct *nc, grid gd, cmdline *cmd, buffers *b, hdf_meta *hm, hid_t *f_id)
{
	int status;
	int nv;
	int ivar;

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
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->xhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->xhid, "axis", strlen("X"), "X");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->yhid, "long_name", strlen("y-coordinate in Cartesian system"), "y-coordinate in Cartesian system");
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->yhid, "units", strlen("km"), "km");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->yhid, "axis", strlen("Y"), "Y");if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->zhid, "long_name", strlen("z-coordinate in Cartesian system"), "z-coordinate in Cartesian system");
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
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
	status = nc_put_att_text(nc->ncid, nc->x0id, "long_name", strlen("westmost grid index from LOFS data"), "westmost grid index from LOFS data");
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->x1id, "long_name", strlen("eastmost grid index from LOFS data"), "eastmost grid index from LOFS data");
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->y0id, "long_name", strlen("southmost grid index from LOFS data"), "southmost grid index from LOFS data");
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->y1id, "long_name", strlen("northmost grid index from LOFS data"), "northmost grid index from LOFS data");
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->z0id, "long_name", strlen("bottom grid index from LOFS data"), "bottom grid index from LOFS data");
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");
	status = nc_put_att_text(nc->ncid, nc->z1id, "long_name", strlen("top grid index from LOFS data"), "top grid index from LOFS data");
	if (status != NC_NOERR) ERROR_STOP("nc_put_att_text failed");

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

//		do this later
//		bufsize = (long) (gd.NX) * (long) (gd.NY) * (long) sizeof(float);
//		if ((b->twodfield = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our 2D swaths buffer array");
/* n2d, twodvarname and twodvarid are global vars declared above to enable the H5Giterate function */
		n2d = 0;
		H5Giterate(*f_id, "/00000/2D/static",NULL,twod_first_pass_hdf2nc,NULL);
		H5Giterate(*f_id, "/00000/2D/swath",NULL,twod_first_pass_hdf2nc,NULL);

		twodvarname = (const char **)malloc(n2d*sizeof(char *));
		twodvarid =   (int *)        malloc(n2d*sizeof(int));

		hm->n2dswaths = n2d;
		printf("There are %i 2D static/swath fields (the former domain of the 2D files).\n",n2d);
//		bufsize = (long) (NX) * (long) (NY) * (long) (n2d) * (long) sizeof(float);
//		if ((twodbuf0 = twodbuffer = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our 3D variable buffer array");

		for (i2d=0; i2d<n2d; i2d++)
		{
			twodvarname[i2d] = (char *)malloc(50*sizeof(char)); // 50 characters per variable
		}

		n2d = 0;
		ncid=nc->ncid; //Needed for iteration func.
		H5Giterate(*f_id, "/00000/2D/static",NULL,twod_second_pass_hdf2nc,NULL);
		H5Giterate(*f_id, "/00000/2D/swath",NULL,twod_second_pass_hdf2nc,NULL);

		for (i2d=0; i2d<n2d; i2d++) free((void *)twodvarname[i2d]); //shaddap compiler
		free(twodvarname);
//		free(twodvarid);

		/* And, like magic, we have populated our netcdf id arrays for all the swath slices */
	}


//OK we are going for simple here. We create a single varname
//array that contains all available variables, plus any
//additional requested variables. We then check for duplicates
//and remove them - this alphabetizes the order of the variables,
//however, which could be annoying and/or useful

	if (!cmd->do_allvars)
	{
		cmd->nvar = cmd->nvar_cmdline;
		nv=0; /* liar! */
	}
	else nv=hm->nvar_available;
//With faking nvar_available to zero this loop now sorts/uniqs
//all possible variables lists (just command line, just allvars,
//or a combination of the two)
	if (cmd->nvar !=0||cmd->do_allvars) 
	{
		char **varname_tmp;
		int ndupes = 0;
		int i,j;

		varname_tmp = (char **)malloc(MAXVARIABLES * sizeof(char *));
		for (i=0; i < MAXVARIABLES; i++) varname_tmp[i] = (char *)(malloc(MAXSTR * sizeof(char)));

		for (i=0; i<nv; i++)
		{
			strcpy(varname_tmp[i],hm->varname_available[i]);
		}
		for (i=0;i<cmd->nvar_cmdline; i++)
		{
			strcpy(varname_tmp[i+nv],cmd->varname_cmdline[i]);
		}

		sortchararray (varname_tmp,nv+cmd->nvar_cmdline);

		strcpy(nc->varname[0],varname_tmp[0]);
		j=1;
		/* Get rid of accidental duplicates */
		for (i=0; i<nv+cmd->nvar_cmdline-1; i++)
		{
			if(strcmp(varname_tmp[i],varname_tmp[i+1]))
			{
				strcpy(nc->varname[j],varname_tmp[i+1]); //THIS IS WHERE VARNAME IS SET!
				j++;
			}
		}
		cmd->nvar = j;
		ndupes = nv+cmd->nvar_cmdline-cmd->nvar;
		if(ndupes!=0)printf("We got rid of %i duplicate requested variables\n",ndupes);
		for (i=0; i < MAXVARIABLES; i++) free(varname_tmp[i]);
		free(varname_tmp);
	}

// Now that we have nvar, allocate the netcdf varnameid array

	nc->varnameid        = (int *) malloc(cmd->nvar*sizeof(int));

// This is our main "loop over all requested variable names" loop that
// sets all the metadata shit.

	for (ivar = 0; ivar < cmd->nvar; ivar++)
	{
		/* u v and w live on their own mesh (Arakawa C grid)*/

		/* Recommend, however, for making netcdf files, to just
		 * request uinterp vinterp winterp which NOW are cacluated
		 * HERE rather than saved in LOFS */

		/* We are going to preserve u v w with the extra point for
		 * saving u v and w easily while facilitating averaging
		 * which requires 1 fewer point */

		if(!strcmp(nc->varname[ivar],"u"))
		{
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = nc->nzh_dimid;
			nc->dims[2] = nc->nyh_dimid;
			nc->dims[3] = nc->nxf_dimid;
		}
		else if (!strcmp(nc->varname[ivar],"v"))
		{
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = nc->nzh_dimid;
			nc->dims[2] = nc->nyf_dimid;
			nc->dims[3] = nc->nxh_dimid;
		}
		else if (!strcmp(nc->varname[ivar],"w"))
		{
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = nc->nzf_dimid;
			nc->dims[2] = nc->nyh_dimid;
			nc->dims[3] = nc->nxh_dimid;
		}
		else /* scalar grid */
		{
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = nc->nzh_dimid;
			nc->dims[2] = nc->nyh_dimid;
			nc->dims[3] = nc->nxh_dimid;
		}

// I'm now going to create truly 2D files when we ask for them rather
//than putting them into a 3D container as I have done in the past

		if(gd.X0==gd.X1)
		{
			nc->twodslice = TRUE;
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = nc->nzh_dimid;
			nc->dims[2] = nc->nyh_dimid;
			nc->start[0] = 0;
			nc->start[1] = 0;
			nc->start[2] = 0;
			nc->edges[0] = 1;
			nc->edges[1] = gd.NZ;
			nc->edges[2] = gd.NY;
			status = nc_def_var (nc->ncid, nc->varname[ivar], NC_FLOAT, 3, nc->dims, &(nc->varnameid[ivar]));
		}
		else if(gd.Y0==gd.Y1)
		{
			nc->twodslice = TRUE;
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = nc->nzh_dimid;
			nc->dims[2] = nc->nxh_dimid;
			nc->start[0] = 0;
			nc->start[1] = 0;
			nc->start[2] = 0;
			nc->edges[0] = 1;
			nc->edges[1] = gd.NZ;
			nc->edges[2] = gd.NX;
			status = nc_def_var (nc->ncid, nc->varname[ivar], NC_FLOAT, 3, nc->dims, &(nc->varnameid[ivar]));
		}
		else if(gd.Z0==gd.Z1)
		{
			nc->twodslice = TRUE;
			nc->dims[0] = nc->time_dimid;
			nc->dims[1] = nc->nyh_dimid;
			nc->dims[2] = nc->nxh_dimid;
			nc->start[0] = 0;
			nc->start[1] = 0;
			nc->start[2] = 0;
			nc->edges[0] = 1;
			nc->edges[1] = gd.NY;
			nc->edges[2] = gd.NX;
			status = nc_def_var (nc->ncid, nc->varname[ivar], NC_FLOAT, 3, nc->dims, &(nc->varnameid[ivar]));
		}
		else status = nc_def_var (nc->ncid, nc->varname[ivar], NC_FLOAT, 4, nc->dims, &(nc->varnameid[ivar]));

//		printf("dims:  %5i %5i %5i %5i %s\n",dims[0],dims[1],dims[2],dims[3],varname[ivar]);
//		printf("start: %5i %5i %5i %5i %s\n",start[0],start[1],start[2],start[3],varname[ivar]);
//		printf("edges: %5i %5i %5i %5i %s\n",edges[0],edges[1],edges[2],edges[3],varname[ivar]);
//		printf("\n");

		if (status != NC_NOERR) 
		{
			printf ("Cannot nc_def_var for var #%i %s, status = %i, message = %s\n", ivar, nc->varname[ivar],status,nc_strerror(status));
			ERROR_STOP("nc_def_var failed");
		}

// We are still before the nc_enddef call, in case you are lost

		if(!strcmp(nc->varname[ivar],"uinterp"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","eastward_wind","m/s");
		else if(!strcmp(nc->varname[ivar],"vinterp"))	set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","northward_wind","m/s");
		else if(!strcmp(nc->varname[ivar],"winterp"))	set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","upward_wind","m/s");
		else if(!strcmp(nc->varname[ivar],"prespert"))	set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","pressure_perturbation","hPa");
		else if(!strcmp(nc->varname[ivar],"thpert"))	set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","potential_temperature_perturbation","K");
		else if(!strcmp(nc->varname[ivar],"rhopert"))	set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","density_perturbation","kg/m^3");
		else if(!strcmp(nc->varname[ivar],"khh"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","horizontal_subgrid_eddy_scalar_diffusivity","m^2/s");
		else if(!strcmp(nc->varname[ivar],"khv"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","vertical_subgrid_eddy_scalar_diffusivity","m^2/s");
		else if(!strcmp(nc->varname[ivar],"kmh"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","horizontal_subgrid_eddy_momentum_viscosity","m^2/s");
		else if(!strcmp(nc->varname[ivar],"kmv"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","vertical_subgrid_eddy_momentum_viscosity","m^2/s");
		else if(!strcmp(nc->varname[ivar],"xvort"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","x_vorticity","s^-1");
		else if(!strcmp(nc->varname[ivar],"yvort"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","y_vorticity","s^-1");
		else if(!strcmp(nc->varname[ivar],"zvort"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","z_vorticity","s^-1");
		else if(!strcmp(nc->varname[ivar],"vortmag"))	set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","vorticity_magnitude","s^-1");
		else if(!strcmp(nc->varname[ivar],"hvort"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","horizontal_vorticity_magnitude","s^-1");
		else if(!strcmp(nc->varname[ivar],"streamvort"))set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","streamwise_vorticity","s^-1");
		else if(!strcmp(nc->varname[ivar],"dbz"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","radar_reflectivity_simulated","dBZ");
		else if(!strcmp(nc->varname[ivar],"qvpert"))	set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","water_vapor_perturbation_mixing_ratio","cm^-3");
		else if(!strcmp(nc->varname[ivar],"qc"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","cloud_water_mixing_ratio","g/kg");
		else if(!strcmp(nc->varname[ivar],"qr"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","rain_water_mixing_ratio","g/kg");
		else if(!strcmp(nc->varname[ivar],"qi"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","cloud_ice_mixing_ratio","g/kg");
		else if(!strcmp(nc->varname[ivar],"qs"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","now_mixing_ratio","g/kg");
		else if(!strcmp(nc->varname[ivar],"qg"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","hail_mixing_ratio","g/kg");
		else if(!strcmp(nc->varname[ivar],"nci"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","ice_number_concenctration","cm^-3");
		else if(!strcmp(nc->varname[ivar],"ncr"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","rain_number_concenctration","cm^-3");
		else if(!strcmp(nc->varname[ivar],"ncs"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","snow_number_concenctration","cm^-3");
		else if(!strcmp(nc->varname[ivar],"ncg"))		set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","hail_number_concenctration","cm^-3");
		else if(!strcmp(nc->varname[ivar],"hwin_sr"))	set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","storm_relative_horizontal_wind_speed","m/s");
		else if(!strcmp(nc->varname[ivar],"windmag_sr"))set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","storm_relative_wind_speed","m/s");
		else if(!strcmp(nc->varname[ivar],"hwin_gr"))	set_nc_meta(nc->ncid,nc->varnameid[ivar],"standard_name","ground_relative_horizontal_wind_speed","m/s");
		if (cmd->gzip) status=nc_def_var_deflate(nc->ncid, nc->varnameid[ivar], 1, 1, 1);
	}

/* Write sounding data */

	status = nc_def_var (nc->ncid, "u0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->u0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta(nc->ncid,nc->u0id,"standard_name","base_state_u","m/s");
	status = nc_def_var (nc->ncid, "v0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->v0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta(nc->ncid,nc->v0id,"standard_name","base_state_v","m/s");
	status = nc_def_var (nc->ncid, "th0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->th0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta(nc->ncid,nc->th0id,"standard_name","base_state_theta","m/s");
	status = nc_def_var (nc->ncid, "pres0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->pres0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta(nc->ncid,nc->pres0id,"standard_name","base_state_pressure","m/s");
	status = nc_def_var (nc->ncid, "pi0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->pi0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta(nc->ncid,nc->pi0id,"standard_name","base_state_pi","m/s");
	status = nc_def_var (nc->ncid, "qv0", NC_FLOAT, 1, &nc->nzh_dimid, &nc->qv0id); if (status != NC_NOERR) ERROR_STOP("nc_def_var failed");
	set_nc_meta(nc->ncid,nc->pi0id,"standard_name","base_state_qv","m/s");

//	expose this in main for clarity since it denotes when we can start actually doing omething
//	status = nc_enddef (nc->ncid); if (status != NC_NOERR) ERROR_STOP("nc_enddef failed");
}

void nc_write_1d_data (ncstruct nc, grid gd, mesh msh, sounding snd, cmdline cmd)
{
	double timearray[1];
	//ORF for writing single time in unlimited time dimension/variable
	const size_t timestart = 0;
	const size_t timecount = 1;
	int status;

	status = nc_put_var_float (nc.ncid,nc.xhid,msh.xhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.yhid,msh.yhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.zhid,msh.zhout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.xfid,msh.xfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.yfid,msh.yfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.zfid,msh.zfout); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.u0id,snd.u0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.v0id,snd.v0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.th0id,snd.th0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.pres0id,snd.pres0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.pi0id,snd.pi0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_float (nc.ncid,nc.qv0id,snd.qv0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_float failed");
	status = nc_put_var_int (nc.ncid,nc.x0id,&gd.X0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	status = nc_put_var_int (nc.ncid,nc.y0id,&gd.Y0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	status = nc_put_var_int (nc.ncid,nc.x1id,&gd.X1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	status = nc_put_var_int (nc.ncid,nc.y1id,&gd.Y1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	status = nc_put_var_int (nc.ncid,nc.z0id,&gd.Z0); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	status = nc_put_var_int (nc.ncid,nc.z1id,&gd.Z1); if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
	timearray[0] = cmd.time;
	status = nc_put_vara_double (nc.ncid,nc.timeid,&timestart,&timecount,timearray);
	if (status != NC_NOERR) ERROR_STOP ("nc_put_var_int failed");
}

void set_readahead(readahead *rh,ncstruct nc,cmdline cmd)
{
	int ivar;

	for (ivar = 0; ivar < cmd.nvar; ivar++)
	{
		if(!strcmp(nc.varname[ivar],"uinterp")&&!cmd.use_interp) {rh->u=1;}
		if(!strcmp(nc.varname[ivar],"vinterp")&&!cmd.use_interp) {rh->v=1;}
		if(!strcmp(nc.varname[ivar],"winterp")&&!cmd.use_interp) {rh->w=1;}
		if(!strcmp(nc.varname[ivar],"hwin_sr")) {rh->u=rh->v=1;}
		if(!strcmp(nc.varname[ivar],"hdiv")) {rh->u=rh->v=1;}
		if(!strcmp(nc.varname[ivar],"windmag_sr")) {rh->u=rh->v=rh->w=1;}
		if(!strcmp(nc.varname[ivar],"zvort")) {rh->u=rh->v=1;}
		if(!strcmp(nc.varname[ivar],"xvort")) {rh->v=rh->w=1;}
		if(!strcmp(nc.varname[ivar],"yvort")) {rh->u=rh->w=1;}
		if(!strcmp(nc.varname[ivar],"hvort")) {rh->u=rh->v=rh->w=1;}
		if(!strcmp(nc.varname[ivar],"vortmag")) {rh->u=rh->v=rh->w=1;}
		if(!strcmp(nc.varname[ivar],"streamvort")) {rh->u=rh->v=rh->w=1;}
		if(!strcmp(nc.varname[ivar],"hwin_gr")) {rh->u=rh->v=rh->w=1;}
//		if(!strcmp(varname[ivar],"streamfrac")) {u_rh=v_rh=w_rh=1;}
//		if(!strcmp(varname[ivar],"xvort_baro")) {thrhopert_rh=1;}
//		if(!strcmp(varname[ivar],"yvort_baro")) {thrhopert_rh=1;}
//		if(!strcmp(varname[ivar],"yvort_stretch")) {u_rh=w_rh=yvort_rh=1;}
//		if(!strcmp(varname[ivar],"yvort_tilt")) {v_rh=xvort_rh=yvort_rh=1;}
//		if(!strcmp(varname[ivar],"xvort_stretch")) {v_rh=w_rh=xvort_rh=1;}
//		if(!strcmp(varname[ivar],"xvort_tilt")) {u_rh=yvort_rh=zvort_rh=1;}
//		if(!strcmp(varname[ivar],"zvort_stretch")) {u_rh=v_rh=zvort_rh=1;}
//		if(!strcmp(varname[ivar],"zvort_tilt")) {w_rh=xvort_rh=zvort_rh=1;}
//Might use these, but I usually only read cloud stuff once, so this
//isn't used... yet
//		if(!strcmp(varname[ivar],"qcloud")) {qc_rh=qi_rh=1;}
//		if(!strcmp(varname[ivar],"qprecip")) {qr_rh=qs_rh=qg_rh=1;}
	}

}

void malloc_3D_arrays (buffers *b, grid gd, readahead rh,cmdline cmd)
{
	if (cmd.nvar>0)
	{
		long bufsize,totbufsize;

		bufsize = (long) (gd.NX+2) * (long) (gd.NY+2) * (long) (gd.NZ+1) * (long) sizeof(float);
		totbufsize = bufsize;

		if ((b->buf0 = b->buf = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our 3D variable buffer array");
		if (rh.u)
		{
			if ((b->ustag = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our ustag buffer array");
			totbufsize+=bufsize;
		}
		if (rh.v)
		{
			if ((b->vstag = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our vstag buffer array");
			totbufsize+=bufsize;
		}
		if (rh.w)
		{
			if ((b->wstag = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our wstag buffer array");
			totbufsize+=bufsize;
		}
		if (rh.u||rh.v||rh.w)
		{
			if ((b->dum00 = b->dum0 = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our first 3D temp calculation array");
			totbufsize+=bufsize;
		}
		if (rh.vortmag||rh.hvort||rh.streamvort)//Not really readahead, but if we calculated these we need another array
		{
			if ((b->dum10 = b->dum1 = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our second 3D temp calculation array");
			totbufsize+=bufsize;
		}
		printf("We allocated a total of %8.5f GB of memory for floating point buffers\n",1.0e-9*totbufsize);
	}
}

void do_the_swaths(hdf_meta hm, ncstruct nc, dir_meta dm, grid gd, cmdline cmd)
{
	int i2d,ix,iy,status;
	float *twodfield,*swathbuf,*writeptr;
	float bufsize;

	printf("Working on 2D static fields and swaths ("); 

	bufsize = (long) (gd.NX) * (long) (gd.NY) * (long) sizeof(float);
	if ((twodfield = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our 2D swath buffer array");
	bufsize = (long) (gd.NX) * (long) (gd.NY) * (long) (hm.n2dswaths) * (long) sizeof(float);
	if ((swathbuf = (float *) malloc ((size_t)bufsize)) == NULL) ERROR_STOP("Cannot allocate our 3D swaths buffer array");

//	read_hdf_mult_md(swathbuf0,topdir,timedir,nodedir,ntimedirs,dn,dirtimes,alltimes,ntottimes,t0,"swaths",X0,Y0,X1,Y1,0,n2d,nx,ny,nz,nodex,nodey);

	read_lofs_buffer(swathbuf,"swaths",dm,hm,gd,cmd);

	for (i2d=0;i2d<hm.n2dswaths;i2d++)
	{
		for (iy=0; iy<gd.NY; iy++)
			for (ix=0; ix<gd.NX; ix++)
				twodfield[P2(ix,iy,gd.NX)] = swathbuf[P3(ix,iy,i2d,gd.NX,gd.NY)];
		writeptr = twodfield;
		status = nc_put_vara_float (nc.ncid, twodvarid[i2d], nc.s2, nc.e2, writeptr);
	}
	free(swathbuf);
	free(twodfield);
	printf(")\n");
}
