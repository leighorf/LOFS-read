#include "../include/lofs-read.h"
#include "../include/lofs-dirstruct.h"
#include "../include/lofs-hdf2nc.h"

void get0dint (hid_t f_id, char *varname, int *var)
{
	hid_t dataset_id;
	int status;

	/* shouldn't need these any more */
	if ((dataset_id = H5Dopen (f_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
//	printf("varname = %s, dataset_id = %i\n",varname,dataset_id);
	if ((status = H5Dread (dataset_id, H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,var)) < 0) ERROR_STOP("Could not H5Dread");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP("Could not H5Dclose");
}

void get0dfloat (hid_t f_id, char *varname, float *var)
{
	hid_t dataset_id;
	int status;

//	printf("Getting %s\n",varname);
	if ((dataset_id = H5Dopen (f_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
	if ((status = H5Dread (dataset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,var)) < 0) ERROR_STOP("Could not H5Dread");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP("Could not H5Dclose");
}

void get1ddouble (hid_t f_id, char *varname, double *var, int p0, int np)
{
	int rank;
	hsize_t count[1], dims[1];
	hsize_t offset_in[1],offset_out[1];
	hid_t dataset_id,dataspace_id,memoryspace_id;
	herr_t status;

	rank = 1;
	offset_in[0] = p0;
	offset_out[0] = 0; //Always
	count[0] = np;
	dims[0] = np;

	if ((dataset_id = H5Dopen (f_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) ERROR_STOP("Could not H5Dget_space");
	if ((memoryspace_id = H5Screate_simple(rank,dims,NULL)) < 0) ERROR_STOP("Could not H5Screate_simple");
	if ((status = H5Sselect_hyperslab (dataspace_id,H5S_SELECT_SET,offset_in,NULL,count,NULL)) < 0) ERROR_STOP("Could not H5Sselect_hyperslab");
	if ((status = H5Sselect_hyperslab (memoryspace_id,H5S_SELECT_SET,offset_out,NULL,count,NULL)) < 0) ERROR_STOP ("Could not H5Sselect_hyperslab");
	if ((status = H5Dread (dataset_id, H5T_NATIVE_DOUBLE,memoryspace_id,dataspace_id,H5P_DEFAULT,var)) < 0) ERROR_STOP ("Could not H5Dread");
	if ((status = H5Sclose(memoryspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Sclose(dataspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP("Could not H5Dclose");

}

void get1dfloat (hid_t f_id, char *varname, float *var, int p0, int np)
{
	int rank;
	hsize_t count[1], dims[1];
	hsize_t offset_in[1],offset_out[1];
	hid_t dataset_id,dataspace_id,memoryspace_id;
	herr_t status;

	rank = 1;
	offset_in[0] = p0;
	offset_out[0] = 0; //Always
	count[0] = np;
	dims[0] = np;

	if ((dataset_id = H5Dopen (f_id, varname,H5P_DEFAULT)) < 0)
	{
		printf("get1dfloat: %s failed\n",varname);
		ERROR_STOP("Could not H5Dopen");
	}
	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) ERROR_STOP("Could not H5Dget_space");
	if ((memoryspace_id = H5Screate_simple(rank,dims,NULL)) < 0) ERROR_STOP("Could not H5Screate_simple");
	if ((status = H5Sselect_hyperslab (dataspace_id,H5S_SELECT_SET,offset_in,NULL,count,NULL)) < 0) ERROR_STOP("Could not H5Sselect_hyperslab");
	if ((status = H5Sselect_hyperslab (memoryspace_id,H5S_SELECT_SET,offset_out,NULL,count,NULL)) < 0) ERROR_STOP ("Could not H5Sselect_hyperslab");
	if ((status = H5Dread (dataset_id, H5T_NATIVE_FLOAT,memoryspace_id,dataspace_id,H5P_DEFAULT,var)) < 0) ERROR_STOP ("Could not H5Dread");
	if ((status = H5Sclose(memoryspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Sclose(dataspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP("Could not H5Dclose");

}

//ORF this will fill our hdf_meta struct and some other things
void get_hdf_metadata(dir_meta dm, hdf_meta *hm, cmdline *cmd, ncstruct *nc, char *argv[], hid_t *f_id, zfpacc *zfpacc)
{
	hid_t g_id,d_id,attr_id,attr_memtype;
	herr_t status;
	H5O_info_t dset_info;
	H5G_info_t group_info;
	int i,j,nattr;
	double *zptr;
	char groupname[MAXSTR];
	char attrname[MAXATTR][MAXSTR]; //ORF FIX TODO
	htri_t existence;
	hid_t lapl_id;

//ORF I (finally) changed CM1's 'nodex' and 'nodey' to 'rankx' and 'ranky'
//because they refer to MPI ranks not 'nodes'
//But make sure we are backward compatible!!
//Note: Watch for leading slashes, H5Lexists is very picky, see docs.

	existence = H5Lexists(*f_id, "grid", H5P_DEFAULT);
//	printf("existence of /grid = %i\n",existence);
	existence = H5Lexists(*f_id, "grid/rankx", H5P_DEFAULT);
//	printf("existence of /grid/rankx = %i\n",existence);
//  If rankx doesn't exist, neither does ranky so that's enough checking
	if (existence > 0)
	{
		get0dint (*f_id, "grid/rankx", &hm->rankx);
		get0dint (*f_id, "grid/ranky", &hm->ranky);
	}
	else //HDF5 will spew out error crap but won't exit
	{
		get0dint (*f_id, "grid/nodex", &hm->rankx);
		get0dint (*f_id, "grid/nodey", &hm->ranky);
	}

	get0dint (*f_id, "grid/nx", &hm->nx);
	get0dint (*f_id, "grid/ny", &hm->ny);
	get0dint (*f_id, "grid/nz", &hm->nz);

	sprintf(groupname,"%05i/3D",0);//All vars in 3D group 00000 are always available
	g_id = H5Gopen(*f_id,groupname,H5P_DEFAULT);
	H5Gget_info(g_id,&group_info);
	hm->nvar_available = group_info.nlinks;

	/* Allocate our varaible name arrays, available in the lofs files, and requested variables */
	hm->varname_available = (char **)malloc(hm->nvar_available * sizeof(char *));
	for (i=0; i < hm->nvar_available; i++) hm->varname_available[i] = (char *)(malloc(MAXSTR * sizeof(char)));
	cmd->varname_cmdline = (char **)malloc(cmd->nvar_cmdline * sizeof(char *));
	for (i=0; i < cmd->nvar_cmdline; i++) cmd->varname_cmdline[i] = (char *)(malloc(MAXSTR * sizeof(char)));

	for (i = 0; i < hm->nvar_available; i++)
	{
	    H5Lget_name_by_idx(g_id,".",H5_INDEX_NAME,H5_ITER_INC,i,hm->varname_available[i],40,H5P_DEFAULT);  //ORF TODO make 40 a constant somewhere
	}
// XXXXXX
	for (i=0; i<cmd->nvar_cmdline; i++)
	{
		strcpy(cmd->varname_cmdline[i],argv[i+cmd->argc_hdf2nc_min+cmd->optcount]);//HERE IS WHERE WE POPULATE VARNAME_CMDLINE
		strcpy(nc->var3d[i].varname,cmd->varname_cmdline[i]);
		if(cmd->verbose)printf("var3d[%i].varname = %s\n",i,nc->var3d[i].varname);
	}
/*
 *
 * ORF 2021-03-26
 * Now we get 3D zfp_accuracy_LOFS data that we need
 * We loop over our commandline requested 3D variables
 * We check if it's native LOFS
 * If so we sweep through each 3D variable's attributes
 * We look for the "zfp_accuracy_LOFS" attribute
 * If it exists, we grab the value and stash it in the appropriate place
 * It will be written to the netCDF file as well as the ZFP value for the netCDF data,
 * which is allowed to be whatever you want... it's up to the user to make sure that they
   aren't choosing stupid accuracy parameters
 *
 */

//		if ((d_id = H5Dopen (g_id,hm->varname_available[i],H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
//
//	for (i = 0; i < hm->nvar_available; i++)
	for (i = 0; i < cmd->nvar_cmdline; i++)
	{
		nc->var3d[i].is_LOFS_var=0;
		for (j=0; j<hm->nvar_available; j++) //check against all available vars to see if we have ZFP data to retrieve
		{
//			printf("%s %s\n",hm->varname_available[j],nc->var3d[i].varname);
			if (same(hm->varname_available[j],nc->var3d[i].varname)) // Is this a native LOFS variable?
			{
				nc->var3d[i].is_LOFS_var=1;
				break;
			}
		}
		if (nc->var3d[i].is_LOFS_var)
		{
			if ((d_id = H5Dopen (g_id,nc->var3d[i].varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
			if ((status = H5Oget_info(d_id,&dset_info,H5O_INFO_NUM_ATTRS)) < 0) ERROR_STOP("Could not H5Oget_info");
			nattr=dset_info.num_attrs;
			for (j=0; j<nattr; j++)
			{
				if ((status=H5Aget_name_by_idx(d_id,".",H5_INDEX_NAME,H5_ITER_INC,j,attrname[j],40,H5P_DEFAULT))<0) ERROR_STOP("Could not H5Aget_name_by_idex");
				if ((attr_id=H5Aopen(d_id,attrname[j],H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Aopen");
//				printf ("%s %s\n",attrname[j],"zfp_accuracy_LOFS");
				if (same(attrname[j],"zfp_accuracy"))//Just called zfp_accuracy in the LOFS HDF5 files
				{
					zptr=&(nc->var3d[i].zfpacc_LOFS);
					if ((attr_memtype=H5Aget_type(attr_id)) < 0) ERROR_STOP("Could not H5Aget_type");
					if ((status=H5Aread(attr_id,attr_memtype,zptr)) < 0) ERROR_STOP("Could not H5Aread"); //grab ZFP accuracy from LOFS 3dvar
					if(cmd->verbose)printf("%s is a native LOFS variable with zfpacc_LOFS = %f\n",nc->var3d[i].varname,nc->var3d[i].zfpacc_LOFS);
					break;
				}
			}
			H5Dclose(d_id);

/* ORF 2022-08-23*/
/* I struggled with this crap. I wrote the above code thinking I could be clever but in
 * reality I have to be dumb. And yes, I do this if else ladder shit too often but... C */

			/* Here is where we populate the zfpacc->lofs structure with the values that
			 * LOFS contains. We can then compare these against the netcdf chosen values
			 * to make sure we are not compressing with more accuracy than was originally
			 * saved, for instance! */

			if      (same(nc->var3d[i].varname,"u"))          zfpacc->lofs->u =           nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"v"))          zfpacc->lofs->v =           nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"w"))          zfpacc->lofs->w =           nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"uinterp"))    zfpacc->lofs->uinterp =     nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vinterp"))    zfpacc->lofs->vinterp =     nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"winterp"))    zfpacc->lofs->winterp =     nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"prespert"))   zfpacc->lofs->prespert =    nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"thrhopert"))  zfpacc->lofs->thrhopert =   nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"dbz"))        zfpacc->lofs->dbz =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"qc"))         zfpacc->lofs->qc =          nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"qi"))         zfpacc->lofs->qi =          nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"qr"))         zfpacc->lofs->qr =          nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"qg"))         zfpacc->lofs->qg =          nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"qs"))         zfpacc->lofs->qs =          nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"nci"))        zfpacc->lofs->nci =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ncg"))        zfpacc->lofs->ncg =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ncr"))        zfpacc->lofs->ncr =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ncs"))        zfpacc->lofs->ncs =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"qhl"))        zfpacc->lofs->qhl =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ccn"))        zfpacc->lofs->ccn =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ccw"))        zfpacc->lofs->ccw =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"crw"))        zfpacc->lofs->ccw =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"cci"))        zfpacc->lofs->cci =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"csw"))        zfpacc->lofs->csw =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"chw"))        zfpacc->lofs->chw =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"chl"))        zfpacc->lofs->chl =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vhw"))        zfpacc->lofs->vhw =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"zhl"))        zfpacc->lofs->zhl =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"zhw"))        zfpacc->lofs->zhw =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"zrw"))        zfpacc->lofs->zrw =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vhl"))        zfpacc->lofs->vhl =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"qv"))         zfpacc->lofs->qv =          nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"qvpert"))     zfpacc->lofs->qvpert =      nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"thpert"))     zfpacc->lofs->thpert =      nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"th"))         zfpacc->lofs->th =          nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"prs"))        zfpacc->lofs->prs =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"pi"))         zfpacc->lofs->pi =          nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"pipert"))     zfpacc->lofs->pipert =      nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"rho"))        zfpacc->lofs->rho =         nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"rhopert"))    zfpacc->lofs->rhopert =     nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"tke_sg"))     zfpacc->lofs->tke_sg =      nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"kmh"))        zfpacc->lofs->km =          nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"kmv"))        zfpacc->lofs->km =          nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"khh"))        zfpacc->lofs->kh =          nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"khv"))        zfpacc->lofs->kh =          nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"xvort"))      zfpacc->lofs->xvort =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"yvort"))      zfpacc->lofs->yvort =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"zvort"))      zfpacc->lofs->zvort =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vortmag"))    zfpacc->lofs->vortmag =       nc->var3d[i].zfpacc_LOFS;
/*
 * But, these should be calculated with LOFT (Kelton's calc code)
 * outside of CM1 like the other stuff... some of it already is...
			else if (same(nc->var3d[i].varname,"ub_cor"))    zfpacc->lofs->ub_cor =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ub_hadv"))    zfpacc->lofs->ub_hadv =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ub_hediff"))    zfpacc->lofs->ub_hediff =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ub_hidiff"))    zfpacc->lofs->ub_hidiff =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ub_hturb"))    zfpacc->lofs->ub_hturb =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ub_pbl"))    zfpacc->lofs->ub_pbl =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ub_pgrad"))    zfpacc->lofs->ub_pgrad =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ub_rdamp"))    zfpacc->lofs->ub_rdamp =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ub_subs"))    zfpacc->lofs->ub_subs =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ub_vadv"))    zfpacc->lofs->ub_vadv =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ub_vediff"))    zfpacc->lofs->ub_vediff =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ub_vidiff"))    zfpacc->lofs->ub_vidiff =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"ub_vturb"))    zfpacc->lofs->ub_vturb =       nc->var3d[i].zfpacc_LOFS;

			else if (same(nc->var3d[i].varname,"vb_cor"))    zfpacc->lofs->vb_cor =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vb_hadv"))    zfpacc->lofs->vb_hadv =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vb_hediff"))    zfpacc->lofs->vb_hediff =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vb_hidiff"))    zfpacc->lofs->vb_hidiff =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vb_hturb"))    zfpacc->lofs->vb_hturb =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vb_pbl"))    zfpacc->lofs->vb_pbl =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vb_pgrad"))    zfpacc->lofs->vb_pgrad =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vb_rdamp"))    zfpacc->lofs->vb_rdamp =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vb_subs"))    zfpacc->lofs->vb_subs =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vb_vadv"))    zfpacc->lofs->vb_vadv =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vb_vediff"))    zfpacc->lofs->vb_vediff =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vb_vidiff"))    zfpacc->lofs->vb_vidiff =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"vb_vturb"))    zfpacc->lofs->vb_vturb =       nc->var3d[i].zfpacc_LOFS;

			else if (same(nc->var3d[i].varname,"wb_buoy"))    zfpacc->lofs->wb_buoy =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"wb_hadv"))    zfpacc->lofs->wb_hadv =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"wb_hediff"))    zfpacc->lofs->wb_hediff =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"wb_hidiff"))    zfpacc->lofs->wb_hidiff =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"wb_hturb"))    zfpacc->lofs->wb_hturb =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"wb_pgrad"))    zfpacc->lofs->wb_pgrad =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"wb_rdamp"))    zfpacc->lofs->wb_rdamp =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"wb_vadv"))    zfpacc->lofs->wb_vadv =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"wb_vediff"))    zfpacc->lofs->wb_vediff =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"wb_vidiff"))    zfpacc->lofs->wb_vidiff =       nc->var3d[i].zfpacc_LOFS;
			else if (same(nc->var3d[i].varname,"wb_vturb"))    zfpacc->lofs->wb_vturb =       nc->var3d[i].zfpacc_LOFS;
*/
			else
			{
				fprintf(stderr,"This cannot happen. Goodbye.\n");
				ERROR_STOP("LOFS variable is not an LOFS variable LOL");
			}
		}
	}
	H5Gclose(g_id);
}

