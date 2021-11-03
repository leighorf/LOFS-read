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
void get_hdf_metadata(dir_meta dm, hdf_meta *hm, cmdline *cmd, ncstruct *nc, char *argv[], hid_t *f_id)
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

	// ORF FIX the order here is not the same order when we do our variables
	// So the indexing wil be wrong
	// We sort everything, this determines the order we do shit in
	// FUCKBAGS1 WHY NOT SORT HERE NOW THEN INDEXING WILL REMAIN GOOD
	for (i=0; i<cmd->nvar_cmdline; i++)
	{
		strcpy(cmd->varname_cmdline[i],argv[i+cmd->argc_hdf2nc_min+cmd->optcount]);//HERE IS WHERE WE POPULATE VARNAME_CMDLINE
		strcpy(nc->var3d[i].varname,cmd->varname_cmdline[i]);
		printf("var3d[%i].varname = %s\n",i,nc->var3d[i].varname);
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
				printf("%s is a native LOFS variable\n",nc->var3d[i].varname);
				break;
			}
		}
		if (nc->var3d[i].is_LOFS_var)
		{
			if ((d_id = H5Dopen (g_id,nc->var3d[i].varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
			if ((status = H5Oget_info(d_id,&dset_info)) < 0) ERROR_STOP("Could not H5Oget_info");
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
					if (1) printf("%10s: %s = %12.7f\n",nc->var3d[i].varname,attrname[j],*zptr);
					break;
				}
			}
			H5Dclose(d_id);
		}
	}
	H5Gclose(g_id);
}

