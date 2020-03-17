#include "../include/lofs-read.h"
#include "../include/dirstruct.h"
#include "../include/hdf2nc.h"

void
get0dint (hid_t f_id, char *varname, int *var)
{
	hid_t dataset_id;
	int status;

	/* shouldn't need these any more */
	if ((dataset_id = H5Dopen (f_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
//	printf("varname = %s, dataset_id = %i\n",varname,dataset_id);
	if ((status = H5Dread (dataset_id, H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,var)) < 0) ERROR_STOP("Could not H5Dread");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP("Could not H5Dclose");
}

void
get0dfloat (hid_t f_id, char *varname, float *var)
{
	hid_t dataset_id;
	int status;

	printf("Getting %s\n",varname);
	if ((dataset_id = H5Dopen (f_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
	if ((status = H5Dread (dataset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,var)) < 0) ERROR_STOP("Could not H5Dread");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP("Could not H5Dclose");
}

void
get1ddouble (hid_t f_id, char *varname, double *var, int p0, int np)
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

void
get1dfloat (hid_t f_id, char *varname, float *var, int p0, int np)
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

//	printf("get1dfloat: %s\n",varname);
	if ((dataset_id = H5Dopen (f_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) ERROR_STOP("Could not H5Dget_space");
	if ((memoryspace_id = H5Screate_simple(rank,dims,NULL)) < 0) ERROR_STOP("Could not H5Screate_simple");
	if ((status = H5Sselect_hyperslab (dataspace_id,H5S_SELECT_SET,offset_in,NULL,count,NULL)) < 0) ERROR_STOP("Could not H5Sselect_hyperslab");
	if ((status = H5Sselect_hyperslab (memoryspace_id,H5S_SELECT_SET,offset_out,NULL,count,NULL)) < 0) ERROR_STOP ("Could not H5Sselect_hyperslab");
	if ((status = H5Dread (dataset_id, H5T_NATIVE_FLOAT,memoryspace_id,dataspace_id,H5P_DEFAULT,var)) < 0) ERROR_STOP ("Could not H5Dread");
	if ((status = H5Sclose(memoryspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Sclose(dataspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP("Could not H5Dclose");

}

//ORF this will fill our hdf_meta struct
void get_hdf_metadata(dir_meta dm, hdf_meta *hm, cmdline *cmd, char *argv[], hid_t *f_id)
{
    hid_t g_id;
    H5G_info_t group_info;
    int i,status;
	char groupname[MAXSTR];

	/* Already open 
	if ((*f_id = H5Fopen (dm.firstfilename, H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
    {
		fprintf(stderr,"\n\nget_hdf_metadata: Unable to read metadata from %s, bailing!\n", dm.firstfilename);
   		exit(0);
    }
	*/
    get0dint (*f_id, "grid/nodex", &hm->nodex);
    get0dint (*f_id, "grid/nodey", &hm->nodey);
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
	    H5Lget_name_by_idx(g_id,".",H5_INDEX_NAME,H5_ITER_INC,i,hm->varname_available[i],40,H5P_DEFAULT); //40 characters per varname
		//ORF TODO make 40 a constant somewhere
	}
	H5Gclose(g_id);

	for (i=0; i<cmd->nvar_cmdline; i++)
	{
		strcpy(cmd->varname_cmdline[i],argv[i+cmd->argc_hdf2nc_min+cmd->optcount]);//HERE IS WHERE WE POPULATE VARNAME_CMDLINE
	} printf("\n");
}

