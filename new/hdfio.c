#include "include/lofs-read.h"
#include "include/dirstruct.h"

void
get0dint (hid_t file_id, char *varname, int *var)
{
	hid_t dataset_id;
	int status;

	/* shouldn't need these any more */
	if ((dataset_id = H5Dopen (file_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
//	printf("varname = %s, dataset_id = %i\n",varname,dataset_id);
	if ((status = H5Dread (dataset_id, H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,var)) < 0) ERROR_STOP("Could not H5Dread");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP("Could not H5Dclose");
}

void
get0dfloat (hid_t file_id, char *varname, float *var)
{
	hid_t dataset_id;
	int status;

	if ((dataset_id = H5Dopen (file_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
	if ((status = H5Dread (dataset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,var)) < 0) ERROR_STOP("Could not H5Dread");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP("Could not H5Dclose");
}

void
get1ddouble (hid_t file_id, char *varname, double *var, int p0, int np)
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

	if ((dataset_id = H5Dopen (file_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
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
get1dfloat (hid_t file_id, char *varname, float *var, int p0, int np)
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
	if ((dataset_id = H5Dopen (file_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
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
void get_hdf_metadata(dir_meta dm, hdf_meta *hm)
{
    hid_t file_id;
    int status;

	if ((file_id = H5Fopen (dm.firstfilename, H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
    {
        fprintf(stderr,"\n\nget_hdf_metadata: Unable to read metadata from %s, bailing!\n", dm.firstfilename);
        exit(0);
    }
    get0dint (file_id, "grid/nodex", &hm->nodex);
    get0dint (file_id, "grid/nodey", &hm->nodey);
    get0dint (file_id, "grid/nx", &hm->nx);
    get0dint (file_id, "grid/ny", &hm->ny);
    get0dint (file_id, "grid/nz", &hm->nz);
    if ((status = H5Fclose (file_id)) < 0)
    {
        fprintf(stderr,"\n\n10900: get_hdf_metadata: OH NO! Can't close hdf file with sd_id = %i\n",file_id);
        fprintf(stderr, "This simply should not ever happen.  Exiting out of fear.\n");
        exit(-1);
    }

}
