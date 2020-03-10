#include "include/lofs-read.h"

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
//
//FORTRAN wrapper
void get0dint_(hid_t *file_id,char *varname, int *var)
{
	get0dint(*file_id,varname,var);
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

//FORTRAN wrapper
void get0dfloat_(hid_t *file_id,char *varname, float *var)
{
	get0dfloat(*file_id,varname,var);
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

void
get1dint (hid_t file_id, char *varname, int *var, int p0, int np)
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
	if ((status = H5Dread (dataset_id, H5T_NATIVE_INT,memoryspace_id,dataspace_id,H5P_DEFAULT,var)) < 0) ERROR_STOP ("Could not H5Dread");
	if ((status = H5Sclose(memoryspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Sclose(dataspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP("Could not H5Dclose");
}

void
get2dfloat (hid_t file_id, char *varname, float *var, int y0, int ny, int x0, int nx)
{
	int rank;
	hsize_t count[2], dims[2];
	hsize_t offset_in[2],offset_out[2];
	hid_t dataset_id,dataspace_id,memoryspace_id;
	int status;

	rank = 2;
	offset_in[0] = y0;
	offset_in[1] = x0;
	offset_out[0] = 0;
	offset_out[1] = 0;
	count[0] = ny;
	count[1] = nx;
	dims[0] = ny-y0;
	dims[1] = nx-x0;

	if ((dataset_id = H5Dopen (file_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) ERROR_STOP ("Could not H5Dget_space");
	if ((memoryspace_id = H5Screate_simple(rank,dims,NULL)) < 0) ERROR_STOP("Could not H5Screate_simple");
	if ((status = H5Sselect_hyperslab (dataspace_id,H5S_SELECT_SET,offset_in,NULL,count,NULL)) < 0) ERROR_STOP ("Could not H5Sselect_hyperslab");
	if ((status = H5Sselect_hyperslab (memoryspace_id,H5S_SELECT_SET,offset_out,NULL,count,NULL)) < 0) ERROR_STOP ("Could not H5Sselect_hyperslab");
	if ((status = H5Dread (dataset_id, H5T_NATIVE_FLOAT,memoryspace_id,dataspace_id,H5P_DEFAULT,var)) < 0) ERROR_STOP ("Could not H5Dread");
	if ((status = H5Sclose(memoryspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Sclose(dataspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP ("Could not H5Dclose");
}

void
get3dfloat (hid_t file_id, char *varname, float *var, int z0, int nz, int y0, int ny, int x0, int nx)
{
	int rank;
	hsize_t count[3], dims[3];
	hsize_t offset_in[3],offset_out[3];
	hid_t dataset_id,dataspace_id,memoryspace_id;
	int status;

	rank = 3;
	offset_in[0] = z0;
	offset_in[1] = y0;
	offset_in[2] = x0;
	offset_out[0] = 0;
	offset_out[1] = 0;
	offset_out[2] = 0;
	count[0] = nz;
	count[1] = ny;
	count[2] = nx;
	dims[0] = nz-z0;
	dims[1] = ny-y0;
	dims[2] = nx-x0;

	if ((dataset_id = H5Dopen (file_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP ("Could not H5Dopen");
	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) ERROR_STOP ("Could not H5Dget_space");
	if ((memoryspace_id = H5Screate_simple(rank,dims,NULL)) < 0) ERROR_STOP("Could not H5Screate_simple");
	if ((status = H5Sselect_hyperslab (dataspace_id,H5S_SELECT_SET,offset_in,NULL,count,NULL)) < 0) ERROR_STOP ("Could not H5Sselect_hyperslab");
	if ((status = H5Sselect_hyperslab (memoryspace_id,H5S_SELECT_SET,offset_out,NULL,count,NULL)) < 0) ERROR_STOP ("Could not H5Sselect_hyperslab");
	if ((status = H5Dread (dataset_id, H5T_NATIVE_FLOAT,memoryspace_id,dataspace_id,H5P_DEFAULT,var)) < 0) ERROR_STOP ("Could not H5Dread");
	if ((status = H5Sclose(memoryspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Sclose(dataspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP ("Could not H5Dclose");
}

/* I have decided to create a new routine which must be called before
 * any 3D data is read with read_hdf_mult. This will reduce the number
 * of redundant reads to the HDF files for each 3D read. While it
 * requires the user to call an extra function (only once mind you) it's just the right way
 * to do things. The pointers nx,ny,nz,nodex,nodey are updated and are
 * expected to be passed as arguments to read_hdf_mult so read_hdf_mult
 * doesn't have to get them every time. nx,ny,nodex,nodey are all that
 * are required to recreate the entire domain decomposition. */

void get_hdf_metadata(char *hdffilename, int *nx, int *ny, int *nz, int *nodex, int *nodey)
{
    hid_t file_id;
    int status;

	if ((file_id = H5Fopen (hdffilename, H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
    {
        fprintf(stderr,"\n\nget_hdf_metadata: Unable to read metadata from %s, bailing!\n", hdffilename);
        exit(0);
    }
    get0dint (file_id, "grid/nodex", nodex);
    get0dint (file_id, "grid/nodey", nodey);
    get0dint (file_id, "grid/nx", nx);
    get0dint (file_id, "grid/ny", ny);
    get0dint (file_id, "grid/nz", nz);
    if ((status = H5Fclose (file_id)) < 0)
    {
        fprintf(stderr,"\n\n10900: get_hdf_metadata: OH NO! Can't close hdf file with sd_id = %i\n",file_id);
        fprintf(stderr, "This simply should not ever happen.  Exiting out of fear.\n");
        exit(-1);
    }

}
