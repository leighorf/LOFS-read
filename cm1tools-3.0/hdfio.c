#include "hdforf.h"
#include "errorf.h"

void
open_cm1_hdf_file (hid_t *file_id, char *base, int itime, int inode)
{
	char fname[200];
	itime = 30; //ORF test 5/17/11
	sprintf (fname, "%s%05i_%06i.cm1hdf5", base, itime, inode);

	if ((*file_id = H5Fopen (fname, H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
      {
          fprintf(stderr,"Could not open file %s for reading\n",fname);
          exit(-1);
      }
}

//FORTRAN wrapper
void open_cm1_hdf_file_(hid_t *file_id, char *base, int *itime, int *inode)
{
	open_cm1_hdf_file(file_id,base,*itime,*inode);
}

void
open_cm1_hdf_file_rw (hid_t *file_id, char *base, int itime, int inode)
{
	char fname[200];
	sprintf (fname, "%s%05i_%04i.hdf", base, itime, inode);
	if ((*file_id = H5Fopen (fname, H5F_ACC_RDWR,H5P_DEFAULT)) < 0)
      {
          fprintf(stderr,"Could not open file %s for reading/writing\n",fname);
          exit(-1);
      }
}

/*
void close_cm1_hdf_file (int32 sd_id)
{
    int status;
    if ((status = SDend (sd_id)) == FAIL)
    {
        fprintf(stderr,"close_cm1_hdf_file: OH NO! Can't close hdf file with sd_id = %i\n",sd_id);
    }
}
*/

void
get0dint (hid_t file_id, char *varname, int *var)
{
	hid_t rank, dims[3], start[3], stride[3];
	hid_t dataset_id;
	int status;

	/* shouldn't need these any more */
	rank = 1;
	dims[0] = 1;
	start[0] = 0;
	stride[0] = 1;
	if ((dataset_id = H5Dopen (file_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
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
	hid_t rank, dims[1], start[1], stride[1];
	hid_t dataset_id;
	int status;

	/* shouldn't need these any more */
	rank = 1;
	dims[0] = 1;
	start[0] = 0;
	stride[0] = 1;
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

//	fprintf(stderr,"get1dfloat %i to %i points, %s\n",p0,np,varname);

	if ((dataset_id = H5Dopen (file_id, varname,H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) ERROR_STOP("Could not H5Dget_space");
	if ((memoryspace_id = H5Screate_simple(rank,dims,NULL)) < 0) ERROR_STOP("Could not H5Screate_simple");
	if ((status = H5Sselect_hyperslab (dataspace_id,H5S_SELECT_SET,offset_in,NULL,count,NULL)) < 0) ERROR_STOP("Could not H5Sselect_hyperslab");
	if ((status = H5Sselect_hyperslab (memoryspace_id,H5S_SELECT_SET,offset_out,NULL,count,NULL)) < 0) ERROR_STOP ("Could not H5Sselect_hyperslab");
	if ((status = H5Dread (dataset_id, H5T_NATIVE_FLOAT,memoryspace_id,dataspace_id,H5P_DEFAULT,var)) < 0) ERROR_STOP ("Could not H5Dread");
	if ((status = H5Sclose(memoryspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Sclose(dataspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP("Could not H5Dclose");

//	fprintf(stderr,"Got %s!\n",varname);
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

#define MAXHDF (5000)
/* If the first file is not less than 5000 then it must be at least 2050 AD */

/* Sometimes we need to get metadata before we read in a 3D variable.
Since all of the nodefiles have the same metadata, we need to pick one
to read from. Since we may have culled the files to save on disk space,
we can't always assume '0'. So this function goes through and finds a
file, stating at zero. It breaks at the first file it finds. */

int first_hdf_index(char *base, int itime)
{
    int i;
    char full_file_name[300];

    for (i = 0; i < MAXHDF; i++)
    {
        sprintf (full_file_name, "%s%05i_%06i.cm1hdf5", base, itime, i);
	  itime = 30;//ORF test 5/17/11
	  fprintf(stderr,"Trying to access metadata from %s...\n",full_file_name);
    
        if (access(full_file_name,F_OK)==0)
        {
//		fprintf(stderr,"Got metadata from %s\n",full_file_name);
            break; /* Current value of i will point to a hdf file from which we can get our metadata */
        }
    }

    if (i == MAXHDF)
    {
        fprintf(stderr,"\n\nfirst_hdf_index: WOOP WOOP where are all the HDF files??\n");
        fprintf(stderr,"Unable to find files with base %s at time %i.\nExiting.\n",base,itime);
        exit(-1);
    }
    return i;
}

//FORTRAN wrapper
void first_hdf_index_(char *base, int *itime, int *index)
{
	*index = first_hdf_index(base,*itime);
}


/* I have decided to create a new routine which must be called before
 * any 3D data is read with read_hdf_mult. This will reduce the number
 * of redundant reads to the HDF files for each 3D read. While it
 * requires the user to call an extra function (only once mind you) it's just the right way
 * to do things. The pointers nx,ny,nz,nodex,nodey are updated and are
 * expected to be passed as arguments to read_hdf_mult so read_hdf_mult
 * doesn't have to get them every time. nx,ny,nodex,nodey are all that
 * are required to recreate the entire domain decomposition. */

void get_hdf_metadata(char *hdffilename, int node, int itime, int *nx, int *ny, int *nz, int *nodex, int *nodey)
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
