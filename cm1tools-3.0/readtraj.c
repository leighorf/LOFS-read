#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "errorf.h"

main (int argc, char **argv)
{

	hid_t file_id;
	int status;
	int prank,trank,ttimes,nparcels,iparcel,tlevel;
	hsize_t count[3];
	hsize_t offset_in[3],offset_out[3];
	hsize_t dims[3], maxdims[3];
	hid_t dataset_id,dataspace_id,memoryspace_id;

	float *time;
	float *zc,*w;
	int i;
	char dname[50];


	if ((file_id = H5Fopen (argv[1], H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
	{
		fprintf(stderr,"\n\nget_hdf_metadata: Unable to read metadata from %s, bailing!\n", argv[0]);
		exit(0);
	}

	if ((dataset_id = H5Dopen(file_id,"/parcel/time",H5P_DEFAULT)) < 0) ERROR_STOP("Cannot H5Dopen");
	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) ERROR_STOP("Cannot H5Dget_space");
	if ((trank = H5Sget_simple_extent_ndims(dataspace_id)) != 1) ERROR_STOP("/parcel/time should have a rank of 1, OH NO!");
	if ((H5Sget_simple_extent_dims(dataspace_id,dims,maxdims)) < 0) ERROR_STOP("Cannot H5Sget_simple_extent_dims"); //dims[0] will equal number of time levels in file
	ttimes = dims[0];
	if ((time = (float *)malloc(ttimes*sizeof(float)))==NULL) ERROR_STOP("Cannot malloc time array");
	if ((status = H5Dread (dataset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,time)) < 0) ERROR_STOP ("Could not H5Dread");
//	for (i=0; i < ttimes-1; i++) printf("%f ",time[i+1]-time[i]);printf("\n");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP ("Could not H5Dclose");
	

//Get number of parcels first

	strcpy(dname,"/parcel/w");
	if ((dataset_id = H5Dopen(file_id,dname,H5P_DEFAULT)) < 0) ERROR_STOP("Cannot H5Dopen");
	if ((dataspace_id = H5Dget_space(dataset_id)) < 0) ERROR_STOP("Cannot H5Dget_space");
	if ((prank = H5Sget_simple_extent_ndims(dataspace_id)) != 2) ERROR_STOP("parcel data should have a rank of 2, OH NO!");
	if ((H5Sget_simple_extent_dims(dataspace_id,dims,maxdims)) < 0) ERROR_STOP("Cannot H5Sget_simple_extent_dims"); //dims[0] will equal number of time levels in file
	if ((status = H5Sclose(dataspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
	if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP ("Could not H5Dclose");

	nparcels = dims[1];
	if ((w = (float *)malloc(nparcels*sizeof(float)))==NULL) ERROR_STOP("Cannot malloc w array");
	if ((zc = (float *)malloc(nparcels*sizeof(float)))==NULL) ERROR_STOP("Cannot malloc zc array");

	for (tlevel = 0; tlevel < dims[0]; tlevel++)
	{
		//printf("Requesting time level %i which is %f seconds\n",tlevel,time[tlevel]);

		strcpy(dname,"/parcel/w");
		if ((dataset_id = H5Dopen(file_id,dname,H5P_DEFAULT)) < 0) ERROR_STOP("Cannot H5Dopen");
		if ((dataspace_id = H5Dget_space(dataset_id)) < 0) ERROR_STOP("Cannot H5Dget_space");
		if ((prank = H5Sget_simple_extent_ndims(dataspace_id)) != 2) ERROR_STOP("parcel data should have a rank of 2, OH NO!");
		if ((H5Sget_simple_extent_dims(dataspace_id,dims,maxdims)) < 0) ERROR_STOP("Cannot H5Sget_simple_extent_dims"); //dims[0] will equal number of time levels in file

	//	printf("%s dims: %i parcels and %i time levels\n",dname,dims[1],dims[0]);


	//tlevel = 210;


		offset_in[0] = tlevel;
		offset_in[1] = 0;
		count[0] = 1;
		count[1] = nparcels;
		dims[0] = 1;
		dims[1] = nparcels;
		offset_out[0] = 0;
		offset_out[0] = 0;

		// Can reuse memory space but can't reuse dataspace as it is linked to floating point data in dataset on disk

		if ((memoryspace_id = H5Screate_simple(prank,dims,NULL)) < 0) ERROR_STOP("Could not H5Screate_simple");
		if ((status = H5Sselect_hyperslab (dataspace_id,H5S_SELECT_SET,offset_in,NULL,count,NULL)) < 0) ERROR_STOP ("Could not H5Sselect_hyperslab");
		if ((status = H5Sselect_hyperslab (memoryspace_id,H5S_SELECT_SET,offset_out,NULL,count,NULL)) < 0) ERROR_STOP ("Could not H5Sselect_hyperslab");
		if ((status = H5Dread (dataset_id, H5T_NATIVE_FLOAT,memoryspace_id,dataspace_id,H5P_DEFAULT,w)) < 0) ERROR_STOP ("Could not H5Dread");
	//	if ((status = H5Sclose(memoryspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
		if ((status = H5Sclose(dataspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
		if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP ("Could not H5Dclose");

		strcpy(dname,"/parcel/z");
		if ((dataset_id = H5Dopen(file_id,dname,H5P_DEFAULT)) < 0) ERROR_STOP("Cannot H5Dopen");
		if ((dataspace_id = H5Dget_space(dataset_id)) < 0) ERROR_STOP("Cannot H5Dget_space");
		if ((prank = H5Sget_simple_extent_ndims(dataspace_id)) != 2) ERROR_STOP("parcel data should have a rank of 2, OH NO!");
		if ((H5Sget_simple_extent_dims(dataspace_id,dims,maxdims)) < 0) ERROR_STOP("Cannot H5Sget_simple_extent_dims"); //dims[0] will equal number of time levels in file
	//	printf("%s dims: %i parcels and %i time levels\n",dname,dims[1],dims[0]);

		nparcels = dims[1];

	//	if ((zc = (float *)malloc(nparcels*sizeof(float)))==NULL) ERROR_STOP("Cannot malloc zc array");

		if ((status = H5Sselect_hyperslab (dataspace_id,H5S_SELECT_SET,offset_in,NULL,count,NULL)) < 0) ERROR_STOP ("Could not H5Sselect_hyperslab");
		if ((status = H5Dread (dataset_id, H5T_NATIVE_FLOAT,memoryspace_id,dataspace_id,H5P_DEFAULT,zc)) < 0) ERROR_STOP ("Could not H5Dread");
		if ((status = H5Sclose(memoryspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
		if ((status = H5Sclose(dataspace_id)) < 0) ERROR_STOP("Could not H5Sclose");
		if ((status = H5Dclose (dataset_id)) < 0) ERROR_STOP ("Could not H5Dclose");


		for (i = 0; i < nparcels; i++)
		{
			if (w[i] > 12.0 || w[i] < -12.0) printf("time = %7.3f\t z = %10.5f w = %6.3f parcel id = %i\n",time[tlevel],zc[i],w[i],i);
		}
	}
}
