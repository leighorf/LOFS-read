#include "../include/lofs-read.h"
#include "../include/lofs-dirstruct.h"
#include "../include/lofs-hdf2nc.h"

// See hdfio.c get_hdf_metadata to finish this, just cycle through and
// list all LOFS zfpacc

void list_LOFS_zfpacc (hdf_meta hm, hid_t *f_id)
{
	hid_t g_id,d_id,attr_id,attr_memtype;
	herr_t status;
	H5O_info_t dset_info;
	H5G_info_t group_info;
	int i,j,nattr;
	double zval;
	char groupname[MAXSTR];
	char attrname[MAXATTR][MAXSTR]; //ORF FIX TODO
	htri_t existence;
	hid_t lapl_id;

	sprintf(groupname,"%05i/3D",0);//All vars in 3D group 00000 are always available
	g_id = H5Gopen(*f_id,groupname,H5P_DEFAULT);
	H5Gget_info(g_id,&group_info);
	hm.nvar_available = group_info.nlinks;
	printf("nvar avail = %i\n",hm.nvar_available);
	printf("g_id = %i varname = %s\n",g_id,hm.varname_available[0]);

	for (i = 0; i < hm.nvar_available; i++)
	{
		if ((d_id = H5Dopen (g_id,hm.varname_available[i],H5P_DEFAULT)) < 0) ERROR_STOP("Could not H5Dopen");
		if ((status = H5Oget_info(d_id,&dset_info)) < 0) ERROR_STOP("Could not H5Oget_info");
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
				printf("LOFS variable %10s: zfpacc_LOFS = %14.7f\n",hm.varname_available[i],zval);
				break;
			}
		}
		H5Dclose(d_id);
	}
}
