#include "../include/lofs-dirstruct.h"
#include "../include/lofs-limits.h"
#include "../include/lofs-read.h"

/* is between or lies on actually */
int
is_between_int (int x, int y, int z)
{
	int retval = 0;
	if (z >= x && z <= y) retval = 1;
	return retval;
}
int
is_between (double x, double y, double z)
{
	int retval = 0;
	if (z >= x && z <= y) retval = 1;
	return retval;
}
int
is_between_l (double x, double y, double z)
{
	int retval = 0;
	double eps = 0.001;
//	printf("is_between_l: %lf %lf %lf\n",x,y,z);
	if (z >= x && z < y) retval = 1;
	return retval;
}

// Need to make this foolproof, may need to go beyond just
// is_between_fuzzy to determine which hdf5 file for 'edge' cases
int
is_between_fuzzy (double x, double y, double z)
{
	int retval = 0;
	double eps = 0.001;
//	printf("is_between_fuzzy: %lf %lf %lf\n",x,y,z);
	if (z > x && z < y) retval = 1;
	if (fabs(z-y)<eps) retval = 0; // on the edge of nowhere
	if (fabs(z-x)<eps) retval = 1; // on the edge of somewhere
	return retval;
}
int
is_not_between_int (int x, int y, int z)
{
	int retval = 0;
	if (z < x || z > y) retval = 1;
	return retval;
}


char *replace_str(char *str, char *orig, char *rep)
{
  static char buffer[4096];
  char *p;

  if(!(p = strstr(str, orig)))  // Is 'orig' even in 'str'?
    return str;

  strncpy(buffer, str, p-str); // Copy characters from 'str' start to 'orig' st$
  buffer[p-str] = '\0';

  sprintf(buffer+(p-str), "%s%s", rep, p+strlen(orig));

  return buffer;
}

//2019-04-30 taken from hdf2d2nc.c - iterate over a HDF group, first
//pass gets the number of items, second pass gets the variable names

int n2d,n2dstatic,n2dswath,i2d;
const char **twodvarname;

herr_t twod_first_pass(hid_t loc_id, const char *name, void *opdata)
{
    H5G_stat_t statbuf;
    H5Gget_objinfo(loc_id, name, FALSE, &statbuf);
    n2d++;
    return 0;
}

herr_t twod_second_pass(hid_t loc_id, const char *name, void *opdata)
{
    H5G_stat_t statbuf;
    H5Gget_objinfo(loc_id, name, FALSE, &statbuf);
    if (!strcmp(name,"prespert_min_sfc_move") || !strcmp(name,"hwin_max_sfc_move"))
	{
		strcpy((char *)twodvarname[i2d],name);
	    i2d++;
	}
    return 0;
}

// We now incorporate 2D swaths from 3D files
//
// Trick: Don't have to change API here. Instead:
//
// Make swaths a "special" varname. If "swaths" is chosen, nz is
// repurposed as being the number of swaths... we should have our
// mallocs done so we don't have to worry about doing that in here. For
// now, it's an all or nothing deal - all the swaths or none of them.
// Could be more fine-grained later.
//
// I feel a refactoring coming on. For one, for a given time we can
// cache all the metadata (everything in the hdf structure) - and
// implicit to LOFS you could chose any time and you'd have all the
// metadata for all the times. So some sort of caching mechanism, either
// using global variables or files.
//
// Also now that I've got the swath stuff working it's time to create a
// data structure for all the metadata stuff that is passed redundantly
// to this routine - reduce the number of arguments by a lot. This
// structure should perhaps be world readable so we don't need to keep
// passing it to routines.
//
// To avoid accidental weirdness, I make a buffer exclusively for the 3D
// swath stack rather than reusing the 3D one
//
//extern int debug;

//ORF 2020-02-07 nodedir not used here!!
//void
//read_hdf_mult_md (float *gf, char *topdir, char **timedir, char **nodedir, int ntimedirs, int dn,
//		double *dirtimes, double *alltimes, int ntottimes, double dtime, char *varname,
//		int gx0, int gy0, int gxf, int gyf, int gz0, int gzf,
//		int nx, int ny, int nz, int nodex, int nodey)
//
//
//

void read_lofs_buffer(float *buf, char *varname, dir_meta dm, hdf_meta hm, requested_cube rc, cmdline cmd)
{
	int i, tb;
	int maxfilelength = 512;
	char **nodefile;
	char datasetname[100];
	int numhdf;

	hsize_t count3[3];
	hsize_t offset_in3[3],offset_out3[3];
	hsize_t count2[2];
	hsize_t offset_in2[2],offset_out2[2];
	hid_t file_id,dataset_id,dataspace_id,memoryspace_id;
	hid_t swath_dataset_id,swath_dataspace_id,swath_memoryspace_id, single_swath_memoryspace_id;
	int status;
	int numi, numj;
	int ihdf;

	int ix,iy;

	int gnx, gny, gnz;
	int fx0=0, fxf=0, fy0=0, fyf=0;	/* file indices */
	int sx0, sxf, sy0, syf;
	int ax0, axf, ay0, ayf;
	int snx, sny, snz, dxleft, dxright, dybot, dytop, ixnode, iynode, k, nxnode;
	int rank;
	double eps;
 // filetimes is the actual times stored in the hdf5 file dirtimes
 // contains the list of times contained in the directory names, which
 // corresponds to the first time in the hdf files.
	int time_is_in_file;

 // Some silly 'heartbeat' output - letter written to screen is
 // proportional to fraction of 2D horizontal extent of file read; z
 // means you read the whole horizontal space in the file, a means you
 // read a tiny subset, etc.

	char *alph="abcdefghijklmnopqrstuvwxyz";
	int letter;
	int retval;
	int ntimes,timeindex;
	hsize_t dims[1],maxdims[1];
	double *filetimes;

	typedef struct hdfstruct
	{
		int x0, xf, y0, yf;
		int sx0, sxf, sy0, syf;
		int ax0, axf, ay0, ayf;
		int myi, myj;
	} HDFstruct;

	HDFstruct **hdf;

	if(cmd.verbose)
	{
		printf("\nrc.X0 = %5i\t rc.X1 = %5i\n",rc.X0,rc.X1);
		printf("rc.Y0 = %5i\t rc.Y1 = %5i\n",rc.Y0,rc.Y1);
		printf("rc.Z0 = %5i\t rc.Z1 = %5i\n",rc.Z0,rc.Z1);
		printf("\nReading %s...\n",varname);
	}

	gnx = rc.X1 - rc.X0 + 1;
	gny = rc.Y1 - rc.Y0 + 1;
	gnz = rc.Z1 - rc.Z0 + 1;

	numhdf = hm.rankx * hm.ranky;

	/* ORF 2021-10-18 Note: numhdf will be larger than the actual number
	 * of saved hdf files if the user did not save the full domain. That
	 * is OK because we only allocate short character arrays 
	 * here for creating the file names, so allocating more than we need
	 * is no big whoop */

	if ((hdf = (HDFstruct **) malloc (numhdf * sizeof (HDFstruct *))) == NULL) ERROR_STOP("Insufficient memory for HDFstruct");
	for (i = 0; i < numhdf; i++)
	{
		if ((hdf[i] = (HDFstruct *) malloc (sizeof (HDFstruct))) == NULL) ERROR_STOP("Insufficient memory for HDFstruct");
	}

	if ((nodefile = (char **) malloc (numhdf * sizeof (char *))) == NULL)  ERROR_STOP("Insufficient memory for nodefile");

	eps=1.0e-3;
	if (!is_between_fuzzy(dm.alltimes[0],dm.alltimes[dm.ntottimes-1],cmd.time))
	{
		// edge case: Asking for last time will fail unless we check
		// this
		if(fabs(dm.alltimes[dm.ntottimes-1]-cmd.time)>eps)
		{
			fprintf(stderr,"Out of range: %f must be within the range %f to %f\n",cmd.time,dm.alltimes[0],dm.alltimes[dm.ntottimes-1]);
			ERROR_STOP("Requested time not within range\n");
		}
	}
	if(cmd.debug)
	{
		for (i=0; i < dm.ntimedirs; i++)
			printf("SANITY CHECK: dirtimes = %lf\n",dm.dirtimes[i]);
	}
	for (i=0; i < dm.ntimedirs-1; i++)
	{
		if(cmd.debug) printf("DEBUG: dirtimes[%i] = %lf, dirtimes[%i] = %lf, time = %lf\n",i,dm.dirtimes[i],i+1,dm.dirtimes[i+1],cmd.time);
		if (is_between_fuzzy(dm.dirtimes[i],dm.dirtimes[i+1],cmd.time)) break;
	}
	tb = i;
	if (dm.ntimedirs == 1) tb = 0;
	for (i = 0; i < numhdf; i++)
	{
		if ((nodefile[i] = (char *) malloc (MAXSTR*sizeof(char))) == NULL) ERROR_STOP("Insufficient memory");
		/* Construct the cm1hdf5 file name and path, which is full of tasty metadata! */
		sprintf (nodefile[i], "%s/%s/%07i/%s_%07i.cm1hdf5", dm.topdir,dm.timedir[tb],(dm.dn!=-1)?((i/dm.dn)*dm.dn):0,dm.timedir[tb],i);
	}

	eps = 1.0e-3;

	time_is_in_file = FALSE;
	for (i=0; i<dm.ntottimes; i++)
	{
		if (fabs(cmd.time-dm.alltimes[i])<eps)
		{
			time_is_in_file = TRUE;
			break;
		}
	}

	if (time_is_in_file == FALSE)
	{
		fprintf(stderr,"Requested time %lf was not saved\n",cmd.time);
		fprintf(stderr,"Available times follow:");
		for (i=0; i<dm.ntottimes; i++)
		{
			fprintf(stderr," %f",dm.alltimes[i]);
		}
		fprintf(stderr,"\n");
	}
	if (time_is_in_file == FALSE)	ERROR_STOP("Invalid time requested");

	/* We build our decomposition from metadata stored in each hdf file */

	numi = hm.nx / hm.rankx;
	numj = hm.ny / hm.ranky;
	for (ihdf = 0; ihdf < numhdf; ihdf++) /* just i not ihdf */
	{
		hdf[ihdf]->myj = ihdf / hm.rankx;
		hdf[ihdf]->myi = ihdf % hm.rankx;
		hdf[ihdf]->x0 = hdf[ihdf]->myi * numi;
		hdf[ihdf]->xf = (hdf[ihdf]->myi + 1) * numi - 1;
		hdf[ihdf]->y0 = hdf[ihdf]->myj * numj; 
		hdf[ihdf]->yf = (hdf[ihdf]->myj + 1) * numj - 1;
		if (cmd.debug)
			fprintf (stderr, "myj = %i myi =%i x0 = %i xf = %i y0 = %i yf = %i\n",
				 hdf[ihdf]->myj, hdf[ihdf]->myi, hdf[ihdf]->x0, hdf[ihdf]->xf, hdf[ihdf]->y0, hdf[ihdf]->yf);
	}

	/* first check if our requested subcube lies within our data */

	if (is_not_between_int (0, hm.nx - 1, rc.X0)) ERROR_STOP("Chosen x data out of range");
	if (is_not_between_int (0, hm.ny - 1, rc.Y0)) ERROR_STOP("Chosen y data out of range");
	if (is_not_between_int (0, hm.nz - 1, rc.Z0)) ERROR_STOP("Chosen z data out of range");

	if (is_not_between_int (0, hm.nx - 1, rc.X1)) ERROR_STOP("Chosen x data out of range");
	if (is_not_between_int (0, hm.ny - 1, rc.Y1)) ERROR_STOP("Chosen y data out of range");
//ORF: we don't do this check for swaths, they are handled differently,
//but we are using the z dimension
	if ((strcmp(varname,"swaths")) && is_not_between_int (0, hm.nz - 1, rc.Z1)) ERROR_STOP("Chosen z data out of range");

	for (i = 0; i < numhdf; i++)
	{
		if (is_between_int (hdf[i]->x0, hdf[i]->xf, rc.X0) && is_between_int (hdf[i]->y0, hdf[i]->yf, rc.Y0))
		{
			fx0 = hdf[i]->myi;
			fy0 = hdf[i]->myj;
			if (cmd.debug) fprintf (stderr, "found fx0,fy0 = %i,%i\n", fx0, fy0);
		}
		if (is_between_int (hdf[i]->x0, hdf[i]->xf, rc.X1) && is_between_int (hdf[i]->y0, hdf[i]->yf, rc.Y1))
		{
			fxf = hdf[i]->myi;
			fyf = hdf[i]->myj;
			if (cmd.debug) fprintf (stderr, "found fxf,fyf = %i,%i\n", fxf, fyf);
		}
	}

	snx = numi;
	sny = numj;
	dxleft = rc.X0 % snx;
	dxright = rc.X1 % snx;
	dybot = rc.Y0 % sny;
	dytop = rc.Y1 % sny;
	nxnode = hm.rankx;

	/*

	   There are four coordinate systems to worry about:

	   1. Node (spans fx0->fxf and fy0->fyf) - indexes the hdf files we read.
	   2. Global (spans 0->nx-1 and 0->ny-1) - the whole model domain.  Calling program provides indices relative to global.
	   3. Requested array (spans 0->gnx-1 and 0->gny-1) - is passed back to calling program.
	   4. Node-relative array (spans 0->snx-1 and 0->sny-1 - what we ultimately must read from each hdf file.

	   Note: It is assumed each hdf file contains all points in z which is
	   the usual way these decompositions are done.

	   First do x pass through indices. Logic is easier with one x pass
	   and one y pass rather than doing both in the same loops

	 */

	/* x pass */
	for (iynode = fy0; iynode <= fyf; iynode++)
	{
		for (ixnode = fx0; ixnode <= fxf; ixnode++)
		{
			k = ixnode + iynode * nxnode;	/* since hdf struct is 1D array */

			if (fx0 == fxf)	/* if we're only on one file */
			{
				sx0 = dxleft;
				sxf = dxright;
				ax0 = 0;
				axf = dxright - dxleft;
			}
			else if (ixnode == fx0)	/* on left edge (x increases to right y increases up) */
			{
				sx0 = dxleft;
				sxf = snx - 1;
				ax0 = 0;
				axf = snx - 1 - dxleft;
			}
			else if (ixnode == fxf)	/* on right edge */
			{
				sx0 = 0;
				sxf = dxright;
				ax0 = (snx - dxleft) + snx * (ixnode - fx0 - 1);
				axf = ax0 + dxright;
			}
			else	/* in middle */
			{
				sx0 = 0;
				sxf = snx - 1;
				ax0 = (snx - dxleft) + snx * (ixnode - fx0 - 1);
				axf = ax0 + snx - 1;
			}
			hdf[k]->sx0 = sx0;
			hdf[k]->sxf = sxf;
			hdf[k]->ax0 = ax0;
			hdf[k]->axf = axf;
		}
	}

	/* y pass */

	for (iynode = fy0; iynode <= fyf; iynode++)
	{
		for (ixnode = fx0; ixnode <= fxf; ixnode++)
		{
			k = ixnode + iynode * nxnode;	/* since hdf struct is 1D array */

			if (fy0 == fyf)	/* if we're only on one file */
			{
				sy0 = dybot;
				syf = dytop;
				ay0 = 0;
				ayf = dytop - dybot;
			}
			else if (iynode == fy0)	/* on bottom edge */
			{
				sy0 = dybot;
				syf = sny - 1;
				ay0 = 0;
				ayf = sny - 1 - dybot;
			}
			else if (iynode == fyf)	/* on top edge */
			{
				sy0 = 0;
				syf = dytop;
				ay0 = (sny - dybot) + sny * (iynode - fy0 - 1);
				ayf = ay0 + dytop;
			}
			else	/* in middle */
			{
				sy0 = 0;
				syf = sny - 1;
				ay0 = (sny - dybot) + sny * (iynode - fy0 - 1);
				ayf = ay0 + sny - 1;
			}
			hdf[k]->sy0 = sy0;
			hdf[k]->syf = syf;
			hdf[k]->ay0 = ay0;
			hdf[k]->ayf = ayf;
		}
	}

	/* We need to split the above part of the code off and only call it
	 * once while using this program. Currently we recalculate this
	 * every time when it's not necessary; further we should shove all
	 * the metadata goodness into a structure... maybe just add to the
	 * current HDFStruct... and make it a global structure (TODO)) */

/* Now we have collected all of our array index values delineating node
relative, domain relative, and retrieved-array relative positions. Need
to read data from each node file, assemble retrieved array and send back
as a 1D buffer. Sine this is C and all, 3D arrays are already 1D arrays,
really. See P3 macro in lofs-read.h */

// ORF 7/7/09
// HDF5 isn't just a data format, it's a way of life.
//
// Save ourselves array assembling issues by creating a memory space with the
// dimensions of the array we are sending back. This is actually pretty cool.
//
// First, initialize ZFP lossy floating point compression

//	H5Z_zfp_initialize();

	rank=3;
	count3[0]=gnz;count3[1]=gny;count3[2]=gnx;
	status=memoryspace_id = H5Screate_simple(rank,count3,NULL);
	if (status < 0) ERROR_STOP("H5Screate_simple failed");
	if(cmd.verbose)printf("\n");
	for (iynode = fy0; iynode <= fyf; iynode++)
	{
		for (ixnode = fx0; ixnode <= fxf; ixnode++)
		{
			k = ixnode + iynode * nxnode;

			if (cmd.debug) printf("Working on %s\n",nodefile[k]);
			if ((file_id = H5Fopen (nodefile[k], H5F_ACC_RDONLY,H5P_DEFAULT)) < 0)
			{
				fprintf (stderr, "\n\nread_hdf_mult: Could not start to read %s!\n", nodefile[k]);
				fprintf (stderr, "Perhaps:\n");
				fprintf (stderr, "\t(a) %s does not exist\n",nodefile[k]);
				fprintf (stderr, "\t(b) you are asking for data outside the range of the model domain\n");
				fprintf (stderr, "\t(c) you are asking for data from hdf files you culled to save disk space\n");
				fprintf (stderr, "\t(d) %s is a corrupt file\n",nodefile[k]);
				fprintf (stderr, "\nExiting.\n");
				exit (-1);
			}
			if (cmd.debug) printf("Varname = %s\n",varname);

//			printf("DEBUG: nodefile = %s\n",nodefile[k]);

	/* ORF 2016-12-21
	 * Because of single precision floating point time variables in
	 * cm1r16, we have junk in the least significant digits of our
	 * time that carry over into our directory names and floating
	 * point links. It occurs to me that we will proablby never save more
	 * frequently than every 0.1 seconds or thereabouts and that we
	 * could probably just remove all those insigificant digtis (store
	 * it as X.YYY). But, first, I am going to try just doing floating
	 * point math with the /times array and pick the right group to
	 * read from, since the groups are zero padded integers (like the
	 * old way) (converted to a character string) starting at zero, and the 
	 * floating point links were done as a test (we don't need to use
	 * them). If performance suffers from having to read the /times
	 * group every time, I can consider changing the significant
	 * digits.
	 * Unfortunately I see no easy way to have this information cached
	 * since each hdf5 dump can contain any amount of times stacked
	 * within it, and up to now just sorting the times has been
	 * sufficient... but there is no way to preserve the group name
	 * and the time while sorting the way it is currently done.
	 * I'm getting closer to just requiring a more advanced "makevisit" being run once
	 * and having all my shit rely on the existence of that file...
	 * but I'm not quite there yet...
	 *
	 * Note from the future (2019-05-01): We haven't had problems with the floating
	 * point stuff and we've been saving data as frequently as every
	 * 1/6 second. So maybe we are cool.
	 *
	 * */


/* Since our data format requires identical temporal data for all HDF
 * files in the same directory, we just do this once, in the first pass
 * through.
 *
 * 2019-04-30 Here we also collect all of the new 2D snapshot/swath names and the total number.
 * This is now stored in the LOFS HDF5 files */

			if (iynode==fy0 && ixnode==fx0)
			{
				if ((dataset_id = H5Dopen(file_id,"times",H5P_DEFAULT)) < 0) ERROR_STOP("Cannot H5Dopen");
				if ((dataspace_id = H5Dget_space(dataset_id)) < 0) ERROR_STOP("Cannot H5Dget_space");
				if ((H5Sget_simple_extent_dims(dataspace_id,dims,maxdims)) < 0)
					ERROR_STOP("Cannot H5Sget_simple_extent_dims"); //dims[0] will equal number of time levels in file
				H5Dclose (dataset_id);
				H5Sclose (dataspace_id);
				if ((filetimes = (double *)malloc(dims[0]*sizeof(double)))==NULL) ERROR_STOP("Cannot malloc filetimes array");
				get1ddouble(file_id,"times",filetimes,0,dims[0]);
				/* Now we have our times array, the size of the times
				 * array (the number of times in the file) and the
				 * requested time. There is probably a faster way to
				 * do this but the truth is we will never have more
				 * than say 500 array elements to consider so we might
				 * as well just rip through one by one*/

				/* The actual groups in the HDF files are zero padded,
				 * starting with zero, incrementing by 1. So getting
				 * to the right time is easy peasy */

				ntimes = dims[0];

				for (i=0; i<ntimes; i++)
				{
//					printf("i = %i \t filetimes[%2i] = %f \t time = %f\n",i,i,filetimes[i],time);
					if (fabs(filetimes[i]-cmd.time)<eps)
						break;
				}

				/* Heh, this never fails, it will always send back the
				 * last time in the file if we don't find our time */

				timeindex = i;

				/* 2019-05-01. New swath code.

				We are still in the "only do once" part of the 2d loop. We must get
				our 2d information (number, names) and memoryspace worked out. THEN
				break from the "only do once" part and do the rest, still will have
				to check if doing 2D
				*/

				if (!strcmp(varname,"swaths")) //"swaths" is a special variable name here, kind of kludgey...
				{
					n2d = 0; // number of 2D swaths
					H5Giterate(file_id, "/00000/2D/static",NULL,twod_first_pass,NULL);
					n2dstatic = n2d; // number of "static" swaths (the stuff I added like dbz_500m)
					H5Giterate(file_id, "/00000/2D/swath",NULL,twod_first_pass,NULL);
/* OVERRIDE for Larry */
					n2d=2;
					n2dstatic=0;
					n2dswath = n2d - n2dstatic; // number of actual swaths (George's, and the bunch I added)


//					printf("n2dstatic = %i n2dswath = %i n2d = %i\n",n2dstatic,n2dswath,n2d);

					twodvarname = (const char **)malloc(n2d*sizeof(char *));

					// Our memoryspace is a 3D array of
					// dimension [n2d][gny][gnx]
					// n2d = number of swaths
					// gny, gnx is as it has always been, see
					// above "four coordinate systems" comment
					rank=3;
					count3[0]=n2d;count3[1]=gny;count3[2]=gnx;
					swath_memoryspace_id = H5Screate_simple(rank,count3,NULL);

					for (i2d=0; i2d<n2d; i2d++)
					{
						twodvarname[i2d] = (char *)malloc(50*sizeof(char)); // 50 characters per variable
					} //ORF TODO MAKE 50 A CONSTANT

					i2d=0;
					H5Giterate(file_id, "/00000/2D/static",NULL,twod_second_pass,NULL);
					H5Giterate(file_id, "/00000/2D/swath",NULL,twod_second_pass,NULL);

					/* Unlike with hdf2.c, we iterate to only get
					 * a list of the swath variable names, not
					 * set up an id array because HDF5 has a
					 * different way... we'll get access to
					 * the data a bit further down */
				}
			}

			/* Now we are in our 2D loop, after 1st pass stuff above */

			if (!strcmp(varname,"swaths"))
			{
				snx = hdf[k]->sxf - hdf[k]->sx0 + 1;
				sny = hdf[k]->syf - hdf[k]->sy0 + 1;
				offset_in2[0] = hdf[k]->sy0;
				offset_in2[1] = hdf[k]->sx0;
				offset_out3[1] = hdf[k]->ay0;
				offset_out3[2] = hdf[k]->ax0;
				count2[0] = sny;
				count2[1] = snx;
				count3[0] = 1;
				count3[1] = sny;
				count3[2] = snx;

				for (i2d = 0; i2d < n2dstatic; i2d++)
				{
					sprintf(datasetname,"/%05i/2D/static/%s",timeindex,twodvarname[i2d]);

					if ((swath_dataset_id = H5Dopen (file_id, datasetname, H5P_DEFAULT)) < 0)
					{
							fprintf(stderr,"Cannot open dataset %s in file %s\n",datasetname,nodefile[k]);
							ERROR_STOP("H5Dopen failed");
					}

					swath_dataspace_id = H5Dget_space(swath_dataset_id);
					rank = H5Sget_simple_extent_ndims(swath_dataspace_id);
					if (rank != 2) ERROR_STOP("Rank has to equal 2 - something hideously fucked up, goodbye!");

					offset_out3[0] = i2d;

					if(i2d == 0 && cmd.verbose) //seriously you want this
					{
						letter = (int)(25.0*(float)(snx*sny)/(float)(numi*numj));
						printf ("%c",alph[letter]); fflush (stdout);
					}

					status=H5Sselect_hyperslab (swath_dataspace_id,H5S_SELECT_SET,offset_in2,NULL,count2,NULL); if (status < 0) ERROR_STOP("select_hyperslab");
					status=H5Sselect_hyperslab (swath_memoryspace_id,H5S_SELECT_SET,offset_out3,NULL,count3,NULL); if (status < 0) ERROR_STOP("select_hyperslab");
					status=H5Dread (swath_dataset_id,H5T_NATIVE_FLOAT,swath_memoryspace_id,swath_dataspace_id,H5P_DEFAULT,buf); if (status < 0) ERROR_STOP("h5dread");
					H5Dclose (swath_dataset_id);
					H5Sclose (swath_dataspace_id);
				}

				for (i2d = n2dstatic; i2d < n2d; i2d++)
				{
					sprintf(datasetname,"/%05i/2D/swath/%s",timeindex,twodvarname[i2d]);

					if ((swath_dataset_id = H5Dopen (file_id, datasetname, H5P_DEFAULT)) < 0)
					{
							fprintf(stderr,"Cannot open dataset %s in file %s\n",datasetname,nodefile[k]);
							ERROR_STOP("H5Dopen failed");
					}

					swath_dataspace_id = H5Dget_space(swath_dataset_id);
					rank = H5Sget_simple_extent_ndims(swath_dataspace_id);
					if (rank != 2) ERROR_STOP("Rank has to equal 2 - something hideously fucked up, goodbye!");

					offset_out3[0] = i2d;

					status=H5Sselect_hyperslab (swath_dataspace_id,H5S_SELECT_SET,offset_in2,NULL,count2,NULL); if (status < 0) ERROR_STOP("select_hyperslab");
					status=H5Sselect_hyperslab (swath_memoryspace_id,H5S_SELECT_SET,offset_out3,NULL,count3,NULL); if (status < 0) ERROR_STOP("select_hyperslab");
					status=H5Dread (swath_dataset_id,H5T_NATIVE_FLOAT,swath_memoryspace_id,swath_dataspace_id,H5P_DEFAULT,buf); if (status < 0) ERROR_STOP("h5dread");
					H5Dclose (swath_dataset_id);
					H5Sclose (swath_dataspace_id);
				}
			}
			else
			{
				sprintf(datasetname,"/%05i/3D/%s",timeindex,varname);
				if ((dataset_id = H5Dopen (file_id, datasetname, H5P_DEFAULT)) < 0)
				{
						fprintf(stderr,"Cannot open dataset %s in file %s\n",datasetname,nodefile[k]);
						ERROR_STOP("H5Dopen failed");
				}

				dataspace_id = H5Dget_space(dataset_id);
				rank = H5Sget_simple_extent_ndims(dataspace_id);
				if (rank == 3)
				{
					snx = hdf[k]->sxf - hdf[k]->sx0 + 1; 
					sny = hdf[k]->syf - hdf[k]->sy0 + 1;
					snz = rc.Z1 - rc.Z0 + 1;
					offset_in3[0] = rc.Z0;
					offset_in3[1] = hdf[k]->sy0;
					offset_in3[2] = hdf[k]->sx0;
					offset_out3[0] = 0; //orf 6/30/10 gz0 was a bug? should always be 0?
					offset_out3[1] = hdf[k]->ay0;
					offset_out3[2] = hdf[k]->ax0;
					count3[0] = snz;
					count3[1] = sny;
					count3[2] = snx;
				}
				else
				{
					ERROR_STOP("Rank should be 3, WTF?");
				}
				if(cmd.verbose)
				{
					letter = (int)(25.0*(float)(snx*sny)/(float)(numi*numj));
					printf ("%c",alph[letter]); fflush (stdout);
				}
				if ((H5Sselect_hyperslab (dataspace_id,H5S_SELECT_SET,offset_in3,NULL,count3,NULL)) < 0)
				{
						fprintf(stderr,"\nCannot select hyperslab dataspace_id %i in file %s\n",dataspace_id,nodefile[k]);
						H5Eprint(H5E_DEFAULT,NULL);
						ERROR_STOP("H5Sselect_hyperslab failed");
				}
				if ((H5Sselect_hyperslab (memoryspace_id,H5S_SELECT_SET,offset_out3,NULL,count3,NULL)) < 0)
				{
						fprintf(stderr,"\nCannot select hyperslab memoryspace_id %i in file %s\n",memoryspace_id,nodefile[k]);
						H5Eprint(H5E_DEFAULT,NULL);
						ERROR_STOP("H5Sselect_hyperslab failed");
				}
				if ((retval=H5Dread (dataset_id,H5T_NATIVE_FLOAT,memoryspace_id,dataspace_id,H5P_DEFAULT,buf)) < 0)
				{
						fprintf(stderr,"\nCannot read hyperslab for %s in %s\n",datasetname,nodefile[k]);
						fprintf(stderr,"return value = %i\n",retval);
						H5Eprint(H5E_DEFAULT,NULL);
						ERROR_STOP("H5Dread failed");
				}
				
				H5Dclose (dataset_id);
				H5Sclose (dataspace_id);
			}
			H5Fclose (file_id);
		}
		if(cmd.verbose)printf("\n");
	}
    if (!strcmp(varname,"swaths")) {
	    H5Sclose (swath_memoryspace_id);
    }
    H5Sclose (memoryspace_id);
	/* free all pointers */
	for (i = 0; i < numhdf; i++)
	{
		free (hdf[i]);
		free (nodefile[i]);
	}
	free (hdf);
	free (nodefile);
	if(cmd.verbose)printf("\n");
}
