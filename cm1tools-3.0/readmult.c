/*
    This file is part of cm1tools, written by Leigh Orf (http://orf5.com)

    cm1tools is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cm1tools is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cm1tools.  If not, see <http://www.gnu.org/licenses/>.

*/


/*

Leigh Orf 8/31/2009

Function:

read_mult_hdf_md is a function which reads in data from multiple 3D
model HDF files as if they were one big file, assembling requested data
into a floating point buffer. Code works with a modified version of
Geroge Bryan's CM1 v16 model (email leigh.orf@cmich.edu for I/O code).

*/

#include "hdforf.h"
#include "errorf.h"
#define P3(x,y,z,mx,my) (((z)*(mx)*(my))+((y)*(mx))+(x))

#define TRUE  (1)
#define FALSE (0)

/* is between or lies on actually */
int
is_between (int x, int y, int z)
{
	int retval = 0;
	if (z >= x && z <= y) retval = 1;
	return retval;
}
int
is_between_l (int x, int y, int z)
{
	int retval = 0;
	if (z >= x && z < y) retval = 1;
	return retval;
}

int
is_not_between (int x, int y, int z)
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

void
link_hdf_files (char *topdir, char **timedir, char **nodedir, int ntimedirs, int dn, int *dirtimes, int *alltimes, int ntottimes,
		int gx0, int gy0, int gxf, int gyf,
		int nx, int ny, int nz, int nodex, int nodey)
{
	int dbg = 0;
	int i, j, tb;
	int maxfilelength = 512;
	char **nodefile;
	char *linkname,*subtopdir;
	char metafile[100];
	char datasetname[100];
	int idum,numhdf;
	char ctime[6];

	int retval;
	hsize_t count[3], dims[3];
	hsize_t offset_in[3],offset_out[3];
	hid_t file_id,dataset_id,dataspace_id,memoryspace_id;
	int status;
	float time;
	int numi, numj;
	int ix, iy, iz;
	int ihdf;

	int gnx, gny, gnz;
	int fx0, fxf, fy0, fyf;	/* file indices */
	int sx0, sxf, sy0, syf;
	int ax0, axf, ay0, ayf;
	int ax, ay, az;
	int snx, sny, snz, dxleft, dxright, dybot, dytop, ixnode, iynode, k, nxnode, nynode;
	short int *bufi, *bufi0;
	float *buff, *buff0;
	float val;
	int rank;
	int ir=0,jr;
 // filetimes is the actual times stored in the hdf5 file dirtimes
 // contains the list of times contained in the directory names, which
 // corresponds to the first time in the hdf files.
	int *filetimes;
	int time_is_in_file;

 // Some silly output - letter is proportional to fraction of 2D
 // horizontal extent of file read; z means you read the whole
 // horizontal space in the file, a means you read a tiny subset.

	char *alph="abcdefghijklmnopqrstuvwxyz";
	int letter;

	char q_varname[100];

	typedef struct hdfstruct
	{
		int x0, xf, y0, yf;
		int sx0, sxf, sy0, syf;
		int ax0, axf, ay0, ayf;
		int myi, myj;
	} HDFstruct;

	HDFstruct **hdf;

	gnx = gxf - gx0 + 1;
	gny = gyf - gy0 + 1;

	numhdf = nodex * nodey;

	linkname = (char *)malloc(512*sizeof(char));
	subtopdir = (char *)malloc(512*sizeof(char));

	if ((hdf = (HDFstruct **) malloc (numhdf * sizeof (HDFstruct *))) == NULL) ERROR_STOP("Insufficient memory");
	for (i = 0; i < numhdf; i++)
	{
		if ((hdf[i] = (HDFstruct *) malloc (sizeof (HDFstruct))) == NULL) ERROR_STOP("Insufficient memory");
	}

	if ((nodefile = (char **) malloc (numhdf * sizeof (char *))) == NULL)  ERROR_STOP("Insufficient memory");


	printf("ntimedirs = %i\n",ntimedirs);
	for (j=0; j < ntimedirs; j++)//ORF changed
	{
	tb = j;
	if (ntimedirs == 1) tb = 0;
	for (i = 0; i < numhdf; i++)
	{
		int fooi;
		if ((nodefile[i] = (char *) malloc (maxfilelength*sizeof(char))) == NULL) ERROR_STOP("Insufficient memory");
		sprintf (nodefile[i], "%s/%s/%06i/%s_%06i.cm1hdf5", topdir,timedir[tb],(dn!=-1)?((i/dn)*dn):0,timedir[tb],i);
	}

	numi = nx / nodex;
	numj = ny / nodey;
	for (ihdf = 0; ihdf < numhdf; ihdf++) /* just i not ihdf */
	{
		hdf[ihdf]->myj = ihdf / nodex;
		hdf[ihdf]->myi = ihdf % nodex;
		hdf[ihdf]->x0 = hdf[ihdf]->myi * numi;
		hdf[ihdf]->xf = (hdf[ihdf]->myi + 1) * numi - 1;
		hdf[ihdf]->y0 = hdf[ihdf]->myj * numj; 
		hdf[ihdf]->yf = (hdf[ihdf]->myj + 1) * numj - 1;
		if (dbg)
			fprintf (stderr, "myj = %i myi =%i x0 = %i xf = %i y0 = %i yf = %i\n",
				 hdf[ihdf]->myj, hdf[ihdf]->myi, hdf[ihdf]->x0, hdf[ihdf]->xf, hdf[ihdf]->y0, hdf[ihdf]->yf);
	}

	/* first check if our requested subcube lies within our data */

	if (is_not_between (0, nx - 1, gx0)) ERROR_STOP("Chosen x data out of range");
	if (is_not_between (0, ny - 1, gy0)) ERROR_STOP("Chosen y data out of range");

	if (is_not_between (0, nx - 1, gxf)) ERROR_STOP("Chosen x data out of range");
	if (is_not_between (0, ny - 1, gyf)) ERROR_STOP("Chosen y data out of range");

	for (i = 0; i < numhdf; i++)
	{
		if (is_between (hdf[i]->x0, hdf[i]->xf, gx0) && is_between (hdf[i]->y0, hdf[i]->yf, gy0))
		{
			fx0 = hdf[i]->myi;
			fy0 = hdf[i]->myj;
			if (dbg) fprintf (stderr, "found fx0,fy0 = %i,%i\n", fx0, fy0);
		}
		if (is_between (hdf[i]->x0, hdf[i]->xf, gxf) && is_between (hdf[i]->y0, hdf[i]->yf, gyf))
		{
			fxf = hdf[i]->myi;
			fyf = hdf[i]->myj;
			if (dbg) fprintf (stderr, "found fxf,fyf = %i,%i\n", fxf, fyf);
		}
	}

	snx = numi;
	sny = numj;
	dxleft = gx0 % snx;
	dxright = gxf % snx;
	dybot = gy0 % sny;
	dytop = gyf % sny;
	nxnode = nodex;
	nynode = nodey;
// ORF AMS2014 help determine subdomain node numbers
	printf("fx0 fy0 fxf fyf: %i %i %i %i\n",fx0,fy0,fxf,fyf);

	printf("ORF: debug: fy0 = %i fyf = %i fx0 = %i fxf = %i\n",fy0,fyf,fx0,fxf);

	subtopdir = replace_str(topdir,"3D","sub3D");
	retval=rename("sub3D",subtopdir);

//	printf("ORF: debug: Just did this: mv(%s,%s)\n","sub3D",subtopdir);
//	if(retval!=0 )
//	{
//		printf("mv failure: errno = %i %s\n",errno,subtopdir);
//		ERROR_STOP("Failed to mv sub3D directory");
//	}

	for (iynode = fy0; iynode <= fyf; iynode++)
	{
		for (ixnode = fx0; ixnode <= fxf; ixnode++)
		{
			k = ixnode + iynode * nxnode;

			linkname = replace_str(nodefile[k],"3D","sub3D");

			retval=link(nodefile[k],linkname);

//			printf("ORF: debug: Just did this: link(%s,%s)\n",nodefile[k],linkname);
			if(retval!=0 && errno!=EEXIST)
			{
				printf("Link failure: errno = %i %s %s\n",errno,nodefile[k],linkname);
				ERROR_STOP("Failed to create hard link");
			}
//			printf ("linking: %s %s\n",retval,nodefile[k],linkname);
		}
	}
	}
	/* free all pointers */
	for (i = 0; i < numhdf; i++)
	{
		free (hdf[i]);
		free (nodefile[i]);
	}
	free (hdf);
	free (nodefile);
}

void
read_hdf_mult_md (float *gf, char *topdir, char **timedir, char **nodedir, int ntimedirs, int dn, int *dirtimes, int *alltimes, int ntottimes, int itime, char *varname,
		int gx0, int gy0, int gxf, int gyf, int gz0, int gzf,
		int nx, int ny, int nz, int nodex, int nodey)
{
	int dbg = 0;
	int i, j, tb;
	int maxfilelength = 512;
	char **nodefile;
	char metafile[100];
	char datasetname[100];
	int idum,numhdf;
	char ctime[6];

	hsize_t count[3], dims[3];
	hsize_t offset_in[3],offset_out[3];
	hid_t file_id,dataset_id,dataspace_id,memoryspace_id;
	int status;
	float time;
	int numi, numj;
	int ix, iy, iz;
	int ihdf;

	int gnx, gny, gnz;
	int fx0, fxf, fy0, fyf;	/* file indices */
	int sx0, sxf, sy0, syf;
	int ax0, axf, ay0, ayf;
	int ax, ay, az;
	int snx, sny, snz, dxleft, dxright, dybot, dytop, ixnode, iynode, k, nxnode, nynode;
	short int *bufi, *bufi0;
	float *buff, *buff0;
	float val;
	int rank;
	int ir=0,jr;
 // filetimes is the actual times stored in the hdf5 file dirtimes
 // contains the list of times contained in the directory names, which
 // corresponds to the first time in the hdf files.
	int *filetimes;
	int time_is_in_file;

 // Some silly output - letter is proportional to fraction of 2D
 // horizontal extent of file read; z means you read the whole
 // horizontal space in the file, a means you read a tiny subset.

	char *alph="abcdefghijklmnopqrstuvwxyz";
	int letter;
	int retval;

	char q_varname[100];

	typedef struct hdfstruct
	{
		int x0, xf, y0, yf;
		int sx0, sxf, sy0, syf;
		int ax0, axf, ay0, ayf;
		int myi, myj;
	} HDFstruct;

	HDFstruct **hdf;

	gnx = gxf - gx0 + 1;
	gny = gyf - gy0 + 1;
	gnz = gzf - gz0 + 1;

	numhdf = nodex * nodey;

	if ((hdf = (HDFstruct **) malloc (numhdf * sizeof (HDFstruct *))) == NULL) ERROR_STOP("Insufficient memory");
	for (i = 0; i < numhdf; i++)
	{
		if ((hdf[i] = (HDFstruct *) malloc (sizeof (HDFstruct))) == NULL) ERROR_STOP("Insufficient memory");
	}

	if ((nodefile = (char **) malloc (numhdf * sizeof (char *))) == NULL)  ERROR_STOP("Insufficient memory");


	// ORF 7/14/11 get last available time. We can only do this by
	// reading the times array from one of the files in the directory
	// with the largest time in its name
	
	// ORF 11/17/12 Wait, we already have 'alltimes' array with all
	// available times! Just use that!

	if (itime < alltimes[0] || itime > alltimes[ntottimes-1])
	{
		fprintf(stderr,"Out of range: %i must be within the range %i to %i\n",itime,alltimes[0],alltimes[ntottimes-1]);
		ERROR_STOP("Requested time not within range\n");
	}

	for (i=0; i < ntimedirs-1; i++)
	{
		if (is_between_l(dirtimes[i],dirtimes[i+1],itime)) break;
	}
	tb = i;
	if (ntimedirs == 1) tb = 0;
// ORF 3/7/13 Hmm, nodefile contains an array with EVERY SINGLE cm1hdf file in
// it - do we really need this?
	for (i = 0; i < numhdf; i++)
	{
		int fooi;
		if ((nodefile[i] = (char *) malloc (maxfilelength*sizeof(char))) == NULL) ERROR_STOP("Insufficient memory");
		sprintf (nodefile[i], "%s/%s/%06i/%s_%06i.cm1hdf5", topdir,timedir[tb],(dn!=-1)?((i/dn)*dn):0,timedir[tb],i);
	}

// ORF 7/12/11 Check if requested time actually exists in the files
// using the times array
//
// ORF 11/17/12 Forget that, we already have alltimes array. Just sweep
// through and check for a match. A binary search algorithm would be
// best, but we are lazy and the array is short

// Note: This check is unnecessary with VisIt reading from the .cm1visit
// file as it's already been read in. We really should just generated
// .cm1visit files on the fly and have hdf2nc and hdf2v5d just go from
// those, that way all checks are done only once.

	time_is_in_file = FALSE;
	for (i=0; i<ntottimes; i++)
	{
		if (itime == alltimes[i])
		{
			time_is_in_file = TRUE;
			break;
		}
	}

	if (time_is_in_file == FALSE)
	{
		fprintf(stderr,"Requested time %i was not saved\n",itime);
		fprintf(stderr,"Available times follow:");
		for (i=0; i<ntottimes; i++)
		{
			fprintf(stderr," %i",alltimes[i]);
		}
		fprintf(stderr,"\n");
	}
	if (time_is_in_file == FALSE)	ERROR_STOP("Invalid time requested");


	/* We build our decomposition from metadata stored in each hdf
	 * file. This replaces the metafile creation etc. Uses same
	 * coordinates as were set in param.F with myi and myj etc.
	 * EXCEPT I start at zero where as the model starts at one (C vs
	 * FORTRAN) */

	numi = nx / nodex;
	numj = ny / nodey;
	for (ihdf = 0; ihdf < numhdf; ihdf++) /* just i not ihdf */
	{
		hdf[ihdf]->myj = ihdf / nodex;
		hdf[ihdf]->myi = ihdf % nodex;
		hdf[ihdf]->x0 = hdf[ihdf]->myi * numi;
		hdf[ihdf]->xf = (hdf[ihdf]->myi + 1) * numi - 1;
		hdf[ihdf]->y0 = hdf[ihdf]->myj * numj; 
		hdf[ihdf]->yf = (hdf[ihdf]->myj + 1) * numj - 1;
		if (dbg)
			fprintf (stderr, "myj = %i myi =%i x0 = %i xf = %i y0 = %i yf = %i\n",
				 hdf[ihdf]->myj, hdf[ihdf]->myi, hdf[ihdf]->x0, hdf[ihdf]->xf, hdf[ihdf]->y0, hdf[ihdf]->yf);
	}

	/* first check if our requested subcube lies within our data */

	if (is_not_between (0, nx - 1, gx0)) ERROR_STOP("Chosen x data out of range");
	if (is_not_between (0, ny - 1, gy0)) ERROR_STOP("Chosen y data out of range");
	if (is_not_between (0, nz - 1, gz0)) ERROR_STOP("Chosen z data out of range");

	if (is_not_between (0, nx - 1, gxf)) ERROR_STOP("Chosen x data out of range");
	if (is_not_between (0, ny - 1, gyf)) ERROR_STOP("Chosen y data out of range");
	if (is_not_between (0, nz - 1, gzf)) ERROR_STOP("Chosen z data out of range");

	for (i = 0; i < numhdf; i++)
	{
		if (is_between (hdf[i]->x0, hdf[i]->xf, gx0) && is_between (hdf[i]->y0, hdf[i]->yf, gy0))
		{
			fx0 = hdf[i]->myi;
			fy0 = hdf[i]->myj;
			if (dbg) fprintf (stderr, "found fx0,fy0 = %i,%i\n", fx0, fy0);
		}
		if (is_between (hdf[i]->x0, hdf[i]->xf, gxf) && is_between (hdf[i]->y0, hdf[i]->yf, gyf))
		{
			fxf = hdf[i]->myi;
			fyf = hdf[i]->myj;
			if (dbg) fprintf (stderr, "found fxf,fyf = %i,%i\n", fxf, fyf);
		}
	}

	snx = numi;
	sny = numj;
	dxleft = gx0 % snx;
	dxright = gxf % snx;
	dybot = gy0 % sny;
	dytop = gyf % sny;
	nxnode = nodex;
	nynode = nodey;

	/*

	   There are four coordinate systems to worry about:

	   1. Node (spans fx0->fxf and fy0->fyf) - indexes the hdf files we read.
	   2. Global (spans 0->nx-1 and 0->ny-1) - the whole model domain.  Calling program provides cartesian coords relative to global.
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

/* Now we have collected all of our array index values delineating node
relative, domain relative, and retrieved-array relative positions. Need
to read data from each node file, assemble retrieved array and send back
as unrolled 1D array which must be reassembled by user */

// ORF 7/7/09
// HDF5 isn't just a data format, it's a way of life.
//
// Save ourselves array assembling issues by creating a memory space with the
// dimensions of the array we are sending back. This is actually pretty cool.

	rank=3;
	count[0]=gnz;count[1]=gny;count[2]=gnx;
	status=memoryspace_id = H5Screate_simple(rank,count,NULL);
	if (status < 0) ERROR_STOP("H5Screate_simple failed");
	printf("\n");
	for (iynode = fy0; iynode <= fyf; iynode++)
	{
		for (ixnode = fx0; ixnode <= fxf; ixnode++)
		{
			k = ixnode + iynode * nxnode;

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
			if (dbg) printf("Varname = %s\n",varname);

			sprintf(datasetname,"%05i/3d/%s",itime,varname);
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
				snz = gzf - gz0 + 1;
				offset_in[0] = gz0;
				offset_in[1] = hdf[k]->sy0;
				offset_in[2] = hdf[k]->sx0;
				offset_out[0] = 0; //orf 6/30/10 gz0 was a bug? should always be 0?
				offset_out[1] = hdf[k]->ay0;
				offset_out[2] = hdf[k]->ax0;
				count[0] = snz;
				count[1] = sny;
				count[2] = snx;
			}
			else if (rank == 2) // I haven't been using this code with 2D data for years...
			{
				snz = 1; //really a 2d array but need snz for malloc
				snx = hdf[k]->sxf - hdf[k]->sx0 + 1;
				sny = hdf[k]->syf - hdf[k]->sy0 + 1;
				offset_in[0] = hdf[k]->sy0;
				offset_in[1] = hdf[k]->sx0;
				offset_out[0] = hdf[k]->ay0;
				offset_out[1] = hdf[k]->ax0;
				count[0] = sny;
				count[1] = snx;
			}
			else 
			{
				ERROR_STOP("Rank must be 2 or 3");
			}
			letter = (int)(25.0*(float)(snx*sny)/(float)(numi*numj));
			printf ("%c",alph[letter]); fflush (stdout);
			if ((H5Sselect_hyperslab (dataspace_id,H5S_SELECT_SET,offset_in,NULL,count,NULL)) < 0)
			{
					fprintf(stderr,"\nCannot select hyperslab dataspace_id %i in file %s\n",dataspace_id,nodefile[k]);
					H5Eprint(H5E_DEFAULT,NULL);
					ERROR_STOP("H5Sselect_hyperslab failed");
			}
			if ((H5Sselect_hyperslab (memoryspace_id,H5S_SELECT_SET,offset_out,NULL,count,NULL)) < 0)
			{
					fprintf(stderr,"\nCannot select hyperslab memoryspace_id %i in file %s\n",memoryspace_id,nodefile[k]);
					H5Eprint(H5E_DEFAULT,NULL);
					ERROR_STOP("H5Sselect_hyperslab failed");
			}
			if ((retval=H5Dread (dataset_id,H5T_NATIVE_FLOAT,memoryspace_id,dataspace_id,H5P_DEFAULT,gf)) < 0)
			{
					fprintf(stderr,"\nCannot read hyperslab for %s in %s\n",datasetname,nodefile[k]);
					fprintf(stderr,"return value = %i\n",retval);
					H5Eprint(H5E_DEFAULT,NULL);
					ERROR_STOP("H5Dread failed");
			}
			
			//ORF TEST
			//printf("%f\n",gf[0]);
			H5Dclose (dataset_id);
			H5Sclose (dataspace_id);
			H5Fclose (file_id);
		}
		printf("\n");
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
}
