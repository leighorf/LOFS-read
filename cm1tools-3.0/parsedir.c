#include <dirent.h>
#include <errno.h>
#include <ctype.h>
#include <limits.h>
#include "lofs-read.h"

#define MAXSTR (512)

extern int regenerate_cache;
DIR *dip;
struct dirent *dit;

void
open_directory(char basedir_full[MAXSTR])
/* This routine must be called before any call to
 * readdir (which iterates over files in the directory) */
{
	if ((dip = opendir (basedir_full)) == NULL)
	{
		perror ("opendir");
		fprintf (stderr, "Directory = %s\n", basedir_full);
		ERROR_STOP ("Can't open directory");
	}
}
void
close_directory()
{
	if (closedir (dip) == -1)
	{
		perror ("closedir");
		ERROR_STOP ("Can't close directory");
	}
}

int
isNumeric (const char *s)
{
	int i;

	if (s == NULL || *s == '\0')
		return 0;

	for (i = 0; i < strlen (s); i++)
		if (!isdigit (s[i]))
			return 0;
	return 1;
}

static int
cmpstringp (const void *p1, const void *p2)
{
	return strcmp (*(char *const *) p1, *(char *const *) p2);
}

static int
cmpintp (const void *p1, const void *p2)
{
	return *(int *) p1 - *(int *) p2;
}

static int
cmpfloatp (const void *p1, const void *p2)
{
	static int retval;
	float fdiff;
	fdiff = *(float *) p1 - *(float *) p2;
	if (fdiff < 0.0) retval= -1;
	if (fdiff > 0.0) retval= 1;
	if (fdiff == 0.0) retval= 0;
	return retval;
}

static int
cmpdoublep (const void *p1, const void *p2)
{
	double fdiff;
	static int retval;
	fdiff = *(double *) p1 - *(double *) p2;
	if (fdiff < 0.0) retval= -1;
	if (fdiff > 0.0) retval= 1;
	if (fdiff == 0.0) retval= 0;
	return retval;
}

void
sortchararray (char **strarray, int nel)
{
	qsort (strarray, nel, sizeof (char *), cmpstringp);
}

void
sortintarray (int *intarray, int nel)
{
	qsort (intarray, nel, sizeof (int), cmpintp);
}

void
sortfloatarray (float *floatarray, int nel)
{
	qsort (floatarray, nel, sizeof (float), cmpfloatp);
}

void
sortdoublearray (double *floatarray, int nel)
{
	qsort (floatarray, nel, sizeof (double), cmpdoublep);
}

void
get_sorted_node_dirs (char *topdir, char *timedir, char **nodedir, int *dn, int nnodedirs,int debug)
{
	int i,j,iret,ns;
	char tmpstr[256]; //size of dit->d_name
	char timedir_full[MAXSTR];

	FILE *fp;

	if (regenerate_cache||(fp = fopen(".cm1hdf5_sorted_node_dirs","r")) == NULL)
	{
		sprintf (timedir_full, "%s/%s", topdir, timedir);

		open_directory (timedir_full);
		j = 0;
		while ((dit = readdir (dip)) != NULL)
		{
			strcpy (tmpstr, dit->d_name);
			ns = strlen (tmpstr);
			if (ns == 7 && isNumeric (tmpstr))
			{
				strcpy (nodedir[j++], tmpstr);
			}
		}

		close_directory();
		sortchararray (nodedir, nnodedirs);
		*dn = (j != 1) ? (atoi (nodedir[1]) - atoi (nodedir[0])) : -1;
		/* What if only one node directory?? In that case, send back -1
		 * and this will tell us to set node directory to 000000 */

//How about regeneratecache		if (!regenerate_cache)
		if ((fp = fopen(".cm1hdf5_sorted_node_dirs","w")) != NULL)
		{
			fprintf(fp,"%i\n",*dn);
			fprintf(fp,"%i\n",nnodedirs);
			for (i=0; i<nnodedirs; i++)
			{
				fprintf(fp,"%s\n",nodedir[i]);
			}
			fclose(fp);
		}
	}
	else
	{
		if (!regenerate_cache)
		if ((fp = fopen(".cm1hdf5_sorted_node_dirs","r")) != NULL)
		{
			if((iret=fscanf(fp,"%i\n",dn))==EOF)ERROR_STOP("fscanf failed");
			if((iret=fscanf(fp,"%i",&nnodedirs))==EOF)ERROR_STOP("fscanf failed");
			for (i=0; i<nnodedirs;i++)
			{
				if((iret=fscanf(fp,"%s",nodedir[i]))==EOF)ERROR_STOP("fscanf failed");
			}
			fclose(fp);
			printf("Read cached nodedir\n");
		}
	}
}

void
get_sorted_time_dirs (char *basedir, char **timedir, double *times, int ntimedirs, char *filebase,int debug)
{
	int i, j, k, ns, iret;
	char tmpstr[256]; // size of dit->d_name
	char firstbase[MAXSTR], base[MAXSTR],dstring[MAXSTR];
	char rhsstr[8];
	char lhsstr[6];

	FILE *fp;

	if (regenerate_cache||(fp = fopen(".cm1hdf5_sorted_time_dirs","r")) == NULL) // First time, haven't created metadata files yet
	{
		printf("Grabbing and caching metadata (only done once):");
		j = 0;
		open_directory (basedir);
		while ((dit = readdir (dip)) != NULL) 
		{
			strcpy (tmpstr, dit->d_name);
			if(debug)printf("%s ",tmpstr);fflush(stdout);
			ns = strlen (tmpstr);
			if (ns > 15)
			{
				for (i = 0; i < 7; i++)
					rhsstr[i] = tmpstr[ns - 7 + i];
				rhsstr[7] = '\0';
				if (isNumeric (rhsstr))
				{
					for (i = 0; i < 5; i++)
						lhsstr[i] = tmpstr[ns - 13 + i]; // 13 = 7+5+1
					lhsstr[5] = '\0';
				}
				else
				{
					ERROR_STOP("Something is wrong with our timedir format");
				}
				if (isNumeric (lhsstr))
				{
					sprintf(dstring,"%s.%s",lhsstr,rhsstr);
//					printf("dstring = %s\n",dstring);
					times[j] = atof(dstring); //ding
					if(debug) printf("j = %i times = %lf\n",j,times[j]);
				}
				else
				{
					ERROR_STOP("Something is wrong with our timedir format");
				}
// Sicne we sometimes use . in the base, we will work backwards and
// choose the stuff *before* the second dot *from the right*
				i = ns - 1;
				while (tmpstr[i--]!='.');
				while (tmpstr[i--]!='.');
				i++;
				if (debug) printf("i = %i, j = %i\n",i,j);
				for (k = 0; k < i; k++)
				{
					if (j == 0) firstbase[k] = tmpstr[k];
					else             base[k] = tmpstr[k];
				}
				firstbase[k] = '\0';
				base[k] = '\0';

				if (j > 0)
				{
					if (debug) printf("base = %s\n",base);
					if (strcmp (firstbase, base) != 0)
					{
						printf ("ACK: we have more than one base!\n");
						printf ("firstbase = %s, base = %s\n", firstbase, base);
						ERROR_STOP ("Only one basename allowed per directory");
					}
				}
				if (j == 0)
				{
					printf ("\nfirstbase = %s\n", firstbase);
					strcpy(filebase,firstbase);
				}

				strcpy (timedir[j], dit->d_name);
				j++;
			}
			else if( !strcmp(tmpstr,".") || !strcmp(tmpstr,".."))
			{
				// do nothing; we have happened upon '.' or '..'
			}
			else ERROR_STOP("Something wrong with file names in timedir");
		}
		close_directory();
//		}
		sortchararray (timedir, ntimedirs);
		sortdoublearray (times, ntimedirs);
//recache		if (!regenerate_cache)
		if ((fp = fopen(".cm1hdf5_sorted_time_dirs","w")) != NULL) 
		{
			fprintf(fp,"%i\n",ntimedirs);
			for (i = 0; i < ntimedirs; i++) fprintf(fp,"%s %lf\n",timedir[i],times[i]);
			fclose(fp);
		}
	}
	else
	{
		if (!regenerate_cache)
		if ((fp = fopen(".cm1hdf5_sorted_time_dirs","r")) != NULL) 
		{
			if((iret=fscanf(fp,"%i\n",&ntimedirs))==EOF)ERROR_STOP("fscanf failed");
			for (i = 0; i < ntimedirs; i++) if((iret=fscanf(fp,"%s %lf",timedir[i],&(times[i])))==EOF)ERROR_STOP("fscanf failed");
			fclose(fp);
			printf("Read cached sorted time dirs\n");
		}
	}
}

int
get_num_time_dirs (char *basedir,int debug)
{
	int i, j, ns, iret;
	char tmpstr[256]; // size of dit->d_name
	char firstbase[MAXSTR], base[MAXSTR];
	char rhsstr[8]; // 7 digits to right hand side of . in time which is SSSSS.FFFFFFF (plus null termination char)
	char lhsstr[6]; // 5 digits to left hand side of . in time (plus null termination char)
	char sd;
	
	FILE *fp;

	if (regenerate_cache||(fp = fopen(".cm1hdf5_num_time_dirs","r")) == NULL) // First time, haven't created metadata files yet
	{
		open_directory (basedir);
		j = 0;
		while ((dit = readdir (dip)) != NULL)
		{
			strcpy (tmpstr, dit->d_name);
			ns = strlen (tmpstr);
			if (debug) printf("ns = %i\n",ns);
			if (ns > 15 )
			{
				for (i = 0; i < 7; i++) rhsstr[i] = tmpstr[ns - 7 + i];
				rhsstr[7] = '\0';
				if(debug) printf("DEBUG: rhsstr = %s\n",rhsstr);
				if (isNumeric (rhsstr))
				{
					for (i = 0; i < 5; i++) lhsstr[i] = tmpstr[ns - 13 + i]; // 13 = 7+5+1
					lhsstr[5] = '\0';

					if (isNumeric (lhsstr))
					{
						i = 0;
						sd = 'a';
						while (sd != '.')
						{
							sd = tmpstr[i];
							if (j == 0)
								firstbase[i] = sd;
							else
								base[i] = sd;
							i++;
							if(debug) printf("DEBUG: firstbase = %s\n",firstbase);
							if(debug) printf("DEBUG: base = %s\n",base);
						}
						if (j == 0)
							firstbase[i - 1] = '\0';
						else
							base[i - 1] = '\0';
						if (j > 0)
						{
							if (strcmp (firstbase, base) != 0)
							{
								printf ("ACK: we have more than one base!\n");
								printf ("firstbase = %s, base = %s\n", firstbase, base);
								ERROR_STOP ("Only one basename allowed per directory");
							}
						}
					}
				}
				else
				{
					ERROR_STOP("Something wrong with file names in timedir");
				}
				j++;
			}
			else if( !strncmp(tmpstr,".",1) || !strncmp(tmpstr,"..",2))
			{
				// do nothing; we have happened upon '.' or '..'
			}
			else ERROR_STOP("Something wrong with file names in timedir");
		}
		close_directory();
//recache		if (!regenerate_cache)
		if ((fp = fopen(".cm1hdf5_num_time_dirs","w")) != NULL)
		{
			fprintf(fp,"%i\n",j);
			fclose(fp);
		}
		else
		{
			printf("CANNOT WRITE .num_time_dirs\n");
		}
	}
	else
	{
		if (!regenerate_cache)
		if ((fp = fopen(".cm1hdf5_num_time_dirs","r")) != NULL)
		{
			if((iret=fscanf(fp,"%i",&j))==EOF)ERROR_STOP("fscanf failed");
			fclose(fp);
			printf("Read cached num_time_dirs of %i\n",j);
		}
	}

	return j;
}

int
get_num_node_dirs (char *topdir, char *timedir,int debug)
{
	int j, ns,iret;
	char timedir_full[MAXSTR];
	char tmpstr[256]; // size of dit->d_name

	FILE *fp;

	if (regenerate_cache||(fp = fopen(".cm1hdf5_num_node_dirs","r")) == NULL)
	{
		sprintf (timedir_full, "%s/%s", topdir, timedir);
		open_directory (timedir_full);
		j = 0;
		while ((dit = readdir (dip)) != NULL)
		{
			strcpy (tmpstr, dit->d_name);
			ns = strlen (tmpstr);
			if (ns == 7 && isNumeric (tmpstr))
				j++;
		}
		close_directory();
//recache		if (!regenerate_cache)
		if ((fp = fopen(".cm1hdf5_num_node_dirs","w")) != NULL)
		{
			fprintf(fp,"%i\n",j);
			fclose(fp);
		}
	}
	else
	{
		if (!regenerate_cache)
		if ((fp = fopen(".cm1hdf5_num_node_dirs","r")) != NULL)
		{
			if((iret=fscanf(fp,"%i",&j))==EOF)ERROR_STOP("fscanf failed");
			fclose(fp);
			printf("Read cached num node dirs\n");
		}
	}

	return j;
}

int get_nodemask(char basedir_full[MAXSTR])
{
	char tmpstr[256]; //size of dit->d_name
	int j,mask;
	char *foochar="a";//set to anything

	open_directory (basedir_full);
	while ((dit = readdir (dip)) != NULL)
	{
		strcpy (tmpstr, dit->d_name);
		foochar = strstr(tmpstr,".cm1hdf5");
		mask = (foochar == NULL) ? 0 : 1;
		if (foochar != NULL) break; /* Once we find a single cm1hdf5 file we are done */
	}
	close_directory();
	return mask;
}
int get_nfiles()
{
	char tmpstr[256]; //size of dit->d_name
	char *foochar="a";//set to anything
	int nfiles;
	nfiles=0;
	while ((dit = readdir (dip)) != NULL)
	{
		strcpy (tmpstr, dit->d_name);
//		printf("tmpstr = %s\n",tmpstr);
		foochar = strstr(tmpstr,".cm1hdf5");
		if (foochar != NULL) nfiles++;
	}
	close_directory();
	return nfiles;
}
void get_unsorted_file_list(char** cm1hdf5file)
{
	char tmpstr[256]; //size of dit->d_name
	char *foochar="a";//set to anything
	int i;
	i=0;
	while ((dit = readdir (dip)) != NULL) /* Populate file list */
	{
		strcpy (tmpstr, dit->d_name);
		foochar = strstr(tmpstr,".cm1hdf5");
		if (foochar != NULL) {strcpy (cm1hdf5file[i],tmpstr);i++;}
	}
	close_directory();
}

double *
get_all_available_times (char *topdir, char **timedir, int ntimedirs, char **nodedir, int nnodedirs, int *ntottimes, char *firstfilename, int *firsttimedirindex,
		int *saved_X0, int *saved_Y0, int *saved_X1, int *saved_Y1, int debug)
{
	int i, j, k, iret, itime;
	char basedir_full[MAXSTR], tmpstr[256]; // size of dit->d_name
	hid_t file_id,dataset_id,dataspace_id;
	hsize_t dims[1],maxdims[1];
	double *filetimes;
	double *alltimes;
	char *foochar="a";//set to anything
	int *nodedirmask;
	int firstnodedir,lastnodedir,nfiles;
	char **cm1hdf5file;
	char *hdf5filename;

	FILE *fp;

	hdf5filename = (char *)malloc(512*sizeof(char));
	alltimes = (double *)malloc(sizeof(double)); //to keep compiler from complaining
	nodedirmask = (int *)malloc(nnodedirs*sizeof(int));
	if (regenerate_cache||(fp = fopen(".cm1hdf5_all_available_times","r")) == NULL)
	{
		*ntottimes = 0;
		k = 0;

/*

3/14/18

New: This routine, called once per data set where data is cached
to ASCII text files, also now retrieves actual saved domain index
(x0,y0,x1,y1) parameters. We cannot assume it's the full domain as I
regularly save subdomains focused on the mesocyclone, and it's possible
that files are culled to save disk space. Without passing values to
X0 etc. these saved bounds are the default. Hence it's a quick way to
construct, for instance, a netCDF file with all the saved LOFS data at
a given time. No longer will I need to hunt in the node directories for
the files containing these values!

This required a lot more code than I thought it would. I have broken out
several subroutines to make the code more clear.

Usually a subdomain is saved in my big CM1 simulations. The actual
indices can only be retrieved by diving into the LOFS directory
structure. This new functionality will mean that if hdf2nc (or
makevisit) is called with no horizontal indices, it will choose the
saved ones rather than throwing an error when it can't find file
00000...

First we evaluate the which node directories contain actual cm1hdf5
files (there can be empty ones), creating a 1d integer mask array.
Then, we go to the 1st of these node directories to get x0,y0 from the
"smallest numbered file" and then we go to the last to get x1 y1 in the
"largest numbered file" since cm1 does a 2D decomposition which can be
expressed in a simple 1d array.

TODO: Save Z0 in cm1hdf5 files so we can retrieve that as well.

*/

		firstnodedir=lastnodedir=*saved_X0=*saved_Y0=*saved_X1=*saved_Y1=0;

		itime = 0; /* Only need one time */
		{
			int found_something=0;
			for (j=0; j < nnodedirs; j++)
			{
				sprintf (basedir_full, "%s/%s/%s", topdir, timedir[itime], nodedir[j]);
				if(debug)printf("get_all_available_times: topdir = %s\t timedir[%i] = %s\t nodedir[%i] = %s\n",topdir,itime,timedir[itime],j,nodedir[j]);

				nodedirmask[j] = get_nodemask(basedir_full);
//				printf("nodedirmask[%i]=%i\n",j,nodedirmask[j]);
				if (nodedirmask[j] != 0) found_something = 1;
			}
//Edge case: No friggin hdf5 files anywhere! THANK YOU SCRUBBER PROCESS YOU HEARTLESS RO-BOT
			if (found_something == 0) ERROR_STOP("We don't have any hdf5 files at all! What is this nonsense??");
		}

/*

Now nodedirmask contains zero or more zeroes, followed by at least
one one, followed by zero or more zeroes. Sweep through and get the
first and last node dirs that contain the metadata we need. Plants
crave electrolytes.

*/

		/* Tested good with big-assed 15 meter data where we only saved
		 * the center of the domain*/
		firstnodedir=0;
		lastnodedir=nnodedirs-1;
		for (j=0; j < nnodedirs-1; j++)
		{
			if((nodedirmask[j+1]-nodedirmask[j]) ==  1) firstnodedir=j+1;
			if((nodedirmask[j+1]-nodedirmask[j]) == -1) lastnodedir=j;
		}
		printf("firstnodedir=%i lastnodedir=%i\n",firstnodedir,lastnodedir);

		j=firstnodedir;
		{
			sprintf (basedir_full, "%s/%s/%s", topdir, timedir[itime], nodedir[j]);
			open_directory(basedir_full);
			nfiles = get_nfiles();
			fprintf(stderr,"Number of cm1hdf5 files in our first (%i) node directory: %i\n",j,nfiles);
			/* Allocate file name array */
			cm1hdf5file = (char **)malloc(nfiles * sizeof(char *));
			for (i=0; i < nfiles; i++) cm1hdf5file[i] = (char *)(malloc(MAXSTR * sizeof(char)));

			open_directory(basedir_full);
			get_unsorted_file_list(cm1hdf5file);
			
			/* Sort 'em */

			sortchararray (cm1hdf5file, nfiles);

			/* Now, finally, pluck our data! */
			sprintf(hdf5filename,"%s/%s",basedir_full,cm1hdf5file[0]); /* <-- index zero is earliest time */
			if ((file_id = H5Fopen (hdf5filename, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
			{
				fprintf (stderr, "Cannot open %s, even though it should exist!\n", hdf5filename);
				ERROR_STOP ("Cannot open hdf file");
			}
			get0dint(file_id,"grid/x0",saved_X0);
			get0dint(file_id,"grid/y0",saved_Y0);
			H5Fclose(file_id);
			fprintf(stderr,"Setting X0 to saved_X0 which is %i\n",*saved_X0);
			fprintf(stderr,"Setting Y0 to saved_Y0 which is %i\n",*saved_Y0);
		}
		j=lastnodedir;
		{
			sprintf (basedir_full, "%s/%s/%s", topdir, timedir[itime], nodedir[j]);
			open_directory(basedir_full);
			nfiles = get_nfiles();
			fprintf(stderr,"Number of cm1hdf5 files in our last (%i) node directory: %i\n",j,nfiles);
			/* Allocate file name array */
			cm1hdf5file = (char **)malloc(nfiles * sizeof(char *));
			for (i=0; i < nfiles; i++) cm1hdf5file[i] = (char *)(malloc(MAXSTR * sizeof(char)));

			open_directory (basedir_full);
			get_unsorted_file_list(cm1hdf5file);
			
			/* Sort 'em */
			sortchararray (cm1hdf5file, nfiles);

			/* Now, finally, pluck our data! */

			sprintf(hdf5filename,"%s/%s",basedir_full,cm1hdf5file[nfiles-1]);/* <-- index nfiles-1 is latest time */
			if ((file_id = H5Fopen (hdf5filename, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
			{
				fprintf (stderr, "Cannot open %s, even though it should exist!\n", firstfilename);
				ERROR_STOP ("Cannot open hdf file");
			}

			get0dint(file_id,"grid/x1",saved_X1);
			get0dint(file_id,"grid/y1",saved_Y1);
			H5Fclose(file_id);
			fprintf(stderr,"Setting X1 to saved_X1 which is %i\n",*saved_X1);
			fprintf(stderr,"Setting Y1 to saved_Y1 which is %i\n",*saved_Y1);
		}

// Thus endeth the code that finds the saved domain bounds in X and Y. Z
// shouldn't be too hard... except for the fact that unless I save Z0
// as metadata, I have to go find the dimensions of one of the saved
// 3D files... I think I'll add a new piece of metadata to the grid
// group...

		for (i = 0; i < ntimedirs; i++)
		{
			for (j=0; j < nnodedirs; j++)
				//WHY ARE WE DOING THIS? Because now that I am doing big assed runs I have situations
				//where when saving subdomains, some of the node
				//directories are completely empty, so I can no longer
				//just choose 'nodedir[0]' to find my first hdf5 file.
			{
				sprintf (basedir_full, "%s/%s/%s", topdir, timedir[i], nodedir[j]);
				if(debug)printf("get_all_available_times: topdir = %s\t timedir[%i] = %s\t nodedir[%i] = %s\n",topdir,i,timedir[i],j,nodedir[j]);

				open_directory(basedir_full);
				while ((dit = readdir (dip)) != NULL)
				{
					strcpy (tmpstr, dit->d_name);
					foochar = strstr(tmpstr,".cm1hdf5");
					if(debug)printf("get_all_available_times: %s\n",basedir_full);
					if (foochar != NULL) break;	// Got a cm1hdf5 file
				}	
				close_directory();
				if (foochar != NULL) break;	// Got a cm1hdf5 file
			}
			*firsttimedirindex = j; //Might use this someday

			if (!foochar)
			{
				/* This should never happen if we have gotten
				 * this far! */
				printf("Argh: something wrong with file names in %s\n",basedir_full);
				ERROR_STOP("Files are messed up\n");
			}

/*

   firstfilename is the first cm1hdf5 file retrieved from a directory
   that contains model data. All cm1hdf5 files contain metadata for the
   full simulation, and also contains metadata for their own time and/or
   location. Here, we grab the first cm1hdf5 file in a given time/node
   directory and extract all the times the file contains, which is
   the same for all of the other files in the directory. We keep
   firstfilename for use later for extracting global metadata such as
   nx, ny, nz etc... for these any cm1hdf5 file will do! In this case,
   firstfilename will be from a file in the last time directory since
   here we are sweeping through all times to construct the all_times
   array.

*/

			sprintf(firstfilename,"%s/%s",basedir_full,tmpstr);
			if(debug)
			{
				printf("get_all_available_times: firstfilename = %s\n",firstfilename);
			}
			else
			{
				printf(".");fflush(stdout);
			}

			if ((file_id = H5Fopen (firstfilename, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
			{
				fprintf (stderr, "Cannot open firstfilename which is %s, even though we have already opened it!\n", firstfilename);
				ERROR_STOP ("Cannot open hdf file");
			}

// ORF 8/12/11 This code was taken from readmult. We have a lot of
// redundant I/O that could be reduced. Should probably have the option
// of making a master file contaning metadata for VisIt. It could be
// the cm1_3D file, perhaps.

// ORF 2/5/13 Concerning above: We did that. This is only called when
// doing the 1 shot grokking
//
// ORF 6/1/16 Actually now I cache everything, so only need to do
// all this shit once per hdf2nc or whatever. This routine is 
// expensive if you have a lot of files.

			if ((dataset_id = H5Dopen(file_id,"times",H5P_DEFAULT)) < 0) ERROR_STOP("Cannot H5Dopen");
			if ((dataspace_id = H5Dget_space(dataset_id)) < 0) ERROR_STOP("Cannot H5Dget_space");
			if ((H5Sget_simple_extent_dims(dataspace_id,dims,maxdims)) < 0) ERROR_STOP("Cannot H5Sget_simple_extent_dims"); //dims[0] will equal number of time levels in file
			if ((filetimes = (double *)malloc(dims[0]*sizeof(double)))==NULL) ERROR_STOP("Cannot malloc filetimes array");
			get1ddouble(file_id,"times",filetimes,0,dims[0]);
			
			(*ntottimes) += dims[0];

			alltimes = (k==0)?(double *)malloc(dims[0]*sizeof(double)):(double *)realloc(alltimes,(*ntottimes)*sizeof(double));

			for (j=0; j<dims[0];j++)alltimes[j+k] = filetimes[j];
			k = *ntottimes;

			free (filetimes);
			H5Fclose(file_id);

		}
		printf("\n");
//recache		if (!regenerate_cache)
		if ((fp = fopen(".cm1hdf5_all_available_times","w")) != NULL)
		{
			fprintf(fp,"%s\n",firstfilename);
			fprintf(fp,"%i %i %i %i\n",*saved_X0,*saved_Y0,*saved_X1,*saved_Y1);
			fprintf(fp,"%i\n",*ntottimes);
			for (i=0; i<*ntottimes; i++)
			{
				fprintf(fp,"%lf\n",alltimes[i]);
			}
			fclose(fp);
		}
	}
	else
	{
		if (!regenerate_cache)
		if ((fp = fopen(".cm1hdf5_all_available_times","r")) != NULL)
		{
			iret=fscanf(fp,"%s",firstfilename);
			if(iret!=EOF) fprintf(stderr,"Cached: firstfilename = %s\n",firstfilename);
				else ERROR_STOP("fscanf firstfilename failed");
			iret=fscanf(fp,"%i %i %i %i",saved_X0,saved_Y0,saved_X1,saved_Y1);
			if(iret!=EOF){
				fprintf(stderr,"Cached: saved_X0  = %6i\n",*saved_X0);
				fprintf(stderr,"Cached: saved_Y0  = %6i\n",*saved_Y0);
				fprintf(stderr,"Cached: saved_X1  = %6i\n",*saved_X1);
				fprintf(stderr,"Cached: saved_Y1  = %6i\n",*saved_Y1);
			}
				else ERROR_STOP("fscanf saved_[XY][01] failed");
			iret=fscanf(fp,"%i",ntottimes);
			if(iret!=EOF) fprintf(stderr,"Cached: ntottimes = %6i\n",*ntottimes);
			alltimes = (double *)malloc(*ntottimes * sizeof(double));
			for (i=0; i<*ntottimes; i++)
			{
				if((iret=fscanf(fp,"%lf",&(alltimes[i])))==EOF)ERROR_STOP("fscanf alltimes failed");

			}
			printf("Read all cached metadata from dot files\n");
		}
	}

	return(alltimes);
}
