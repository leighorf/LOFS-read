#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <limits.h>
#include "hdforf.h"
#include "errorf.h"

#define MAXSTR (512)

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
get_sorted_node_dirs (char *topdir, char *timedir, char **nodedir, int *dn, int nnodedirs)
{
	DIR *dip;
	struct dirent *dit;
	int i,j,iret,ns;
	char tmpstr[256]; //size of dit->d_name
	char timedir_full[MAXSTR];

	FILE *fp;

	if ((fp = fopen(".cm1hdf5_sorted_node_dirs","r")) == NULL)
	{
		sprintf (timedir_full, "%s/%s", topdir, timedir);

		if ((dip = opendir (timedir_full)) == NULL)
		{
			perror ("opendir");
			fprintf (stderr, "Directory = %s\n", timedir_full);
			ERROR_STOP ("Can't open directory");
		}

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

		if (closedir (dip) == -1)
		{
			perror ("closedir");
			ERROR_STOP ("Can't close directory");
		}
		sortchararray (nodedir, nnodedirs);
		*dn = (j != 1) ? (atoi (nodedir[1]) - atoi (nodedir[0])) : -1;
		/* What if only one node directory?? In that case, send back -1
		 * and this will tell us to set node directory to 000000 */

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
get_sorted_time_dirs (char *basedir, char **timedir, double *times, int ntimedirs, char *filebase)
{
	DIR *dip;
	struct dirent *dit;
	int i, j, k, ns, iret;
	char tmpstr[256]; // size of dit->d_name
	char firstbase[MAXSTR], base[MAXSTR],dstring[MAXSTR];
	char rhsstr[8];
	char lhsstr[6];

	FILE *fp;

	if ((fp = fopen(".cm1hdf5_sorted_time_dirs","r")) == NULL) // First time, haven't created metadata files yet
	{
		j = 0;
		if ((dip = opendir (basedir)) == NULL)
		{
			perror ("opendir");
			fprintf (stderr, "Directory = %s\n", basedir);
			ERROR_STOP ("Can't open directory");
		}
		while ((dit = readdir (dip)) != NULL) //ORF the problem was times was passed as int * so we are past the "hang" problem at least
		{
			printf("doing strcpy to tmpstr...");
			strcpy (tmpstr, dit->d_name);
			printf("...done\n");
			printf("Working on %s\n",tmpstr);
			ns = strlen (tmpstr);
			printf("ns = %i\n",ns);
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
					printf("dstring = %s\n",dstring);
					times[j] = atof(dstring); //ding
					printf("j = %i times = %lf\n",j,times[j]);
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
				printf("i = %i, j = %i\n",i,j);
				for (k = 0; k < i; k++)
				{
					if (j == 0) firstbase[k] = tmpstr[k];
					else             base[k] = tmpstr[k];
				}
				firstbase[k] = '\0';
				base[k] = '\0';

				if (j > 0)
				{
					printf("base = %s\n",base);
					if (strcmp (firstbase, base) != 0)
					{
						printf ("ACK: we have more than one base!\n");
						printf ("firstbase = %s, base = %s\n", firstbase, base);
						ERROR_STOP ("Only one basename allowed per directory");
					}
				}
				if (j == 0)
				{
					printf ("firstbase = %s\n", firstbase);
					strcpy(filebase,firstbase);
				}

				printf("Doing strcpy...");
				strcpy (timedir[j], dit->d_name);
				printf("...done\n");
				j++;
			}
			else if( !strcmp(tmpstr,".") || !strcmp(tmpstr,".."))
			{
				// do nothing; we have happened upon '.' or '..'
			}
			else ERROR_STOP("Something wrong with file names in timedir");
			printf("Made it to end of while dit loop\n");
		}
		if (closedir (dip) == -1)
		{
			perror ("closedir");
			ERROR_STOP ("Can't close directory");
		}
		sortchararray (timedir, ntimedirs);
		sortdoublearray (times, ntimedirs);
		if ((fp = fopen(".cm1hdf5_sorted_time_dirs","w")) != NULL) 
		{
			fprintf(fp,"%i\n",ntimedirs);
			for (i = 0; i < ntimedirs; i++) fprintf(fp,"%s %lf\n",timedir[i],times[i]);
			fclose(fp);
		}
	}
	else
	{
		if ((fp = fopen(".cm1hdf5_sorted_time_dirs","r")) != NULL) 
		{
			if((iret=fscanf(fp,"%i\n",&ntimedirs))==EOF)ERROR_STOP("fscanf failed");
			for (i = 0; i < ntimedirs; i++) if((iret=fscanf(fp,"%s %lf",timedir[i],&(times[i])))==EOF)ERROR_STOP("fscanf failed");
			fclose(fp);
			printf("Read cached sorted time dirs\n");
//			for (i = 0; i < ntimedirs; i++) printf("FART: %s %lf\n",timedir[i],times[i]);
		}
	}
}

int
get_num_time_dirs (char *basedir)
{
	DIR *dip;
	struct dirent *dit;
	int i, j, ns, iret;
	char tmpstr[256]; // size of dit->d_name
	char firstbase[MAXSTR], base[MAXSTR];
	char rhsstr[8]; // 7 digits to right hand side of . in time which is SSSSS.FFFFFFF (plus null termination char)
	char lhsstr[6]; // 5 digits to left hand side of . in time (plus null termination char)
	char sd;
	
	FILE *fp;

	if ((fp = fopen(".cm1hdf5_num_time_dirs","r")) == NULL) // First time, haven't created metadata files yet
	{
		printf("Caching num_time_dirs...\n");
		if ((dip = opendir (basedir)) == NULL)
		{
			perror ("opendir");
			fprintf (stderr, "Directory = %s\n", basedir);
			ERROR_STOP ("Can't open directory");
		}
		j = 0;
		while ((dit = readdir (dip)) != NULL)
		{
			strcpy (tmpstr, dit->d_name);
	//printf("ORF: debug: tmpstr = %s\n",tmpstr);
			ns = strlen (tmpstr);
			printf("ns = %i\n",ns);
			if (ns > 15 )
			{
				for (i = 0; i < 7; i++) rhsstr[i] = tmpstr[ns - 7 + i];
				rhsstr[7] = '\0';
	//printf("ORF: debug: rhsstr = %s\n",rhsstr);
				if (isNumeric (rhsstr))
				{
					for (i = 0; i < 5; i++) lhsstr[i] = tmpstr[ns - 13 + i]; // 13 = 7+5+1
					lhsstr[5] = '\0';

					if (isNumeric (lhsstr))
					{
//						itime = strtol (rhsstr, NULL, 10); if (itime==LONG_MIN||itime==LONG_MAX) ERROR_STOP("strtol failed"); //does nothing
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
		//printf("ORF: debug: firstbase = %s\n",firstbase);
		//printf("ORF: debug: base = %s\n",base);
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
		if (closedir (dip) == -1)
		{
			perror ("closedir");
			ERROR_STOP ("Can't close directory");
		}
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
get_num_node_dirs (char *topdir, char *timedir)
{
	DIR *dip;
	struct dirent *dit;
	int j, ns,iret;
	char timedir_full[MAXSTR];
	char tmpstr[256]; // size of dit->d_name

	FILE *fp;

	if ((fp = fopen(".cm1hdf5_num_node_dirs","r")) == NULL)
	{
		sprintf (timedir_full, "%s/%s", topdir, timedir);
		if ((dip = opendir (timedir_full)) == NULL)
		{
			perror ("opendir");
	//		printf ("ORF DEBUG: timedir_full = %s\n", timedir_full);
			ERROR_STOP ("Can't open directory");
			return 0;
		}

		j = 0;
		while ((dit = readdir (dip)) != NULL)
		{
			strcpy (tmpstr, dit->d_name);
			ns = strlen (tmpstr);
			if (ns == 7 && isNumeric (tmpstr))
				j++;
		}
		if (closedir (dip) == -1)
		{
			perror ("closedir");
			ERROR_STOP ("Can't close directory");
		}
		if ((fp = fopen(".cm1hdf5_num_node_dirs","w")) != NULL)
		{
			fprintf(fp,"%i\n",j);
			fclose(fp);
		}
	}
	else
	{
		if ((fp = fopen(".cm1hdf5_num_node_dirs","r")) != NULL)
		{
			if((iret=fscanf(fp,"%i",&j))==EOF)ERROR_STOP("fscanf failed");
			fclose(fp);
			printf("Read cached num node dirs\n");
		}
	}

	return j;
}

//ORF WE NO LONGER CALL THIS ROUTINE
//We were basically calling it several times in the get_all_times
//routine. So now it is sent back from get_all_available_times below.
//Note that "firstfilename" will be in the *LAST* available time
//directory. Who cares, right? But we should probably clean this shit
//up. 2015-10-16
void
get_first_hdf_file_name (char *topdir, char *timedir, char *nodedir, char *filename)
{
	DIR *dip;
	struct dirent *dit;
	char basedir_full[MAXSTR], tmpstr[256]; // size of dit->d_name

	sprintf (basedir_full, "%s/%s/%s", topdir, timedir, nodedir);

	if ((dip = opendir (basedir_full)) == NULL)
	{
		perror ("opendir");
		fprintf (stderr, "Directory = %s\n", basedir_full);
		ERROR_STOP ("Can't open directory");
	}
//	dit = readdir (dip);
	while ((dit = readdir (dip)) != NULL)
	{
		strcpy (tmpstr, dit->d_name);
//	      printf("get_first_hdf_filename: tmpstr = %s\n",tmpstr);
		if (strstr(tmpstr,".cm1hdf5")!=NULL) break;
	}
	if (closedir (dip) == -1)
	{
		perror ("closedir");
		ERROR_STOP ("Can't close directory");
	}
	sprintf (filename, "%s/%s", basedir_full, tmpstr);
}

double *
get_all_available_times (char *topdir, char **timedir, int ntimedirs, char **nodedir, int nnodedirs, int *ntottimes, char *firstfilename, int *firsttimedirindex)
{
	DIR *dip;
	struct dirent *dit;
	int i, j, k,iret;
	char basedir_full[MAXSTR], tmpstr[256]; // size of dit->d_name
	hid_t file_id,dataset_id,dataspace_id;
	hsize_t dims[1],maxdims[1];
	double *filetimes;
	double *alltimes;
	char *foochar="a";//set to anything

	FILE *fp;

//	firstfilename = (char *)malloc(512*sizeof(char));
//
//

	alltimes = (double *)malloc(sizeof(double)); //to keep compiler from complaining
	if ((fp = fopen(".cm1hdf5_all_available_times","r")) == NULL)
	{
		*ntottimes = 0;
		k = 0;
		for (i = 0; i < ntimedirs; i++)
		{
			for (j=0; j < nnodedirs; j++) //WHY ARE WE DOING THIS? Because now that I am doing big assed runs I have situations
				//where when saving subdomains, some of the node
				//directories are completely empty, so I can no longer
				//just choose 'nodedir[0]' to find my first hdf5 file.
			{
				sprintf (basedir_full, "%s/%s/%s", topdir, timedir[i], nodedir[j]);
				printf("get_all_available_times: topdir = %s\t timedir[%i] = %s\t nodedir[%i] = %s\n",topdir,i,timedir[i],j,nodedir[j]);

				if ((dip = opendir (basedir_full)) == NULL)
				{
					perror ("opendir");
					fprintf (stderr, "Directory = %s\n", basedir_full);
					ERROR_STOP ("Can't open directory");
				}
				while ((dit = readdir (dip)) != NULL)
				{
					strcpy (tmpstr, dit->d_name);
					foochar = strstr(tmpstr,".cm1hdf5");
					printf("get_all_available_times: %s\n",basedir_full);
					if(foochar != NULL) printf("get_all_available_times: foochar = %c\n", *foochar); else printf("NULL\n");
					printf("get_all_available_times: tmpstr = %s\n",tmpstr);
					if (foochar != NULL) break;	// Got one
				}	
				if (foochar != NULL) break;	// Got one
			}
			*firsttimedirindex = j; //Might use this someday

			if (!foochar)
			{
				printf("Argh: something wrong with file names in %s\n",basedir_full);
				ERROR_STOP("Files are messed up\n");
			}

			sprintf(firstfilename,"%s/%s",basedir_full,tmpstr);
			printf("get_all_available_times: firstfilename = %s\n",firstfilename);

			if ((file_id = H5Fopen (firstfilename, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
			{
				fprintf (stderr, "Cannot open firstfilename which is %s, even though we have already opened it!\n", firstfilename);
				ERROR_STOP ("Cannot open hdf file");
			}

	// ORF 8/12/11 This code was taken from readmult. We have a lot of
	// redundant I/O that could be reduced. Should probably have the option
	// of making a master file contaning metadata for VisIt. It could be the
	// cm1_3D file, perhaps.

	// ORF 2/5/13 Concerning above: We did that. This is only called when
	// doing the 1 shot grokking
	//
	// ORF 6/1/16 Actually now I cache everything, so only need to do all this shit once per hdf2nc or whatever. This routine is very expensive.

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

			if (closedir (dip) == -1)
			{
				perror ("closedir");
				ERROR_STOP ("Can't close directory");
			}
		}
		if ((fp = fopen(".cm1hdf5_all_available_times","w")) != NULL)
		{
			fprintf(fp,"%s\n",firstfilename);
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
		if ((fp = fopen(".cm1hdf5_all_available_times","r")) != NULL)
		{
			if((iret=fscanf(fp,"%s",firstfilename))==EOF)ERROR_STOP("fscanf failed");
			if((iret=fscanf(fp,"%i",ntottimes))==EOF)ERROR_STOP("fscanf failed");
			alltimes = (double *)malloc(*ntottimes * sizeof(double));
			for (i=0; i<*ntottimes; i++)
			{
				if((iret=fscanf(fp,"%lf",&(alltimes[i])))==EOF)ERROR_STOP("fscanf failed");
				printf("Reading from cached file: alltimes = %f\n",alltimes[i]);
			}
			printf("Read cached firstfilename and all times\n");
		}
	}

//	free (firstfilename);
	return(alltimes);
}
