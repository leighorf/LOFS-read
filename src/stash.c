/*
 *
 * This code was written by Leigh Orf <leigh.orf@.wisc.edu> 24-26 October 2016
 *
 * Its purpose is to stash useful data in the "user block" of hdf5
 * files. It does not create the user block; this must be done in the
 * code that wrote the file with h5pset_userblock.

 * This code inserts namelist.input data, sounding data, and the first
 * few lines of the output file, and some other things, including an
 * optional stash block. The code first gets the userblock size from
 * the hdf5 file and only writes data up to this size so as to avoid
 * clobbering the HDF5 data, which would suck.

 * I've decided to create a virtual file for the user block data using
 * fmemopen which allows us to treat a block of memory as a stream
 * (file) that can be rewound etc. This means we will not be constantly
 * reading and writing this file to the FS, just to memory, which should
 * speed things up.

 * The idea is that this code will be called serially but repeatedly,
 * iterating over hundreds or thousands of HDF5 files that have user
 * blocks already created (if they don't, the program doesn't modify the
 * hdf5 files). I'm thinking I'll be calling it like:

 *  for i in (find . -name \*hdf5); do
 *      stash --namelist=namelist.input --output=cm1.out \\
 *     	 --sounding=input_sounding --pbs=cm1.pbs --stash=README.stash --hdf=$i
 *  done

 * This results in no separate userblock.txt file since it's destroyed
 * after the code runs (my first version of this code did a regular
 * fopen on the userblock file).

 * I have stuck an EOF character at the end of our virtual file so we
 * don't end up writing the entire buffer to the HDF5 file (unless it
 * is completely filled of course) since there is no EOF character by
 * default in the virtual file.

 * Regarding the use of a strcpy for time information, the following is
 * from the ctime man page:

 * "The return value points to a statically allocated string which
 * might be overwritten by subsequent calls to any of the date and time
 * functions."
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <getopt.h>
#include <hdf5.h>

/* MAXCOLS is the number of characters in a line of text to scan; assumes
 * formated output with line breaks */

/* Currently the user block size I use, minus one, for null byte */
#define MAXFILEBYTES (32767)
#define MAXFILENAME    (200)
#define MAXCOLS        (200)
#define MAXHOSTNAME    (100)
#define MAXDATE        (100)

void exit_malloc(int size);
void exit_fopen(char *filename);
void exit_hdf(char *routine, char *filename);
void exit_fmemopen(int size);
void parse_cmdline(int argc, char *argv[],
		char *namelist_filename,int *got_namelist,char *output_filename,
		int *got_output,char *sounding_filename,int *got_sounding,
		char *hdf_filename, int *got_hdf,char *pbs_filename,int *got_pbs,
		char *stash_filename, int *got_stash);

int main (int argc, char *argv[])
{
	FILE *fp_namelist,*fp_output,*fp_sounding,*fp_userblock,*fp_userblock2;
	FILE *fp_hdf,*fp_pbs;
	char *namelist_filename,*output_filename,*sounding_filename,*hdf_filename,*pbs_filename;
	char c;
	int n_output_header_lines = 10;
	int n_orftxt_lines;
	float pct;

	char line[MAXCOLS];
	char *my_hostname;
	char *current_date_and_time;
	char *hdf_atime,*hdf_mtime,*hdf_ctime;
	char *rval;
	char *cwd;
	void *vfile,*vfile2;

	hid_t f_id, plist_id;
	hsize_t userblock_size;

	time_t rawtime;
	struct stat statbuffer;

	int i,j;
	int chrsize;

	int nchar_userblock;

	char *stash_filename;
	FILE *fp_stash;
	int got_stash=0,got_namelist=0,got_sounding=0,got_output=0,got_pbs=0,got_hdf=0;

	char *alph="abcdefghijklmnopqrstuvwxyz";

	/* A friendly header */

	const char *orftxt[] = {
		"********************************************************************\n",
		"*******************Welcome to the HDF5 userblock********************\n",
		"********************************************************************\n",
		"This file was created with CM1 with I/O modifications\n",
		"by Leigh Orf <leigh.orf@wisc.edu>. Useful information is stashed here\n",
		"including the namelist.input file and the sounding file.\n",
		"Have a nice day!\n"
	};

	const char *divider =
		"--------------------------------------------------------------------\n";

	/* allocate space for strings */

	chrsize=MAXDATE*sizeof(char *);
	if ((current_date_and_time = (char *) malloc(chrsize)) == NULL) exit_malloc(chrsize);
	if ((hdf_atime =             (char *) malloc(chrsize)) == NULL) exit_malloc(chrsize);
	if ((hdf_mtime =             (char *) malloc(chrsize)) == NULL) exit_malloc(chrsize);
	if ((hdf_ctime =             (char *) malloc(chrsize)) == NULL) exit_malloc(chrsize);
	chrsize=MAXHOSTNAME*sizeof(char *);
	if ((my_hostname =           (char *) malloc(chrsize)) == NULL) exit_malloc(chrsize);
	chrsize=MAXFILENAME*sizeof(char *);
	if ((namelist_filename =     (char *) malloc(chrsize)) == NULL) exit_malloc(chrsize);
	if ((output_filename =       (char *) malloc(chrsize)) == NULL) exit_malloc(chrsize);
	if ((sounding_filename =     (char *) malloc(chrsize)) == NULL) exit_malloc(chrsize);
	if ((hdf_filename =          (char *) malloc(chrsize)) == NULL) exit_malloc(chrsize);
	if ((cwd =                   (char *) malloc(chrsize)) == NULL) exit_malloc(chrsize);
	if ((pbs_filename =          (char *) malloc(chrsize)) == NULL) exit_malloc(chrsize);
	if ((stash_filename =        (char *) malloc(chrsize)) == NULL) exit_malloc(chrsize);

	/* get command line options */

	parse_cmdline(argc,argv,
			namelist_filename,&got_namelist,output_filename,&got_output,sounding_filename,&got_sounding,
			hdf_filename,&got_hdf,pbs_filename,&got_pbs,stash_filename,&got_stash);

	/* turn off buffering to stdout */
	setbuf(stdout,NULL);

	/* Get the machine name and current working directory */
	if ((gethostname(my_hostname,(size_t)MAXHOSTNAME))<0) strcpy(my_hostname,"UNKNOWN");
	if ((getcwd(cwd,chrsize))==NULL) strcpy(cwd,"UNKNOWN");

	/* Current date and time */
	if ((rawtime=time(NULL))<(time_t)0) strcpy(current_date_and_time,"UNKNOWN");
	else  strcpy(current_date_and_time,ctime(&rawtime));

	/* Get the stat timestamp info on the hdf5 file */
	if ((stat(hdf_filename,&statbuffer))<0)
	{
		strcpy(hdf_atime,"UNKNOWN");
		strcpy(hdf_mtime,"UNKNOWN");
		strcpy(hdf_ctime,"UNKNOWN");
	}
	else
	{
		strcpy(hdf_atime,ctime(&(statbuffer.st_atime)));
		strcpy(hdf_mtime,ctime(&(statbuffer.st_mtime)));
		strcpy(hdf_ctime,ctime(&(statbuffer.st_ctime)));
	}

	chrsize=MAXCOLS*sizeof(char *);
	if ((fp_namelist  =  fopen(namelist_filename,"r"))   == NULL) exit_fopen(namelist_filename);
	if ((fp_output    =  fopen(output_filename,"r"))     == NULL) exit_fopen(output_filename);
	if ((fp_sounding  =  fopen(sounding_filename,"r"))   == NULL) exit_fopen(sounding_filename);

	chrsize=MAXFILEBYTES*sizeof(char);
	if ((vfile  =  (void *) calloc(MAXFILEBYTES,sizeof(char))) == NULL) exit_malloc(chrsize);
	if ((vfile2 =  (void *) calloc(MAXFILEBYTES,sizeof(char))) == NULL) exit_malloc(chrsize);
	if ((fp_userblock  =  fmemopen(vfile,chrsize,"r+"))   == NULL) exit_fmemopen(chrsize);
	if ((fp_userblock2 =  fmemopen(vfile2,chrsize,"r+"))  == NULL) exit_fmemopen(chrsize);

	n_orftxt_lines = (sizeof(orftxt)/sizeof(const char *));

   /* 
	* First we construct our userblock. Then we go through and count the
	* number of (new)lines. The number of lines in the userblock is then
	* the first line in the userblock data that is written to the hdf5
	* file. This makes it easy to just grab it if we need it, i.e.,
	*      head -$(head -1 blah.hdf) blah.hdf > userblock.txt
    */

	/* Print header to userblock file */
	for (i=0; i<n_orftxt_lines; i++) fprintf(fp_userblock,"%s",orftxt[i]);
	fprintf(fp_userblock,"%s",divider);
	/* Print hostname and timestamp info to userblock file */
	fprintf(fp_userblock,"The following information was collected at the time the userblock was modified,\n");
	fprintf(fp_userblock,"not when the hdf5 file was written: \n");
	fprintf(fp_userblock,"\tHostname = %s\n",my_hostname);
	fprintf(fp_userblock,"\tCurrent working directory = %s\n",cwd);
	fprintf(fp_userblock,"\tLocal time = %s",current_date_and_time);
	fprintf(fp_userblock,"\tHDF file atime = %s",hdf_atime);
	fprintf(fp_userblock,"\tHDF file mtime = %s",hdf_mtime);
	fprintf(fp_userblock,"\tHDF file ctime = %s",hdf_ctime);
	fprintf(fp_userblock,"%s",divider);

	if(got_stash) /* Custom stash block file */
	{
		if ((fp_stash = fopen(stash_filename,"r"))==NULL) exit_fopen(stash_filename);
		fprintf(fp_userblock,"User-selected optional stash data from file %s follows:\n",stash_filename);
		while((c = fgetc(fp_stash))!=EOF) fputc(c,fp_userblock);
		fprintf(fp_userblock,"%s",divider);
		fclose(fp_stash);
	}

	if(got_pbs)
	{
		if ((fp_pbs = fopen(pbs_filename,"r"))==NULL) exit_fopen(pbs_filename);
		fprintf(fp_userblock,"PBS file %s follows:\n",pbs_filename);
		while((c = fgetc(fp_pbs))!=EOF) fputc(c,fp_userblock);
		fprintf(fp_userblock,"%s",divider);
		fclose(fp_pbs);
	}

	/* Print head of model output file to userblock file */
	fprintf(fp_userblock,"First %i lines of model output file %s:\n",n_output_header_lines,output_filename);
	for (i=0; i<n_output_header_lines; i++)
	{
		rval = fgets(line,sizeof(line),fp_output);
		fputs(line,fp_userblock);
	}
	fprintf(fp_userblock,"%s",divider);

	/* Print namelist file to userblock file */
	fprintf(fp_userblock,"Namelist file %s:\n",namelist_filename);
	while((rval = fgets(line,sizeof(line),fp_namelist))!=NULL) fputs(line,fp_userblock);

	fprintf(fp_userblock,"%s",divider);

	/* Print sounding file to userblock file */
	fprintf(fp_userblock,"Sounding file %s:\n",sounding_filename);
	while((rval = fgets(line,sizeof(line),fp_sounding))!=NULL) fputs(line,fp_userblock);

	fprintf(fp_userblock,"%s",divider);
	fputc(EOF,fp_userblock); // ORF stick our own EOF character here

	fclose(fp_namelist);
	fclose(fp_output);
	fclose(fp_sounding);
	rewind(fp_userblock);
/*
 *  We have created the file that we'll put in the userblock of the hdf5
 *  file. Now, get the userblock size from the hdf5 file to make damned
 *  sure we don't accidentally clobber HDF5 data in the file
*/

	if ((f_id = H5Fopen(hdf_filename,H5F_ACC_RDONLY,H5P_DEFAULT)) < 0) exit_hdf("H5Fopen",hdf_filename);
	if ((plist_id = H5Fget_create_plist(f_id)) < 0) exit_hdf("H5Fget_create_plist",hdf_filename);
	if ((H5Pget_userblock(plist_id,&userblock_size)) < 0) exit_hdf("H5Pget_userblock",hdf_filename);
	if ((H5Pclose(plist_id)) < 0) exit_hdf("H5Pclose",hdf_filename);
	if ((H5Fclose(f_id)) < 0) exit_hdf("H5Fclose",hdf_filename);

	/* Check that we have a userblock at all */

	if(userblock_size == 0) // There is no userblock in this file
	{
		printf("*");
		exit (-1);
	}

	i=0;

	/* get the # of chars in userblock data */

	while((c = fgetc(fp_userblock))!=EOF) i++;

	nchar_userblock = i;

	/* OK I really want to have the 1st line of the header be the number
	 * of lines of ascii text so I can head it to a file easily so..... */

	rewind(fp_userblock);

	/* calculate the number of lines by counting the newline characters */

	j=0;
	for(i=0;i<nchar_userblock;i++)if(fgetc(fp_userblock)=='\n')j++;

	fprintf(fp_userblock2,"%i\n",j+1); /* Add 1 for the line we just added! */

	rewind(fp_userblock);
	for(i=0;i<nchar_userblock;i++)fputc((char)fgetc(fp_userblock),fp_userblock2);
	fputc(EOF,fp_userblock2); // ORF stick our own EOF character here
	fclose(fp_userblock);

	/* Finally, do a byte by byte copy to the user block region of the
	 * HDF5 file which is at the very beginning */

	if ((fp_hdf = fopen(hdf_filename,"r+"))==NULL) exit_fopen(hdf_filename);

	i=0;
	rewind(fp_userblock2);
	while((c = fgetc(fp_userblock2))!=EOF)
	{
		i++;
		if(i==userblock_size) break; //We have entirely filled the user block, do not write more
		fputc(c,fp_hdf);
	}

	fclose(fp_userblock2);
	fclose(fp_hdf);

	pct = (float)i/(float)userblock_size;
	j = (int)(pct*25);
	pct*=100.0;
	printf("%c",alph[j]);

	free(current_date_and_time);
	free(my_hostname);
	free(hdf_atime);
	free(hdf_mtime);
	free(hdf_ctime);
	free(cwd);
	free(namelist_filename);
	free(output_filename);
	free(sounding_filename);
	free(hdf_filename);
	free(pbs_filename);
	free(stash_filename);

	exit(0);

}

void exit_malloc(int size)
{
	printf("You have got to be shitting me, we can't malloc %i friggin bytes?\n",size);
	printf("I refuse to believe this will ever happen. Antwerp.\n");
	exit(-1);
}

void exit_fopen(char *filename)
{
	printf("Well, that sucks, we can't open %s, bailing.\n",filename);
	exit(-1);
}

void exit_hdf(char *routine, char *filename)
{
	printf("HDF routine %s failed for %s, bailing\n",routine,filename);
	exit(-1);
}


void exit_fmemopen(int size)
{
	printf("fmemopen failed with a requested size of %i bytes, dammit\n",size);
	exit(-1);
}

void parse_cmdline(int argc, char *argv[],
	 char *namelist_filename, int *got_namelist, char *output_filename, int *got_output,
	 char *sounding_filename, int *got_sounding,char *hdf_filename, int *got_hdf,
	 char *pbs_filename, int *got_pbs, char *stash_filename, int *got_stash)
{
	if(argc==1)
	{
		printf("Usage: stash --namelist [namelist_file] --sounding [sounding_file] --output [cm1out_file] --hdf [hdf5_file] --pbs [opt_pbs_file] --stash [opt_stash-file]\n");
		exit(0);
	}
	while (1)
	{
		static struct option long_options[] = {
			{"namelist", required_argument, 0 , 'n'},
			{"output",   required_argument, 0,  'o'},
			{"sounding", required_argument, 0,  's'},
			{"hdf",      required_argument, 0,  'h'},
			{"pbs",      optional_argument, 0,  'p'},
			{"stash",    optional_argument, 0,  'z'}
		};

		int r;
		int option_index = 0;
		r = getopt_long (argc, argv,"n:o:s:h:z:",long_options,&option_index);
		if (r == -1) break;

		switch(r)
		{
			case 'n':
				strcpy(namelist_filename,optarg);
				*got_namelist=1;
				break;
			case 'o':
				strcpy(output_filename,optarg);
				*got_output=1;
				break;
			case 's':
				strcpy(sounding_filename,optarg);
				*got_sounding=1;
				break;
			case 'h':
				strcpy(hdf_filename,optarg);
				*got_hdf=1;
				break;
			case 'p':
				strcpy(pbs_filename,optarg);
				*got_pbs=1;
				break;
			case 'z':
				strcpy(stash_filename,optarg);
				*got_stash = 1;
				break;
		}
	}

	if(*got_namelist==0)
	{
		printf("No namelist file specified, exiting\n");
		exit(0);
	}
	if(*got_sounding==0)
	{
		printf("No sounding file specified, exiting\n");
		exit(0);
	}
	if(*got_output==0)
	{
		printf("No output file specified, exiting\n");
		exit(0);
	}
	if(*got_hdf==0)
	{
		printf("No hdf file specified, exiting\n");
		exit(0);
	}
}
