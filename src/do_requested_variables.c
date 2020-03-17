#include <omp.h>
#include "../include/lofs-read.h"
#include "../include/dirstruct.h"
#include "../include/hdf2nc.h"
#include "../include/limits.h"


void do_requested_variables(buffers *b, ncstruct nc, grid gd,readahead rh,dir_meta dm,hdf_meta hm,cmdline cmd)
{
	int ivar,status;
	requested_cube rc;
	char *var;

 //For flexibility we always set rc in case we need to read outside of what we
 //requested at the command line such as what occurs with staggered variables. The
 //copy_grid_to_requested_cube just duplicates the grid structure stuff to rc (requested
 //cube). rc is a required argument to read_lofs_buffer.

 //For doing calculations that involve derivatives (usually of u, v, or w) we have to read
 //data in slightly differently (see for instance do_readahead where we read in ustag
 //vstag wstag on a bigger mesh than the scalar variables).

	copy_grid_to_requested_cube(&rc,gd);

	var = (char *) malloc (MAXSTR * sizeof(char));

	for (ivar = 0; ivar < cmd.nvar; ivar++)
	{
		var=nc.varname[ivar];

		if(same(var,"uinterp")) {;}
		if(same(var,"vinterp")) {;}
		if(same(var,"winterp")) {;}
		if(same(var,"hwin_sr")) {;}
		if(same(var,"hwin_gr")) {;}
		if(same(var,"windmag_sr")) {;}
		if(same(var,"hdiv")) {;}
		if(same(var,"xvort")) {;}
		if(same(var,"yvort")) {;}
		if(same(var,"zvort")) {;}
		if(same(var,"vortmag")) {;}
		if(same(var,"streamvort")) {;}
		else
		{
			read_lofs_buffer(b->buf,nc.varname[ivar],dm,hm,rc,cmd);
		}
		status = nc_put_vara_float (nc.ncid, nc.varnameid[ivar], nc.start, nc.edges, b->buf);
	}
}
