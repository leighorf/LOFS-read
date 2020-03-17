#include <omp.h>
#include "../include/lofs-read.h"
#include "../include/dirstruct.h"
#include "../include/hdf2nc.h"
#include "../include/limits.h"

void calc_uinterp()
{
}
void calc_vinterp()
{
}
void calc_winterp()
{
}
void calc_hwin_sr()
{
}
void calc_hwin_gr()
{
}
void calc_windmag_sr()
{
}
void calc_hdiv()
{
}
void calc_xvort()
{
}
void calc_yvort()
{
}
void calc_zvort()
{
}
void calc_vortmag()
{
}
void calc_streamvort()
{
}

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

		if(same(var,"uinterp")) {calc_uinterp();}
		if(same(var,"vinterp")) {calc_vinterp();}
		if(same(var,"winterp")) {calc_winterp();}
		if(same(var,"hwin_sr")) {calc_hwin_sr();}
		if(same(var,"hwin_gr")) {calc_hwin_gr();}
		if(same(var,"windmag_sr")) {calc_windmag_sr();}
		if(same(var,"hdiv")) {calc_hdiv();}
		if(same(var,"xvort")) {calc_xvort();}
		if(same(var,"yvort")) {calc_yvort();}
		if(same(var,"zvort")) {calc_zvort();}
		if(same(var,"vortmag")) {calc_vortmag();}
		if(same(var,"streamvort")) {calc_streamvort();}
		else
		{
			read_lofs_buffer(b->buf,nc.varname[ivar],dm,hm,rc,cmd);
		}
		status = nc_put_vara_float (nc.ncid, nc.varnameid[ivar], nc.start, nc.edges, b->buf);
	}
}

