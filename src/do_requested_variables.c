#include <omp.h>
#include "../include/lofs-read.h"
#include "../include/dirstruct.h"
#include "../include/hdf2nc.h"
#include "../include/limits.h"
#include "../include/macros.h"

#define UINTERP BUF
void calc_uinterp(buffers *b, grid gd, cmdline cmd)
{
	int ix,iy,iz,nx,ny,nz;
	nx=gd.NX;ny=gd.NY;nz=gd.NZ;

	if(cmd.verbose) printf("Calculating uinterp...\n");
#pragma omp parallel for private(ix,iy,iz)
	for(iz=0; iz<nz; iz++)
	for(iy=0; iy<ny; iy++)
	for(ix=0; ix<nx; ix++)
		UINTERP(ix,iy,iz) = 0.5*(UA(ix,iy,iz)+UA(ix+1,iy,iz));
}

#define VINTERP BUF
void calc_vinterp(buffers *b, grid gd, cmdline cmd)
{
	int ix,iy,iz,nx,ny,nz;
	nx=gd.NX;ny=gd.NY;nz=gd.NZ;

	if(cmd.verbose) printf("Calculating vinterp...\n");

#pragma omp parallel for private(ix,iy,iz)
	for(iz=0; iz<nz; iz++)
	for(iy=0; iy<ny; iy++)
	for(ix=0; ix<nx; ix++)
		VINTERP(ix,iy,iz) = 0.5*(VA(ix,iy,iz)+VA(ix,iy+1,iz));
}

#define WINTERP BUF
void calc_winterp(buffers *b, grid gd, cmdline cmd)
{
	int ix,iy,iz,nx,ny,nz;
	nx=gd.NX;ny=gd.NY;nz=gd.NZ;

	if(cmd.verbose) printf("Calculating winterp...\n");

#pragma omp parallel for private(ix,iy,iz)
	for(iz=0; iz<nz; iz++)
	for(iy=0; iy<ny; iy++)
	for(ix=0; ix<nx; ix++)
		WINTERP(ix,iy,iz) = 0.5*(WA(ix,iy,iz)+WA(ix,iy,iz+1));
}

#define HWIN_SR BUF
void calc_hwin_sr(buffers *b, grid gd, cmdline cmd)
{
	int ix,iy,iz,nx,ny,nz;
	float usr,vsr;
	nx=gd.NX;ny=gd.NY;nz=gd.NZ;

	if(cmd.verbose) printf("Calculating hwin_sr...\n");

#pragma omp parallel for private(ix,iy,iz,usr,vsr)
	for(iz=0; iz<nz; iz++)
	for(iy=0; iy<ny; iy++)
	for(ix=0; ix<nx; ix++)
	{
		usr = 0.5*(UA(ix,iy,iz)+UA(ix+1,iy,iz));
		vsr = 0.5*(VA(ix,iy,iz)+VA(ix,iy+1,iz));
		HWIN_SR(ix,iy,iz) = sqrt(usr*usr+vsr*vsr);
	}
}
#define HWIN_SR BUF
void calc_hwin_gr(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int ix,iy,iz,nx,ny,nz;
	float usr,vsr;
	nx=gd.NX;ny=gd.NY;nz=gd.NZ;

	if(cmd.verbose) printf("Calculating hwin_gr...\n");

#pragma omp parallel for private(ix,iy,iz,usr,vsr)
	for(iz=0; iz<nz; iz++)
	for(iy=0; iy<ny; iy++)
	for(ix=0; ix<nx; ix++)
	{
		usr = 0.5*(UA(ix,iy,iz)+UA(ix+1,iy,iz)) + msh.umove;
		vsr = 0.5*(VA(ix,iy,iz)+VA(ix,iy+1,iz)) + msh.vmove;
		HWIN_SR(ix,iy,iz) = sqrt(usr*usr+vsr*vsr);
	}
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

void do_requested_variables(buffers *b, ncstruct nc, grid gd, mesh msh, readahead rh,dir_meta dm,hdf_meta hm,cmdline cmd)
{
	int ivar,status;
	requested_cube rc;
	char *var;
	float *wp;

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

		if(same(var,"uinterp"))			calc_uinterp(b,gd,cmd);
		else if(same(var,"vinterp"))	calc_vinterp(b,gd,cmd);
		else if(same(var,"winterp"))	calc_winterp(b,gd,cmd);
		else if(same(var,"hwin_sr"))	calc_hwin_sr(b,gd,cmd);
		else if(same(var,"hwin_gr"))	calc_hwin_gr(b,gd,msh,cmd);
		else if(same(var,"windmag_sr"))	{calc_windmag_sr();}
		else if(same(var,"hdiv")) {calc_hdiv();}
		else if(same(var,"xvort")) {calc_xvort();}
		else if(same(var,"yvort")) {calc_yvort();}
		else if(same(var,"zvort")) {calc_zvort();}
		else if(same(var,"vortmag")) {calc_vortmag();}
		else if(same(var,"streamvort")) {calc_streamvort();}
		else
		{
			read_lofs_buffer(b->buf,nc.varname[ivar],dm,hm,rc,cmd);
		}
		if(cmd.verbose)printf("Writing %s...\n",var);
		status = nc_put_vara_float (nc.ncid, nc.varnameid[ivar], nc.start, nc.edges, b->buf);
	}
}

