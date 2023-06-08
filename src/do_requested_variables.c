#include <omp.h>
#include "../include/lofs-read.h"
#include "../include/lofs-dirstruct.h"
#include "../include/lofs-hdf2nc.h"
#include "../include/lofs-limits.h"
#include "../include/lofs-macros.h"
#include "../include/lofs-constants.h"
#include "./calc/calcvort.c"
#include "./calc/calcmomentum.c"

/*******************************************************************************/

void do_pipert(buffers *b, grid gd, sounding *snd, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	#pragma omp parallel for private(i,j,k) 
	// We're going to calculate over the whole
	// domain size read because often times
	// derivatives operate on this data. The
	// index macros ensure that -1 is considered in-bounds.
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
    	calc_pipert(b->ppert, snd->pres0, b->buf0, i, j, k, ni, nj);
	}
	}
	}

}

/*******************************************************************************/
//NEED INTERPOLATED OPTION FOR THIS
#define WBUOY BUFp
void do_wbuoy(buffers *b, grid gd, sounding *snd, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	#pragma omp parallel for private(i,j,k) 
	for(k=1; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
    	calc_buoyancy(b->thrhopert, snd->th0, b->buf0, i, j, k, ni, nj);
	}
	}
	}

	// lower boundary condition
	#pragma omp parallel for private(i, j)
	for (j=-1; j<nj+1; j++) {
	for (i=-1; i<ni+1; i++) {
		WBUOY(i, j, 0) = 0.0;
	}
	}

}

/*
#define WBUOY BUFp
#define WBUOY_INTERP BUFp
void do_wbuoy_interp(buffers *b, grid gd, sounding *snd, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	#pragma omp parallel for private(i,j,k) 
	for(k=1; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
    	calc_buoyancy(b->thrhopert, snd->th0, b->buf0, i, j, k, ni, nj);
	}
	}
	}

	// lower boundary condition
	#pragma omp parallel for private(i, j)
	for (j=-1; j<nj+1; j++) {
	for (i=-1; i<ni+1; i++) {
		WBUOY(i, j, 0) = 0.0;
	}
	}
	//
//Interpolate to scalar mesh
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
		WBUOY_INTERP(i, j, k) = 0.5*(WBUOY(i, j, k) + WBUOY(i, j, k+1));
	}
	}
	}

}
*/

/*******************************************************************************/

#define WBUOY BUFp
#define WBUOY_INTERP BUFp
void do_wbuoy_interp(buffers *b, grid gd, sounding *snd, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
    float dz,val;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	// need to calculate pipert first
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
    	calc_buoyancy(b->thrhopert, snd->th0, b->buf0, i, j, k, ni, nj);
	}
	}
	}

	// lower boundary condition
	#pragma omp parallel for private(i, j)
	for (j=-1; j<nj+1; j++) {
	for (i=-1; i<ni+1; i++) {
		WBUOY(i, j, 0) = 0.0;
	}
	}
//Interpolate to scalar mesh
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
		WBUOY_INTERP(i, j, k) = 0.5*(WBUOY(i, j, k) + WBUOY(i, j, k+1));
	}
	}
	}

}

/*******************************************************************************/

#define WPGRAD BUFp
void do_wpgrad(buffers *b, grid gd, sounding *snd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
    float dz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	// need to calculate pipert first
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
    	calc_pipert(b->ppert, snd->pres0, b->dum0, i, j, k, ni, nj);
	}
	}
	}

	#pragma omp parallel for private(i,j,k,dz) 
	for(k=1; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
		dz = 1./(msh.rdz * MF(k)); 
    	calc_pgrad_w(b->dum0, b->thrhopert, snd->qv0, snd->th0, b->buf0, dz, i, j, k, ni, nj);
	}
	}
	}

	// lower boundary condition
	#pragma omp parallel for private(i, j)
	for (j=-1; j<nj+1; j++) {
	for (i=-1; i<ni+1; i++) {
		WPGRAD(i, j, 0) = 0.0;
	}
	}

}

/*******************************************************************************/

#define UPGRAD BUFp
void do_upgrad(buffers *b, grid gd, sounding *snd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
    float dx;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	// need to calculate pipert first
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
    	calc_pipert(b->ppert, snd->pres0, b->dum0, i, j, k, ni, nj);
	}
	}
	}

	#pragma omp parallel for private(i,j,k,dx) 
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=0; i<ni+1; i++) {
		dx = 1./(msh.rdx * UF(i)); 
    	calc_pgrad_u(b->dum0, b->thrhopert, snd->qv0, snd->th0, b->buf0, dx, i, j, k, ni, nj);
	}
	}
	}

}

/*******************************************************************************/

#define VPGRAD BUFp
void do_vpgrad(buffers *b, grid gd, sounding *snd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
    float dy;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	// need to calculate pipert first
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
    	calc_pipert(b->ppert, snd->pres0, b->dum0, i, j, k, ni, nj);
	}
	}
	}

	#pragma omp parallel for private(i,j,k,dy) 
	for(k=0; k<nk+1; k++) {
	for(j=0; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
		dy = 1./(msh.rdy * VF(j)); 
    	calc_pgrad_v(b->dum0, b->thrhopert, snd->qv0, snd->th0, b->buf0, dy, i, j, k, ni, nj);
	}
	}
	}
}

/*******************************************************************************/

#define WPGRADINTERP BUFp
#define WPGRAD BUFp
void do_wpgrad_interp(buffers *b, grid gd, sounding *snd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
    float dz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	// need to calculate pipert first
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
    	calc_pipert(b->ppert, snd->pres0, b->dum0, i, j, k, ni, nj);
	}
	}
	}

	#pragma omp parallel for private(i,j,k,dz) 
	for(k=1; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
		dz = 1./(msh.rdz * MF(k)); 
    	calc_pgrad_w(b->dum0, b->thrhopert, snd->qv0, snd->th0, b->buf0, dz, i, j, k, ni, nj);
	}
	}
	}

	// lower boundary condition
	#pragma omp parallel for private(i, j)
	for (j=-1; j<nj+1; j++) {
	for (i=-1; i<ni+1; i++) {
		WPGRAD(i, j, 0) = 0.0;
	}
	}
//Interpolate to scalar mesh
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
		WPGRADINTERP(i, j, k) = 0.5*(WPGRAD(i, j, k) + WPGRAD(i, j, k+1));
	}
	}
	}

}

/*******************************************************************************/

#define UPGRAD BUFp
#define UPGRADINTERP BUFp
void do_upgrad_interp(buffers *b, grid gd, sounding *snd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
    float dx;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	// need to calculate pipert first
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
    	calc_pipert(b->ppert, snd->pres0, b->dum0, i, j, k, ni, nj);
	}
	}
	}

	#pragma omp parallel for private(i,j,k,dx) 
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=0; i<ni+1; i++) {
		dx = 1./(msh.rdx * UF(i)); 
    	calc_pgrad_u(b->dum0, b->thrhopert, snd->qv0, snd->th0, b->buf0, dx, i, j, k, ni, nj);
	}
	}
	}

//Interpolate to scalar mesh
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=0; i<ni; i++) {
		UPGRADINTERP(i, j, k) = 0.5*(UPGRAD(i, j, k) + UPGRAD(i+1, j, k));
	}
	}
	}

}

/*******************************************************************************/

#define VPGRADINTERP BUFp
#define VPGRAD BUFp
void do_vpgrad_interp(buffers *b, grid gd, sounding *snd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
    float dy;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	// need to calculate pipert first
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
    	calc_pipert(b->ppert, snd->pres0, b->dum0, i, j, k, ni, nj);
	}
	}
	}

	#pragma omp parallel for private(i,j,k,dy) 
	for(k=0; k<nk+1; k++) {
	for(j=0; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
		dy = 1./(msh.rdy * VF(j)); 
    	calc_pgrad_v(b->dum0, b->thrhopert, snd->qv0, snd->th0, b->buf0, dy, i, j, k, ni, nj);
	}
	}
	}

//Interpolate to scalar mesh
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk+1; k++) {
	for(j=0; j<nj; j++) {
	for(i=-1; i<ni+1; i++) {
		VPGRADINTERP(i, j, k) = 0.5*(VPGRAD(i, j, k) + VPGRAD(i, j+1, k));
	}
	}
	}

}

/*
#define RHO BUFp

void calc_rho(buffers *b, grid gd, mesh msh, cmdline cmd,dir_meta dm,hdf_meta hm, sounding *snd, requested_cube rc)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float usr,vsr;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;
	float pi;
	float foo;
	int ifoo;

rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1;
rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;
//	How Kelton does it (pointwise) in calcvort.c
//	for(k=0; k<nk+1; k++) {
//	for(j=-1; j<nj+1; j++) {
//	for(i=-1; i<ni+1; i++) {
 //   	calc_pipert(b->ppert, snd->pres0, b->buf0, i, j, k, ni, nj);
//	}
//	}
//	}

	read_lofs_buffer(b->buf0,"rhopert",dm,hm,rc,cmd);
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk+1; k++)
	for(j=-1; j<nj+1; j++)
	for(i=-1; i<ni+1; i++)
	{
		RHO(i,j,k) = BUFp(i,j,k) + snd->rho0[k];
	}

}
*/

/*******************************************************************************/

void buf_u(buffers *b,grid gd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		BUFp(i,j,k)=UAp(i,j,k);
}

void buf_v(buffers *b,grid gd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		BUFp(i,j,k)=VAp(i,j,k);
}

void buf_w(buffers *b,grid gd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		BUFp(i,j,k)=WAp(i,j,k);
}

void buf_prespert(buffers *b,grid gd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		BUFp(i,j,k)=PRESPERTp(i,j,k)
}

void buf_thrhopert(buffers *b,grid gd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		BUFp(i,j,k)=THRHOPERTp(i,j,k)
}

/*******************************************************************************/

/* More Memory Management Madness!!!!!
 * Check if we've already buffered u,v,w before doing interps */

#define UINTERP BUFp
#define U BUFp
void calc_uinterp(buffers *b, grid gd, mesh msh, cmdline cmd,dir_meta dm,hdf_meta hm,requested_cube rc,readahead rh)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
	rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1;
	rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;

	if(rh.u)
	{
		buf_u(b,gd); 
	}
	else
	{
		printf("reading...");FL;
		read_lofs_buffer(b->buf0,"u",dm,hm,rc,cmd);
	}

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		UINTERP(i,j,k) = 0.5*(U(i,j,k)+U(i+1,j,k));
}

/*******************************************************************************/

#define VINTERP BUFp
#define V BUFp
void calc_vinterp(buffers *b, grid gd, mesh msh, cmdline cmd,dir_meta dm,hdf_meta hm,requested_cube rc,readahead rh)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
	rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1;
	rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;

	if(rh.v)
	{
		buf_v(b,gd); 
	}
	else
	{
		printf("reading...");FL;
		read_lofs_buffer(b->buf0,"v",dm,hm,rc,cmd);
	}

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		VINTERP(i,j,k) = 0.5*(V(i,j,k)+V(i,j+1,k));
}

/*******************************************************************************/

#define WINTERP BUFp
#define W BUFp
void calc_winterp(buffers *b, grid gd, mesh msh, cmdline cmd,dir_meta dm,hdf_meta hm,requested_cube rc,readahead rh)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
	rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1+1;
	rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;

	if(rh.w)
	{
		buf_w(b,gd); 
	}
	else
	{
		printf("reading...");FL;
		read_lofs_buffer(b->buf0,"w",dm,hm,rc,cmd);
	}

#pragma omp parallel for private(i,j,k)
	for(k=0;k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		WINTERP(i,j,k) = 0.5*(W(i,j,k)+W(i,j,k+1));
}

/*******************************************************************************/
#define KHH_INTERP BUFp
#define KHH BUFp
void calc_khh_interp(buffers *b, grid gd, mesh msh, cmdline cmd,dir_meta dm,hdf_meta hm,requested_cube rc)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
	rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1;
	rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;

	read_lofs_buffer(b->buf0,"khh",dm,hm,rc,cmd);

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		KHH_INTERP(i,j,k) = 0.5*(KHH(i,j,k)+KHH(i,j,k+1));
}

/*******************************************************************************/
#define KMH_INTERP BUFp
#define KMH BUFp
void calc_kmh_interp(buffers *b, grid gd, mesh msh, cmdline cmd,dir_meta dm,hdf_meta hm,requested_cube rc)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
	rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1;
	rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;

	read_lofs_buffer(b->buf0,"kmh",dm,hm,rc,cmd);

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		KMH_INTERP(i,j,k) = 0.5*(KMH(i,j,k)+KMH(i,j,k+1));
}

/*******************************************************************************/
#define KHV_INTERP BUFp
#define KHV BUFp
void calc_khv_interp(buffers *b, grid gd, mesh msh, cmdline cmd,dir_meta dm,hdf_meta hm,requested_cube rc)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
	rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1;
	rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;

	read_lofs_buffer(b->buf0,"khv",dm,hm,rc,cmd);

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		KHV_INTERP(i,j,k) = 0.5*(KHV(i,j,k)+KHV(i,j,k+1));
}

/*******************************************************************************/
#define KMV_INTERP BUFp
#define KMV BUFp
void calc_kmv_interp(buffers *b, grid gd, mesh msh, cmdline cmd,dir_meta dm,hdf_meta hm,requested_cube rc)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
	rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1;
	rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;

	read_lofs_buffer(b->buf0,"kmv",dm,hm,rc,cmd);

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		KMV_INTERP(i,j,k) = 0.5*(KMV(i,j,k)+KMV(i,j,k+1));
}

/*******************************************************************************/
//else if(same(var,"qiqvpert"))	   {CL;calc_qiqvpert(b,gd,msh,cmd,dm,hm,rc);}
//void do_requested_variables(buffers *b, ncstruct nc, grid gd, mesh msh, sounding *snd, readahead rh,dir_meta dm,hdf_meta hm,cmdline cmd)
//Just read in sequentially and add, don't bother with readahead (except
//the use of readahead flag to allocate an additional buffer in
//do_readahead)
//
#define QIQVPERT BUFp
void calc_qiqvpert(buffers *b, grid gd, mesh msh, cmdline cmd,dir_meta dm,hdf_meta hm,requested_cube rc)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float usr,vsr;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
	rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1;
	rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;

//Just read two fields and add them

	read_lofs_buffer(b->buf0,"qvpert",dm,hm,rc,cmd);

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) = BUFp(i,j,k);

	read_lofs_buffer(b->buf0,"qi",dm,hm,rc,cmd);

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		QIQVPERT(i,j,k) += TEM1p(i,j,k);
}

#define RHO BUFp

void calc_rho(buffers *b, grid gd, mesh msh, cmdline cmd,dir_meta dm,hdf_meta hm, sounding *snd, requested_cube rc)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float usr,vsr;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;
	float pi;
	float foo;
	int ifoo;

rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1;
rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;
//	How Kelton does it (pointwise) in calcvort.c
//	for(k=0; k<nk+1; k++) {
//	for(j=-1; j<nj+1; j++) {
//	for(i=-1; i<ni+1; i++) {
 //   	calc_pipert(b->ppert, snd->pres0, b->buf0, i, j, k, ni, nj);
//	}
//	}
//	}

	read_lofs_buffer(b->buf0,"rhopert",dm,hm,rc,cmd);
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk+1; k++)
	for(j=-1; j<nj+1; j++)
	for(i=-1; i<ni+1; i++)
	{
		RHO(i,j,k) = BUFp(i,j,k) + snd->rho0[k];
	}

}

#define QV BUFp

void calc_qv(buffers *b, grid gd, mesh msh, cmdline cmd,dir_meta dm,hdf_meta hm, sounding *snd, requested_cube rc)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float usr,vsr;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;
	float pi;
	float foo;
	int ifoo;

rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1;
rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;

	read_lofs_buffer(b->buf0,"qvpert",dm,hm,rc,cmd);
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk+1; k++)
	for(j=-1; j<nj+1; j++)
	for(i=-1; i<ni+1; i++)
	{
		QV(i,j,k) = BUFp(i,j,k) + 1000.0*snd->qv0[k];
	}

}

#define TEMPC BUFp

void calc_tempC(buffers *b, grid gd, mesh msh, cmdline cmd,dir_meta dm,hdf_meta hm, sounding *snd, requested_cube rc)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float usr,vsr;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;
	float pi;
	float foo;
	int ifoo;

rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1;
rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;
//	How Kelton does it (pointwise) in calcvort.c
//	for(k=0; k<nk+1; k++) {
//	for(j=-1; j<nj+1; j++) {
//	for(i=-1; i<ni+1; i++) {
 //   	calc_pipert(b->ppert, snd->pres0, b->buf0, i, j, k, ni, nj);
//	}
//	}
//	}

	read_lofs_buffer(b->buf0,"prespert",dm,hm,rc,cmd);
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk+1; k++)
	for(j=-1; j<nj+1; j++)
	for(i=-1; i<ni+1; i++)
	{
		TEM1p(i,j,k) = BUFp(i,j,k) + snd->pres0[k];
	}

	read_lofs_buffer(b->buf0,"thpert",dm,hm,rc,cmd);
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk+1; k++)
	for(j=-1; j<nj+1; j++)
	for(i=-1; i<ni+1; i++)
	{
		pi = pow((TEM1p(i,j,k)/100000.0),rovcp);
		BUFp(i,j,k) = (BUFp(i,j,k)+snd->th0[k])*pi;
		TEMPC(i,j,k) -= 273.15;
	}
}

#define QTOT BUFp
//qtot=qvpert+qc+qi+qs+qr+qg
void calc_qtot(buffers *b, grid gd, mesh msh, cmdline cmd,dir_meta dm,hdf_meta hm,requested_cube rc)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float usr,vsr;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	rc.X0=gd.X0-1; rc.Y0=gd.Y0-1; rc.Z0=gd.Z0;
	rc.X1=gd.X1+1; rc.Y1=gd.Y1+1; rc.Z1=gd.Z1;
	rc.NX=gd.X1-gd.X0+1; rc.NY=gd.Y1-gd.Y0+1; rc.NZ=gd.Z1-gd.Z0+1;

	read_lofs_buffer(b->buf0,"qvpert",dm,hm,rc,cmd);

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) = BUFp(i,j,k);

	read_lofs_buffer(b->buf0,"qc",dm,hm,rc,cmd);

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) += BUFp(i,j,k);

	read_lofs_buffer(b->buf0,"qi",dm,hm,rc,cmd);

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) += BUFp(i,j,k);

	read_lofs_buffer(b->buf0,"qr",dm,hm,rc,cmd);

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) += BUFp(i,j,k);

	read_lofs_buffer(b->buf0,"qs",dm,hm,rc,cmd);

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) += BUFp(i,j,k);

	read_lofs_buffer(b->buf0,"qg",dm,hm,rc,cmd);

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) += BUFp(i,j,k);

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		QTOT(i,j,k) = TEM1p(i,j,k);
}

/*******************************************************************************/

#define HWIN_SR BUFp
void calc_hwin_sr(buffers *b, grid gd, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float usr,vsr;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k,usr,vsr)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
	{
		usr = 0.5*(UAp(i,j,k)+UAp(i+1,j,k));
		vsr = 0.5*(VAp(i,j,k)+VAp(i,j+1,k));
		HWIN_SR(i,j,k) = sqrt(usr*usr+vsr*vsr);
	}
}

/*******************************************************************************/

#define HWIN_GR BUFp
void calc_hwin_gr(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float ugr,vgr;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k,ugr,vgr)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
	{
		ugr = 0.5*(UAp(i,j,k)+UAp(i+1,j,k)) + msh.umove;
		vgr = 0.5*(VAp(i,j,k)+VAp(i,j+1,k)) + msh.vmove;
		HWIN_GR(i,j,k) = sqrt(ugr*ugr+vgr*vgr);
	}
}

/*******************************************************************************/

#define U_GR BUFp
void calc_u_gr(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float ugr,vgr;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k,ugr)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
	{
		ugr = 0.5*(UAp(i,j,k)+UAp(i+1,j,k)) + msh.umove;
		U_GR(i,j,k) = ugr;
	}
}
/*******************************************************************************/

#define V_GR BUFp
void calc_v_gr(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float ugr,vgr;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k,ugr)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
	{
		vgr = 0.5*(VAp(i,j,k)+VAp(i+1,j,k)) + msh.vmove;
		V_GR(i,j,k) = vgr;
	}
}

/*******************************************************************************/

#define WINDMAG_SR BUFp
void calc_windmag_sr(buffers *b, grid gd, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float usr,vsr,w;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k,usr,vsr,w)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
	{
		usr = 0.5*(UAp(i,j,k)+UAp(i+1,j,k));
		vsr = 0.5*(VAp(i,j,k)+VAp(i,j+1,k));
		  w = 0.5*(WAp(i,j,k)+WAp(i,j,k+1));
		WINDMAG_SR(i,j,k) = sqrt(usr*usr+vsr*vsr+w*w);
	}
}

/*******************************************************************************/

#define HDIV BUFp
void calc_hdiv(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dudx,dvdy,rdx,rdy;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;
	rdx=msh.rdx; rdy=msh.rdy;

#pragma omp parallel for private(i,j,k,dudx,dvdy)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
	{
		dudx = (UAp(i+1,j,k)-UAp(i,j,k))*rdx*UH(i);
		dvdy = (VAp(i,j+1,k)-VAp(i,j,k))*rdy*VH(j);
		HDIV(i,j,k) = dudx + dvdy;
	}
}

/*******************************************************************************/

/*
 *typedef struct buffers
 *{
 *    float *ustag, *vstag, *wstag;
 *    float *buf0, *buf, *dum0, *dum1;
 *    float *threedbuf;
 *} buffers;
 */
#define XVORT BUFp
void do_xvort(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dy,dz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	#pragma omp parallel for private(i,j,k,dy,dz) 
	for(k=1; k<nk+1; k++) {
	for(j=0; j<nj+1; j++) {
	for(i=0; i<ni+1; i++) {
		dy = 1./(msh.rdy * VF(j)); 
		dz = 1./(msh.rdz * MF(k)); 
		calc_xvort(b->vstag, b->wstag, b->buf0, dy, dz, i, j, k, ni, nj);
	}
	}
	}

//This is dependent upon our current free slip bc, see CM1 for other decisions
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni+1; i++)
	{
		BUFp(i,j,0)=BUFp(i,j,1);
		BUFp(i,j,nk)=BUFp(i,j,nk-1);
	}

	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		XVORT(i,j,k) = 0.25 * (BUFp(i,j,k)+BUFp(i,j+1,k)+BUFp(i,j,k+1)+BUFp(i,j+1,k+1));
}

/*******************************************************************************/

#define YVORT BUFp
void do_yvort(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dx,dz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	#pragma omp parallel for private(i,j,k,dx,dz) 
	for(k=1; k<nk+1; k++) {
	for(j=0; j<nj+1; j++) {
	for(i=0; i<ni+1; i++) {
		dz = 1./(msh.rdz * MF(k)); 
		dx = 1./(msh.rdx * UF(i)); 
		calc_yvort(b->ustag, b->wstag, b->buf0, dx, dz, i, j, k, ni, nj);
	}
	}
	}

//This is dependent upon our current free slip bc, see CM1 for other decisions
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni+1; i++)
	{
		BUFp(i,j,0)=BUFp(i,j,1);
		BUFp(i,j,nk)=BUFp(i,j,nk-1);
	}

	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		YVORT(i,j,k) = 0.25 * (YVORT(i,j,k)+YVORT(i+1,j,k)+YVORT(i,j,k+1)+YVORT(i+1,j,k+1));
}

/*******************************************************************************/

#define ZVORT BUFp
void do_zvort(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dx,dy;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	#pragma omp parallel for private(i,j,k,dx,dy) 
	for(k=0; k<nk+1; k++) {
	for(j=0; j<nj+1; j++) {
	for(i=0; i<ni+1; i++) {
		dy = 1./(msh.rdy * VF(j));
		dx = 1./(msh.rdx * UF(i)); 
    	calc_zvort(b->ustag, b->vstag, b->buf0, dx, dy, i, j, k, ni, nj);
	}
	}
	}

	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		ZVORT(i,j,k) = 0.25 * (ZVORT(i,j,k)+ZVORT(i+1,j,k)+ZVORT(i,j+1,k)+ZVORT(i+1,j+1,k));
}

/*******************************************************************************/

#define XVSTRETCH BUFp
void do_xvort_stretch(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dy, dz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	// calculate xvort and average it to the scalar mesh
    do_xvort(b, gd, msh, cmd);
	// xvort will be in the buf0 array, and we want/need
	// to copy it to dum0 so that the stretching rate
	// can be stored in buf0 for writing.
	long size = (nk+1)*(nj+2)*(ni+2);
	for (i=0; i<size; i++) {
		b->dum0[i] = b->buf0[i];
	}

	#pragma omp parallel for private(i,j,k,dy,dz) 
	for (k=0; k<nk; k++) {
	for (j=0; j<nj; j++) {
	for (i=0; i<ni; i++) {
		// This routine uses a centerd difference - so dy 
		// needs to be set appropriately for the increased distance
		//
		dy = 1./(msh.rdy * VF(j));
		dz = 1./(msh.rdz * MF(k));
		calc_xvort_stretch(b->vstag, b->wstag, b->dum0, b->buf0, dy, dz, i, j, k, nx, ny);
	}
	}
	}

}

/*******************************************************************************/

#define YVSTRETCH BUFp
void do_yvort_stretch(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dx, dz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	// calculate yvort and average it to the scalar mesh
    do_yvort(b, gd, msh, cmd);
	// yvort will be in the buf0 array, and we want/need
	// to copy it to dum0 so that the stretching rate
	// can be stored in buf0 for writing.
	long size = (nk+1)*(nj+2)*(ni+2);
	for (i=0; i<size; i++) {
		b->dum0[i] = b->buf0[i];
	}

	#pragma omp parallel for private(i,j,k,dx,dz) 
	for (k=0; k<nk; k++) {
	for (j=0; j<nj; j++) {
	for (i=0; i<ni; i++) {
		// This routine uses a centerd difference - so dy 
		// needs to be set appropriately for the increased distance
		//
		dx = 1./(msh.rdx * UF(i));
		dz = 1./(msh.rdz * MF(k));
		calc_yvort_stretch(b->ustag, b->wstag, b->dum0, b->buf0, dx, dz, i, j, k, nx, ny);
	}
	}
	}

}

/*******************************************************************************/

#define ZVSTRETCH BUFp
void do_zvort_stretch(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dx, dy;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	// calculate zvort
    do_zvort(b, gd, msh, cmd);
	// zvort will be in the buf0 array, and we want/need
	// to copy it to dum0 so that the stretching rate
	// can be stored in buf0 for writing.
	long size = (nk+1)*(nj+2)*(ni+2);
	for (i=0; i<size; i++) {
		b->dum0[i] = b->buf0[i];
	}

	#pragma omp parallel for private(i,j,k,dx,dy) 
	for (k=0; k<nk; k++) {
	for (j=0; j<nj; j++) {
	for (i=0; i<ni; i++) {
		// This routine uses a centerd difference - so dy 
		// needs to be set appropriately for the increased distance
		//
		dx = 1./(msh.rdx * UF(i));
		dy = 1./(msh.rdy * VF(j));
		calc_zvort_stretch(b->ustag, b->vstag, b->dum0, b->buf0, dx, dy, i, j, k, nx, ny);
	}
	}
	}

}
/*******************************************************************************/

//ORF 2023-06-08
//We do not use Kelton's code for tilting terms because he pre-allocates dwdx and dwdy as full
//3D arrays - these are temporary and can be calculated locally here

#define ZVTILT BUFp
void do_zvort_tilt(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dx, dy;
	float dwdx,dwdy,dwdxc[4],dwdyc[4],xvort,yvort,zvort_tilt;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	long size = (nk+1)*(nj+2)*(ni+2);

    do_xvort(b, gd, msh, cmd);
	// xvort will be in the buf0 array, and we want to copy it to dum0
#pragma omp parallel for private(i) 
	for (i=0; i<size; i++) {
		b->dum0[i] = b->buf0[i];
	}
    do_yvort(b, gd, msh, cmd);
	// yvort will be in the buf0 array, and we want to copy it to dum1
#pragma omp parallel for private(i) 
	for (i=0; i<size; i++) {
		b->dum1[i] = b->buf0[i];
	}

#pragma omp parallel for private(i,j,k,dx,dy,dwdxc,dwdyc,dwdx,dwdy,xvort,yvort,zvort_tilt) 
	for (k=0; k<nk; k++) {
	for (j=0; j<nj; j++) {
	for (i=0; i<ni; i++) {
		dx = 1./(msh.rdx * UF(i));
		dy = 1./(msh.rdy * VF(j));
		dwdxc[0] = ( WAp(i, j, k) - WAp(i-1, j, k) ) / dx;
		dwdxc[1] = ( WAp(i+1, j, k) - WAp(i, j, k) ) / dx;
		dwdxc[2] = ( WAp(i, j+1, k) - WAp(i-1, j+1, k) ) / dx;
		dwdxc[3] = ( WAp(i+1, j+1, k) - WAp(i, j+1, k) ) / dx;
		dwdyc[0] = ( WAp(i, j, k) - WAp(i, j-1, k) ) / dy;
		dwdyc[1] = ( WAp(i, j+1, k) - WAp(i, j, k) ) / dy;
		dwdyc[2] = ( WAp(i, j, k+1) - WAp(i, j-1, k+1) ) / dy;
		dwdyc[3] = ( WAp(i, j+1, k+1) - WAp(i, j, k+1) ) / dy;
		dwdx = 0.25*(dwdxc[0]+dwdxc[1]+dwdxc[2]+dwdxc[3]);
		dwdy = 0.25*(dwdyc[0]+dwdyc[1]+dwdyc[2]+dwdyc[3]);
		xvort=TEMp(i,j,k);
		yvort=TEM1p(i,j,k);
		zvort_tilt=xvort*dwdx+yvort*dwdy;
		ZVTILT(i,j,k) = zvort_tilt;
	}
	}
	}
}

/*******************************************************************************/

#define XVTILT BUFp
void do_xvort_tilt(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dy, dz;
	float dudy,dudz,dudyc[4],dudzc[4],yvort,zvort,xvort_tilt;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	long size = (nk+1)*(nj+2)*(ni+2);

    do_zvort(b, gd, msh, cmd);
	// zvort will be in the buf0 array, and we want to copy it to dum0
#pragma omp parallel for private(i) 
	for (i=0; i<size; i++) {
		b->dum0[i] = b->buf0[i];
	}
    do_yvort(b, gd, msh, cmd);
	// yvort will be in the buf0 array, and we want to copy it to dum1
#pragma omp parallel for private(i) 
	for (i=0; i<size; i++) {
		b->dum1[i] = b->buf0[i];
	}

#pragma omp parallel for private(i,j) 
	for (j=0; j<nj; j++)
	for (i=0; i<ni; i++)
		XVTILT(i,j,0) = 0.0;

#pragma omp parallel for private(i,j,k,dy,dz,dudyc,dudzc,dudy,dudz,yvort,zvort,xvort_tilt) 
	for (k=1; k<nk; k++) {
	for (j=0; j<nj; j++) {
	for (i=0; i<ni; i++) {
		dy = 1./(msh.rdy * VF(j));
		dz = 1./(msh.rdz * MF(k));
		dudyc[0] = ( UAp(i, j, k) - UAp(i, j-1, k) ) / dy;
		dudyc[1] = ( UAp(i+1, j, k) - UAp(i+1, j-1, k) ) / dy;
		dudyc[2] = ( UAp(i, j+1, k) - UAp(i, j, k) ) / dy;
		dudyc[3] = ( UAp(i+1, j+1, k) - UAp(i+1, j, k) ) / dy;
		dudzc[0] = ( UAp(i, j, k) - UAp(i, j, k-1) ) / dz;
		dudzc[1] = ( UAp(i+1, j, k) - UAp(i+1, j, k-1) ) / dz;
		dudzc[2] = ( UAp(i, j, k+1) - UAp(i, j, k) ) / dz;
		dudzc[0] = ( UAp(i+1, j, k+1) - UAp(i+1, j, k) ) / dz;
		dudy = 0.25*(dudyc[0]+dudyc[1]+dudyc[2]+dudyc[3]);
		dudz = 0.25*(dudzc[0]+dudzc[1]+dudzc[2]+dudzc[3]);
		zvort=TEMp(i,j,k);
		yvort=TEM1p(i,j,k);
		xvort_tilt=zvort*dudz+yvort*dudy;
		XVTILT(i,j,k) = xvort_tilt;
	}
	}
	}
}

/*******************************************************************************/

#define YVTILT BUFp
void do_yvort_tilt(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dz, dx;
	float dvdx,dvdz,dvdxc[4],dvdzc[4],zvort,xvort,yvort_tilt;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	long size = (nk+1)*(nj+2)*(ni+2);

    do_zvort(b, gd, msh, cmd);
	// zvort will be in the buf0 array, and we want to copy it to dum0
#pragma omp parallel for private(i) 
	for (i=0; i<size; i++) {
		b->dum0[i] = b->buf0[i];
	}
    do_xvort(b, gd, msh, cmd);
	// xvort will be in the buf0 array, and we want to copy it to dum1
#pragma omp parallel for private(i) 
	for (i=0; i<size; i++) {
		b->dum1[i] = b->buf0[i];
	}

#pragma omp parallel for private(i,j) 
	for (j=0; j<nj; j++)
	for (i=0; i<ni; i++)
		YVTILT(i,j,0) = 0.0;

#pragma omp parallel for private(i,j,k,dz,dx,dvdzc,dvdxc,dvdz,dvdx,xvort,zvort,yvort_tilt) 
	for (k=1; k<nk; k++) {
	for (j=0; j<nj; j++) {
	for (i=0; i<ni; i++) {
		dx = 1./(msh.rdx * UF(i));
		dz = 1./(msh.rdz * MF(k));
		dvdxc[0] = ( VAp(i, j, k) - VAp(i-1, j, k) ) / dx;
		dvdxc[1] = ( VAp(i+1, j, k) - VAp(i, j, k) ) / dx;
		dvdxc[2] = ( VAp(i, j+1, k) - VAp(i-1, j+1, k) ) / dx;
		dvdxc[3] = ( VAp(i+1, j+1, k) - VAp(i, j+1, k) ) / dx;
		dvdzc[0] = ( VAp(i, j, k) - VAp(i, j, k-1) ) / dz;
		dvdzc[1] = ( VAp(i, j+1, k) - VAp(i, j+1, k-1) ) / dz;
		dvdzc[2] = ( VAp(i, j, k+1) - VAp(i, j, k) ) / dz;
		dvdzc[3] = ( VAp(i, j+1, k+1) - VAp(i, j+1, k) ) / dz;
		dvdx = 0.25*(dvdxc[0]+dvdxc[1]+dvdxc[2]+dvdxc[3]);
		dvdz = 0.25*(dvdzc[0]+dvdzc[1]+dvdzc[2]+dvdzc[3]);
		zvort=TEMp(i,j,k);
		xvort=TEM1p(i,j,k);
		yvort_tilt=zvort*dvdz+xvort*dvdx;
		YVTILT(i,j,k) = yvort_tilt;
	}
	}
	}
}

/*******************************************************************************/

#define XVBARO BUFp
void do_xvort_baro(buffers *b, grid gd, sounding *snd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dy;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	#pragma omp parallel for private(i,j,k,dy) 
	for (k=0; k<nk; k++) {
	for (j=0; j<nj; j++) {
	for (i=0; i<ni; i++) {
		// This routine uses a centerd difference - so dy 
		// needs to be set appropriately for the increased distance
		dy = 1./(msh.rdy * VH(j-1)) + 1./(msh.rdy * VH(j+1));
		calc_xvort_baro(b->thrhopert, snd->th0, snd->qv0, b->buf0, dy, i, j, k, nx, ny);
	}
	}
	}

}

/*******************************************************************************/

#define YVBARO BUFp
void do_yvort_baro(buffers *b, grid gd, sounding *snd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dx;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	#pragma omp parallel for private(i,j,k,dx) 
	for (k=0; k<nk; k++) {
	for (j=0; j<nj; j++) {
	for (i=0; i<ni; i++) {
		// This routine uses a centerd difference - so dy 
		// needs to be set appropriately for the increased distance
		dx = 1./(msh.rdx * UH(i-1)) + 1./(msh.rdx * UH(i+1));
		calc_yvort_baro(b->thrhopert, snd->th0, snd->qv0, b->buf0, dx, i, j, k, nx, ny);
	}
	}
	}

}

/*******************************************************************************/

#define XVSOL BUFp
void do_xvort_solenoid(buffers *b, grid gd, sounding *snd, mesh msh, cmdline cmd) 
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dy, dz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	// need to calculate pipert first
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
    	calc_pipert(b->ppert, snd->pres0, b->dum0, i, j, k, ni, nj);
	}
	}
	}

	#pragma omp parallel for private(i,j,k,dy,dz) 
	for (k=1; k<nk; k++) {
	for (j=0; j<nj; j++) {
	for (i=0; i<ni; i++) {
		// This routine uses a centerd difference - so dy 
		// needs to be set appropriately for the increased distance
		dy = 1./(msh.rdy * VH(j-1)) + 1./(msh.rdy * VH(j+1));
		dz = 1./(msh.rdz * MH(k-1)) + 1./(msh.rdy * MH(k+1));
		calc_xvort_solenoid(b->dum0, b->thrhopert, snd->th0, snd->qv0, b->buf0, dy, dz, i, j, k, nx, ny);
	}
	}
	}

	// lower boundary condition
	#pragma omp parallel for private(i, j)
	for (j=-1; j<nj+1; j++) {
	for (i=-1; i<ni+1; i++) {
		XVSOL(i, j, 0) = 0.0;
	}
	}

}
/*******************************************************************************/

#define YVSOL BUFp
void do_yvort_solenoid(buffers *b, grid gd, sounding *snd, mesh msh, cmdline cmd) 
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dx, dz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	// need to calculate pipert first
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
    	calc_pipert(b->ppert, snd->pres0, b->dum0, i, j, k, ni, nj);
	}
	}
	}

	#pragma omp parallel for private(i,j,k,dx,dz) 
	for (k=1; k<nk; k++) {
	for (j=0; j<nj; j++) {
	for (i=0; i<ni; i++) {
		// This routine uses a centerd difference - so dy 
		// needs to be set appropriately for the increased distance
		dx = 1./(msh.rdx * UH(i-1)) + 1./(msh.rdx * UH(i+1));
		dz = 1./(msh.rdz * MH(k-1)) + 1./(msh.rdz * MH(k+1));
		calc_yvort_solenoid(b->dum0, b->thrhopert, snd->th0, snd->qv0, b->buf0, dx, dz, i, j, k, nx, ny);
	}
	}
	}

	// lower boundary condition
	#pragma omp parallel for private(i, j)
	for (j=-1; j<nj+1; j++) {
	for (i=-1; i<ni+1; i++) {
		YVSOL(i, j, 0) = 0.0;
	}
	}

}

/*******************************************************************************/

#define ZVSOL BUFp
void do_zvort_solenoid(buffers *b, grid gd, sounding *snd, mesh msh, cmdline cmd) 
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dx, dy;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	// need to calculate pipert first
	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk+1; k++) {
	for(j=-1; j<nj+1; j++) {
	for(i=-1; i<ni+1; i++) {
    	calc_pipert(b->ppert, snd->pres0, b->dum0, i, j, k, ni, nj);
	}
	}
	}

	#pragma omp parallel for private(i,j,k,dx,dy) 
	for (k=0; k<nk+1; k++) {
	for (j=0; j<nj; j++) {
	for (i=0; i<ni; i++) {
		// This routine uses a centerd difference - so dy 
		// needs to be set appropriately for the increased distance
		dx = 1./(msh.rdx * UH(i-1)) + 1./(msh.rdx * UH(i+1));
		dy = 1./(msh.rdy * VH(j-1)) + 1./(msh.rdy * VH(j+1));
		calc_zvort_solenoid(b->dum0, b->thrhopert, b->buf0, dx, dy, i, j, k, nx, ny);
	}
	}
	}

}

/*******************************************************************************/

void calc_hvort(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dx,dy,dz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	#pragma omp parallel for private(i,j,k) 
	for(k=1; k<nk; k++) {
		dz = 1./(msh.rdz * MF(k)); 
		// If we parallelize over j and give each thread a row to
		// work on, then we take advantage of caching and memory
		// linearity. 
		for(j=0; j<nj+1; j++) {
			dy = 1./(msh.rdy * VF(j)); 
			for(i=0; i<ni; i++) {
    			calc_xvort(b->vstag, b->wstag, b->dum0, dy, dz, i, j, k, ni, nj);
			}
		}
	}
//This is dependent upon our current free slip bc, see CM1 for other decisions
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni; i++)
	{
		TEMp(i,j,0)=TEMp(i,j,1);
		TEMp(i,j,nk)=TEMp(i,j,nk-1);
	}

#define XVORT BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		XVORT(i,j,k) = 0.25 * (TEMp(i,j,k)+TEMp(i,j+1,k)+TEMp(i,j,k+1)+TEMp(i,j+1,k+1));


#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) = XVORT(i,j,k)*XVORT(i,j,k);

#define YVORT BUFp
	#pragma omp parallel for private(i,j,k) 
	for(k=1; k<nk; k++) {
		dz = 1./(msh.rdz * MF(k)); 
		for(j=0; j<nj; j++) {
			for(i=0; i<ni+1; i++) {
				dx = 1./(msh.rdx * UF(i)); 
    			calc_yvort(b->ustag, b->wstag, b->dum0, dx, dz, i, j, k, ni, nj);
			}
		}
	}
//This is dependent upon our current free slip bc, see CM1 for other
//decisions
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		TEMp(i,j,0)=TEMp(i,j,1);
		TEMp(i,j,nk)=TEMp(i,j,nk-1);
	}
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		YVORT(i,j,k) = 0.25 * (TEMp(i,j,k)+TEMp(i+1,j,k)+TEMp(i,j,k+1)+TEMp(i+1,j,k+1));

	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) += YVORT(i,j,k)*YVORT(i,j,k);

#define HVORT BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		HVORT(i,j,k) = sqrt(TEM1p(i,j,k));
}

/*******************************************************************************/

void calc_vortmag(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dx,dy,dz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	#pragma omp parallel for private(i,j,k) 
	for(k=1; k<nk; k++) {
		dz = 1./(msh.rdz * MF(k)); 
		// If we parallelize over j and give each thread a row to
		// work on, then we take advantage of caching and memory
		// linearity. 
		for(j=0; j<nj+1; j++) {
			dy = 1./(msh.rdy * VF(j)); 
			for(i=0; i<ni; i++) {
    			calc_xvort(b->vstag, b->wstag, b->buf0, dy, dz, i, j, k, ni, nj);
			}
		}
	}
//This is dependent upon our current free slip bc, see CM1 for other decisions
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni; i++)
	{
		BUFp(i,j,0)=BUFp(i,j,1);
		BUFp(i,j,nk)=BUFp(i,j,nk-1);
	}

#define XVORT BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		XVORT(i,j,k) = 0.25 * (BUFp(i,j,k)+BUFp(i,j+1,k)+BUFp(i,j,k+1)+BUFp(i,j+1,k+1));


#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) = XVORT(i,j,k)*XVORT(i,j,k);

	#pragma omp parallel for private(i,j,k) 
	for(k=1; k<nk; k++) {
		dz = 1./(msh.rdz * MF(k)); 
		for(j=0; j<nj; j++) {
			for(i=0; i<ni+1; i++) {
				dx = 1./(msh.rdx * UF(i)); 
    			calc_yvort(b->ustag, b->wstag, b->buf0, dx, dz, i, j, k, ni, nj);
			}
		}
	}
//This is dependent upon our current free slip bc, see CM1 for other
//decisions
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		BUFp(i,j,0)=BUFp(i,j,1);
		BUFp(i,j,nk)=BUFp(i,j,nk-1);
	}
#define YVORT BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		YVORT(i,j,k) = 0.25 * (BUFp(i,j,k)+BUFp(i+1,j,k)+BUFp(i,j,k+1)+BUFp(i+1,j,k+1));

	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) += YVORT(i,j,k)*YVORT(i,j,k);

	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk; k++) {
		for(j=0; j<nj+1; j++) {
			dy = 1./(msh.rdy * VF(j)); 
			for(i=0; i<ni+1; i++) {
				dx = 1./(msh.rdx * UF(i)); 
    			calc_zvort(b->ustag, b->vstag, b->buf0, dx, dy, i, j, k, ni, nj);
			}
		}
	}
#define ZVORT BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		ZVORT(i,j,k) = 0.25 * (BUFp(i,j,k)+BUFp(i+1,j,k)+BUFp(i,j+1,k)+BUFp(i+1,j+1,k));

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) += ZVORT(i,j,k)*ZVORT(i,j,k);

#define VORTMAG BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		VORTMAG(i,j,k) = sqrt(TEM1p(i,j,k));
}

/*******************************************************************************/

#define STREAMVORT BUFp
void calc_streamvort(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dx,dy,dz;
	float uinterp,vinterp,winterp;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	#pragma omp parallel for private(i,j,k) 
	for(k=1; k<nk; k++) {
		dz = 1./(msh.rdz * MF(k)); 
		// If we parallelize over j and give each thread a row to
		// work on, then we take advantage of caching and memory
		// linearity. 
		for(j=0; j<nj+1; j++) {
			dy = 1./(msh.rdy * VF(j)); 
			for(i=0; i<ni; i++) {
    			calc_xvort(b->vstag, b->wstag, b->buf0, dy, dz, i, j, k, ni, nj);
			}
		}
	}
//This is dependent upon our current free slip bc, see CM1 for other decisions
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni; i++)
	{
		BUFp(i,j,0)=BUFp(i,j,1);
		BUFp(i,j,nk)=BUFp(i,j,nk-1);
	}

#define XVORT BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		XVORT(i,j,k) = 0.25 * (XVORT(i,j,k)+XVORT(i,j+1,k)+XVORT(i,j,k+1)+XVORT(i,j+1,k+1));


#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) =  XVORT(i,j,k)*0.5*(UAp(i,j,k)+UAp(i+1,j,k));

	#pragma omp parallel for private(i,j,k) 
	for(k=1; k<nk; k++) {
		dz = 1./(msh.rdz * MF(k)); 
		for(j=0; j<nj; j++) {
			for(i=0; i<ni+1; i++) {
				dx = 1./(msh.rdx * UF(i)); 
    			calc_yvort(b->ustag, b->wstag, b->buf0, dx, dz, i, j, k, ni, nj);
			}
		}
	}
//This is dependent upon our current free slip bc, see CM1 for other
//decisions
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		BUFp(i,j,0)=BUFp(i,j,1);
		BUFp(i,j,nk)=BUFp(i,j,nk-1);
	}
#define YVORT BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		YVORT(i,j,k) = 0.25 * (YVORT(i,j,k)+YVORT(i+1,j,k)+YVORT(i,j,k+1)+YVORT(i+1,j,k+1));

	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) += YVORT(i,j,k)*0.5*(VAp(i,j,k)+VAp(i,j+1,k));

	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk; k++) {
		for(j=0; j<nj+1; j++) {
			dy = 1./(msh.rdy * VF(j)); 
			for(i=0; i<ni+1; i++) {
				dx = 1./(msh.rdx * UF(i)); 
    			calc_zvort(b->ustag, b->vstag, b->buf0, dx, dy, i, j, k, ni, nj);
			}
		}
	}
#define ZVORT BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		ZVORT(i,j,k) = 0.25 * (ZVORT(i,j,k)+ZVORT(i+1,j,k)+ZVORT(i,j+1,k)+ZVORT(i+1,j+1,k));

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1p(i,j,k) += ZVORT(i,j,k)*0.5*(WAp(i,j,k)+WAp(i,j,k+1));

#define STREAMVORT BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
	{
		uinterp=0.5*(UAp(i,j,k)+UAp(i+1,j,k));
		vinterp=0.5*(VAp(i,j,k)+VAp(i,j+1,k));
		winterp=0.5*(WAp(i,j,k)+WAp(i,j,k+1));
		STREAMVORT(i,j,k) = TEM1p(i,j,k)/(sqrt(uinterp*uinterp+vinterp*vinterp+winterp*winterp));
	}
}

/*******************************************************************************/

void z_progress_bar(int iz, int nz)
{
	int index;
	char *alph="abcdefghijklmnopqrstuvwxyz";
	index=(int)((float)iz*26.0/(float)nz);
	printf("%c",alph[index]); fflush (stdout);
}

/*******************************************************************************/

//ORF TODO: since we only read from snd, don't pass a pointer (in
//keeping with my approach of only passing a pointer when we're going to
//change information - if we are just going to read it pass values)

void do_requested_variables(buffers *b, ncstruct nc, grid gd, mesh msh, sounding *snd, readahead rh,dir_meta dm,hdf_meta hm,cmdline cmd)
{
	int ix,iy,iz,nx,ny,buf0nx,buf0ny,i,ixoff,iyoff,ivar,status;
	long int bufsize;
	requested_cube rc;
	char *var;
	float *twodbuf;
	size_t writestart[4],writeedges[4];

// For flexibility we always set rc in case we need to read outside of what we
// requested at the command line such as what occurs with staggered variables. The
// copy_grid_to_requested_cube just duplicates the grid structure stuff to rc (requested
// cube). rc is a required argument to read_lofs_buffer.

// For doing calculations that involve derivatives (usually of u, v, or w) we have to read
// data in slightly differently (see for instance do_readahead where we read in ustag
// vstag wstag on a bigger mesh than the scalar variables).

	copy_grid_to_requested_cube(&rc,gd);

	var = (char *) malloc (MAXSTR * sizeof(char));
	if(cmd.twodwrite)
	{
		bufsize=gd.NX*gd.NY*sizeof(float);
		twodbuf = (float *)malloc(bufsize);
	}
//	else
//	{
//		bufsize=(long)gd.NX*(long)gd.NY*(long)gd.NZ*(long)sizeof(float);
//		threedbuf = (float *)malloc(bufsize*sizeof(float));
//	}

	for (ivar = 0; ivar < cmd.nvar; ivar++)
	{
		buf0nx=gd.NX+2;ixoff=1;
		buf0ny=gd.NY+2;iyoff=1;

		var=nc.var3d[ivar].varname;
		printf("%s: ",var);FL;

		if(same(var,"u"))
		{
			if(!rh.u)
			{
				buf0nx=gd.NX;ixoff=0;
				buf0ny=gd.NY;iyoff=0;
				read_lofs_buffer(b->buf,var,dm,hm,rc,cmd);
			}
			else
			{
				buf_u(b,gd); 
			}
		}
		else if(same(var,"v"))
		{
			if(!rh.v)
			{
				buf0nx=gd.NX;ixoff=0;
				buf0ny=gd.NY;iyoff=0;
				read_lofs_buffer(b->buf,var,dm,hm,rc,cmd);
			}
			else
			{
				buf_v(b,gd);
			}
		}
		else if(same(var,"w"))
		{
			if(!rh.w)
			{
				buf0nx=gd.NX;ixoff=0;
				buf0ny=gd.NY;iyoff=0;
				read_lofs_buffer(b->buf,var,dm,hm,rc,cmd);
			}
			else
			{
				buf_w(b,gd);
			}
		}
		else if(same(var,"prespert"))
		{
			if(!rh.ppert)
			{
				buf0nx=gd.NX;ixoff=0;
				buf0ny=gd.NY;iyoff=0;
				read_lofs_buffer(b->buf,var,dm,hm,rc,cmd);
			}
			else
			{
				buf_prespert(b,gd);
			}
		}
		else if(same(var,"thrhopert"))
		{
			if(!rh.thrhopert)
			{
				buf0nx=gd.NX;ixoff=0;
				buf0ny=gd.NY;iyoff=0;
				read_lofs_buffer(b->buf,var,dm,hm,rc,cmd);
			}
			else
			{
				buf_thrhopert(b,gd);
			}
		}
		else if(same(var,"pipert"))        {CL;do_pipert(b,gd,snd,cmd);}
		else if(same(var,"wb_buoy"))       {CL;do_wbuoy(b,gd,snd,cmd);}
		else if(same(var,"wb_buoy_interp"))  {CL;do_wbuoy_interp(b,gd,snd,cmd);}
		else if(same(var,"ub_pgrad"))      {CL;do_upgrad(b,gd,snd,msh,cmd);}
		else if(same(var,"vb_pgrad"))      {CL;do_vpgrad(b,gd,snd,msh,cmd);}
		else if(same(var,"wb_pgrad"))      {CL;do_wpgrad(b,gd,snd,msh,cmd);}
		else if(same(var,"ub_pgrad_interp"))      {CL;do_upgrad_interp(b,gd,snd,msh,cmd);}
		else if(same(var,"vb_pgrad_interp"))      {CL;do_vpgrad_interp(b,gd,snd,msh,cmd);}
		else if(same(var,"wb_pgrad_interp"))      {CL;do_wpgrad_interp(b,gd,snd,msh,cmd);}
		else if(same(var,"uinterp"))	   {CL;calc_uinterp(b,gd,msh,cmd,dm,hm,rc,rh);}
		else if(same(var,"vinterp"))	   {CL;calc_vinterp(b,gd,msh,cmd,dm,hm,rc,rh);}
		else if(same(var,"winterp"))	   {CL;calc_winterp(b,gd,msh,cmd,dm,hm,rc,rh);}
		else if(same(var,"hwin_sr"))	   {CL;calc_hwin_sr(b,gd,cmd);}
		else if(same(var,"hwin_gr"))	   {CL;calc_hwin_gr(b,gd,msh,cmd);}
		else if(same(var,"windmag_sr"))	   {CL;calc_windmag_sr(b,gd,cmd);}
		else if(same(var,"u_gr"))		   {CL;calc_u_gr(b,gd,msh,cmd);}
		else if(same(var,"v_gr"))		   {CL;calc_v_gr(b,gd,msh,cmd);}
		else if(same(var,"hdiv")) 		   {CL;calc_hdiv(b,gd,msh,cmd);}
		else if(same(var,"xvort"))		   {CL;do_xvort(b,gd,msh,cmd);}
		else if(same(var,"yvort"))		   {CL;do_yvort(b,gd,msh,cmd);}
		else if(same(var,"zvort"))		   {CL;do_zvort(b,gd,msh,cmd);}
		else if(same(var,"xvort_stretch")) {CL;do_xvort_stretch(b,gd,msh,cmd);} 
		else if(same(var,"xvort_tilt")) {CL;do_xvort_tilt(b,gd,msh,cmd);} 
	    else if(same(var,"xvort_baro"))    {CL;do_xvort_baro(b,gd,snd,msh,cmd);}
		else if(same(var,"yvort_stretch")) {CL;do_yvort_stretch(b,gd,msh,cmd);} 
		else if(same(var,"yvort_tilt")) {CL;do_yvort_tilt(b,gd,msh,cmd);} 
	    else if(same(var,"yvort_baro"))    {CL;do_yvort_baro(b,gd,snd,msh,cmd);}
		else if(same(var,"zvort_stretch")) {CL;do_zvort_stretch(b,gd,msh,cmd);} 
		else if(same(var,"zvort_tilt")) {CL;do_zvort_tilt(b,gd,msh,cmd);} 
	    else if(same(var,"xvort_solenoid")){CL;do_xvort_solenoid(b,gd,snd,msh,cmd);}
	    else if(same(var,"yvort_solenoid")){CL;do_yvort_solenoid(b,gd,snd,msh,cmd);}
	    else if(same(var,"zvort_solenoid")){CL;do_zvort_solenoid(b,gd,snd,msh,cmd);}
		else if(same(var,"hvort"))		   {CL;calc_hvort(b,gd,msh,cmd);}
		else if(same(var,"vortmag"))	   {CL;calc_vortmag(b,gd,msh,cmd);}
		else if(same(var,"streamvort"))	   {CL;calc_streamvort(b,gd,msh,cmd);}
		else if(same(var,"qiqvpert"))	   {CL;calc_qiqvpert(b,gd,msh,cmd,dm,hm,rc);}
		else if(same(var,"qtot"))	   {CL;calc_qtot(b,gd,msh,cmd,dm,hm,rc);}
		else if(same(var,"tempC"))	   {CL;calc_tempC(b,gd,msh,cmd,dm,hm,snd,rc);}
		else if(same(var,"rho"))	   {CL;calc_rho(b,gd,msh,cmd,dm,hm,snd,rc);}
		else if(same(var,"qv"))	   {CL;calc_qv(b,gd,msh,cmd,dm,hm,snd,rc);}
		else if(same(var,"kmh_interp"))	   {CL;calc_kmh_interp(b,gd,msh,cmd,dm,hm,rc);}
		else if(same(var,"kmv_interp"))	   {CL;calc_kmv_interp(b,gd,msh,cmd,dm,hm,rc);}
		else if(same(var,"khh_interp"))	   {CL;calc_khh_interp(b,gd,msh,cmd,dm,hm,rc);}
		else if(same(var,"khv_interp"))	   {CL;calc_khv_interp(b,gd,msh,cmd,dm,hm,rc);}
//void do_requested_variables(buffers *b, ncstruct nc, grid gd, mesh msh, sounding *snd, readahead rh,dir_meta dm,hdf_meta hm,cmdline cmd)
//			read_lofs_buffer(b->buf,nc.var3d[ivar].varname,dm,hm,rc,cmd);
		else
		{
//		buf0nx=gd.NX+2;ixoff=1;
//		buf0ny=gd.NY+2;iyoff=1;
//						b->threedbuf[P3(ix,iy,iz,gd.NX,gd.NY)] = b->buf[P3(ix+ixoff,iy+iyoff,iz,buf0nx,buf0ny)];
			printf("reading...");FL;
			//			!!!!ORF!!!! just commented out ... is set above
			//			differently???
			buf0nx=gd.NX;ixoff=0;
			buf0ny=gd.NY;iyoff=0;
			read_lofs_buffer(b->buf,nc.var3d[ivar].varname,dm,hm,rc,cmd);
		}
		printf("writing...");FL;

// ORF I tried this and it kind of sucked performance wise for large-ish
// data, but I leave the option available. This method writes our data
// in 2D XY slices rather than one big 3D chunk (saves a bit of memory).
// Newsflash: 1 big 3D chunk writes a shit-ton faster.

		if(cmd.twodwrite)
		{
			for (iz=0; iz<gd.NZ; iz++)
			{
//				if(cmd.verbose&&cmd.gzip) z_progress_bar(iz,gd.NZ);
#pragma omp parallel for private(iy,ix)
				for(iy=0;iy<gd.NY;iy++)
				{
					for(ix=0;ix<gd.NX;ix++)
					{
						twodbuf[P2(ix,iy,gd.NX)] = b->buf[P3(ix+ixoff,iy+iyoff,iz,buf0nx,buf0ny)];
					}
				}
				if(gd.X0==gd.X1)      //YZ slice
				{
					writestart[0]=0;  writeedges[0]=1;     //time
					writestart[1]=iz; writeedges[1]=1;     //z
					writestart[2]=0;  writeedges[2]=gd.NY; //y
				}
				else if(gd.Y0==gd.Y1) //XZ slice
				{
					writestart[0]=0;  writeedges[0]=1;     //time
					writestart[1]=iz; writeedges[1]=1;     //z
					writestart[2]=0;  writeedges[2]=gd.NX; //x
				}
				else if(gd.Z0==gd.Z1) //XY slice
				{
					writestart[0]=0;  writeedges[0]=1;     //time
					writestart[1]=0;  writeedges[1]=gd.NY; //y
					writestart[2]=0;  writeedges[2]=gd.NX; //x
				}
				else                  //XYZ rectangular prism (cubey thing)
				{
					writestart[0]=0;  writeedges[0]=1;     //time
					writestart[1]=iz; writeedges[1]=1;     //z
					writestart[2]=0;  writeedges[2]=gd.NY; //y
					writestart[3]=0;  writeedges[3]=gd.NX; //x
				}
				status = nc_put_vara_float (nc.ncid, nc.var3d[ivar].varnameid, writestart, writeedges, twodbuf);
			}
		}
		else//This is the default
		{
#pragma omp parallel for private(ix,iy,iz)
			for (iz=0; iz<gd.NZ; iz++)
			{
				for(iy=0;iy<gd.NY;iy++)
				{
					for(ix=0;ix<gd.NX;ix++)
					{
						b->threedbuf[P3(ix,iy,iz,gd.NX,gd.NY)] = b->buf[P3(ix+ixoff,iy+iyoff,iz,buf0nx,buf0ny)];
//						if(iz==gd.NZ-1)printf("%f\n",b->threedbuf[P3(ix,iy,iz,gd.NX,gd.NY)]);
					}
				}
			}

		status = nc_put_vara_float (nc.ncid, nc.var3d[ivar].varnameid, nc.start, nc.edges, b->threedbuf);
		}
		BL;
	}
	if(cmd.twodwrite)free(twodbuf);
}
