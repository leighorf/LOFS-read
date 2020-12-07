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

#define UINTERP BUFp
void calc_uinterp(buffers *b, grid gd, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		UINTERP(i,j,k) = 0.5*(UAp(i,j,k)+UAp(i+1,j,k));
}

/*******************************************************************************/

#define VINTERP BUFp
void calc_vinterp(buffers *b, grid gd, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		VINTERP(i,j,k) = 0.5*(VAp(i,j,k)+VAp(i,j+1,k));
}

/*******************************************************************************/

#define WINTERP BUFp
void calc_winterp(buffers *b, grid gd, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		WINTERP(i,j,k) = 0.5*(WAp(i,j,k)+WAp(i,j,k+1));
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

	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		XVORT(i,j,k) = 0.25 * (TEMp(i,j,k)+TEMp(i,j+1,k)+TEMp(i,j,k+1)+TEMp(i,j+1,k+1));
}

/*******************************************************************************/

#define YVORT BUFp
void do_yvort(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dx,dz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

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
}

/*******************************************************************************/

#define ZVORT BUFp
void do_zvort(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dx,dy;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk; k++) {
		for(j=0; j<nj+1; j++) {
			dy = 1./(msh.rdy * VF(j)); //We might want to just do rdy, much cleaner... and dx,dy,dz are always in the denominator
			for(i=0; i<ni+1; i++) {
				dx = 1./(msh.rdx * UF(i)); 
    			calc_zvort(b->ustag, b->vstag, b->dum0, dx, dy, i, j, k, ni, nj);
			}
		}
	}

	#pragma omp parallel for private(i,j,k) 
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		ZVORT(i,j,k) = 0.25 * (TEMp(i,j,k)+TEMp(i+1,j,k)+TEMp(i,j+1,k)+TEMp(i+1,j+1,k));
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
#define YVORT BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		YVORT(i,j,k) = 0.25 * (TEMp(i,j,k)+TEMp(i+1,j,k)+TEMp(i,j,k+1)+TEMp(i+1,j,k+1));

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
    			calc_zvort(b->ustag, b->vstag, b->dum0, dx, dy, i, j, k, ni, nj);
			}
		}
	}
#define ZVORT BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		ZVORT(i,j,k) = 0.25 * (TEMp(i,j,k)+TEMp(i+1,j,k)+TEMp(i,j+1,k)+TEMp(i+1,j+1,k));

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
		TEM1p(i,j,k) =  XVORT(i,j,k)*0.5*(UAp(i,j,k)+UAp(i+1,j,k));

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
#define YVORT BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		YVORT(i,j,k) = 0.25 * (TEMp(i,j,k)+TEMp(i+1,j,k)+TEMp(i,j,k+1)+TEMp(i+1,j,k+1));

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
    			calc_zvort(b->ustag, b->vstag, b->dum0, dx, dy, i, j, k, ni, nj);
			}
		}
	}
#define ZVORT BUFp
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		ZVORT(i,j,k) = 0.25 * (TEMp(i,j,k)+TEMp(i+1,j,k)+TEMp(i,j+1,k)+TEMp(i+1,j+1,k));

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

/*******************************************************************************/

void z_progress_bar(int iz, int nz)
{
	int index;
	char *alph="abcdefghijklmnopqrstuvwxyz";
	index=(int)((float)iz*26.0/(float)nz);
	printf("%c",alph[index]); fflush (stdout);
}

/*******************************************************************************/

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

		var=nc.varname[ivar];
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
		else if(same(var,"pipert"))     {CL;do_pipert(b,gd,snd,cmd);}
		else if(same(var,"wb_buoy"))    {CL;do_wbuoy(b,gd,snd,cmd);}
		else if(same(var,"ub_pgrad"))   {CL;do_upgrad(b,gd,snd,msh,cmd);}
		else if(same(var,"vb_pgrad"))   {CL;do_vpgrad(b,gd,snd,msh,cmd);}
		else if(same(var,"wb_pgrad"))   {CL;do_wpgrad(b,gd,snd,msh,cmd);}
		else if(same(var,"uinterp"))	{CL;calc_uinterp(b,gd,cmd);}
		else if(same(var,"vinterp"))	{CL;calc_vinterp(b,gd,cmd);}
		else if(same(var,"winterp"))	{CL;calc_winterp(b,gd,cmd);}
		else if(same(var,"hwin_sr"))	{CL;calc_hwin_sr(b,gd,cmd);}
		else if(same(var,"hwin_gr"))	{CL;calc_hwin_gr(b,gd,msh,cmd);}
		else if(same(var,"windmag_sr"))	{CL;calc_windmag_sr(b,gd,cmd);}
		else if(same(var,"hdiv")) 		{CL;calc_hdiv(b,gd,msh,cmd);}
		else if(same(var,"xvort"))		{CL;do_xvort(b,gd,msh,cmd);}
		else if(same(var,"yvort"))		{CL;do_yvort(b,gd,msh,cmd);}
		else if(same(var,"zvort"))		{CL;do_zvort(b,gd,msh,cmd);}
		else if(same(var,"hvort"))		{CL;calc_hvort(b,gd,msh,cmd);}
		else if(same(var,"vortmag"))	{CL;calc_vortmag(b,gd,msh,cmd);}
		else if(same(var,"streamvort"))	{CL;calc_streamvort(b,gd,msh,cmd);}
// At some point after we've calculated all the stuff that requires
// buffered u v w stuff we need to repurpose some of those buffers for
// calculating other things - like temperature, which requires pressure
// and density since I don't save pipert
//		else if(same(var,"tempk"))	{CL;calc_tempk(b,gd,msh,cmd);}
		else
		{
			printf("reading...");FL;
			buf0nx=gd.NX;ixoff=0;
			buf0ny=gd.NY;iyoff=0;
			read_lofs_buffer(b->buf,nc.varname[ivar],dm,hm,rc,cmd);
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
				else                  //XYZ tetrahedra (cubey thing)
				{
					writestart[0]=0;  writeedges[0]=1;     //time
					writestart[1]=iz; writeedges[1]=1;     //z
					writestart[2]=0;  writeedges[2]=gd.NY; //y
					writestart[3]=0;  writeedges[3]=gd.NX; //x
				}
				status = nc_put_vara_float (nc.ncid, nc.varnameid[ivar], writestart, writeedges, twodbuf);
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
					}
				}
			}

		status = nc_put_vara_float (nc.ncid, nc.varnameid[ivar], nc.start, nc.edges, b->threedbuf);
		}
		BL;
	}
	if(cmd.twodwrite)free(twodbuf);
}
