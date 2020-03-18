#include <omp.h>
#include "../include/lofs-read.h"
#include "../include/dirstruct.h"
#include "../include/hdf2nc.h"
#include "../include/limits.h"
#include "../include/macros.h"

/* Note, we use George's i,j,k and ni,nj,nk approach although we personally prefer ix,iy,iz and
 * nx,ny,nz. Because some of our macros use the nx,ny,nz approach we copy ni,nj,nk to a local nx,ny,nz
 * in some of the functions that is then used in the macro */

/*******************************************************************************/

#define UINTERP BUF
void calc_uinterp(buffers *b, grid gd, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		UINTERP(i,j,k) = 0.5*(UA(i,j,k)+UA(i+1,j,k));
}

/*******************************************************************************/

#define VINTERP BUF
void calc_vinterp(buffers *b, grid gd, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		VINTERP(i,j,k) = 0.5*(VA(i,j,k)+VA(i,j+1,k));
}

/*******************************************************************************/

#define WINTERP BUF
void calc_winterp(buffers *b, grid gd, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		WINTERP(i,j,k) = 0.5*(WA(i,j,k)+WA(i,j,k+1));
}

/*******************************************************************************/

#define HWIN_SR BUF
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
		usr = 0.5*(UA(i,j,k)+UA(i+1,j,k));
		vsr = 0.5*(VA(i,j,k)+VA(i,j+1,k));
		HWIN_SR(i,j,k) = sqrt(usr*usr+vsr*vsr);
	}
}

/*******************************************************************************/

#define HWIN_GR BUF
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
		ugr = 0.5*(UA(i,j,k)+UA(i+1,j,k)) + msh.umove;
		vgr = 0.5*(VA(i,j,k)+VA(i,j+1,k)) + msh.vmove;
		HWIN_GR(i,j,k) = sqrt(ugr*ugr+vgr*vgr);
	}
}

/*******************************************************************************/

#define WINDMAG_SR BUF
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
		usr = 0.5*(UA(i,j,k)+UA(i+1,j,k));
		vsr = 0.5*(VA(i,j,k)+VA(i,j+1,k));
		  w = 0.5*(WA(i,j,k)+WA(i,j,k+1));
		WINDMAG_SR(i,j,k) = sqrt(usr*usr+vsr*vsr+w*w);
	}
}

/*******************************************************************************/

#define HDIV BUF
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
		dudx = (UA(i+1,j,k)-UA(i,j,k))*rdx*UH(i);
		dvdy = (VA(i,j+1,k)-VA(i,j,k))*rdy*VH(j);
		HDIV(i,j,k) = dudx + dvdy;
	}
}

/*******************************************************************************/

#define XVORT BUF
void calc_xvort(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dwdy,dvdz,rdy,rdz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;
	rdy=msh.rdy; rdz=msh.rdz;

#pragma omp parallel for private(i,j,k,dwdy,dvdz)
	for(k=1; k<nk; k++)
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni; i++)
	{
		dwdy = (WA(i,j,k)-WA(i,j-1,k))*rdy*VF(j);//was i
		dvdz = (VA(i,j,k)-VA(i,j,k-1))*rdz*MF(k);
		TEM(i,j,k) = dwdy - dvdz;
	}
//This is dependent upon our current free slip bc, see CM1 for other decisions
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni; i++)
	{
		TEM(i,j,0)=TEM(i,j,1);
		TEM(i,j,nk)=TEM(i,j,nk-1);
	}

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		XVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i,j+1,k)+TEM(i,j,k+1)+TEM(i,j+1,k+1));
}

/*******************************************************************************/

#define YVORT BUF
void calc_yvort(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dudz,dwdx,rdx,rdz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;
	rdx=msh.rdx; rdz=msh.rdz;

#pragma omp parallel for private(i,j,k,dudz,dwdx)
	for(k=1; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		dudz = (UA(i,j,k)-UA(i,j,k-1))*rdz*MF(k);
		dwdx = (WA(i,j,k)-WA(i-1,j,k))*rdx*UF(i);
		TEM(i,j,k) = dudz - dwdx;
	}
//This is dependent upon our current free slip bc, see CM1 for other
//decisions
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		TEM(i,j,0)=TEM(i,j,1);
		TEM(i,j,nk)=TEM(i,j,nk-1);
	}
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		YVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i+1,j,k)+TEM(i,j,k+1)+TEM(i+1,j,k+1));
}

/*******************************************************************************/

#define ZVORT BUF
void calc_zvort(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dvdx,dudy,rdx,rdy;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;
	rdx=msh.rdx; rdy=msh.rdy;

#pragma omp parallel for private(i,j,k,dvdx,dudy)
	for(k=0; k<nk; k++)
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni+1; i++)
	{
		dvdx = (VA(i,j,k)-VA(i-1,j,k))*rdx*UF(i);
		dudy = (UA(i,j,k)-UA(i,j-1,k))*rdy*VF(j);
		TEM(i,j,k) = dvdx - dudy;
	}
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		ZVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i+1,j,k)+TEM(i,j+1,k)+TEM(i+1,j+1,k));
}

/*******************************************************************************/

void calc_hvort(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dwdy,dvdz,dudz,dwdx,rdx,rdy,rdz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;
	rdx=msh.rdx; rdy=msh.rdy; rdz=msh.rdz;

#pragma omp parallel for private(i,j,k,dwdy,dvdz)
	for(k=1; k<nk; k++)
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni; i++)
	{
		dwdy = (WA(i,j,k)-WA(i,j-1,k))*rdy*VF(j);//was i
		dvdz = (VA(i,j,k)-VA(i,j,k-1))*rdz*MF(k);
		TEM(i,j,k) = dwdy - dvdz;
	}
//This is dependent upon our current free slip bc, see CM1 for other decisions
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni; i++)
	{
		TEM(i,j,0)=TEM(i,j,1);
		TEM(i,j,nk)=TEM(i,j,nk-1);
	}

#define XVORT BUF
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		XVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i,j+1,k)+TEM(i,j,k+1)+TEM(i,j+1,k+1));


#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1(i,j,k) = XVORT(i,j,k)*XVORT(i,j,k);

#define YVORT BUF
#pragma omp parallel for private(i,j,k,dudz,dwdx)
	for(k=1; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		dudz = (UA(i,j,k)-UA(i,j,k-1))*rdz*MF(k);
		dwdx = (WA(i,j,k)-WA(i-1,j,k))*rdx*UF(i);
		TEM(i,j,k) = dudz - dwdx;
	}
//This is dependent upon our current free slip bc, see CM1 for other
//decisions
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		TEM(i,j,0)=TEM(i,j,1);
		TEM(i,j,nk)=TEM(i,j,nk-1);
	}
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		YVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i+1,j,k)+TEM(i,j,k+1)+TEM(i+1,j,k+1));

	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1(i,j,k) += YVORT(i,j,k)*YVORT(i,j,k);

#define HVORT BUF
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		HVORT(i,j,k) = sqrt(TEM1(i,j,k));
}

/*******************************************************************************/

void calc_vortmag(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dwdy,dvdz,dudz,dwdx,dvdx,dudy,rdx,rdy,rdz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;
	rdx=msh.rdx; rdy=msh.rdy; rdz=msh.rdz;

#pragma omp parallel for private(i,j,k,dwdy,dvdz)
	for(k=1; k<nk; k++)
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni; i++)
	{
		dwdy = (WA(i,j,k)-WA(i,j-1,k))*rdy*VF(j);//was i
		dvdz = (VA(i,j,k)-VA(i,j,k-1))*rdz*MF(k);
		TEM(i,j,k) = dwdy - dvdz;
	}
//This is dependent upon our current free slip bc, see CM1 for other decisions
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni; i++)
	{
		TEM(i,j,0)=TEM(i,j,1);
		TEM(i,j,nk)=TEM(i,j,nk-1);
	}

#define XVORT BUF
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		XVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i,j+1,k)+TEM(i,j,k+1)+TEM(i,j+1,k+1));


#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1(i,j,k) = XVORT(i,j,k)*XVORT(i,j,k);

#pragma omp parallel for private(i,j,k,dudz,dwdx)
	for(k=1; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		dudz = (UA(i,j,k)-UA(i,j,k-1))*rdz*MF(k);
		dwdx = (WA(i,j,k)-WA(i-1,j,k))*rdx*UF(i);
		TEM(i,j,k) = dudz - dwdx;
	}
//This is dependent upon our current free slip bc, see CM1 for other
//decisions
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		TEM(i,j,0)=TEM(i,j,1);
		TEM(i,j,nk)=TEM(i,j,nk-1);
	}
#define YVORT BUF
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		YVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i+1,j,k)+TEM(i,j,k+1)+TEM(i+1,j,k+1));

	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1(i,j,k) += YVORT(i,j,k)*YVORT(i,j,k);

#pragma omp parallel for private(i,j,k,dvdx,dudy)
	for(k=0; k<nk; k++)
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni+1; i++)
	{
		dvdx = (VA(i,j,k)-VA(i-1,j,k))*rdx*UF(i);
		dudy = (UA(i,j,k)-UA(i,j-1,k))*rdy*VF(j);
		TEM(i,j,k) = dvdx - dudy;
	}
#define ZVORT BUF
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		ZVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i+1,j,k)+TEM(i,j+1,k)+TEM(i+1,j+1,k));

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1(i,j,k) += ZVORT(i,j,k)*ZVORT(i,j,k);

#define VORTMAG BUF
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		VORTMAG(i,j,k) = sqrt(TEM1(i,j,k));
}

/*******************************************************************************/

#define STREAMVORT BUF
void calc_streamvort(buffers *b, grid gd, mesh msh, cmdline cmd)
{
	int i,j,k,ni,nj,nk,nx,ny,nz;
	float dwdy,dvdz,dudz,dwdx,dvdx,dudy,rdx,rdy,rdz;
	float uinterp,vinterp,winterp;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;
	rdx=msh.rdx; rdy=msh.rdy; rdz=msh.rdz;

#pragma omp parallel for private(i,j,k,dwdy,dvdz)
	for(k=1; k<nk; k++)
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni; i++)
	{
		dwdy = (WA(i,j,k)-WA(i,j-1,k))*rdy*VF(j);//was i
		dvdz = (VA(i,j,k)-VA(i,j,k-1))*rdz*MF(k);
		TEM(i,j,k) = dwdy - dvdz;
	}
//This is dependent upon our current free slip bc, see CM1 for other decisions
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni; i++)
	{
		TEM(i,j,0)=TEM(i,j,1);
		TEM(i,j,nk)=TEM(i,j,nk-1);
	}

#define XVORT BUF
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		XVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i,j+1,k)+TEM(i,j,k+1)+TEM(i,j+1,k+1));


#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1(i,j,k) =  XVORT(i,j,k)*0.5*(UA(i,j,k)+UA(i+1,j,k));

#pragma omp parallel for private(i,j,k,dudz,dwdx)
	for(k=1; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		dudz = (UA(i,j,k)-UA(i,j,k-1))*rdz*MF(k);
		dwdx = (WA(i,j,k)-WA(i-1,j,k))*rdx*UF(i);
		TEM(i,j,k) = dudz - dwdx;
	}
//This is dependent upon our current free slip bc, see CM1 for other
//decisions
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		TEM(i,j,0)=TEM(i,j,1);
		TEM(i,j,nk)=TEM(i,j,nk-1);
	}
#define YVORT BUF
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		YVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i+1,j,k)+TEM(i,j,k+1)+TEM(i+1,j,k+1));

	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1(i,j,k) += YVORT(i,j,k)*0.5*(VA(i,j,k)+VA(i,j+1,k));

#pragma omp parallel for private(i,j,k,dvdx,dudy)
	for(k=0; k<nk; k++)
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni+1; i++)
	{
		dvdx = (VA(i,j,k)-VA(i-1,j,k))*rdx*UF(i);
		dudy = (UA(i,j,k)-UA(i,j-1,k))*rdy*VF(j);
		TEM(i,j,k) = dvdx - dudy;
	}
#define ZVORT BUF
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		ZVORT(i,j,k) = 0.25 * (TEM(i,j,k)+TEM(i+1,j,k)+TEM(i,j+1,k)+TEM(i+1,j+1,k));

#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
		TEM1(i,j,k) += ZVORT(i,j,k)*0.5*(WA(i,j,k)+WA(i,j,k+1));

#define STREAMVORT BUF
#pragma omp parallel for private(i,j,k)
	for(k=0; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni; i++)
	{
		uinterp=0.5*(UA(i,j,k)+UA(i+1,j,k));
		vinterp=0.5*(VA(i,j,k)+VA(i,j+1,k));
		winterp=0.5*(WA(i,j,k)+WA(i,j,k+1));
		STREAMVORT(i,j,k) = TEM1(i,j,k)/(sqrt(uinterp*uinterp+vinterp*vinterp+winterp*winterp));
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
		BUF(i,j,k)=UA(i,j,k);
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
		BUF(i,j,k)=VA(i,j,k);
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
		BUF(i,j,k)=WA(i,j,k);
}

/*******************************************************************************/


void do_requested_variables(buffers *b, ncstruct nc, grid gd, mesh msh, readahead rh,dir_meta dm,hdf_meta hm,cmdline cmd)
{
	int ivar,status;
	requested_cube rc;
	char *var;
	float *wp;

// For flexibility we always set rc in case we need to read outside of what we
// requested at the command line such as what occurs with staggered variables. The
// copy_grid_to_requested_cube just duplicates the grid structure stuff to rc (requested
// cube). rc is a required argument to read_lofs_buffer.

// For doing calculations that involve derivatives (usually of u, v, or w) we have to read
// data in slightly differently (see for instance do_readahead where we read in ustag
// vstag wstag on a bigger mesh than the scalar variables).

	copy_grid_to_requested_cube(&rc,gd);

	var = (char *) malloc (MAXSTR * sizeof(char));

	for (ivar = 0; ivar < cmd.nvar; ivar++)
	{
		var=nc.varname[ivar];
		printf("%s: ",var);FL;

		if(same(var,"u"))
		{
			if(!rh.u)
			{
				read_lofs_buffer(b->buf,var,dm,hm,rc,cmd); wp=b->buf;
			}
			else
			{
				buf_u(b,gd); wp=b->buf;
			}
		}
		else if(same(var,"v"))
		{
			if(!rh.v)
			{
				read_lofs_buffer(b->buf,var,dm,hm,rc,cmd); wp=b->buf;
			}
			else
			{
				buf_v(b,gd); wp=b->buf;
			}
		}
		else if(same(var,"w"))
		{
			if(!rh.w)
			{
				read_lofs_buffer(b->buf,var,dm,hm,rc,cmd); wp=b->buf;
			}
			else
			{
				buf_w(b,gd); wp=b->buf;
			}
		}
		else if(same(var,"uinterp"))	{CL;calc_uinterp(b,gd,cmd);wp=b->buf;}
		else if(same(var,"vinterp"))	{CL;calc_vinterp(b,gd,cmd);wp=b->buf;}
		else if(same(var,"winterp"))	{CL;calc_winterp(b,gd,cmd);wp=b->buf;}
		else if(same(var,"hwin_sr"))	{CL;calc_hwin_sr(b,gd,cmd);wp=b->buf;}
		else if(same(var,"hwin_gr"))	{CL;calc_hwin_gr(b,gd,msh,cmd);wp=b->buf;}
		else if(same(var,"windmag_sr"))	{CL;calc_windmag_sr(b,gd,cmd);wp=b->buf;}
		else if(same(var,"hdiv")) 		{CL;calc_hdiv(b,gd,msh,cmd);wp=b->buf;}
		else if(same(var,"xvort"))		{CL;calc_xvort(b,gd,msh,cmd);wp=b->buf;}
		else if(same(var,"yvort"))		{CL;calc_yvort(b,gd,msh,cmd);wp=b->buf;}
		else if(same(var,"zvort"))		{CL;calc_zvort(b,gd,msh,cmd);wp=b->buf;}
		else if(same(var,"hvort"))		{CL;calc_hvort(b,gd,msh,cmd);wp=b->buf;}
		else if(same(var,"vortmag"))	{CL;calc_vortmag(b,gd,msh,cmd);wp=b->buf;}
		else if(same(var,"streamvort"))	{CL;calc_streamvort(b,gd,msh,cmd);wp=b->buf;}
		else
		{
			printf("reading...");FL;
			read_lofs_buffer(b->buf,nc.varname[ivar],dm,hm,rc,cmd);
			wp=b->buf;
			BL;
		}
		printf("writing...");FL;
		status = nc_put_vara_float (nc.ncid, nc.varnameid[ivar], nc.start, nc.edges, wp);
		BL;
	}
}

