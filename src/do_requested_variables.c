#include <omp.h>
#include "../include/lofs-read.h"
#include "../include/lofs-dirstruct.h"
#include "../include/lofs-hdf2nc.h"
#include "../include/lofs-limits.h"
#include "../include/lofs-macros.h"

/* Note, we use George's i,j,k and ni,nj,nk approach although we personally prefer ix,iy,iz and
 * nx,ny,nz. Because some of our macros use the nx,ny,nz approach we copy ni,nj,nk to a local nx,ny,nz
 * in some of the functions that is then used in the macro */

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

#define XVORT BUFp
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
		dwdy = (WAp(i,j,k)-WAp(i,j-1,k))*rdy*VF(j);//was i
		dvdz = (VAp(i,j,k)-VAp(i,j,k-1))*rdz*MF(k);
		TEMp(i,j,k) = dwdy - dvdz;
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
		dudz = (UAp(i,j,k)-UAp(i,j,k-1))*rdz*MF(k);
		dwdx = (WAp(i,j,k)-WAp(i-1,j,k))*rdx*UF(i);
		TEMp(i,j,k) = dudz - dwdx;
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
		dvdx = (VAp(i,j,k)-VAp(i-1,j,k))*rdx*UF(i);
		dudy = (UAp(i,j,k)-UAp(i,j-1,k))*rdy*VF(j);
		TEMp(i,j,k) = dvdx - dudy;
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
	float dwdy,dvdz,dudz,dwdx,rdx,rdy,rdz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;
	rdx=msh.rdx; rdy=msh.rdy; rdz=msh.rdz;

#pragma omp parallel for private(i,j,k,dwdy,dvdz)
	for(k=1; k<nk; k++)
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni; i++)
	{
		dwdy = (WAp(i,j,k)-WAp(i,j-1,k))*rdy*VF(j);//was i
		dvdz = (VAp(i,j,k)-VAp(i,j,k-1))*rdz*MF(k);
		TEMp(i,j,k) = dwdy - dvdz;
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
#pragma omp parallel for private(i,j,k,dudz,dwdx)
	for(k=1; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		dudz = (UAp(i,j,k)-UAp(i,j,k-1))*rdz*MF(k);
		dwdx = (WAp(i,j,k)-WAp(i-1,j,k))*rdx*UF(i);
		TEMp(i,j,k) = dudz - dwdx;
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
	float dwdy,dvdz,dudz,dwdx,dvdx,dudy,rdx,rdy,rdz;

	ni=gd.NX;nj=gd.NY;nk=gd.NZ;
	nx=ni; ny=nj; nz=nk;
	rdx=msh.rdx; rdy=msh.rdy; rdz=msh.rdz;

#pragma omp parallel for private(i,j,k,dwdy,dvdz)
	for(k=1; k<nk; k++)
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni; i++)
	{
		dwdy = (WAp(i,j,k)-WAp(i,j-1,k))*rdy*VF(j);//was i
		dvdz = (VAp(i,j,k)-VAp(i,j,k-1))*rdz*MF(k);
		TEMp(i,j,k) = dwdy - dvdz;
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

#pragma omp parallel for private(i,j,k,dudz,dwdx)
	for(k=1; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		dudz = (UAp(i,j,k)-UAp(i,j,k-1))*rdz*MF(k);
		dwdx = (WAp(i,j,k)-WAp(i-1,j,k))*rdx*UF(i);
		TEMp(i,j,k) = dudz - dwdx;
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

#pragma omp parallel for private(i,j,k,dvdx,dudy)
	for(k=0; k<nk; k++)
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni+1; i++)
	{
		dvdx = (VAp(i,j,k)-VAp(i-1,j,k))*rdx*UF(i);
		dudy = (UAp(i,j,k)-UAp(i,j-1,k))*rdy*VF(j);
		TEMp(i,j,k) = dvdx - dudy;
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
		dwdy = (WAp(i,j,k)-WAp(i,j-1,k))*rdy*VF(j);//was i
		dvdz = (VAp(i,j,k)-VAp(i,j,k-1))*rdz*MF(k);
		TEMp(i,j,k) = dwdy - dvdz;
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

#pragma omp parallel for private(i,j,k,dudz,dwdx)
	for(k=1; k<nk; k++)
	for(j=0; j<nj; j++)
	for(i=0; i<ni+1; i++)
	{
		dudz = (UAp(i,j,k)-UAp(i,j,k-1))*rdz*MF(k);
		dwdx = (WAp(i,j,k)-WAp(i-1,j,k))*rdx*UF(i);
		TEMp(i,j,k) = dudz - dwdx;
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

#pragma omp parallel for private(i,j,k,dvdx,dudy)
	for(k=0; k<nk; k++)
	for(j=0; j<nj+1; j++)
	for(i=0; i<ni+1; i++)
	{
		dvdx = (VAp(i,j,k)-VAp(i-1,j,k))*rdx*UF(i);
		dudy = (UAp(i,j,k)-UAp(i,j-1,k))*rdy*VF(j);
		TEMp(i,j,k) = dvdx - dudy;
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

void do_requested_variables(buffers *b, ncstruct nc, grid gd, mesh msh, readahead rh,dir_meta dm,hdf_meta hm,cmdline cmd)
{
	int ix,iy,iz,nx,ny,buf0nx,buf0ny,i,ixoff,iyoff,ivar,status;
	long int bufsize;
	requested_cube rc;
	char *var;
	float *twodbuf,*threedbuf;
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
		else if(same(var,"uinterp"))	{CL;calc_uinterp(b,gd,cmd);}
		else if(same(var,"vinterp"))	{CL;calc_vinterp(b,gd,cmd);}
		else if(same(var,"winterp"))	{CL;calc_winterp(b,gd,cmd);}
		else if(same(var,"hwin_sr"))	{CL;calc_hwin_sr(b,gd,cmd);}
		else if(same(var,"hwin_gr"))	{CL;calc_hwin_gr(b,gd,msh,cmd);}
		else if(same(var,"windmag_sr"))	{CL;calc_windmag_sr(b,gd,cmd);}
		else if(same(var,"hdiv")) 		{CL;calc_hdiv(b,gd,msh,cmd);}
		else if(same(var,"xvort"))		{CL;calc_xvort(b,gd,msh,cmd);}
		else if(same(var,"yvort"))		{CL;calc_yvort(b,gd,msh,cmd);}
		else if(same(var,"zvort"))		{CL;calc_zvort(b,gd,msh,cmd);}
		else if(same(var,"hvort"))		{CL;calc_hvort(b,gd,msh,cmd);}
		else if(same(var,"vortmag"))	{CL;calc_vortmag(b,gd,msh,cmd);}
		else if(same(var,"streamvort"))	{CL;calc_streamvort(b,gd,msh,cmd);}
		else
		{
			printf("reading...");FL;
			buf0nx=gd.NX;ixoff=0;
			buf0ny=gd.NY;iyoff=0;
			read_lofs_buffer(b->buf,nc.varname[ivar],dm,hm,rc,cmd);
			BL;
		}
		printf("writing...");FL;

// ORF I tried this and it kind of sucked performance wise for large-ish
// data, but I leave the option available. This method writes our data
// in 2D XY slices rather than one big 3D chunk (saves a bit of data).
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
						threedbuf[P3(ix,iy,iz,gd.NX,gd.NY)] = b->buf[P3(ix+ixoff,iy+iyoff,iz,buf0nx,buf0ny)];
					}
				}
			}

		status = nc_put_vara_float (nc.ncid, nc.varnameid[ivar], nc.start, nc.edges, threedbuf);
		}
		BL;
	}
	if(cmd.twodwrite)free(twodbuf);
}
