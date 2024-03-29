#ifndef MACROS
#define MACROS
#define ERROR_STOP(string) { fprintf(stderr," *** Fatal error in %s line %i: %s\n",__FILE__,__LINE__,string); exit(0); }
#define ERROR_WARN(string) fprintf(stderr," *** Warning: %s line %i: %s\n",__FILE__,__LINE__,string);

#define P2(x,y,mx) (((y)*(mx))+(x))
#define P3(x,y,z,mx,my) (long)((long)((long)((long)z)*((long)mx)*((long)my))+(long)((long)((long)y)*((long)mx))+((long)x))
//ORF TODO: cast all P4 to long as well?
#define P4(x, y, z, t, mx, my, mz) (((t)*(mx)*(my)*(mz))+((z)*(mx)*(my))+((y)*(mx))+(x))

#define TRUE (1)
#define FALSE (0)

#define same(a,b) !strcmp(a,b)

//NOTE: these macros must be carefully used.
//You must name your mesh structure mesh
//And you have to pick whether you're in pointer land or not

#define UHp(ix) msh->uh[ix]
#define UFp(ix) msh->uf[ix]
#define VHp(iy) msh->vh[iy]
#define VFp(iy) msh->vf[iy]
#define MHp(iz) msh->mh[iz]
#define MFp(iz) msh->mf[iz]

#define UH(ix) msh.uh[ix]
#define UF(ix) msh.uf[ix]
#define VH(iy) msh.vh[iy]
#define VF(iy) msh.vf[iy]
#define MH(iz) msh.mh[iz]
#define MF(iz) msh.mf[iz]

#define xh(ix) msh->xhout[ix+1]
#define xf(ix) msh->xfout[ix+1]
#define yh(iy) msh->yhout[iy+1]
#define yf(iy) msh->yfout[iy+1]
#define zh(iz) msh->zhout[iz]
#define zf(iz) msh->zfout[iz]

// ORF: OK being a bit clever here ... fun with macros. This will make the code a lot
// easier to compare to native CM1 Fortran90 code that we are copying anyway. I adopt TEM
// for his tem array, UA for ua etc.
//
// These 4D ones are currently not used
#define BUF4D(x,y,z,t) buf0[P4(x+1, y+1, z  ,t, nx+2, ny+2, nz+1)]
#define TEM4D(x,y,z,t) dum0[P4(x  , y  , z  ,t, nx+2, ny+2, nz+1)]
#define  UA4D(x,y,z,t) ustag[P4(x+1, y+1, z, t, nx+2, ny+2, nz+1)]
#define  VA4D(x,y,z,t) vstag[P4(x+1, y+1, z, t, nx+2, ny+2, nz+1)]
#define  WA4D(x,y,z,t) wstag[P4(x+1, y+1, z, t, nx+2, ny+2, nz+1)]
#define  KM4D(x,y,z,t) kmstag[P4(x+1, y+1, z, t, nx+2, ny+2, nz+1)]

#define BUFp(x,y,z)  b->buf0[P3(x+1, y+1, z, nx+2, ny+2)]
#define TEMp(x,y,z)  b->dum0[P3(x+1, y+1, z, nx+2, ny+2)]
#define TEM1p(x,y,z) b->dum1[P3(x+1, y+1, z, nx+2, ny+2)]

#define  UAp(x,y,z) b->ustag[P3(x+1, y+1, z, nx+2, ny+2)]
#define  VAp(x,y,z) b->vstag[P3(x+1, y+1, z, nx+2, ny+2)]
#define  WAp(x,y,z) b->wstag[P3(x+1, y+1, z, nx+2, ny+2)]

#define  PRESPERTp(x,y,z) b->ppert[P3(x+1, y+1, z, nx+2, ny+2)]
#define  THRHOPERTp(x,y,z) b->thrhopert[P3(x+1, y+1, z, nx+2, ny+2)]

#define BUF(x,y,z)  buf0[P3(x+1, y+1, z, nx+2, ny+2)]
#define TEM(x,y,z)  dum0[P3(x+1, y+1, z, nx+2, ny+2)]
#define TEM1(x,y,z) dum1[P3(x+1, y+1, z, nx+2, ny+2)]

#define  UA(x,y,z) ustag[P3(x+1, y+1, z, nx+2, ny+2)]
#define  VA(x,y,z) vstag[P3(x+1, y+1, z, nx+2, ny+2)]
#define  WA(x,y,z) wstag[P3(x+1, y+1, z, nx+2, ny+2)]

#define PCL(t,p,mt) (((p)*(mt))+(t))

#define FL fflush(stdout);
#define CL {printf("calculating...");FL}
#define BL {printf("\n");FL}

#endif
