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

#define UH(ix) uh[ix+1]
#define UF(ix) uf[ix+1]
#define VH(iy) vh[iy+1]
#define VF(iy) vf[iy+1]
#define MH(iz) mh[iz]
#define MF(iz) mf[iz]

#define xh(ix) xh[ix+1]
#define xf(ix) xf[ix+1]
#define yh(iy) yh[iy+1]
#define yf(iy) yf[iy+1]
#define zh(iz) zh[iz]
#define zf(iz) zf[iz]

// ORF: OK being a bit clever here ... fun with macros. This will make the code a lot
// easier to compare to native CM1 Fortran90 code that we are copying anyway. I adopt TEM
// for his tem array, UA for ua etc.

#define BUF(x,y,z) buf0[P3(x,y,z,NX,NY)]
#define TEM(x,y,z) dum0[P3(x,y,z,NX+1,NY+1)]
#define TEM1(x,y,z) dum1[P3(x,y,z,NX+1,NY+1)]
#define  UA(x,y,z) ustag[P3(x+1,y+1,z,NX+2,NY+2)]
#define  VA(x,y,z) vstag[P3(x+1,y+1,z,NX+2,NY+2)]
#define  WA(x,y,z) wstag[P3(x+1,y+1,z,NX+2,NY+2)]

#define PCL(t,p,mt) (((p)*(mt))+(t))
#endif
