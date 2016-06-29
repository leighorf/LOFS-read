/* config.h.  Generated automatically by configure.  */
/* config.h.in.  Generated automatically from configure.in by autoheader.  */
/* Here we define the default values for any preprocessor symbols used
   by autoconf but not known by autoheader.  This file is used to
   create config.h.in by running autoheader.  In most cases, they will
   default to being undefined, which is specified here by #undef.
   Remember, the actual values of these symbols that we use will be
   set by the configure script, and output into config.h. */

#ifndef CONFIG_H
#define CONFIG_H 1

/* Define if your processor stores words with the most significant
   byte first (like Motorola and SPARC, unlike Intel and VAX).  */
/*
#define WORDS_BIGENDIAN
*/

/* Define if the X Window System is missing or not being used.  */
/* #undef X_DISPLAY_MISSING */

/* define if we have POSIX threads */
/* #undef HAVE_PTHREADS */

/* define if we have SunOS threads */
/* #undef HAVE_SUNOS_THREADS */

/* define if we have SGI sproc forking */
/* #undef HAVE_SGI_SPROC */

/* define if we have OpenGL or Mesa */
#define HAVE_OPENGL 1

/* define if we have SGI IrisGL */
/* #undef HAVE_SGI_GL */

/* define if we have PEX */
/* #undef HAVE_PEX */

/* define if we have Tcl */
#define HAVE_LIBTCL 1

/* define if we have NetCDF */
#define HAVE_LIBNETCDF 1

/* define if we have Fortran idate function: */
#define HAVE_IDATE 1

/* define if we have the McIDAS library */
/* #undef HAVE_MCIDAS */

/* define if we are linking to a newer version of the McIDAS library
   than the one distributed with Vis5d: */
/* #undef MCIDAS_SIDECAR */

/* define if we are single-threading: */
#define SINGLE_TASK 1

/* define to installation location of data files (e.g. EARTH.TOPO);
   typically /usr/local/share/vis5d */
#define DATA_PREFIX "/usr/local/share/vis5d+/"

/* define to specify how to mangle identifiers for the Fortran linker: */
#define FORTRANIZE_LOWERCASE 0
#define FORTRANIZE_UPPERCASE 0
#define FORTRANIZE_LOWERCASE_UNDERSCORE 1 /* the default */
#define FORTRANIZE_UPPERCASE_UNDERSCORE 0
#define FORTRANIZE_EXTRA_UNDERSCORE 1

#if FORTRANIZE_LOWERCASE
#  define F77_FUNC(x,X) x
#elif FORTRANIZE_UPPERCASE
#  define F77_FUNC(x,X) X
#elif FORTRANIZE_LOWERCASE_UNDERSCORE
#  define F77_FUNC(x,X) x##_
#elif FORTRANIZE_UPPERCASE_UNDERSCORE
#  define F77_FUNC(x,X) X##_
#endif

/* The number of bytes in a float.  */
#define SIZEOF_FLOAT 4

/* The number of bytes in a int.  */
#define SIZEOF_INT 4

/* The number of bytes in a signed char.  */
#define SIZEOF_SIGNED_CHAR 1

/* Define if you have the XMesaGetBackBuffer function.  */
#define HAVE_XMESAGETBACKBUFFER 1

/* Define if you have the setrlimit function.  */
#define HAVE_SETRLIMIT 1

/* Define if you have the <X11/Xm/MwmUtil.h> header file.  */
/* #undef HAVE_X11_XM_MWMUTIL_H */

/* Define if you have the <netcdf.h> header file.  */
#define HAVE_NETCDF_H 1

/* Define if you have the <sys/lock.h> header file.  */
/* #undef HAVE_SYS_LOCK_H */

/* Define if you have the <sys/prctl.h> header file.  */
#define HAVE_SYS_PRCTL_H 1

/* Define if you have the <sys/sysmp.h> header file.  */
/* #undef HAVE_SYS_SYSMP_H */

/* Define if you have the <sys/types.h> header file.  */
#define HAVE_SYS_TYPES_H 1

/* Define if you have the <sysmp.h> header file.  */
/* #undef HAVE_SYSMP_H */

/* Define if you have the m library (-lm).  */
#define HAVE_LIBM 1

/* Name of package */
#define PACKAGE "vis5d+"

/* Version number of package */
#define VERSION "1.0.2"

/* max. memory to use (MB), 0 for no maximum */
#define VIS5D_MAX_MEM 32

/* Define to the necessary symbol if this constant
                           uses a non-standard name on your system. */
/* #undef PTHREAD_CREATE_JOINABLE */

#endif /* CONFIG_H */
