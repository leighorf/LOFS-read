#define ERROR_STOP(string) { fprintf(stderr," *** Fatal error in %s line %i: %s\n",__FILE__,__LINE__,string); exit(0); }
#define ERROR_WARN(string) fprintf(stderr," *** Warning: %s line %i: %s\n",__FILE__,__LINE__,string);

#define P3(x,y,z,mx,my) (long)((long)((long)((long)z)*((long)mx)*((long)my))+(long)((long)((long)y)*((long)mx))+((long)x))
#define P2(x,y,mx) (((y)*(mx))+(x))
#define TRUE (1)
#define FALSE (0)
