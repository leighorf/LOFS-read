#define ERROR_STOP(string) { fprintf(stderr," *** Fatal error in %s line %i: %s\n",__FILE__,__LINE__,string); exit(0); }
#define ERROR_WARN(string) fprintf(stderr," *** Warning: %s line %i: %s\n",__FILE__,__LINE__,string);
