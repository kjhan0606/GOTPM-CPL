/* These are common variables used by all of the PM functions */
#ifdef DEFINE_PM_COMMON
#define EXTERN
#else
#define EXTERN extern
#endif

/* Particle and work arrays */
EXTERN long np;
EXTERN float zheight,zstart;
EXTERN pmparticletype *pmparticles;
EXTERN treeparticletype *treeparticles;
EXTERN float *work;
EXTERN float *den;
EXTERN int worksize;
EXTERN int densize;

/*   used for correlation functions */
EXTERN float *pk,*npk,*cor,*corw,*tpk,*tnpk,*tcor,*tcorw;

/*   Table for periodic boundary conditions */
EXTERN int *xabp0,*xabp1,*xabm1,*yabp0,*yabp1,*yabm1,*zabp0,*zabp1,*zabm1;
EXTERN int *rcounts,*displs;

/* Processor identity and number */
EXTERN int myrank, nrank;

#undef EXTERN
