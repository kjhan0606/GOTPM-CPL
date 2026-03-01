#define FILEHEADER "p1024.z"
#define indextype long  long
#define MPI_INDEX MPI_LONG_LONG
typedef struct pmparticletype {
	float x,y,z,vx,vy,vz;
#ifdef INDEX
	indextype indx;
#endif
} pmparticletype;
typedef struct pmvrparticletype {
	float vx,vy,vz;
} pmvrparticletype;
typedef struct treeparticletype {
	float x,y,z,vx,vy,vz;
#ifdef INDEX
	indextype indx;
#endif
	struct treeparticletype *next;
} treeparticletype;
float gettime();
#define TIMER_START(A) cputime0[A] = gettime();
#define TIMER_STOP(A)  cputime1[A] = gettime();
#define ELAPSED_TIME(A) (cputime1[A] - cputime0[A])

typedef struct simparameters{
	int myid,nid;
	float omep,omeplam,omepb,omei;
	float hubble,npow;
	float boxsize,bias;
	float amax,anow;
	int nx,ny,nz,nspace;
	double lnx,lny,lnz;
	indextype mx,my,mz,mxmy;
	float zinit,astep,rsmooth,rth;
	float theta,sphere_radius,particle_radius;
	float ken,ktot,const0,poten0;
	float fact1,fact2,pfact;
	int nstep,stepcount,stepnum,iseed,nskip;
	int powreadflag;
	char rvfilename[190],rvprefix[190];
	char powfilename[190],inpapkfilename[190];
	float zmin,zmax;
	int local_nz;
	long long np;
	int xyzshiftflag, ptypesize, nviewer;
	char Viewerfile[190];
} SimParameters;
#ifdef DEFINE_SIM_PARA
#define EXTERN
#else
#define EXTERN extern
#endif
EXTERN SimParameters simpar;
EXTERN float cputime0[100],cputime1[100];
#undef EXTERN



#ifdef XYZDBL
#	ifdef INDEX
#		define XofP(p) fmod((double)((p)->x)+(double)(((p)->indx%simpar.mx)*simpar.nspace)+simpar.nx,simpar.nx)
#		define YofP(p) fmod((double)((p)->y)+(double)((((p)->indx%simpar.mxmy)/simpar.mx)*simpar.nspace)+simpar.ny,simpar.ny)
#		define ZofP(p) fmod((double)((p)->z)+(double)(((p)->indx/simpar.mxmy)*simpar.nspace)+simpar.nz,simpar.nz)
#	else
#		error the pair definition of XYZDBL and !INDEX  is not allowed in the current version.
#	endif
#else
#	define XofP(p) ((p)->x)
#	define YofP(p) ((p)->y)
#	define ZofP(p) ((p)->z)
#endif
