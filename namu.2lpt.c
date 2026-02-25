/* My GOTPM code using the tree corrections developed by Juhan Kim
 * =======================================================================
 * It can now read the external power spectrum. Apr, 02, 2009
 * =======================================================================
 ***************************************************************************
 * In mytreecorrection, there is no preprocessor for the sake of simplicity.
 * parpmseed_force.F is not yet implemented into namu.c
 ***************************************************************************
 * Complete run on the Sun Blade system (CentOS, tachyon.ksc.re.kr) is done 
    with a success (Aug. 17th 2008)
 * Comparison to the original code has been checked with a success. 
   (Jan. 17th 2006)
 * Time evolution has been checked ? (Jan. 18th 2006)
 * FoF group test has been passed ? (Jan. 19th 2006)
 * FoF mass function has been checked when z=0 (Mar. 12th 2006)
 * Simulated matter power spectrum has been checked when z=0 & z=47 
 * (Mar. 13th 2006)
 */
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<mpi.h>
#include<stddef.h>
#include <stdlib.h>
#include<sys/times.h>
#include<omp.h>
#	include "pmheader.h"
#	include "common_simparam.h"
#	include "pm_common.h"
#	include "flow.h"
/*
#define DEFINE_SIM_PARA
#	include "pmheader.h"
#undef DEFINE_SIM_PARA

#define DEFINE_COMMON_SIMPARAM
#	include "common_simparam.h"
#undef DEFINE_COMMON_SIMPARAM

#define DEFINE_PM_COMMON
#	include "pm_common.h"
#undef DEFINE_PM_COMMON
#define DEFINE_SIM_FLOW
#	include "flow.h"
#undef DEFINE_SIM_FLOW
*/

#include "Memory.h"
#include "varpm.h"
//#include <srfftw_mpi.h>
#include "fortran.h"


#define MAXINHOMO 1. /* in % */
#define getsign(i) ((i)>0? 1:-1)
#define max(A,B) ( (A) > (B) ? (A) : (B) )
#define min(A,B) ( (A) < (B) ? (A) : (B) )
#define pow2(A) ( (A)*(A))
#define pow3(A) ( (A)*(A)*(A))
void migrate(size_t *,int,float,float);
int flagpsmeasure(float , float , float ,int );
static int iflagPreFoF;
int flagPreFoF(float, float, float, int);
static int iflagsyncpdata;
int flagsyncpdata(float , float , float );
static int iflagwholeden;
int flagwholeden(float , float , float );
static int iflagbinnedden;
int flagBinnedDen(float , float , float );
static float elapsedtime,treetime;

void pm(int, int, int, int, int, int);
static int nxpos,nsizepmparticletype;
static int nlzfloor,nlzceil,local_nz_varpm;
static float lzfloor,lzceil;
static int now,nowflag;
static int myid,nid;
/*
float omep,omeplam,size,astep,zi,omei,hubble;
*/
static float kenold;
static float poten;
static int iseed,local_z_start,local_ny_after_transpose;
static int local_y_start_after_transpose,total_local_size;
static int stepcount;
static int nforcepos;
static float pmas;

static float nextastep;
void jwrite(char *, int );
static int syncflag;

#define DomainDecomp {\
	int npmin,npmax,npav;\
	long nptot,lnp,tnp;\
	float ratio1,ratio2;\
	lnp = np;\
	MPI_Reduce(&lnp,&nptot,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);\
	if(myid==0) npav = nptot/nid;\
	MPI_Bcast(&npav,1,MPI_INT,0,MPI_COMM_WORLD);\
	MPI_Reduce(&np,&npmin,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);\
	MPI_Reduce(&np,&npmax,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);\
	MPI_Bcast(&npmin,1,MPI_INT,0,MPI_COMM_WORLD);\
	MPI_Bcast(&npmax,1,MPI_INT,0,MPI_COMM_WORLD);\
	ratio1 = (float)(npav-npmin)/(float)npav*100.;\
	ratio2 = (float)(npmax-npav)/(float)npav*100.;\
	if(myid==0)printf("The inhomogeneities are %g %g \n",ratio1,ratio2);\
	if(ratio1 > MAXINHOMO || ratio2 > MAXINHOMO) {\
		float zmin,zmax;\
		int mydomaindecomposition(pmparticletype *, int ,long ,\
				        float ,float , float *, float *,int,int,int,int );\
		TIMER_START(39)\
		lnp = np;\
		MPI_Reduce(&lnp,&tnp,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);\
		MPI_Bcast(&tnp,1,MPI_LONG,0,MPI_COMM_WORLD);\
		zmin = zstart; zmax = zstart+zheight;\
		printf("%g %g %d %d %d %d\n",zmin,zmax,simpar.nx,simpar.ny,simpar.nz,simpar.nspace);\
		np= mydomaindecomposition(pmparticles, np,tnp,zmin,zmax, &zmin, &zmax,simpar.nx,simpar.ny,\
				simpar.nz,simpar.nspace);\
		zstart = zmin; zheight = zmax-zmin;\
		simpar.zmin=zmin;simpar.zmax=zmax;\
		simpar.np = np;\
		TIMER_STOP(39);\
       	if(myid==0) fprintf(stdout,"DomainDecomp. CPU= %f \n", ELAPSED_TIME(39));\
	}\
}


#ifdef INCLUDE_TREE_FORCE
    static float efold=0.9;
#else
    static float efold = 0.5;
#endif
void Savingalldata();
void SaveWholeDen(float *);
void SaveBinnedDen(float *, int );
void SaveBinnedData(float *);
void tree2pmtype(long );
void Savexyslice(float *,int ,int ,int ,int ,int );
void Savexzslice(float *);
/*
void MAIN_(int argc, char *argv[]){
int MAIN__(int argc, char *argv[]){
int MAIN_(int argc, char *argv[]){
int main(int argc, char *argv[]){
#ifdef INTEL
int main(int argc, char *argv[]){
#elif PGCC
int MAIN_(int argc, char *argv[]){
#else
#error You must type the compiler for the fortran and C compatibility setting for main.
#endif
*/
int SecondOrderLPT(int argc, char *argv[]){
	int i,j,k;
	int nx, ny, nz, local_nz;
	int nspace;
	int ierror=1, iwrite=0, icont;
	char c;
	int memflag, InitializeBigMemory();
	FILE *simfile;
//	rfftwnd_mpi_plan plan,iplan;


	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	if(myid==0){
		printf("############################################################################\n\n\n");
		printf("Now generating Initial conditions using 2nd order Linear Perturbation Theory\n\n\n");
		printf("############################################################################\n\n\n");
	}


	if(myid==0){
		printf(" MPI SETTING \n");
		printf(" Total %d nodes are allocated\n",nid);
	}

	if( argc == 2 ) {
		int sflag=1;
		simfile = fopen(argv[1],"r");
		if(!simfile) sflag = 0;
		MPI_Bcast(&sflag,1,MPI_INT,0,MPI_COMM_WORLD);
		if(sflag ==0){
		    if(myid==0) fprintf(stderr,"Can't open %s - exiting\n",argv[1]);
			MPI_Finalize();
			exit(0);
		}
   	}
	else {
		fprintf(stderr,"usage: pyul paramsfile \n");
		MPI_Finalize();
		exit(0);
	}
	{
		/* reading parameters */
		int ReadSimulationParameters(FILE *, int *);
		ReadSimulationParameters(simfile,&icont);
	}

/* THIS IS FOR DEFINITION OF THE OPENMP */
#ifdef _OPENMP
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0){
		printf("##################### OpenMP Setting ##################\n");
		fflush(stdout);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#pragma omp parallel
	{
		if(omp_get_thread_num()==0 && myid ==0) 
			printf("P%d has %d threads\n",myid,omp_get_num_threads());
		fflush(stdout);
	}
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0){
		printf("#######################################################\n");
		fflush(stdout);
	}
#endif
/* THIS IS FOR DEFINITION OF THE OPENMP */

	if(Make_Total_Memory() == 0 ){
		fprintf(stderr,"Can't initialize memory - aborting job\n");
		MPI_Finalize();
		exit(0);
	}

	{
		/* determining local_nz for fftw */
		void set3dFftwInfo(int , int , int , MPI_Comm , int *, int *);
        printf("P%d has input domain (%d : %d: %d) \n",myid,simpar.nx,simpar.ny,simpar.nz);
        set3dFftwInfo(simpar.nx,simpar.ny,simpar.nz,MPI_COMM_WORLD, &local_nz,&local_z_start);
        printf("P%d has z local_domain (%d : %d) \n",myid,local_z_start,local_nz);
//		fftwinit_(&simpar.nx,&simpar.ny,&simpar.nz,&local_nz,&local_z_start);
		printf("P%d has z local_domain (%d : %d) \n",myid,local_z_start,local_nz);
		/* This is the default boundary constrained by the fftw */
		simpar.local_nz = local_nz;
		simpar.zmin = local_z_start;
		simpar.local_z_start = local_z_start;
		simpar.zmax = simpar.zmin + simpar.local_nz;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if(myid==0) {
		fprintf(stderr,"\n######################################\n\n");
		fprintf(stderr,"icont %d\n",icont);
		fprintf(stderr,"\n######################################\n");
	}
	iwrite = 0;
	if( icont == 0 ) { iwrite = 1; }

	nxpos = (int)offsetof(pmparticletype,x) /sizeof(float);
	nsizepmparticletype = (int)sizeof(pmparticletype) /sizeof(float);
	
	printf("1");

	zheight = (float) simpar.local_nz;
	zstart = (float) local_z_start;
	printf("2");
#ifdef OLD_FDA4
	fda4sup_(&simpar.nx,&simpar.ny,&simpar.nz, &simpar.local_nz);
        printf("3");
#endif
	/* PMSEED */
	if(icont ==0){
		float *den1,*den2;
		float *fx,*fy,*denpmseedpad;
		float vamp1,vamp2,damp1,damp2;
		{
			long densize,fftwsize;
			int izstart,izwidth,nadd,nzparticlewidth,nzintoffset;
			float *pvmesh;
			float zfinal,vamp,zstartpos;
			zheight = (float)simpar.nz/(float)nid; /*equal slab width */
			zstart = zheight*myid;
			zfinal = zheight*(myid+1);
			zheight = zfinal-zstart;

			izstart = (int)zstart;
			izwidth = ceil(zstart+zheight)-izstart;

			nzparticlewidth = 0;
			zstartpos = 2*simpar.nz;
			for(i=0;i<simpar.nz;i+=simpar.nspace){
				if((float)i>=zstart && (float)i<(zstart+zheight)) {
					nzparticlewidth ++;
					zstartpos=min(zstartpos,(float)i);
				}
			}
			nzintoffset = zstartpos-izstart;


			densize = 2*(simpar.nx/2+1)*(long) simpar.ny*(long) (simpar.local_nz+0);
			fftwsize = 2*(simpar.nx/2+1)*(long) simpar.ny*(ceil(zstart+zheight)-(int) zstart +1);
			if(densize==0) densize=10;
			den1=(float*) Malloc(sizeof(float)*densize, PPTR(den1));
			den2 = (float*) Malloc(sizeof(float)*densize,PPTR(den2));
			fx = (float*) Malloc(sizeof(float)*densize,PPTR(fx));
			fy = (float*) Malloc(sizeof(float)*densize,PPTR(fy));
			denpmseedpad = (float *) Malloc(2*simpar.ny*simpar.nz*sizeof(float),PPTR(denpmseedpad));
			void twolptpmseedforce_(float *, float *, float *, float*,
					float *,float *, float *, float *, float *, float *, 
					float *, float *, float *, float *, float *, float *, 
					float*, float *, float *, float *);
			twolptpmseedforce_(den1, &simpar.omep,&simpar.omeplam,&simpar.wlam0,
					&simpar.wlam1,&simpar.CsDE, &simpar.ken,
					denpmseedpad,&zstart,&zheight,&damp1, &vamp1, &damp2, &vamp2, &simpar.omei, fx,fy,den2,
					&simpar.fNL,&simpar.gNL);
			Free(denpmseedpad);
			Free(fy);
			Free(fx);


			zfinal = zstart+zheight;

			{
				long nmax;
				nmax = freespace()/sizeof(pmparticletype) - 100L 
					- (densize + fftwsize + nid*6+3*2*(simpar.nx/2+1)*simpar.ny)*sizeof(float)/sizeof(pmparticletype);
				pmparticles=(pmparticletype*) Malloc(sizeof(pmparticletype)*nmax,
						PPTR(pmparticles));
				if(myid==0) printf("Now P%d has maximum %ld particles\n",myid,nmax);
			}
			void TwoLPT_setboundaryparticlesforce(long *, pmparticletype *,
					int, float, int, float, float, float *, float *, float, float, float, float);
			TwoLPT_setboundaryparticlesforce(&np,pmparticles,
					nzparticlewidth,zstartpos,nzintoffset,
					zstart,zheight, den1, den2, damp1, damp2, vamp1, vamp2);
			if(myid==0) printf("Now passing boundaryparitlces\n");

			pmparticles=(pmparticletype*) Realloc(pmparticles,
						sizeof(pmparticletype)*np);
#ifdef XYZDBL
			printf("P%d %g %g %g %g %g %g\n",simpar.myid,XofP(pmparticles),YofP(pmparticles),ZofP(pmparticles),
					XofP(pmparticles+np-1),YofP(pmparticles+np-1),ZofP(pmparticles+np-1));
#else
			printf("P%d %g %g %g %g %g %g\n",myid,pmparticles[0].x,pmparticles[0].y,pmparticles[0].z, 
					pmparticles[np-1].x,pmparticles[np-1].y,pmparticles[np-1].z);
#endif
			{ size_t sizetnp; sizetnp=np;
				/*
				printf("P%d has zheight and zstart %g %g\n",myid,zheight,zstart);
				*/
				migrate(&sizetnp,simpar.nz,zheight,zstart);
				np = sizetnp;
			}
			{/* Maximize the number of pmparticles */
				size_t nmax;
				nmax = freespace()/sizeof(pmparticletype) - 100L 
					- (densize + fftwsize + nid*6+3*2*(simpar.nx/2+1)*simpar.ny)*sizeof(float)/sizeof(pmparticletype);
				pmparticles=(pmparticletype*) Realloc(pmparticles,
						sizeof(pmparticletype)*(np+nmax));
				printf("P%d has nmax = %ld\n",simpar.myid,nmax);
			}
			{
				pmvrparticletype *pvr;
				long nnp;
				pvr = (pmvrparticletype *)((pmparticletype*)(pmparticles+np));
				nnp = 0;
				MPI_Barrier(MPI_COMM_WORLD);
				if(myid==0) printf("Before getinnerparticlesforce");
				void TwoLPT_getinnerparticlesforce(long *, pmvrparticletype *,int, int, float, float,
						float *, float *, float, float);
				TwoLPT_getinnerparticlesforce(&nnp,pvr, nzparticlewidth,nzintoffset,
						zstart, zheight,den1, den2,vamp1, vamp2);
				if(myid==0) printf("P0 completed ");
				MPI_Barrier(MPI_COMM_WORLD);
				if(myid==0) printf("after getinnerparticlesforce\n");
				pmparticles=(pmparticletype*) Realloc(pmparticles,
						sizeof(pmparticletype)*(np+nnp));
				/* This is needed because the address of pvr has changed due to the
				 * Free(ff?) */
				pvr = (pmvrparticletype *)((pmparticletype*)(pmparticles+np));
				/* Some compilers are reported to fail in this loop. 
				 * In that case, one should change the for-loop to the while-loop
				 * to force the in-order loop. */
#ifdef INTEL
				i = nnp-1;
				while(i>=0){ /* copy in reverse order to prevent overwriting */
					pmparticles[i+np].vz = pvr[i].vz;
					pmparticles[i+np].vy = pvr[i].vy;
					pmparticles[i+np].vx = pvr[i].vx;
					i--;
				}
#else
				for(i=nnp-1;i>=0;i--){ /* copy in reverse order to prevent overwriting */
					pmparticles[i+np].vz = pvr[i].vz;
					pmparticles[i+np].vy = pvr[i].vy;
					pmparticles[i+np].vx = pvr[i].vx;
				}
#endif
				if(myid==0) printf("Before setinnerparticlesforce");
				void TwoLPT_setinnerparticlesforce(long *, pmparticletype *, int,int, float, float *, 
						 float, float, float, float *, float *, float, float);
				TwoLPT_setinnerparticlesforce(&np,pmparticles,
						nzparticlewidth,nzintoffset,zstartpos, &simpar.ken,zstart,zheight,simpar.omei,
						den1, den2, damp1, damp2);
				if(myid==0) printf("P0 completed ");
				MPI_Barrier(MPI_COMM_WORLD);
				if(myid==0) printf("after setinnerparticlesforce\n");
			}
			pmparticles=(pmparticletype*)Realloc(pmparticles,
					sizeof(pmparticletype)*np);
		}
		/*
		lzfloor = local_z_start; lzceil = local_z_start + simpar.local_nz;
		*/
		/* Now the C rouintes are taking the code */
		simpar.zmin = zstart;
		simpar.zmax = zstart+zheight;
		now = 0;
		nowflag = 0;
		{ size_t sizetnp; sizetnp=np;
			migrate(&sizetnp,simpar.nz,zheight,zstart);
			np = sizetnp;
		}
        if( myid == 0 ) {
			pmparticletype *bp;
			bp = pmparticles;
			fprintf(stderr,"after migrate\n");
#ifdef XYZDBL
			fprintf(stderr,"P%d particle0 is %g %g %g %g %g %g\n", myid,
					XofP(bp), YofP(bp), ZofP(bp), bp->vx, bp->vy, bp->vz);
#else
			fprintf(stderr,"P%d particle0 is %g %g %g %g %g %g\n", myid,
				bp->x, bp->y, bp->z, bp->vx, bp->vy, bp->vz);
#endif
		}
		{
			long lnp,tlnp;
			lnp = (long ) np;
			MPI_Reduce(&lnp,&tlnp,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
			if(myid == 0) fprintf(stderr,"Total particle number : %ld\n",tlnp);
		}
	}
	else if(icont ==1){
		int nmax;
		void jread(char *);
		nmax = freespace()/sizeof(pmparticletype)-1000L;
		pmparticles = (pmparticletype *)Malloc(sizeof(pmparticletype)*nmax,
				PPTR(pmparticles));
		jread(simpar.rvfilename);
		{ /* Overwright the fftw information */
			simpar.local_nz = local_nz;
			simpar.local_z_start = local_z_start;
		}

		pmparticles = (pmparticletype *)Realloc(pmparticles,np*sizeof(pmparticletype));
		if(0){
			double zmin,zmax;
			for(i=0;i<np;i++){
				zmin = min(zmin,ZofP(pmparticles+i));
				zmax = max(zmax,ZofP(pmparticles+i));
			}
			printf("P%d has data %g < z < %g with np=%ld\n",myid,zmin,zmax,np);
			printf("P%d has pdata %g %g %g,,, %g %g %g\n",myid,XofP(pmparticles),
					YofP(pmparticles),ZofP(pmparticles),pmparticles->vx,pmparticles->vy,
					pmparticles->vz);
			MPI_Finalize();exit(99);
		}
		{ size_t sizetnp; sizetnp=np;
			migrate(&sizetnp,simpar.nz,zheight,zstart);
			np = sizetnp;
		}
	}
	{
		long tlnp,lnp;
		float rngc;
		rngc = (float)simpar.nx*(float)simpar.ny*(float)simpar.nz;
		lnp = (long ) np;
		MPI_Reduce(&lnp,&tlnp,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Bcast(&tlnp,1,MPI_LONG,0,MPI_COMM_WORLD);
		pmas = rngc/(float)tlnp;

	}
	/* MAIN TIME LOOP FOR SIMULATION */
	/* nzheight : variable width == local_nz_varpm
	 * simpar.local_nz : fftw width */
	for(simpar.stepcount=simpar.stepnum; simpar.stepcount<=simpar.nstep; simpar.stepcount++) {
		if(DEMODEL==0){
			float poissoncor_(float *, float*, float *, float *, float *, float *, float *);
			simpar.pcorr = poissoncor_(&simpar.anow, &simpar.amax, &simpar.omep, &simpar.omeplam, 
						&simpar.wlam0, &simpar.wlam1, &simpar.CsDE);
			if(myid==0){
				printf("Poisson Correction factor %g at anow = %g\n",simpar.pcorr, simpar.anow);
			}
		}
		else if(DEMODEL==1){
			void InitializePoissonCorrection(); /* @ pcorr.c */
			if(simpar.myid==0) InitializePoissonCorrection();
			printf("Now exit the testing for pcorr \n");
			MPI_Finalize();
			exit(99);

		}
		int nzheight;
		{ if(icont !=1) Savingalldata(); icont = 0; }
		TIMER_START(51);
		if(0){ /* This is to make projected maps for simulated observers */
			int mkanimate(pmparticletype *, int);
			mkanimate(pmparticles,np);

		}
		/*
		if(0){// extract lightcone data if needed 
			void observerextmain();
			observerextmain();
		}
		*/
		if(simpar.stepcount == 1){
			halfstep.first = 1; halfstep.second = 0;
			void dumpmain(int); dumpmain(1);
		}
		{
			halfstep.first = halfstep.second=0;
			if((iflagsyncpdata=flagsyncpdata(simpar.amax,simpar.anow,simpar.astep))){
				halfstep.first = 1; halfstep.second = 0;
			}
			if( (iflagPreFoF=flagPreFoF(simpar.amax,simpar.anow,simpar.astep, simpar.stepcount))){
				halfstep.first = 1; halfstep.second = 0;
			}
			if(!iflagsyncpdata && !iflagPreFoF){
				halfstep.first = halfstep.second = 0;
			}
		}
halfevolution:

		{/* Tree correction */
			void TreeAndPreFoFmain(int );
			TreeAndPreFoFmain(iflagPreFoF);
		}
                
		{/* Dump data if needed */
			void dumpmain(int);
			dumpmain(iflagsyncpdata);
		}

		if(halfstep.first==1 && halfstep.second==0) {
			halfstep.second = 1;
			goto halfevolution;
		}
		/*
		if(0){// Dump lightcone data if needed
			void observersavemain();
			observersavemain();
		}
		*/

		/* One step forward with particle positions */
		{
			void onestepforwardposition(pmparticletype*,int);
			onestepforwardposition(pmparticles,np);
		}
		TIMER_STOP(51);
		{ size_t sizetnp; sizetnp=np;
			migrate(&sizetnp,simpar.nz,(simpar.zmax-simpar.zmin),simpar.zmin);
			np = sizetnp;
		}
		if(myid==0){
			pmparticletype *bp; bp = pmparticles;
			printf("+P%d %g %g %g %g %g %g\n",myid,bp->x,bp->y,bp->z, bp->vx,bp->vy,bp->vz);
			fprintf(stdout,"Step %d CPU= %f \n", simpar.stepcount,ELAPSED_TIME(51));
			elapsedtime = gettime();
			fprintf(stderr,"Wallclock time %d = %g\n",simpar.stepcount,elapsedtime);
		}
		MPI_Bcast(&elapsedtime, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		simpar.anow += simpar.astep;

		DomainDecomp;
	}



	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	exit(0);
}

