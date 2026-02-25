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
#define DEFINE_SIM_PARA
#	include "pmheader.h"
#undef DEFINE_SIM_PARA

#define DEFINE_COMMON_SIMPARAM
#	include "common_simparam.h"
#undef DEFINE_COMMON_SIMPARAM

#define DEFINE_PM_COMMON
#	include "pm_common.h"
#undef DEFINE_PM_COMMON

#include "Memory.h"
#include "varpm.h"
#include <srfftw_mpi.h>
#include "fortran.h"

#define DEFINE_SIM_FLOW
#	include "flow.h"
#undef DEFINE_SIM_FLOW

#define MAXINHOMO 1. /* in % */
#define getsign(i) ((i)>0? 1:-1)
#define max(A,B) ( (A) > (B) ? (A) : (B) )
#define min(A,B) ( (A) < (B) ? (A) : (B) )
#define pow2(A) ( (A)*(A))
#define pow3(A) ( (A)*(A)*(A))
void migrate(size_t *,int,float,float);
int flagpsmeasure(float , float , float ,int );
int iflagPreFoF, flagPreFoF(float, float, float);
int iflagsyncpdata,flagsyncpdata(float , float , float );
int iflagwholeden,flagwholeden(float , float , float );
int iflagbinnedden,flagBinnedDen(float , float , float );
float elapsedtime,treetime;

void pm(int, int, int, int, int, int);
int nxpos,nsizepmparticletype;
int nlzfloor,nlzceil,local_nz_varpm;
float lzfloor,lzceil;
int now,nowflag;
int myid,nid;
/*
float omep,omeplam,size,astep,zi,omei,hubble;
*/
float kenold;
float poten;
float *den,*fftwmesh;
int iseed,local_z_start,local_ny_after_transpose;
int local_y_start_after_transpose,total_local_size;
int stepcount;
int nforcepos;
float pmas;

float nextastep;
void jwrite(char *, int );
int syncflag;

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


void Savingalldata(){
	int suddenstop;
	int flagsuddenstop(int);
	suddenstop = flagsuddenstop(simpar.stepcount);
	if(suddenstop == 3) {
		MPI_Finalize();exit(0);
	}
	else if( suddenstop) {
		delta_a = simpar.astep;
		jwrite(simpar.rvprefix,simpar.stepcount);
		if(suddenstop == 1) {
			MPI_Finalize();exit(0);
		}
	}
}
void SaveWholeDen(float *fftwmesh){
	int i,j,k;
	int ii,niter,slicepixel;
	FILE *fp;
	char filename[100];
	MPI_Status status;
	float *tmpden;
	sprintf(filename,"WholeDen.%.5d",simpar.stepcount);
	slicepixel = 2*(simpar.nx/2+1)*simpar.ny;
	if(simpar.myid==0){
		int isizeof;
		isizeof = sizeof(SimParameters);
		fp = fopen(filename,"w");
		fwrite(&isizeof,sizeof(int),1,fp);
		fwrite(&simpar,sizeof(SimParameters),1,fp);
		fwrite(fftwmesh,sizeof(float),2*(simpar.nx/2+1)*simpar.ny*simpar.local_nz,fp);
		tmpden = (float*)Malloc(sizeof(float)*slicepixel,PPTR(tmpden));
		for(i=1;i<simpar.nid;i++){
			MPI_Recv(&niter,1,MPI_INT,i,i,MPI_COMM_WORLD,&status);
			for(j=0;j<niter;j++){
				MPI_Recv(tmpden,slicepixel,MPI_FLOAT,i,j,MPI_COMM_WORLD,
						&status);
				fwrite(tmpden,sizeof(float),slicepixel,fp);
			}
		}
		fclose(fp);
		Free(tmpden);
	}
	else {
		MPI_Send(&(simpar.local_nz),1,MPI_INT,0,simpar.myid,MPI_COMM_WORLD);
		for(i=0;i<simpar.local_nz;i++)
			MPI_Send(fftwmesh+i*slicepixel,slicepixel,MPI_FLOAT,0,i,
					MPI_COMM_WORLD);		
	}
}
void SaveBinnedDen(float *fftwmesh, int mode){
	int i,j,k;
	int ii,niter;
	long long slicepixel,bslicepixel;
	int nnx,nny,nnz;
	FILE *fp;
	char filename[100];
	MPI_Status status;
	float *tmpden;
	if(mode==0) sprintf(filename,"BinnedDen.%.5d",simpar.stepcount);
	else if(mode==1) sprintf(filename,"BinnedVx.%.5d",simpar.stepcount);
	else if(mode==2) sprintf(filename,"BinnedVy.%.5d",simpar.stepcount);
	else if(mode==3) sprintf(filename,"BinnedVz.%.5d",simpar.stepcount);
	else {
		fprintf(stderr,"Error in the mode of the binned output \n");
		exit(999);
	}
	slicepixel = 2*(simpar.nx/2+1)*simpar.ny;
	bslicepixel = (simpar.nx/4)*(simpar.ny/4);
	if(simpar.myid==0){
		int isizeof;
		float *binnedden;
		int icount,jj,binslice;
		binslice = (simpar.nx/4) *(simpar.ny/4);
		binnedden = (float*)Malloc(sizeof(float)*binslice,PPTR(binnedden));
		isizeof = sizeof(SimParameters);
		fp = fopen(filename,"w");
		/*
		fwrite(&isizeof,sizeof(int),1,fp);
		fwrite(&simpar,sizeof(SimParameters),1,fp);
		*/
		nnx = simpar.nx/4;
		nny = simpar.ny/4;
		nnz = simpar.nz/4;
		fwrite(&(nnx),sizeof(int),1,fp);
		fwrite(&(nny),sizeof(int),1,fp);
		fwrite(&(nnz),sizeof(int),1,fp);
		icount = ii = 0;
		for(i=0;i<binslice;i++) binnedden[i] = 0;
		for(k=0;k<simpar.local_nz;k++){
			for(j=0;j<simpar.ny;j++) for(i=0;i<simpar.nx;i++){
				binnedden[(i/4)+(simpar.nx/4)*(j/4)] += 
					fftwmesh[i+(2*(simpar.nx/2+1))*j+slicepixel*k]/64.;
			}
			icount ++;
			if(icount%4 ==0){
				fwrite(binnedden,sizeof(float),binslice,fp);
				for(i=0;i<binslice;i++) binnedden[i]= 0;
				icount = 0;
			}
		}
		tmpden = (float*)Malloc(sizeof(float)*slicepixel,PPTR(tmpden));
		for(ii=1;ii<simpar.nid;ii++){
			MPI_Recv(&niter,1,MPI_INT,ii,ii,MPI_COMM_WORLD,&status);
			for(jj=0;jj<niter;jj++){
				MPI_Recv(tmpden,slicepixel,MPI_FLOAT,ii,jj,MPI_COMM_WORLD,
						&status);
				for(j=0;j<simpar.ny;j++) for(i=0;i<simpar.nx;i++){
					binnedden[(i/4)+(simpar.nx/4)*(j/4)] += tmpden[i+(2*(simpar.nx/2+1))*j]/64.;
				}
				icount ++;
				if(icount%4 ==0){
					fwrite(binnedden,sizeof(float),binslice,fp);
					for(i=0;i<binslice;i++) binnedden[i] = 0;
					icount = 0;
				}
			}
		}
		fclose(fp);
		Free(tmpden);
		Free(binnedden);
	}
	else {
		MPI_Send(&(simpar.local_nz),1,MPI_INT,0,simpar.myid,MPI_COMM_WORLD);
		for(i=0;i<simpar.local_nz;i++)
			MPI_Send(fftwmesh+i*slicepixel,slicepixel,MPI_FLOAT,0,i,
					MPI_COMM_WORLD);		
	}
}

void SaveBinnedData(float *fftwmesh){
	float *vel,*fftwvel;
	float den;
	int i,mode;
	int nzheight;
	nzheight = ceil(lzceil) - (int)lzfloor;

	/*
	SaveBinnedDen(fftwmesh,0);
	*/
	{
		vel = (float *) Malloc(((nzheight+4)*(2*(simpar.nx/2+1))*simpar.ny)*
				sizeof(float), PPTR(vel));
		if(simpar.local_nz > 0) fftwvel = (float *) Malloc(
				((simpar.local_nz+0)*(2*(simpar.nx/2+1))*simpar.ny)*
				sizeof(float),PPTR(fftwvel));
		else fftwvel = (float *) Malloc(sizeof(float),PPTR(fftwvel));
	}
	for(mode=1;mode<=3;mode++){
		VarPM_veltsc2(vel,simpar.nx,simpar.ny,simpar.nz, lzfloor,lzceil, pmparticles,
				nsizepmparticletype,nxpos, np,pmas,
                fftwvel, simpar.local_nz,local_z_start,mode);
		for(i=0;i<simpar.local_nz*(2*(simpar.nx/2+1))*simpar.ny;i++){
			den = fftwmesh[i] + 1; /* back to the unnormalized density */
			if(den >0) fftwvel[i] = fftwvel[i]/den; 
			else fftwvel[i] = 0;
		}
		SaveBinnedDen(fftwvel,mode);
	}
	Free(fftwvel);
	Free(vel);
}


void tree2pmtype(long nplong){
	long i;
	for(i=0;i<nplong;i++){
#ifdef MULTIMASS
		pmparticles[i].mass=treeparticles[i].mass ; 
#endif
#ifdef INDTIME
		pmparticles[i].tflag=treeparticles[i].tflag ;
#endif
		pmparticles[i].x = treeparticles[i].x;
		pmparticles[i].y = treeparticles[i].y;
		pmparticles[i].z = treeparticles[i].z;
		pmparticles[i].vx = treeparticles[i].vx;
		pmparticles[i].vy = treeparticles[i].vy;
		pmparticles[i].vz = treeparticles[i].vz;
#ifdef INDEX
		pmparticles[i].indx = treeparticles[i].indx;
#endif
	}
	pmparticles = (pmparticletype*)Realloc(pmparticles,sizeof(pmparticletype)*nplong);
}
void pm2treetype(long nplong){
	long i;
	pmparticles = (pmparticletype*)Realloc(pmparticles,sizeof(treeparticletype)*nplong);
	for(i=nplong-1;i>=0;i--){
#ifdef INDEX
		treeparticles[i].indx = pmparticles[i].indx;
#endif
		treeparticles[i].vz = pmparticles[i].vz;
		treeparticles[i].vy = pmparticles[i].vy;
		treeparticles[i].vx = pmparticles[i].vx;
		treeparticles[i].z = pmparticles[i].z;
		treeparticles[i].y = pmparticles[i].y;
		treeparticles[i].x = pmparticles[i].x;
#ifdef INDTIME
		treeparticles[i].tflag = pmparticles[i].tflag;
#endif
#ifdef MULTIMASS
		treeparticles[i].mass = pmparticles[i].mass;
#endif

	}
}
void Savexzslice(float *sden){
	char filename[100];
	FILE *fp;
	float *tden,*allden;
	int i,j,k,ii,mx,slocal_nz;
	int offset;
	MPI_Status status;
	sprintf(filename,"xzslice.%.5d",simpar.stepcount);
	mx = 2*(simpar.nx/2+1);
	{
		if(simpar.local_nz)
			tden = (float*)Malloc(sizeof(float)*simpar.nx*simpar.local_nz,PPTR(tden));
		for(i=0;i<simpar.nx*simpar.local_nz;i++) tden[i] = 0.;
		for(k=0;k<simpar.local_nz;k++) for(j=0;j<3;j++)for(i=0;i<simpar.nx;i++)
			tden[i+simpar.nx*k] += sden[i+mx*(j+simpar.ny*k)];

	}
	if(simpar.myid==0){
		allden=(float*) Malloc(sizeof(float)*simpar.nx*simpar.nz,PPTR(allden));
		for(i=0;i<simpar.nx*simpar.local_nz;i++) allden[i] = tden[i];
		offset=simpar.local_nz*simpar.nx;
		for(i=1;i<simpar.nid;i++){
			MPI_Recv(&slocal_nz,1,MPI_INT,i,0,MPI_COMM_WORLD,&status);
			MPI_Recv(allden+offset,
					slocal_nz*simpar.nx,MPI_FLOAT,i,1,MPI_COMM_WORLD,&status);
			offset += slocal_nz*simpar.nx;
		}
		fp=fopen(filename,"w");
		fwrite(&simpar.nx,sizeof(int),1,fp);
		fwrite(&simpar.nz,sizeof(int),1,fp);
		fwrite(allden,sizeof(float),simpar.nx*simpar.nz,fp);
		fclose(fp);
		Free(allden);
	}
	else {
		MPI_Send(&(simpar.local_nz),1,MPI_INT,0,0,MPI_COMM_WORLD);
		MPI_Send(tden,simpar.local_nz*simpar.nx,MPI_FLOAT,0,1,MPI_COMM_WORLD);
	}
	if(simpar.local_nz) Free(tden);
}
void Savexyslice(float *sden,int nx,int ny,int nz,int myid,int stepcount){
	char filename[100];
	float *tden;
	FILE *fp;
	int i,j,k,mx;
	int npix;
	if(myid == 0){
		mx = 2*(nx/2+1);
		sprintf(filename,"xyslice.%.5d",stepcount);
		fp = fopen(filename,"w");
		fwrite(&nx,sizeof(int),1,fp);
		fwrite(&ny,sizeof(int),1,fp);
		tden = (float *)Malloc(sizeof(float)*nx*ny,PPTR(tden));
		/*
		for(i=0;i<nx*ny;i++) tden[i] = 0.;
		for(k=0;k<3;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++)
				tden[i+nx*j] += sden[i+mx*(j+ny*k)];
		fwrite(tden,sizeof(float),nx*ny,fp);
		*/
		for(k=0;k<4;k++) {
			for(j=0;j<ny;j++) for(i=0;i<nx;i++){
				npix = i+mx*(j+ny*k);
				tden[i+nx*j] = sden[npix];
			}
			fwrite(tden,sizeof(float),nx*ny,fp);
		}
		fflush(fp);
		fclose(fp);
		Free(tden);
	}
}
/*
void MAIN_(int argc, char *argv[]){
int MAIN__(int argc, char *argv[]){
int MAIN_(int argc, char *argv[]){
int main(int argc, char *argv[]){
*/
#ifdef INTEL
int main(int argc, char *argv[]){
#elif PGCC
int MAIN_(int argc, char *argv[]){
#else
#error You must type the compiler for the fortran and C compatibility setting for main.
#endif
	int i,j,k;
	int nx, ny, nz, local_nz;
	int nspace;
	int ierror=1, iwrite=0, icont;
	char c;
	int memflag, InitializeBigMemory();
	FILE *simfile;
	rfftwnd_mpi_plan plan,iplan;

#ifdef LAM
	printf("%d\n",MPI_Init(NULL,NULL));
#else
	MPI_Init(&argc,&argv);
#endif
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);


	if(myid==0){
		printf(" MPI SETTING \n");
		printf(" Total %d nodes are allocated\n",nid);
	}

	if( argc == 2 ) {
		int sflag=1;
		if(myid==0) simfile = fopen(argv[1],"r");
		if(myid==0 && !simfile) sflag = 0;
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
		fftwinit_(&simpar.nx,&simpar.ny,&simpar.nz,&local_nz,&local_z_start);
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

	zheight = (float) simpar.local_nz;
	zstart = (float) local_z_start;
#ifdef OLD_FDA4
	fda4sup_(&simpar.nx,&simpar.ny,&simpar.nz, &simpar.local_nz);
#endif
	/* PMSEED */
	if(icont ==0){
		float *fx,*fy,*fz,*denpmseedpad;
		float *ffx,*ffy,*ffz;
		{
			long densize;
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

			printf("P%d: %g %g nzzzz: %d %d %g\n",myid,zstart,zstart+zheight,
					nzparticlewidth,nzintoffset,zstartpos);

			densize = 2*(simpar.nx/2+1)*simpar.ny*(simpar.local_nz+0);
			if(densize==0) densize=10;
#ifdef PMSEEDFORCE
			fx = (float*) Malloc(sizeof(float)*densize,PPTR(fx));
			fy = (float*) Malloc(sizeof(float)*densize,PPTR(fy));
			fz = (float*) Malloc(sizeof(float)*densize,PPTR(fz));
			denpmseedpad = (float *) Malloc(2*simpar.ny*simpar.nz*sizeof(float),PPTR(denpmseedpad));
			den=(float*) Malloc(sizeof(float)*densize, PPTR(den));
			pmseedforce_(&np,pmparticles,&nsizepmparticletype,&nxpos,den,
					&simpar.omep,&simpar.omeplam,&simpar.ken,
					denpmseedpad,&zstart,&zheight,&vamp,&simpar.omei, fx,fy,fz);
			Free(den); /* Now free den in the pmseedforce */
#else
			den=(float*) Malloc(sizeof(float)*densize, PPTR(den));
			denpmseedpad = (float *) Malloc(2*simpar.ny*simpar.nz*sizeof(float),PPTR(denpmseedpad));
			pmseed_(&np,pmparticles,&nsizepmparticletype,&nxpos,den,
					&simpar.omep,&simpar.omeplam,&simpar.ken,
					denpmseedpad,&zstart,&zheight,&vamp,&simpar.omei);
#endif
			Free(denpmseedpad);

			MPI_Barrier(MPI_COMM_WORLD);
			if(0){
#define fx(i,j,k) fx[(long)(i-1)+mx*((long)(j-1)*(long)ny*(k-1))]
#define fy(i,j,k) fy[(long)(i-1)+mx*((long)(j-1)*(long)ny*(k-1))]
#define fz(i,j,k) fz[(long)(i-1)+mx*((long)(j-1)*(long)ny*(k-1))]
					long mx;
					mx = 2*(simpar.nx/2+1);
					printf("+P%d has fx,fy,fz %g %g %g\n",myid,fx(simpar.nx,simpar.ny,5),fy(simpar.nx,simpar.ny,5),fz(simpar.nx,simpar.ny,5));
#undef fx
#undef fy
#undef fz
			}

			zfinal = zstart+zheight;
#ifdef PMSEEDFORCE
			ffx = (float*)Malloc(sizeof(float)*2*(simpar.nx/2+1)*(long)simpar.ny*
					(long)(ceil(zstart+zheight)-(int)zstart+1),PPTR(ffx));
			MPI_Barrier(MPI_COMM_WORLD);
			if(myid==0) printf("Now passing fda4 mesh for fx\n");
			adVarPM_fda4_mesh(fx,simpar.local_nz,local_z_start,ffx,zstart,zfinal,simpar.nx,simpar.ny,simpar.nz);
			Free(fx);
			ffy = (float*)Malloc(sizeof(float)*2*(simpar.nx/2+1)*(long)simpar.ny*
					(long)(ceil(zstart+zheight)-(int)zstart+1),PPTR(ffy));
			MPI_Barrier(MPI_COMM_WORLD);
			if(myid==0) printf("Now passing fda4 mesh for fy\n");
			adVarPM_fda4_mesh(fy,simpar.local_nz,local_z_start,ffy,zstart,zfinal,simpar.nx,simpar.ny,simpar.nz);
			Free(fy);
			ffz = (float*)Malloc(sizeof(float)*2*(simpar.nx/2+1)*(long)simpar.ny*
					(long)(ceil(zstart+zheight)-(int)zstart+1),PPTR(ffz));
			MPI_Barrier(MPI_COMM_WORLD);
			if(myid==0) printf("Now passing fda4 mesh for fz\n");
			adVarPM_fda4_mesh(fz,simpar.local_nz,local_z_start,ffz,zstart,zfinal,simpar.nx,simpar.ny,simpar.nz);
			Free(fz);
#else
			pvmesh = (float *) Malloc(sizeof(float)*2*(simpar.nx/2+1)*
							(long)simpar.ny*(ceil(zstart+zheight)-(int)zstart+4), PPTR(pvmesh));
			adVarPM_fda4_mesh(den,simpar.local_nz,local_z_start,pvmesh,zstart,zfinal,simpar.nx,simpar.ny,simpar.nz);
			Free(den);
			paddingpotentialmesh(pvmesh,izstart,izwidth,izwidth,simpar.nx,simpar.ny,simpar.nz,-2,2);
#endif
			MPI_Barrier(MPI_COMM_WORLD);
			if(0){
#define ffx(i,j,k) ffx[(long)(i-1)+mx*((long)(j-1)*(long)ny*(k-1))]
#define ffy(i,j,k) ffy[(long)(i-1)+mx*((long)(j-1)*(long)ny*(k-1))]
#define ffz(i,j,k) ffz[(long)(i-1)+mx*((long)(j-1)*(long)ny*(k-1))]
					long mx;
					mx = 2*(simpar.nx/2+1);
					printf("*P%d has fx,fy,fz %g %g %g\n",myid,ffx(simpar.nx,simpar.ny,5),ffy(simpar.nx,simpar.ny,5),ffz(simpar.nx,simpar.ny,5));
#undef ffx
#undef ffy
#undef ffz
			}
			{
				long nmax;
				nmax = freespace()/sizeof(pmparticletype) - 32L;
				pmparticles=(pmparticletype*) Malloc(sizeof(pmparticletype)*nmax,
						PPTR(pmparticles));
				if(myid==0) printf("Now P%d has maximum %ld particles\n",myid,nmax);
			}
#ifdef PMSEEDFORCE
			setboundaryparticlesforce_(&np,pmparticles,&nsizepmparticletype,&nxpos,
					&nzparticlewidth,&zstartpos,&nzintoffset,&simpar.omep,&simpar.omeplam,
					&simpar.ken,&zstart,&zheight,&vamp,
					&simpar.nspace,&simpar.omei,&simpar.amax,ffx,ffy,ffz);
			if(myid==0) printf("Now passing boundaryparitlces\n");
#else
			setboundaryparticles_(&np,pmparticles,&nsizepmparticletype,&nxpos,
					&nzparticlewidth,&zstartpos,&nzintoffset,&simpar.omep,&simpar.omeplam,
					&simpar.ken,
					pvmesh,&zstart,&zheight,&vamp,&simpar.nspace,&simpar.omei,&simpar.amax);
#endif

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
				nmax = freespace()/sizeof(pmparticletype) - 32L;
				pmparticles=(pmparticletype*) Realloc(pmparticles,
						sizeof(pmparticletype)*(np+nmax));
				printf("P%d has namx = %ld\n",simpar.myid,nmax);
			}
#ifdef PMSEEDFORCE
			{
				pmvrparticletype *pvr;
				long nnp;
				pvr = (pmvrparticletype *)((pmparticletype*)(pmparticles+np));
				nnp = 0;
				MPI_Barrier(MPI_COMM_WORLD);
				if(myid==0) printf("Before getinnerparticlesforce");
				getinnerparticlesforce_(&nnp,pvr,&nsizepmparticletype,&nxpos,
						&nzparticlewidth,&zstartpos,&nzintoffset,&simpar.omep,&simpar.omeplam,
						&simpar.ken,&zstart,&zheight,&vamp,&simpar.nspace,
						&simpar.omei,&simpar.amax, ffx,ffy,ffz);
				if(myid==0) printf("P0 completed ");
				MPI_Barrier(MPI_COMM_WORLD);
				if(myid==0) printf("after getinnerparticlesforce\n");
				Free(ffz);Free(ffy);Free(ffx);
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
				setinnerparticlesforce_(&np,pmparticles,&nsizepmparticletype,&nxpos,
						&nzparticlewidth,&zstartpos,&nzintoffset,&simpar.omep,&simpar.omeplam,
						&simpar.ken,&zstart,&zheight,&vamp,&simpar.nspace,
						&simpar.omei,&simpar.amax);
				if(myid==0) printf("P0 completed ");
				MPI_Barrier(MPI_COMM_WORLD);
				if(myid==0) printf("after setinnerparticlesforce\n");
			}
#else
			setinnerparticles_(&np,pmparticles,&nsizepmparticletype,&nxpos,
					&nzparticlewidth,&zstartpos,&nzintoffset,&simpar.omep,&simpar.omeplam,
					&simpar.ken,pvmesh,&zstart,&zheight,&vamp,
					&simpar.nspace,&simpar.omei,&simpar.amax);
			Free(pvmesh);
#endif
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
		int nzheight;
		{ if(icont !=1) Savingalldata(); icont = 0; }
		TIMER_START(51);
		if(0){ /* This is to make projected maps for simulated observers */
			int mkanimate(pmparticletype *, int);
			mkanimate(pmparticles,np);

		}
		if(0){/* extract lightcone data if needed */
			void observerextmain();
			observerextmain();
		}
		{
			halfstep.first = halfstep.second=0;
			if((iflagsyncpdata=flagsyncpdata(simpar.amax,simpar.anow,simpar.astep))){
				halfstep.first = 1; halfstep.second = 0;
			}
			if( (iflagPreFoF=flagPreFoF(simpar.amax,simpar.anow,simpar.astep))){
				halfstep.first = 1; halfstep.second = 0;
			}
			if(!iflagsyncpdata && !iflagPreFoF){
				halfstep.first = halfstep.second = 0;
			}
		}
halfevolution:
		{/* This is PM */
			void pmmain(int *, int *);
			pmmain(&nowflag,&now);
		}

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
		if(0){/* Dump lightcone data if needed */
			void observersavemain();
			observersavemain();
		}

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

int ReadSimulationParameters(FILE *simfile,int *icont) {
        char filename[180];
	FILE *rvfile;
	int abortflag;

	if(myid==0){
		/*
		fscanf(simfile,"%f %f\n",&(simpar.boxsize),&(simpar.hubble));
		fscanf(simfile,"%f %f %f %f %f %f\n",
   	        &(simpar.npow),&(simpar.omep),&(simpar.omepb), &(simpar.omeplam),
			   &(simpar.bias),&(simpar.rsmooth));
		fscanf(simfile,"%d %d %d %d\n",&(simpar.nx), &(simpar.ny), &(simpar.nz),&(simpar.nspace));
		simpar.mx = simpar.nx/simpar.nspace;
		simpar.my = simpar.ny/simpar.nspace;
		simpar.mz = simpar.nz/simpar.nspace;
		simpar.mxmy = simpar.mx*simpar.my;
		simpar.lnx = (double)simpar.nx;
		simpar.lny = (double)simpar.ny;
		simpar.lnz = (double)simpar.nz;
		fscanf(simfile,"%f %f %f\n",&(simpar.sphere_radius),&(simpar.particle_radius),&(simpar.theta));
		L = simpar.nx;
		L2 = L*L;
		fscanf(simfile,"%f %f %f\n",&(simpar.zinit),&(simpar.astep),&(simpar.anow));
		simpar.amax = 1.0+simpar.zinit;
		fscanf(simfile,"%d %d %d\n",&(simpar.nstep),&(simpar.stepnum),&(simpar.nskip));
   		fscanf(simfile,"%d\n",&(simpar.iseed));
		simpar.rth = simpar.nx/(simpar.boxsize)*8.0;
		fscanf(simfile,"%s %s\n",simpar.rvfilename, simpar.rvprefix);
		fscanf(simfile,"%f %f %f %f\n",
			&(simpar.ken), &(simpar.ktot), &(simpar.const0), &(simpar.poten0));
		if(fscanf(simfile,"%d\n",&(simpar.powreadflag))){
			if(simpar.powreadflag) fscanf(simfile,"%s",simpar.powfilename);
			if(simpar.powreadflag==2) fscanf(simfile,"%s\n",simpar.inpapkfilename);
		}
		else {
			simpar.powreadflag = 0;
			sprintf(simpar.powfilename,"NULL");
		}
   		fclose(simfile);
		simpar.omei = simpar.omep*pow(1.+simpar.zinit,3.)/
			(simpar.omep*pow(1.+simpar.zinit,3.)+simpar.omeplam+
			 (1.-simpar.omep-simpar.omeplam)*pow(1+simpar.zinit,2.));
			 */
		{
			SimParameters  read_sim_parameter_file(FILE *);
			void determine_mpi_misc_param(SimParameters *);
			void mk_default_param(SimParameters *, char *);
			char cosmology[190]="WMAP5";
			mk_default_param(&simpar,cosmology);
			simpar = read_sim_parameter_file(simfile);
			determine_mpi_misc_param(&simpar);
			simpar.stepnum = simpar.stepcount;
		}
	}

	MPI_Bcast(&simpar,sizeof(SimParameters),MPI_BYTE,0,MPI_COMM_WORLD);
	simpar.myid = myid; simpar.nid = nid;


/* Ok now output to stdout to verify all is well */
	if( myid ==0 ) {
		fprintf(stdout,"*** INPUT PARAMETER LIST ***\n\n");
		fprintf(stdout,"*** Cosmological Parameters ***\n");
		fprintf(stdout,"Box size = %g hinv_Mpc   Hubble_parameter = %g\n",
			(simpar.boxsize), (simpar.hubble));
		fprintf(stdout,"npow = %g Omega = %g OmegaBaryon = %g Lambda = %g\nBias_8Mpc = %g Rsmooth = %g\n\n", 
				(simpar.npow),simpar.omep,simpar.omepb,simpar.omeplam,simpar.bias,simpar.rsmooth);
		fprintf(stdout,"\n*** Mesh Parameters ***\n");
		fprintf(stdout,"Nxmesh = %d Nymesh = %d Nzmesh = %d Spacing = %d\n\n",
				simpar.nx,simpar.ny,simpar.nz,simpar.nspace);
		fprintf(stdout,"\n*** Tree Parameters ***\n");
		fprintf(stdout,"SPHERE_RADIUS = %g Theta = %g\n\n",simpar.sphere_radius,simpar.theta);
		fprintf(stdout,"\n*** Time Parameters ***\n");
		fprintf(stdout,"Initial Redshift = %g Delta_a = %g Current a = %g\n\n",simpar.zinit,simpar.astep,
				simpar.anow);
		fprintf(stdout,"\n*** Stepping Parameters ***\n");
		fprintf(stdout,"Nstep = %d Stepnum = %d Nskip = %d\n", 
			simpar.nstep,simpar.stepnum,simpar.nskip);
		fprintf(stdout,"Random Seed = %d\n\n",simpar.iseed);
		fprintf(stdout,"Initial RV File %s  RV File Prefix %s\n\n",
			simpar.rvfilename,simpar.rvprefix);
		if(simpar.powreadflag) fprintf(stdout,"Initial Power File %s \n\n", simpar.powfilename);
		if(simpar.powreadflag==2) fprintf(stdout,"Initial Power Ascii File %s \n\n", simpar.inpapkfilename);
	}
	fflush(stdout);

/* Check to see if the initial conditions file exists - if not then default to g   generating your own initial conditions */

	abortflag = 0;
	if( myid == 0 ) {
/* if strcmp() == 0 then string comparison is successful */
		if( strcmp(simpar.rvfilename,"INITIAL") == 0 ) {
			fprintf(stderr,"Generating Initial Conditions for run....\n");
			*icont = 0;
			simpar.anow = 1.0;
		}
		else {
            sprintf(filename,"%s%.5d",simpar.rvfilename,0);
			rvfile = fopen(filename,"r");
			if( rvfile == NULL ) {
				fprintf(stderr,"Can't open initial rvfile %s - exiting\n",	
					simpar.rvfilename);
				abortflag = 1;
			}	
			else { 
				fprintf(stderr,"Successfully opened rvfile %s - starting simulation\n",simpar.rvfilename);
				fclose(rvfile); 
				*icont = 1;
			}
		}
	}
	MPI_Bcast(&abortflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if( abortflag ) {
		fprintf(stderr,"Aborting can't find initial rvfile\n");
		MPI_Finalize();
		exit(0);
	}

	MPI_Bcast(icont, 1, MPI_INT, 0, MPI_COMM_WORLD);
	return 0;
		 
}
			 
/* Initialization for power spectrum and initial condition generator */
			 
void getpowerparams_(nspace,size,amax,astep,npow,smooth,h,rth,bias,iseed,
		omepb,powflag,powfile,inpapkfile)
int *nspace;
float *size, *amax, *astep;
float *npow;
int *iseed;
float *smooth, *h, *rth, *bias;
float *omepb;
char *powfile,*inpapkfile;
int *powflag;
{
	int i;		 
     *nspace = simpar.nspace;
     *size = simpar.boxsize;
     *amax = simpar.amax;
     *astep = simpar.astep;
     *npow = simpar.npow;
     *iseed = simpar.iseed;
     *smooth = simpar.rsmooth;
     *h = simpar.hubble;
     *rth = simpar.rth;
     *bias = simpar.bias;
	 *omepb = simpar.omepb;
	 *powflag = simpar.powreadflag;
	 for(i=0;i<190;i++) powfile[i] = '\0';
	 if(simpar.powreadflag > 0) sprintf(powfile,"%s",simpar.powfilename);
	 if(simpar.powreadflag ==2) sprintf(inpapkfile,"%s",simpar.inpapkfilename);

}
