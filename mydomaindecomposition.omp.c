#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<mpi.h>
#include<omp.h>

#include "pmheader.h"
#include "Memory.h"
#ifndef _OPENMP
int nullfct0(), nullfct1();
#define omp_get_thread_num() nullfct0()
#define omp_get_num_threads() nullfct1()
#endif


#define NBIN 1000
#define NBINP 1001

void migrate(size_t *,int ,int,int,int,float,float);
int *BIN, *TBIN;
int *OBIN;
static int myid,nid;
float *tzmin,*tzmax;
float *zbound;
typedef struct Binbound{
	float zmin,zmax;
} Binbound;
Binbound *binbound;
int onpmin,onpmax;
int npmin,npmax;
extern float maintreetime;

int mydomaindecomposition(pmparticletype *pmparticles, int np,long tnp,
		float zmin,float zmax, float *newzmin, float *newzmax,int nx,int ny,int nz,int nspace){
	int i,j,k;
	int nbin;
	int npav;
	float frac;
	float zwidth;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);

	zwidth = zmax-zmin;
	frac = (float)(NBIN)/zwidth;
	npav = tnp/nid;
	MPI_Reduce(&np,&onpmin,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
	MPI_Reduce(&np,&onpmax,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);

	int nthreads=1;
#ifdef _OPENMP
#pragma omp parallel
	{
#pragma omp single
		{
			nthreads = omp_get_num_threads();
		}
	}
#endif
	TBIN = (int *)Malloc(sizeof(int)*NBINP*nthreads,PPTR(TBIN));
	for(i=0;i<NBINP*nthreads;i++) TBIN[i] = 0;

#ifdef _OPENMP
#pragma omp parallel private(BIN,i,nbin)
#endif
	{
		int ithread = omp_get_thread_num();
		BIN = TBIN + ithread*NBINP;
		long  istart, ifinal, iwidth;
		iwidth = (np+nthreads-1)/nthreads;
		istart = iwidth *ithread;
		ifinal = istart + iwidth;
		if(ifinal > np) ifinal = np;

		for(i=istart;i<ifinal;i++){
			nbin = (ZofP((pmparticles+i))-zmin)*frac;
			BIN[nbin]++;
		}
	}
	for(j=0;j<NBINP;j++){
		for(i=1;i<nthreads;i++){
			TBIN[j] += TBIN[j+i*NBINP];
		}
	}

	BIN = TBIN;


	binbound = (Binbound *)Malloc(sizeof(Binbound)*NBIN*nid,PPTR(binbound));
	for(i=0;i<NBIN;i++){
		binbound[i].zmin = zmin + i/frac;
		binbound[i].zmax = zmin + (i+1.)/frac;
	}
	binbound[NBIN-1].zmax = zmax; /* confirming the upper boundary value */
//	printf("P%d : np=%d npav=%d zmin=%g zmax=%g\n",myid,np,npav,zmin,zmax);
	if(myid==0){
		int ncount,mid;
		OBIN = (int *)Malloc(sizeof(int)*NBIN*nid,PPTR(OBIN));
		zbound = (float*)Malloc(sizeof(float)*nid,PPTR(zbound));
		for(i=0;i<NBIN;i++) OBIN[i] = BIN[i];
		for(i=1;i<nid;i++){
			MPI_Recv(binbound+i*NBIN,NBIN*sizeof(Binbound),MPI_BYTE,i,0,
					MPI_COMM_WORLD,&status);
			MPI_Recv(OBIN+i*NBIN,NBIN,MPI_INT,i,2,MPI_COMM_WORLD,&status);
		}
		ncount =  0;
		mid = 0;
		for(i=0;i<nid*NBIN;i++){
			ncount += OBIN[i];
			if(ncount >= npav){
				zbound[mid] = (float)(npav-(ncount-OBIN[i]))/(float)OBIN[i]*
						(binbound[i].zmax-binbound[i].zmin)+binbound[i].zmin;
				mid++;
				ncount = OBIN[i]-(npav-(ncount-OBIN[i]));
			}
		}
		zbound[nid-1] = nz;
		for(i=1;i<nid;i++){
			zmin = zbound[i-1];
			zmax = zbound[i];
			MPI_Send(&zmin,1,MPI_FLOAT,i,3,MPI_COMM_WORLD);
			MPI_Send(&zmax,1,MPI_FLOAT,i,4,MPI_COMM_WORLD);
		}
		zmin = 0; zmax = zbound[0];
		Free(zbound);
		Free(OBIN);
	}
	else {
		MPI_Send(binbound,NBIN*sizeof(Binbound),MPI_BYTE,0,0,MPI_COMM_WORLD);
		MPI_Send(BIN,NBIN,MPI_INT,0,2,MPI_COMM_WORLD);
		MPI_Recv(&zmin,1,MPI_FLOAT,0,3,MPI_COMM_WORLD,&status);
		MPI_Recv(&zmax,1,MPI_FLOAT,0,4,MPI_COMM_WORLD,&status);
	}
	Free(binbound);
	Free(TBIN);
	*newzmin = zmin;
	*newzmax = zmax;
//	if(myid==0) printf("P%d has new zmin=%g zmax=%g\n",myid,zmin,zmax);
	{ 
		size_t sizetnp; 
		float zheight,zstart;
		sizetnp=np;
		zheight = zmax-zmin;
		zstart = zmin;
		migrate(&sizetnp,nx,ny,nz,nspace,zheight,zstart);
		np = sizetnp;
	}
	MPI_Reduce(&np,&npmin,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
	MPI_Reduce(&np,&npmax,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
	if(myid==0){
		float oratio1,oratio2;
		float ratio1,ratio2;
		oratio1 = (float)(npav-onpmin)/(float)npav*100.;
		oratio2 = (float)(onpmax-npav)/(float)npav*100.;
		ratio1 = (float)(npav-npmin)/(float)npav*100.;
		ratio2 = (float)(npmax-npav)/(float)npav*100.;
		printf("Before D.D. inhomogeneity %d %d %g  %g \n",onpmin,onpmax,
				oratio1,oratio2);
		printf("After  D.D. inhomogeneity %d %d %g  %g \n",npmin,npmax,
				ratio1,ratio2);
	}
	return np;
}
