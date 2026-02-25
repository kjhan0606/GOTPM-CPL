#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "Memory.h"
#include "pmheader.h"
#include "lightcone.boss.h"
#define MIN(A,B) ((A) <(B) ? (A):(B))
#define PI 3.1415926535L
double Theta,Phi,DTheta,DPhi;
MPI_Status status;
char surveyfilename[100];
int *freework;
enum boolean { NO = 00,YES = 01};
typedef struct Spherical{
	double r,theta,phi;
} Spherical;
double r1, r2;
double theta,phi,theta1,phi1,theta2;
float omega0,H,lambda0;
float qsimp(float (*)(float),float,float);
float trapzd(float (*)(float), float , float , int );
int ssorttype(const void *, const void *);

typedef struct CHECK{
	int nobs;
	int index;
} CHECK;

typedef struct gottgridtype{
	indextype index;
	int iobs,type;
	float potent;
} gottgridtype;

typedef struct savetype{
	indextype index;
	float potent;
} savetype;

int gottsorttype(const void *a,const void *b){
	gottgridtype *aa,*bb;
	aa = (gottgridtype *)a;
	bb = (gottgridtype *)b;
	if(aa->iobs < bb->iobs) return -1;
	else if(aa->iobs > bb->iobs) return +1;
	else {
		if(aa->type < bb->type) return -1;
		else if(aa->type > bb->type) return +1;
		return 0;
	}
}
#define NPSTEP 1000000
int nnp[3];
/* Particle data type for 108 BOSS survey */
static double rmin2,rmax2;
static gottgridtype *gottgrids;

static int snlcp,maxlcp,tsnlcp,mbox;



double SX0,SY0,SZ0;
void a2comovingpixel(double *, double *, double *,int , float ,float , float ,float , 
		float ,float ,float , float );
void FindcomovingpixelWidth(double *, double *, double *,int , float ,float , float ,float , 
		float ,float ,float , float );
void GottPMeshDump(float *den,int tnstep,int obsidstart, 
				double *sx0, double *sy0, double *sz0,int nobs, float redi){
	int nx,ny,nz, nspace;
	float amax;
	float boxsize,hboxsize;
	float omep,omeplam,hubble,omepb;
	int myid,nid;
	pmparticletype *bp;
	float redshift;
	double r,dr1,dr2,distsq,dist,rot;
	double r1sq,r2sq;
	double xp,yp,zp;
	double xpp,ypp,zpp;
	double rminsq,rmaxsq;
	int nbox,ibox;
	long i,j,k;
	long ix,iy,iz;
	float astep,a,preastep;
	int local_z_start,local_nz;
	int *iobs,mobs;

	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;
	nspace = simpar.nspace;
	boxsize = simpar.boxsize;
	hboxsize = boxsize*0.5;
	omep = simpar.omep;
	omeplam = simpar.omeplam;
	omepb = simpar.omepb;
	hubble = simpar.hubble;
	myid = simpar.myid;
	nid = simpar.nid;
	amax = simpar.amax;
	a = simpar.anow;
	local_z_start = simpar.local_z_start;
	local_nz = simpar.local_nz;
	{
		int mx,my,mz,nnx;
		int jj,kk,ii;
		int nbin = 8;
		int *Local_z_start,*nzrank,*Local_nz;
		long ishift,jkshift,kshift;
		float *bden,*rbden,*saveslice;
		MPI_Status mpistatus,rstatus;


		mx = nx/nbin; my = ny/nbin; mz = nz/nbin;

		if(nz%nbin !=0) {
			fprintf(stderr,"Error: indivisible size of mesh %d : %d\n",nz,nbin);
			return;
		}
		nnx = 2*(nx/2+1);
		if(myid==0) {
			FILE *wp;
			char outfile[190];
			rbden = (float *)Malloc(sizeof(float)*mx*my,PPTR(rbden));
			sprintf(outfile,"BinnedPotent.%.5d",simpar.stepcount);
			wp = fopen(outfile,"w");
			fwrite(&mx,sizeof(int),1,wp);
			fwrite(&my,sizeof(int),1,wp);
			fwrite(&mz,sizeof(int),1,wp);
			printf("P%d has outfilename %s with nz=%d\n",myid,outfile,nz);

			for(i=0;i<mx*my;i++) rbden[i]  = 0;
			for(k=0;k<local_nz;k++){
				int ix,iy;
				long jkshift,kshift,mxiy,mxmyk;
				kshift = nnx*ny*k;
				for(j=0;j<ny;j++){
					iy = j/nbin;
					mxiy = mx*iy;
					jkshift = nnx*j+kshift;
					for(i=0;i<nx;i++){
						ix = i/nbin;
						rbden[ix+mxiy] += den[i+jkshift];
					}
				}
			}
			for(i=0;i<mx*my;i++) rbden[i]  = rbden[i]/(nbin*nbin*nbin);
			fwrite(rbden,sizeof(float),mx*my,wp);

			for(k=1;k<nid;k++){
				MPI_Recv(rbden,mx*my,MPI_FLOAT,k,0,MPI_COMM_WORLD,
						&mpistatus);
				fwrite(rbden,sizeof(float),mx*my,wp);
			}
			Free(rbden);
			fclose(wp);
		}
		else {
			if(local_nz != 5) rbden = (float *)Malloc(sizeof(float)*mx*my,PPTR(rbden));
			for(i=0;i<mx*my;i++) rbden[i]  = 0;
			for(k=0;k<local_nz;k++){
				int ix,iy;
				long jkshift,kshift,mxiy,mxmyk;
				kshift = nnx*ny*k;
				for(j=0;j<ny;j++){
					iy = j/nbin;
					mxiy = mx*iy;
					jkshift = nnx*j+kshift;
					for(i=0;i<nx;i++){
						ix = i/nbin;
						rbden[ix+mxiy] += den[i+jkshift];
					}
				}
			}
			for(i=0;i<mx*my;i++) rbden[i]  = rbden[i]/(nbin*nbin*nbin);
			MPI_Send(rbden,mx*my,MPI_FLOAT,0,0,MPI_COMM_WORLD);
			if(local_nz != 5) Free(rbden);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0) printf("Dump binned Grid Survey\n");
}
