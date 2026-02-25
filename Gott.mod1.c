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
		int nbin = 1;
		int *Local_z_start,*nzrank,*Local_nz;
		long ishift,jkshift,kshift;
		float *bden,*rbden,*saveslice;
		MPI_Status mpistatus,rstatus;


		Local_z_start = (int *)Malloc(sizeof(int)*nid,PPTR(Local_z_start));
		Local_nz = (int *)Malloc(sizeof(int)*nid,PPTR(Local_nz));
		nzrank = (int *)Malloc(sizeof(int)*nz,PPTR(nzrank));
		if(myid ==0){
			Local_z_start[0] = local_z_start;
			Local_nz[0] = local_nz;
			for(i=1;i<nid;i++){
				/*
				MPI_Probe(MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&rstatus);
				MPI_Recv(&(Local_z_start[rstatus.MPI_SOURCE]),
							1,MPI_INT,rstatus.MPI_SOURCE,0,
							MPI_COMM_WORLD,&mpistatus);
				*/
				MPI_Recv(Local_z_start+i, 1,MPI_INT,i,0, MPI_COMM_WORLD,&mpistatus);
				MPI_Recv(Local_nz+i, 1,MPI_INT,i,1, MPI_COMM_WORLD,&mpistatus);
			}
		}
		else {
			/* if local_nz=0 then local-z-start=0 in fftw*/ 			
			if(local_nz ==0 ) local_z_start = nz;
			MPI_Send(&local_z_start,1,MPI_INT,0,0,MPI_COMM_WORLD);
			MPI_Send(&local_nz,1,MPI_INT,0,1,MPI_COMM_WORLD);
		}
		if(myid==0){
			for(i=0;i<nz;i++){
				for(j=0;j<nid;j++){
					if(i >= Local_z_start[j]  && i< Local_z_start[j] + Local_nz[j] ){
						nzrank[i] = j;
						break;
					}
				}
			}
		}
		/*
		if(myid==0) for (i=0;i<nid;i++) printf("P%d has %d local_z_start %d : %d\n",myid,i,
						Local_z_start[i], Local_nz[i]);
		MPI_Finalize();exit(99);
//		MPI_Bcast(Local_z_start,nid,MPI_INT,0,MPI_COMM_WORLD);
//		*/

		mx = nx/nbin; my = ny/nbin; mz = nz/nbin;

		if(nz%nbin !=0) {
			fprintf(stderr,"Error: indivisible size of mesh %d : %d\n",nz,nbin);
			return;
		}
		nnx = 2*(nx/2+1);
		saveslice = (float *)Malloc(sizeof(float)*mx*my,PPTR(saveslice));
		/*
		MPI_Barrier(MPI_COMM_WORLD); printf("P%d has reached here %d %d %d\n",myid,mx,my,mz);
		*/
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
			for(k=0;k<nz;k++){
				if(k%nbin ==0) for(i=0;i<mx*my;i++) saveslice[i] = 0;
				if(nzrank[k] != myid){
					MPI_Recv(rbden,mx*my,MPI_FLOAT,nzrank[k],0,MPI_COMM_WORLD,
							&mpistatus);
					for(i=0;i<mx*my;i++) saveslice[i] += rbden[i];
				}
				else {
					int ix,iy;
					long jkshift,kshift,mxiy;
					kshift = nnx*ny*k;
					for(j=0;j<ny;j++){
						iy = j/nbin;
						mxiy = mx*iy;
						jkshift = nnx*j + kshift;
						for(i=0;i<nx;i++){
							ix = i/nbin;
							saveslice[ix+mxiy] += den[i+jkshift];
						}
					}
				}

				if(k%nbin == nbin-1) {
					float inv3nbin = 1./(nbin*nbin*nbin);
					for(i=0;i<mx*my;i++) saveslice[i] = saveslice[i]*inv3nbin;
					fwrite(saveslice,sizeof(float),mx*my,wp);
				}
			}
			fclose(wp);
			Free(rbden);
		}
		else {
			if(local_nz>0) rbden = (float *)Malloc(sizeof(float)*mx*my*local_nz,PPTR(rbden));
			for(i=0;i<mx*my*local_nz;i++) rbden[i]  = 0;
			for(k=0;k<local_nz;k++){
				int ix,iy;
				long jkshift,kshift,mxiy,mxmyk;
				mxmyk = mx*my*k;
				kshift = nnx*ny*k;
				for(j=0;j<ny;j++){
					iy = j/nbin;
					mxiy = mx*iy;
					jkshift = nnx*j+kshift;
					for(i=0;i<nx;i++){
						ix = i/nbin;
						rbden[ix+mxiy+mxmyk] += den[i+jkshift];
					}
				}
			}
			for(k=0;k<local_nz;k++){
				MPI_Send(rbden+mx*my*k,mx*my,MPI_FLOAT,0,0,MPI_COMM_WORLD);
			}
			if(local_nz>0) Free(rbden);
		}
		Free(saveslice);
		Free(nzrank);
		Free(Local_nz);
		Free(Local_z_start);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0) printf("Dump binned Grid Survey\n");
}
