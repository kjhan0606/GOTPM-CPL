#include<stdio.h>
#include<stdlib.h>
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
	long i,j,k,jj;
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

	if(1){
		float redshift;
		redshift = amax/a-1.;
		if(redshift > redi) return;
	}

	preastep = simpar.astep;
	astep = simpar.astep;

	{
		double R1,R2,R3,R4,R;
		double dR1,dR2;
		a2comovingpixel(&R,&dR1,&dR2,nx,amax,a,astep,astep,omep,omeplam,boxsize,hubble);
		r  = R;
		r1 = R-dR1;
		r2 = R+dR2;
		dr1 = dR1;
		dr2 = dR2;
		a2comovingpixel(&R,&dR1,&dR2,nx,amax,a-astep,astep,astep,omep,omeplam,boxsize,hubble);
		rmaxsq = (R+dR2)*(R+dR2);
		a2comovingpixel(&R,&dR1,&dR2,nx,amax,a+astep,astep,astep,omep,omeplam,boxsize,hubble);
		rminsq = (R-dR1)*(R-dR1);
		
	}
	if(r1 <0.) r1 = 0;
	r1sq = r1*r1;
	r2sq = r2*r2;
	rmin2 = r1sq;
	rmax2 = r2sq;

	if(myid==0) printf("Spherical Boss Grid Survey %g %g %g with buffer %g : %g for nobs=%d\n",
					r,dr1,dr2, sqrt(rminsq), sqrt(rmaxsq),nobs);
	snlcp = 0;
	maxlcp = NPSTEP;
	{
		iobs = (int *)Malloc(sizeof(int)*nobs,PPTR(iobs));
		mobs = 0;
		for(i=0;i<nobs;i++){
			float zlow,zhigh;
			zlow = local_z_start-sqrt(rmaxsq);
			zhigh = local_z_start+local_nz+sqrt(rmaxsq);
			if(sz0[i] > zlow && sz0[i] < zhigh) iobs[mobs++] = i;
			else if(sz0[i]+nz > zlow && sz0[i]+nz< zhigh) iobs[mobs++] = i;
			else if(sz0[i]-nz > zlow && sz0[i]-nz< zhigh) iobs[mobs++] = i;
		}		
	}
#ifdef _OPENMP
	{
		CHECK *check;
		int ncheck,nxy;
		long ii;
		long mx,jdx,kdx,jkdx;
		long zshift;
		double zbp;
		mx = 2*(nx/2+1);
		zshift = mx*ny*local_z_start;
		check = (CHECK *)Malloc(sizeof(CHECK)*local_nz*ny*nx,PPTR(check));
		for(ii=0;ii<local_nz*ny*nx;ii++) check[ii].nobs = -1;
		for(k=0;k<local_nz;k++) { 
			zbp = k+local_z_start;
			kdx = mx*ny*k;
#pragma omp parallel for private(jdx,jkdx,i,j,jj,xp,yp,zp,distsq,SX0,SY0,SZ0)
			for(j=0;j<ny;j++){
				jdx = mx*j;
				jkdx = jdx+kdx;
				for(i=0;i<nx;i++){
					for(ii=0;ii<mobs;ii++){ jj = iobs[ii];
						SX0 = sx0[jj]; SY0 = sy0[jj]; SZ0 = sz0[jj];
						xp = i - SX0;
						yp = j - SY0;
						zp = zbp - SZ0;
						BOSSBCondition(xp,hboxsize,boxsize);
						BOSSBCondition(yp,hboxsize,boxsize);
						BOSSBCondition(zp,hboxsize,boxsize);
						distsq = xp*xp+yp*yp+zp*zp;
						if(distsq >= rminsq && distsq < rmaxsq){
							check[kdx+jdx+i].nobs = jj;
							check[kdx+jdx+i].index = kdx+jdx+i;
							goto exclusive_exit2;
						}
					}
exclusive_exit2:;
				}
			}
		}
		ncheck = 0;
		for(ii=0;ii<local_nz*ny*nx;ii++) if(check[ii].nobs>=0) {
				check[ncheck] = check[ii];
				ncheck ++;
		}
		nxy = mx *ny;
		check = Realloc(check,sizeof(CHECK)*ncheck);
		if(ncheck >0) gottgrids = (gottgridtype *)Malloc(sizeof(gottgridtype)*ncheck, PPTR(gottgrids));
#pragma omp parallel for private(ii,SX0,SY0,SZ0,xp,yp,zp,distsq)
		for(ii=0;ii<ncheck;ii++){
			int index;
			index = check[ii].index;
			SX0 = sx0[check[ii].nobs]; SY0 = sy0[check[ii].nobs]; SZ0 = sz0[check[ii].nobs];
			xp = index%mx - SX0;
			yp = (index%(nxy)/mx) - SY0;
			zp = index/(nxy) - SZ0;
			BOSSBCondition(xp,hboxsize,boxsize);
			BOSSBCondition(yp,hboxsize,boxsize);
			BOSSBCondition(zp,hboxsize,boxsize);
			distsq = xp*xp+yp*yp+zp*zp;
			gottgrids[ii].index = index;
			gottgrids[ii].potent = den[index];
			if(distsq>=r1sq && distsq < r2sq) gottgrids[ii].type = 0;
			else if(distsq >=r2sq) gottgrids[ii].type = 1;
			else gottgrids[ii].type = 2;
			gottgrids[ii].iobs = check[ii].nobs;
		}
		Free(check);
		snlcp = ncheck;
	}
#else
	{
		long mx,idx,jdx,kdx,jkdx;
		long zshift;
		double zbp;
		gottgrids = (gottgridtype *)Malloc(sizeof(gottgridtype)*maxlcp, PPTR(gottgrids));
		mx = 2*(nx/2+1);
		zshift = mx*ny*local_z_start;
		for(k=0;k<local_nz;k++) { 
			zbp = k+local_z_start;
			kdx = mx*ny*k;
			for(j=0;j<ny;j++){
				jdx = mx*j;
				jkdx = jdx+kdx;
				for(i=0;i<nx;i++){
					for(jj=0;jj<nobs;jj++){
						SX0 = sx0[jj]; SY0 = sy0[jj]; SZ0 = sz0[jj];
						xp = i - SX0;
						yp = j - SY0;
						zp = zbp - SZ0;
						BOSSBCondition(xp,hboxsize,boxsize);
						BOSSBCondition(yp,hboxsize,boxsize);
						BOSSBCondition(zp,hboxsize,boxsize);
						distsq = xp*xp+yp*yp+zp*zp;
						if(distsq >= rminsq && distsq < rmaxsq){
							indextype iidx;
							iidx = i+jkdx;
							gottgrids[snlcp].index = iidx + zshift;
							gottgrids[snlcp].potent = den[iidx];
							if(distsq>=r1sq && distsq < r2sq) gottgrids[snlcp].type = 0;
							else if(distsq >=r2sq) gottgrids[snlcp].type = 1;
							else gottgrids[snlcp].type = 2;
							gottgrids[snlcp].iobs = jj;
							snlcp ++; 
							if(snlcp >= maxlcp){
								maxlcp += NPSTEP;
								gottgrids = (gottgridtype*)Realloc(gottgrids,
										sizeof(gottgridtype)*maxlcp);
							}
							goto exclusive_exit;
						}
					}
exclusive_exit: ;
				}
			}
		}
	}
#endif

	if(snlcp >0) gottgrids = (gottgridtype*)Realloc(gottgrids, sizeof(gottgridtype)*snlcp);
	else 
		Free(gottgrids);
	{
		savetype *saves;
		int ioffset=0;
		if(snlcp>0) qsort(gottgrids,snlcp,sizeof(gottgridtype),gottsorttype);
		for(j=0;j<nobs;j++){
			int npiobs=0;

			SX0 = sx0[j]; SY0 = sy0[j]; SZ0 = sz0[j];

			nnp[0] = nnp[1] = nnp[2] = 0;

			for(i=ioffset;i<snlcp;i++){
				if(gottgrids[i].iobs == j) npiobs ++;
			}
			for(i=ioffset;i<ioffset+npiobs;i++){
				nnp[gottgrids[i].type]++;
			}
			if(npiobs > 0) saves = (savetype *) Malloc(sizeof(savetype)*npiobs,PPTR(saves));
			for(i=0,k=ioffset;i<npiobs;i++,k++){
				saves[i].potent = gottgrids[k].potent;
				saves[i].index = gottgrids[k].index;
			}
			{
				savetype *pp;
				sprintf(surveyfilename,"Gridlightcone.%.3d.main.%.5d",j+obsidstart,tnstep);
				pp = saves;
				writedowngriddata(surveyfilename,pp,nnp[0],sx0[j],sy0[j],sz0[j]);
		
				sprintf(surveyfilename,"Gridlightcone.%.3d.outer.%.5d",j+obsidstart,tnstep);
				pp += nnp[0];
				writedowngriddata(surveyfilename,pp,nnp[1],sx0[j],sy0[j],sz0[j]);
		
				sprintf(surveyfilename,"Gridlightcone.%.3d.inner.%.5d",j+obsidstart,tnstep);
				pp += nnp[1];
				writedowngriddata(surveyfilename,pp,nnp[2],sx0[j],sy0[j],sz0[j]);
				MPI_Barrier(MPI_COMM_WORLD);
			}
			if(npiobs > 0 ) Free(saves);
			ioffset += npiobs;
		}
	}
	if(snlcp >0) Free(gottgrids);
	Free(iobs);
	if(myid==0) printf("Grid SURVEY %d detected\n",snlcp);
}
