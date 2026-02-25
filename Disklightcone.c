#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <mpi.h>
#include "Memory.h"
#include "pmheader.h"
#include "lightcone.boss.h"
#define MIN(A,B) ((A) <(B) ? (A):(B))
#define PI (3.14159265358979323846L)
MPI_Status status;
int nnp[3];
double Theta,Phi,DTheta,DPhi;
char surveyfilename[100];
int *freework;
enum boolean { NO = 00,YES = 01};
typedef struct Spherical{
	double r,theta,phi;
} Spherical;
typedef struct CHECK{
	int ibox;
	int ip;
} CHECK;
/*
typedef struct InBox{
	int nx,ny,nz;
} InBox;
*/
InBox *box;
double r1, r2;
double theta,phi,theta1,phi1,theta2,phi2;
float omega0,H,lambda0;
float qsimp(float (*)(float),float,float);
float trapzd(float (*)(float), float , float , int );
/* ov: velocity data at half-step before the current step */
/* nv: velocity data at half-step after the current step */

int pslicesorttype(const void *a,const void *b){
	slcparticletype *aa,*bb;
	aa = (slcparticletype *)a;
	bb = (slcparticletype *)b;
	if(aa->type < bb->type) return -1;
	else if(aa->type > bb-> type) return +1;
	else return 0;
}

#define NPSTEP 1000000

static double rmin2,rmax2;
static slcparticletype *slcparticles;
static int snlcp,maxlcp,tsnlcp,mbox;
/* if dist >= r1-maxd && dist < r1 --> located in inner buffer zone 
   if dist >= r1 && dist < r2 --> located in the main zone
   if dist >= r2 && dist < r2+maxd --> located in the outter buffer zone
   */
#define ZWIDTH (20.)
float SX0,SY0,SZ0;
#define twoPI (6.28318530717958647692L)
#define deg2rad (.01745329251994329576L)
void a2comovingpixel(double *, double *, double *,int , float ,float , float ,float ,
		float ,float ,float , float );
void DiskESLightConeData(pmparticletype *pmparticles,long np,
		int tnstep, float anow,float preastep, float astep,float maxd, int flagsync){
	int nx,ny,nz, nspace;
	float  amax, a,  boxsize;
	float omep,omeplam,hubble,omepb;
	int myid,nid;
	pmparticletype *bp;
	float redshift;
	double r,dr1,dr2,distsq,dist,rot;
	double xp,yp,zp;
	double xpp,ypp,zpp;
	double rminsq,rmaxsq;
	double xrot,yrot,zrot,zwidth,mzwidth,xppp,yppp,zppp;
	int nbox,ibox;
	long i,j,k;
	long ix,iy,iz,nxy;
	InBox box[9];

	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;
	nspace = simpar.nspace;
	amax = simpar.amax;
	a = anow;
	boxsize = simpar.boxsize;
	omep = simpar.omep;
	omeplam = simpar.omeplam;
	omepb = simpar.omeplam;
	hubble = simpar.hubble;
	myid = simpar.myid;
	nid = simpar.nid;


	SX0 = nx;
	SY0 = ny;
	SZ0 = nz/2.;

	nbox = 4;

	xrot = atan(0.1L);
	yrot = atan(0.23333L);
	ibox = 0;
	for(j=0;j<2;j++) for(i=0;i<2;i++){
		box[ibox].nx = nx*i-SX0;
		box[ibox].ny = ny*j-SY0;
		box[ibox].nz =     -SZ0;
		ibox ++;
	}

	zwidth = ZWIDTH;
	mzwidth = -1.*zwidth;

	a2comovingpixel(&r,&dr1,&dr2,nx,amax,a,preastep,astep,
			omep,omeplam,boxsize,hubble);
	r1 = (r-dr1);
	if(r1 <0.) r1 = 0;
	r2 = (r+dr2);
	rmin2 = r1*r1;
	rmax2 = r2*r2;
	rminsq = (r1-maxd)*(r1-maxd);
	if(r1<maxd) rminsq = 0.;
	rmaxsq = (r2+maxd)*(r2+maxd);
	if(myid==0) printf("Disk SURVEY %g %g %g at observer=%g %g %g ....",r,dr1,dr2,
			SX0,SY0,SZ0);
	maxlcp = NPSTEP;

#ifdef _OPENMP
	{
		CHECK *check;
		int ii,jj;
		check = (CHECK *)Malloc(sizeof(CHECK)*np,PPTR(check));
#		pragma omp parallel for private(i)
		for(i=0;i<np;i++) check[i].ibox = -1;
#		pragma omp parallel  for private(i,j,xp,yp,zp,distsq,xpp,ypp,zpp,xppp,yppp,zppp)
		for(i=0;i<np;i++){
			for(j=0;j<nbox;j++){
				xp = (double)XofP(pmparticles+i) + (double)box[j].nx;
				yp = (double)YofP(pmparticles+i) + (double)box[j].ny;
				zp = (double)ZofP(pmparticles+i) + (double)box[j].nz;
				distsq = xp*xp + yp*yp + zp*zp;
				if(distsq >= rminsq && distsq < rmaxsq){
					xpp =  xp;
					ypp =  cos(xrot)*yp + sin(xrot)*zp;
					zpp = -sin(xrot)*yp + cos(xrot)*zp; 
		
					yppp = ypp;
					zppp =  cos(yrot)*zpp + sin(yrot)*xpp;
					xppp = -sin(yrot)*zpp + cos(yrot)*xpp;
		
					if(zppp >=mzwidth && zppp < zwidth){
							check[i].ibox = j;
							check[i].ip = i;
					}
				}
			}
		}
		snlcp = 0;
		for(i=0;i<np;i++){
			if(check[i].ibox >=0){
				check[snlcp] = check[i];
				snlcp++;
			}
		}
		check = Realloc(check,sizeof(CHECK)*snlcp);
		if(snlcp >0) slcparticles = (slcparticletype *)Malloc(sizeof(slcparticletype)*snlcp,
				PPTR(slcparticles));
#		pragma omp parallel for private(ii,i,j,xp,yp,zp,distsq)
		for(ii=0;ii<snlcp;ii++){ i = check[ii].ip; j = check[ii].ibox;
			xp = (double)XofP(pmparticles+i) + (double)box[j].nx;
			yp = (double)YofP(pmparticles+i) + (double)box[j].ny;
			zp = (double)ZofP(pmparticles+i) + (double)box[j].nz;
			distsq = xp*xp + yp*yp + zp*zp;
			slcparticles[ii].x = xp;
			slcparticles[ii].y = yp;
			slcparticles[ii].z = zp;
			slcparticles[ii].vx = pmparticles[i].vx;
			slcparticles[ii].vy = pmparticles[i].vy;
			slcparticles[ii].vz = pmparticles[i].vz;
			slcparticles[ii].bp = pmparticles+i;
			if(distsq>=rmin2&& distsq < rmax2) slcparticles[ii].type = 0;
			else if(distsq >=rmax2) slcparticles[ii].type = 1;
			else slcparticles[ii].type = 2;
		}
		Free(check);
	}
#else
	{
		snlcp = 0;
		slcparticles = (slcparticletype *)Malloc(sizeof(slcparticletype)*maxlcp,
				PPTR(slcparticles));
		for(i=0;i<np;i++){
			for(j=0;j<nbox;j++){
				xp = (double)XofP(pmparticles+i) + (double)box[j].nx;
				yp = (double)YofP(pmparticles+i) + (double)box[j].ny;
				zp = (double)ZofP(pmparticles+i) + (double)box[j].nz;
				distsq = xp*xp + yp*yp + zp*zp;
				if(distsq >= rminsq && distsq < rmaxsq){
					xpp =  xp;
					ypp =  cos(xrot)*yp + sin(xrot)*zp;
					zpp = -sin(xrot)*yp + cos(xrot)*zp; 
		
		
					yppp = ypp;
					zppp =  cos(yrot)*zpp + sin(yrot)*xpp;
					xppp = -sin(yrot)*zpp + cos(yrot)*xpp;
		
					if(zppp >=mzwidth && zppp < zwidth){
						{
							slcparticles[snlcp].x = xp;
							slcparticles[snlcp].y = yp;
							slcparticles[snlcp].z = zp;
							slcparticles[snlcp].vx = pmparticles[i].vx;
							slcparticles[snlcp].vy = pmparticles[i].vy;
							slcparticles[snlcp].vz = pmparticles[i].vz;
							slcparticles[snlcp].bp = pmparticles+i;
							if(distsq>=rmin2&& distsq < rmax2) slcparticles[snlcp].type = 0;
							else if(distsq >=rmax2) slcparticles[snlcp].type = 1;
							else slcparticles[snlcp].type = 2;
							snlcp ++;
							if(snlcp >= maxlcp){
								maxlcp += NPSTEP;
								slcparticles = (slcparticletype*)Realloc(slcparticles,
										sizeof(slcparticletype)*maxlcp);
		
							}
						}
					}
				}
			}
		}
	}
	if(snlcp>0) slcparticles = (slcparticletype*)Realloc(slcparticles,
		sizeof(slcparticletype)*snlcp);
	else
		Free(slcparticles);
#endif



	SLCP_OBSTEMWRITE(slcparticles,snlcp,myid,nid,1000);
	if(snlcp >0) Free(slcparticles);
	MPI_Barrier(MPI_COMM_WORLD);
	{
		int tsnlcp;
		MPI_Reduce(&snlcp,&tsnlcp,1,MPI_INT,MPI_SUM, 0,MPI_COMM_WORLD);
		if(myid==0) printf("Disk SURVEY %d detected\n",tsnlcp);
	}
}


void DiskSavingLightConeData(pmparticletype *pmparticles,long np,
		int tnstep, float amax, float a, float astep,float nextastep){
	float boxsize, hboxsize;
	float omei;
	float omep,omepb,omeplam,hubble;
	int myid,nid;
	int nx, ny, nz;
	pmparticletype *bp;
	float redshift;
	double r,dr1,dr2;
	double xp,yp,zp,distsq,dist;
	float rminsq,rmaxsq;
	int nbox;
	long i,j,k;
	float ax,ay,az;
	slcparticletype *tp;
	lightconetype *lightcones;
	float ainv,app,apsq,afact,bfact,vfact1h,vfact2h,rng;
	float asteph;
	float vfact1,vfact2;
	float gv1,gv2;
	float wlam0;
	char nowtsubpower;
	MPI_Status status;

	boxsize = simpar.boxsize;
	hboxsize = simpar.boxsize*0.5; 
	omei = simpar.omei; 
	omep = simpar.omep; 
	omepb = simpar.omepb; 
	omeplam = simpar.omeplam; 
	hubble = simpar.hubble; 
	wlam0 = simpar.wlam0;
	myid = simpar.myid; 
	nid = simpar.nid; 
	nx = simpar.nx; 
	ny = simpar.ny; 
	nz = simpar.nz;
	rng = nx;


	SLCP_OBSTEMREAD(slcparticles,snlcp,myid,nid,1000);
	MPI_Reduce(&snlcp,&tsnlcp,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&tsnlcp,1,MPI_INT,0,MPI_COMM_WORLD);
	if(tsnlcp==0) return;


	asteph = astep*0.5;
	ainv = 1./a;
	app = -4.*pi/3.*(ainv*ainv+(1.+3*wlam0)*omeplam/omep/pow(a,2.+3*wlam0)*pow(amax, 3*wlam0));
	apsq = 8.*pi/3.*(ainv+1./omei-1.+omeplam/omep*(pow(a,-1.-3*wlam0)-1.)*pow(amax,3*wlam0));
	afact = (2.+app*a/apsq)/2.*ainv;
	bfact = ainv*ainv*ainv/apsq;
	vfact1 = (1.-afact*astep)/(1.+afact*astep);
	vfact2 = bfact*astep/(1.+afact*astep)*rng;
	vfact1h = (1.-afact*asteph)/(1.+afact*asteph);
	vfact2h = bfact*asteph/(1.+afact*asteph)*rng;

	vfact1h = (1.-afact*astep);
	vfact2h = bfact*astep*0.5*rng;

	gv1 = vfact1/vfact1h;
	gv2 = (vfact2-vfact1*vfact2h)/vfact2;


#ifdef _OPENMP
#pragma omp parallel for private(bp,i)
#endif
	for(i=0;i<snlcp;i++){
		bp = slcparticles[i].bp;
		slcparticles[i].vx = gv1*slcparticles[i].vx + gv2*(bp->vx-vfact1*slcparticles[i].vx);
		slcparticles[i].vy = gv1*slcparticles[i].vy + gv2*(bp->vy-vfact1*slcparticles[i].vy);
		slcparticles[i].vz = gv1*slcparticles[i].vz + gv2*(bp->vz-vfact1*slcparticles[i].vz);
	}

	nnp[0] = nnp[1] = nnp[2] =0;
	for(i=0;i<snlcp;i++){
		nnp[slcparticles[i].type]++;
	}
	if(snlcp>0) qsort(slcparticles,snlcp,sizeof(slcparticletype),pslicesorttype);
	if(snlcp>0) lightcones = (lightconetype *)Malloc(sizeof(lightconetype)*snlcp,
			PPTR(lightcones));
#ifdef _OPENMP
#pragma omp parallel for private(bp,i)
#endif
	for(i=0;i<snlcp;i++){
		lightcones[i].x = slcparticles[i].x;
		lightcones[i].y = slcparticles[i].y;
		lightcones[i].z = slcparticles[i].z;
		lightcones[i].vx = slcparticles[i].vx;
		lightcones[i].vy = slcparticles[i].vy;
		lightcones[i].vz = slcparticles[i].vz;
		bp = slcparticles[i].bp;
#ifdef INDEX
		lightcones[i].indx = bp->indx;
#endif
	}
	{
		lightconetype *pp;

		sprintf(surveyfilename,"Disklightcone.main.%.5d",tnstep);
		pp = lightcones;
		writedownlightconedata(surveyfilename,pp,nnp[0],SX0,SY0,SZ0);

		sprintf(surveyfilename,"Disklightcone.outer.%.5d",tnstep);
		pp += nnp[0];
		writedownlightconedata(surveyfilename,pp,nnp[1],SX0,SY0,SZ0);

		sprintf(surveyfilename,"Disklightcone.inner.%.5d",tnstep);
		pp += nnp[1];
		writedownlightconedata(surveyfilename,pp,nnp[2],SX0,SY0,SZ0);
	}
	if(snlcp>0) Free(lightcones);
	Free(slcparticles);
}

#undef npstep
