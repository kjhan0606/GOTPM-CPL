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
int bosssorttype(const void *a,const void *b){
	bossparticletype *aa,*bb;
	aa = (bossparticletype *)a;
	bb = (bossparticletype *)b;
	if(aa->iobs < bb->iobs) return -1;
	else if(aa->iobs > bb->iobs) return +1;
	else {
		if(aa->type < bb->type) return -1;
		else if(aa->type > bb->type) return +1;
		return 0;
	}
}
typedef struct CHECK{
	int iobs;
	int ip;
}CHECK;
#define NPSTEP 1000000
int nnp[3];
/* Particle data type for 108 BOSS survey */
static double rmin2,rmax2;
static bossparticletype *bossparticles;

static int snlcp,maxlcp,tsnlcp,mbox;
/* if dist >= r1-maxd && dist < r1 --> located in inner buffer zone 
   if dist >= r1 && dist < r2 --> located in the main zone
   if dist >= r2 && dist < r2+maxd --> located in the outter buffer zone
   */
#define twoPI (6.283185306L)
#define deg2rad (0.017453292516L)


double SX0,SY0,SZ0;
void a2comovingpixel(double *, double *, double *,int , float ,float , float ,float , 
		float ,float ,float , float );
void BossESLightConeData(pmparticletype *pmparticles,long np,int tnstep,float a,
		float preastep, float astep,float maxd,int flagsync,int obsidstart,
		double *sx0, double *sy0, double *sz0,int nobs, float redi){
	int nx,ny,nz, nspace;
	int mobs; 
	float amax;
	float boxsize,hboxsize;
	float omep,omeplam,hubble,omepb;
	int myid,nid;
	pmparticletype *bp;
	float redshift;
	double r,dr1,dr2,distsq,dist,rot;
	double xp,yp,zp;
	double xpp,ypp,zpp;
	double rminsq,rmaxsq;
	int nbox,ibox;
	long i,j,k;
	long ix,iy,iz;
	lightconetype *lightcones;
	int *iobs;
	int jobs,jj;

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

	if(1){
		float redshift;
		redshift = amax/a-1.;
		if(redshift > redi) return;
	}
	/*
	nbox = 3*3*2;
	box = (InBox *)Malloc(sizeof(InBox)*nbox,PPTR(box));
	theta1 = 48.L*deg2rad;
	theta2 = 132.L*deg2rad;
	phi1 = 65.L*deg2rad;
	rot = 0.L;
	*/
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
	if(myid==0) printf("Spherical Boss Survey %g %g %g\n",r,dr1,dr2);
	snlcp = 0;
	maxlcp = NPSTEP;
	{
		iobs = (int *)Malloc(sizeof(int)*nobs,PPTR(iobs));
		mobs = 0;
		for(i=0;i<nobs;i++){
			float zlow,zhigh;
			zlow = simpar.zmin-sqrt(rmaxsq);
			zhigh = simpar.zmax+sqrt(rmaxsq);
			if(sz0[i] > zlow && sz0[i] < zhigh) iobs[mobs++] = i;
			else if(sz0[i]+nz > zlow && sz0[i]+nz< zhigh) iobs[mobs++] = i;
			else if(sz0[i]-nz > zlow && sz0[i]-nz< zhigh) iobs[mobs++] = i;
		}		
	}
#ifdef _OPENMP
	{
		CHECK *check;
		check = (CHECK *)Malloc(sizeof(CHECK)*np,PPTR(check));
#		pragma omp parallel for 
		for(i=0;i<np;i++) check[i].iobs = -1;
#		pragma omp parallel for  private(SX0,SY0,SZ0,xp,yp,zp,i,distsq)
		for(i=0;i<np;i++){
			for(jj=0;jj<mobs;jj++){ j = iobs[jj];
				SX0 = sx0[j]; SY0 = sy0[j]; SZ0 = sz0[j];
				xp = XofP(pmparticles+i) - SX0;
				yp = YofP(pmparticles+i) - SY0;
				zp = ZofP(pmparticles+i) - SZ0;
				/* This should be dealt with much care. */
				BOSSBCondition(xp,hboxsize,boxsize);
				BOSSBCondition(yp,hboxsize,boxsize);
				BOSSBCondition(zp,hboxsize,boxsize);

				distsq = xp*xp+yp*yp+zp*zp;
				if(distsq >= rminsq && distsq < rmaxsq){
					check[i].iobs = j;
					check[i].ip = i;
					goto exclusive_exit;
				}
			}
exclusive_exit:;
		}
		snlcp = 0;
		for(i=0;i<np;i++){
			if(check[i].iobs>=0) {
				check[snlcp] = check[i];
				snlcp ++;
			}
		}
		check = Realloc(check,sizeof(int)*snlcp);
		if(snlcp > 0) bossparticles = (bossparticletype *)Malloc( 
				sizeof(bossparticletype)*snlcp, PPTR(bossparticles));
#		pragma omp parallel for private(i,jobs,jj,SX0,SY0,SZ0,xp,yp,zp,distsq)
		for(i=0;i<snlcp;i++){
			jobs = check[i].iobs;
			jj = check[i].ip;

			SX0 = sx0[jobs]; SY0 = sy0[jobs]; SZ0 = sz0[jobs];
			xp = XofP(pmparticles+jj) - SX0;
			yp = YofP(pmparticles+jj) - SY0;
			zp = ZofP(pmparticles+jj) - SZ0;
			BOSSBCondition(xp,hboxsize,boxsize);
			BOSSBCondition(yp,hboxsize,boxsize);
			BOSSBCondition(zp,hboxsize,boxsize);
			distsq = xp*xp + yp*yp +zp*zp;

			if(distsq>=rmin2&& distsq < rmax2) bossparticles[i].type = 0;
			else if(distsq >=rmax2) bossparticles[i].type = 1;
			else bossparticles[i].type = 2;
			bossparticles[i].iobs = jobs;

			bossparticles[i].vx = pmparticles[jj].vx;
			bossparticles[i].vy = pmparticles[jj].vy;
			bossparticles[i].vz = pmparticles[jj].vz;
			bossparticles[i].bp = pmparticles+jj;
		}

		Free(check);
	}
#else
	{
		bossparticles = (bossparticletype *)Malloc(
						sizeof(bossparticletype)*maxlcp, PPTR(bossparticles));
		for(i=0;i<np;i++){
			for(jj=0;jj<mobs;jj++){ j = iobs[jj];
				SX0 = sx0[j]; SY0 = sy0[j]; SZ0 = sz0[j];
				xp = XofP(pmparticles+i) - SX0;
				yp = YofP(pmparticles+i) - SY0;
				zp = ZofP(pmparticles+i) - SZ0;
				/* This should be dealt with much care. */
				BOSSBCondition(xp,hboxsize,boxsize);
				BOSSBCondition(yp,hboxsize,boxsize);
				BOSSBCondition(zp,hboxsize,boxsize);

				distsq = xp*xp+yp*yp+zp*zp;
				if(distsq >= rminsq && distsq < rmaxsq){
					bossparticles[snlcp].vx = pmparticles[i].vx;
					bossparticles[snlcp].vy = pmparticles[i].vy;
					bossparticles[snlcp].vz = pmparticles[i].vz;
					bossparticles[snlcp].bp = pmparticles+i;
					if(distsq>=rmin2&& distsq < rmax2) bossparticles[snlcp].type = 0;
					else if(distsq >=rmax2) bossparticles[snlcp].type = 1;
					else bossparticles[snlcp].type = 2;
					bossparticles[snlcp].iobs = j;
					snlcp ++; 
					if(snlcp >= maxlcp){
						maxlcp += NPSTEP;
						bossparticles = (bossparticletype*)Realloc(bossparticles,
								sizeof(bossparticletype)*maxlcp);
					}
					goto exclusive_exit;
				}
			}
exclusive_exit:;
		}
	}
	if(snlcp >0) bossparticles = (bossparticletype*)Realloc(bossparticles, sizeof(bossparticletype)*snlcp);
	else 
		Free(bossparticles);

#endif

	if(flagsync){
		int ioffset=0;
		if(snlcp>0) qsort(bossparticles,snlcp,sizeof(bossparticletype),bosssorttype);
		for(j=0;j<nobs;j++){
			int npiobs=0;

			SX0 = sx0[j]; SY0 = sy0[j]; SZ0 = sz0[j];

			nnp[0] = nnp[1] = nnp[2] = 0;
			for(i=ioffset;i<snlcp;i++){
				if(bossparticles[i].iobs == j) npiobs ++;
			}
			for(i=ioffset;i<ioffset+npiobs;i++){
				nnp[bossparticles[i].type]++;
			}
			if(npiobs > 0) lightcones = (lightconetype *)
				Malloc(sizeof(lightconetype)*npiobs,PPTR(lightcones));
			k=ioffset;
			for(i=0;i<npiobs;i++,k++){

				bp = bossparticles[k].bp;
				xp = XofP(bp) - SX0;
				yp = YofP(bp) - SY0;
				zp = ZofP(bp) - SZ0;
				/* This should be dealt with much care. */
				BOSSBCondition(xp,hboxsize,boxsize);
				BOSSBCondition(yp,hboxsize,boxsize);
				BOSSBCondition(zp,hboxsize,boxsize);

				lightcones[i].x = xp;
				lightcones[i].y = yp;
				lightcones[i].z = zp;
				lightcones[i].vx = bossparticles[k].vx;
				lightcones[i].vy = bossparticles[k].vy;
				lightcones[i].vz = bossparticles[k].vz;
#ifdef INDEX
				lightcones[i].indx = bp->indx;
#endif
			}
			{
				lightconetype *pp;
		
				sprintf(surveyfilename,"Bosslightcone.%.3d.main.%.5d",(int)(j+obsidstart),tnstep);
				pp = lightcones;
				writedownlightconedata(surveyfilename,pp,nnp[0],sx0[j],sy0[j],sz0[j]);
		
				sprintf(surveyfilename,"Bosslightcone.%.3d.outer.%.5d",(int)(j+obsidstart),tnstep);
				pp += nnp[0];
				writedownlightconedata(surveyfilename,pp,nnp[1],sx0[j],sy0[j],sz0[j]);
		
				sprintf(surveyfilename,"Bosslightcone.%.3d.inner.%.5d",(int)(j+obsidstart),tnstep);
				pp += nnp[1];
				writedownlightconedata(surveyfilename,pp,nnp[2],sx0[j],sy0[j],sz0[j]);
				MPI_Barrier(MPI_COMM_WORLD);
			}
			if(npiobs > 0 ) Free(lightcones);
			ioffset += npiobs;
		}
	}
	else {
		OBSTEMWRITE(bossparticles,snlcp,myid,nid,obsidstart);
	}
	if(snlcp >0) Free(bossparticles);
	if(myid==0) printf("Boss SURVEY %d detected\n",snlcp);

	Free(iobs);
}
void BossSavingLightConeData(pmparticletype *pmparticles,long np,int tnstep, 
		float amax, float a, float preastep,float astep,
		int obsidstart, double *sx0, double *sy0, double *sz0,int nobs,float redi){
	pmparticletype *bp;
	int nx,ny,nz;
	float redshift;
	double r,dr1,dr2;
#ifdef XYZDBL
	double xp,yp,zp,distsq,dist;
#else
	float xp,yp,zp,distsq,dist;
#endif
	float rminsq,rmaxsq;
	long i,j,k;
	float ax,ay,az;
	bossparticletype *tp;
	lightconetype *lightcones;
	float ainv,app,apsq,afact,bfact,vfact1h,vfact2h,rng;
	float asteph;
	float vfact1,vfact2;
	float gv1,gv2;
	float omei,omep,omepb,omeplam,hubble;
	int myid,nid;
	char nowtsubpower;
	float boxsize,hboxsize;
	float wlam0, wlam1;

	boxsize = simpar.boxsize;
	hboxsize = simpar.boxsize*0.5;
	omei = simpar.omei;
	omep = simpar.omep;
	omepb = simpar.omepb;
	omeplam = simpar.omeplam;
	wlam0 = simpar.wlam0;
	wlam1 = simpar.wlam1;
	hubble = simpar.hubble;
	myid = simpar.myid;
	nid = simpar.nid;
	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;

	if(1){
		float redshift;
		redshift = amax/a-1.;
		if(redshift > redi) return;
	}

	OBSTEMREAD(bossparticles,snlcp,myid,nid,obsidstart);

	MPI_Reduce(&snlcp,&tsnlcp,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&tsnlcp,1,MPI_INT,0,MPI_COMM_WORLD);
	if(tsnlcp==0) return;


	rng = nx;
	asteph = astep*0.5;
	ainv = 1./a;
	/*
	app = -4.*pi/3.*(ainv*ainv+(1.+3*wlam0)*omeplam/omep/pow(a,2.+3*wlam0)*pow(amax, 3*wlam0));
	apsq = 8.*pi/3.*(ainv+1./omei-1.+omeplam/omep*(pow(a,-1.-3*wlam0)-1.)*pow(amax,3*wlam0));
	*/
	appapsq_(&omep, &omeplam, &wlam0, &wlam1, &amax, &a, &app, &apsq);
	afact = (2.+app*a/apsq)/2.*ainv;
	bfact = ainv*ainv*ainv/apsq;
	vfact1 = (1.-afact*astep)/(1.+afact*astep);
	vfact2 = bfact*astep/(1.+afact*astep)*rng;
	/*
	vfact1h = (1.-afact*asteph)/(1.+afact*asteph);
	vfact2h = bfact*asteph/(1.+afact*asteph)*rng;
	*/
	vfact1h = (1.-afact*astep);
	vfact2h = bfact*astep*0.5*rng;

	gv1 = vfact1/vfact1h;
	gv2 = (vfact2-vfact1*vfact2h)/vfact2;



	/* recover the synchronized velocity */
#ifdef _OPENMP
#pragma omp parallel for private(bp,i)
#endif
	for(i=0;i<snlcp;i++){
		bp = bossparticles[i].bp;
		bossparticles[i].vx = gv1*bossparticles[i].vx + gv2*(bp->vx-vfact1*bossparticles[i].vx);
		bossparticles[i].vy = gv1*bossparticles[i].vy + gv2*(bp->vy-vfact1*bossparticles[i].vy);
		bossparticles[i].vz = gv1*bossparticles[i].vz + gv2*(bp->vz-vfact1*bossparticles[i].vz);
	}

	{
		int ioffset=0;
		if(snlcp>0) qsort(bossparticles,snlcp,sizeof(bossparticletype),bosssorttype);
		for(j=0;j<nobs;j++){
			int npiobs = 0;
			nnp[0] = nnp[1] = nnp[2] = 0;

			SX0 = sx0[j]; SY0 = sy0[j]; SZ0 = sz0[j];

			for(i=ioffset;i<snlcp;i++){
				if(bossparticles[i].iobs == j) npiobs ++;
			}
			for(i=ioffset;i<ioffset+npiobs;i++){
				nnp[bossparticles[i].type]++;
			}
			if(npiobs > 0) lightcones = (lightconetype *) 
				Malloc(sizeof(lightconetype)*npiobs,PPTR(lightcones));
			k=ioffset;
			for(i=0;i<npiobs;i++,k++){

				bp = bossparticles[k].bp;
				xp = XofP(bp) - SX0;
				yp = YofP(bp) - SY0;
				zp = ZofP(bp) - SZ0;
				/* This should be dealt with much care. */
				BOSSBCondition(xp,hboxsize,boxsize);
				BOSSBCondition(yp,hboxsize,boxsize);
				BOSSBCondition(zp,hboxsize,boxsize);

				lightcones[i].x = xp;
				lightcones[i].y = yp;
				lightcones[i].z = zp;
				lightcones[i].vx = bossparticles[k].vx;
				lightcones[i].vy = bossparticles[k].vy;
				lightcones[i].vz = bossparticles[k].vz;
#ifdef INDEX
				lightcones[i].indx = bp->indx;
#endif
			}
			{
				lightconetype *pp;
		
				sprintf(surveyfilename,"Bosslightcone.%.3d.main.%.5d",(int)(j+obsidstart),tnstep);
				pp = lightcones;
				writedownlightconedata(surveyfilename,pp,nnp[0],sx0[j],sy0[j],sz0[j]);
		
				sprintf(surveyfilename,"Bosslightcone.%.3d.outer.%.5d",(int)(j+obsidstart),tnstep);
				pp += nnp[0];
				writedownlightconedata(surveyfilename,pp,nnp[1],sx0[j],sy0[j],sz0[j]);
		
				sprintf(surveyfilename,"Bosslightcone.%.3d.inner.%.5d",(int)(j+obsidstart),tnstep);
				pp += nnp[1];
				writedownlightconedata(surveyfilename,pp,nnp[2],sx0[j],sy0[j],sz0[j]);
			}
			if(npiobs > 0 ) Free(lightcones);
			ioffset += npiobs;
		}
		Free(bossparticles);
	}
}
#undef npstep
