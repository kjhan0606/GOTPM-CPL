/*
 * ?? ?Úµ??? ?????? X0,Y0,Z0 ?? À§Ä¡?? ?????Ú°? ?Ã¹Ä·??Ì¼? ?Ú½??? ????À¸??
 * ?Ã¹Ä·??Ì¼? ???Úµ?À» ??????À»??  ?????? ?????? ??Å¸????.
 * ???? ?????? comoving space ?????? ?Ã¹Ä·??Ì¼Ç¿??? ?????? ?Å¸???À§?? ???Â´?.
 * ?????? theta, phi  ?? ???? 1???Ð¸é¿¡???? ?????? ?????Ï´?.
 * ???? 2,3,4 ???Ð¸é¿¡???? ?????Ï°? ?Ò·Á¸?, min,max ?? ????À» ??????
 * ?Ù²Ù¾??? ?Ñ´?. ?×¸??? ?????? ?ÎºÐ¿? ?Ø´??? ???Ú¸? Ã£À» ?? phi ??À» 
 * ?????Ï´Âµ?,
 * ???â¿¡???? ?×³? 1,4???Ð¸??Ì¶??? ??Á¤À» ?ß´?. 
 * wflag ?? ??Àº 2????À¸?? Ç¥???Ø¾? ?Ñ´?.
 * 09/08/2002 ?????? 
 */
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <mpi.h>
#include "Memory.h"
#include "pmheader.h"
#include "lightcone.h"
#define MIN(A,B) ((A) <(B) ? (A):(B))
#define PI (3.14159265358979323846L)
double Theta,Phi,DTheta,DPhi;
char surveyfilename[100];
int *freework;
enum boolean { NO = 00,YES = 01};
typedef struct Spherical{
	double r,theta,phi;
} Spherical;
InBox *box;
double r1, r2;
double theta,phi,theta1,phi1,theta2,phi2;
float omega0,H,lambda0;
float qsimp(float (*)(float),float,float);
float trapzd(float (*)(float), float , float , int );
/* ov: velocity data at half-step before the current step */
/* nv: velocity data at half-step after the current step */
typedef struct slcparticletype{
#ifdef XYZDBL
	double x,y,z;
#else
	float x,y,z;
#endif
	float vx,vy,vz;
	float ovx,ovy,ovz,nvx,nvy,nvz;
	pmparticletype *bp;
	/*
	float dist;
	*/
#ifdef INDEX
	indextype indx;
#endif
	int type;
} slcparticletype;
typedef struct lightconetype{
#ifdef XYZDBL
	double x,y,z;
#else
	float x,y,z;
#endif
	float vx,vy,vz;
	float ovx,ovy,ovz,nvx,nvy,nvz;
#ifdef INDEX
	indextype indx;
#endif
} lightconetype;

int ssorttype(const void *a,const void *b){
	slcparticletype *aa,*bb;
	aa = (slcparticletype *)a;
	bb = (slcparticletype *)b;
	if(aa->type < bb->type) return -1;
	else if(aa->type > bb-> type) return +1;
	else return 0;
}

#define NPSTEP 10000

static double rmin2,rmax2;
static slcparticletype *slcparticles;
static int snlcp,maxlcp,tsnlcp,mbox;
/* if dist >= r1-maxd && dist < r1 --> located in inner buffer zone 
   if dist >= r1 && dist < r2 --> located in the main zone
   if dist >= r2 && dist < r2+maxd --> located in the outter buffer zone
   */
#define ZWIDTH (20.)
#define SX0 ZWIDTH
#define SY0 (ny*0.5)
#define SZ0 ZWIDTH
#define twoPI (6.28318530717958647692L)
#define deg2rad (.01745329251994329576L)
void a2comovingpixel(double *, double *, double *,int , float ,float , float ,float ,
		float ,float ,float , float );
void PSliceExtractingLightConeParticles(pmparticletype *pmparticles,int np,float nextastep,float maxd){
	int nx,ny,nz, nspace;
	int mstep, amax, a, astep, boxsize;
	float omep,omeplam,hubble;
	int myid,nid;
	pmparticletype *bp;
	float redshift;
	double r,dr1,dr2,distsq,dist,rot;
	double xp,yp,zp;
	double xpp,ypp,zpp;
	double rminsq,rmaxsq;
	double yrot,phimin,zwidth,mzwidth,xppp,yppp,zppp;
	int nbox,ibox;
	int i,j,k;
	long ix,iy,iz,nxy;

	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;
	nspace = simpar.nspace;
	mstep = simpar.nstep;
	a = simpar.anow;
	amax = simpar.amax;
	astep = simpar.astep;
	boxsize = simpar.boxsize;
	omep = simpar.omep;
	omeplam = simpar.omeplam;
	hubble = simpar.hubble;
	myid = simpar.myid;
	nid = simpar.nid;

	nbox = 4;
	box = (InBox *)Malloc(sizeof(InBox)*nbox,PPTR(box));
	/*
	theta1 = 75.L*deg2rad;
	theta2 = 105.L*deg2rad;
	phi1 = 30.L*deg2rad;
	phi2 = 330.L*deg2rad;
	*/

	/*
	zrot = 45.L*deg2rad;
	yrot = acos(sqrt(2.)/sqrt(3.));
	*/
	yrot = atan(0.5L);
	phimin = 33.L*deg2rad;
	/*
	nbox = 0;
	for(k=-1;k<=0;k++){
		box[ibox].nx = 0.;
		box[ibox].ny = 0.;
		box[ibox].nz = k*nz-SZ0;
		nbox ++;
	}
	*/
	box[0].nx = -SX0;
	box[0].ny = -SY0;
	box[0].nz = -SZ0;

	box[1].nx = nx-SX0;
	box[1].ny =   -SY0;
	box[1].nz =   -SZ0;

	box[2].nx = nx-SX0;
	box[2].ny = ny-SY0;
	box[2].nz =    -SZ0;

	box[3].nx =  nx-SX0;
	box[3].ny = -ny-SY0;
	box[3].nz =    -SZ0;

	nbox = 4;

	zwidth = ZWIDTH;
	mzwidth = -1.*zwidth;

	a2comovingpixel(&r,&dr1,&dr2,nx,amax,a,astep,nextastep,
			omep,omeplam,boxsize,hubble);
	r1 = (r-dr1);
	if(r1 <0.) r1 = 0;
	r2 = (r+dr2);
	rmin2 = r1*r1;
	rmax2 = r2*r2;
	rminsq = (r1-maxd)*(r1-maxd);
	if(r1<maxd) rminsq = 0.;
	rmaxsq = (r2+maxd)*(r2+maxd);
	if(myid==0) printf("Slice SURVEY %g %g %g at observer=%g %g %g ....",r,dr1,dr2,
			SX0,SY0,SZ0);
	snlcp = 0;
	maxlcp = NPSTEP;
	slcparticles = (slcparticletype *)Malloc(sizeof(slcparticletype)*maxlcp,
			PPTR(slcparticles));
	for(j=0;j<nbox;j++){
		bp = pmparticles;
		for(i=0;i<np;i++){
			xp = (double)XofP(bp) + (double)box[j].nx;
			yp = (double)YofP(bp) + (double)box[j].ny;
			zp = (double)ZofP(bp) + (double)box[j].nz;

			/*
			xpp =  cos(zrot)*xp + sin(zrot)*yp;
			ypp = -sin(zrot)*xp + cos(zrot)*yp;
			*/
			xpp =  xp;
			ypp =  yp;

			zpp = zp;
			yppp = ypp;
			xppp =  cos(yrot)*xpp + sin(yrot)*zpp;
			zppp = -sin(yrot)*xpp + cos(yrot)*zpp;

			distsq = xppp*xppp+yppp*yppp+zppp*zppp;
			if(distsq >= rminsq && distsq < rmaxsq && distsq <51000000.L){
				if(zppp >=mzwidth && zppp < zwidth){
					double distrho;
					distrho = sqrt(xppp*xppp+yppp*yppp);
					phi = acos(xppp/distrho);
					if(phi <= phimin){
						dist = sqrt(distsq);
						slcparticles[snlcp].x = xp;
						slcparticles[snlcp].y = yp;
						slcparticles[snlcp].z = zp;
						/* Write down the velocity data at half-step before for
						 * identifying particles */
						slcparticles[snlcp].ovx = bp->vx;
						slcparticles[snlcp].ovy = bp->vy;
						slcparticles[snlcp].ovz = bp->vz;
						slcparticles[snlcp].bp = bp;
						if(dist>=r1&& dist < r2) slcparticles[snlcp].type = 0;
						else if(dist >=r2 && dist < r2+maxd) slcparticles[snlcp].type = 1;
						else slcparticles[snlcp].type = 2;
						/*
						slcparticles[snlcp].dist = dist;
						*/
#ifdef INDEX
						slcparticles[snlcp].indx = bp->indx;
#endif
						snlcp ++;
						if(snlcp >= maxlcp){
							maxlcp += NPSTEP;
							slcparticles = (slcparticletype*)Realloc(slcparticles,
									sizeof(slcparticletype)*maxlcp);
	
						}
					}
				}
			}
			bp++;
		}
	}
	slcparticles = (slcparticletype*)Realloc(slcparticles,
		sizeof(slcparticletype)*snlcp);
	/*
	for(i=0;i<snlcp;i++){
		dist = slcparticles[i].dist;
		if(dist >=r1 && dist < r2) slcparticles[i].type = 0;
		else if(dist >=r2 && dist < r2+maxd) slcparticles[i].type = 1;
		else slcparticles[i].type = 2;
	}
	*/
	OBSTEMWRITE(slcparticles,snlcp,myid,nid,10);
	if(myid==0) printf("Slice SURVEY %d detected\n",snlcp);
	Free(box);
}

#define writedownlightconedata(A,B,C) {\
	FILE *wp;\
	float x0,y0,z0;\
	int nmlcp;\
	lightconetype *ttmp;\
	x0 = SX0;\
	y0 = SY0;\
	z0 = SZ0;\
	{\
		if(myid==0) {\
			wp=fopen(A,"w");\
			fwrite(&nx,sizeof(int),1,wp);\
			fwrite(&ny,sizeof(int),1,wp);\
			fwrite(&nz,sizeof(int),1,wp);\
			fwrite(&omep,sizeof(float),1,wp);\
			fwrite(&omeplam,sizeof(float),1,wp);\
			fwrite(&omepb,sizeof(float),1,wp);\
			fwrite(&hubble,sizeof(float),1,wp);\
			fwrite(&boxsize,sizeof(float),1,wp);\
			fwrite(&amax,sizeof(float),1,wp);\
			fwrite(&a,sizeof(float),1,wp);\
			fwrite(&astep,sizeof(float),1,wp);\
			fwrite(&x0,sizeof(float),1,wp);\
			fwrite(&y0,sizeof(float),1,wp);\
			fwrite(&z0,sizeof(float),1,wp);\
			fwrite(&rmin2,sizeof(double),1,wp);\
			fwrite(&rmax2,sizeof(double),1,wp);\
			fwrite(&C,sizeof(int),1,wp);\
			fwrite(B,sizeof(lightconetype),C,wp);\
			for(i=1;i<nid;i++){	\
				MPI_Recv(&nmlcp,1,MPI_INT,i,i,MPI_COMM_WORLD,&status);\
				if(nmlcp>0) ttmp=(lightconetype*)Malloc(sizeof(lightconetype)*nmlcp,PPTR(ttmp));\
				MPI_Recv(ttmp,sizeof(lightconetype)*nmlcp,MPI_BYTE,i,i,\
						MPI_COMM_WORLD,&status);\
				fwrite(&nmlcp,sizeof(int),1,wp);\
				fwrite(ttmp,sizeof(lightconetype),nmlcp,wp);\
				if(nmlcp>0) Free(ttmp);\
			}\
			fclose(wp);\
		}\
		else  {\
			MPI_Send(&C,1,MPI_INT,0,myid,MPI_COMM_WORLD);\
			MPI_Send(B,sizeof(lightconetype)*C,MPI_BYTE,0,myid,MPI_COMM_WORLD);\
		}\
	}\
}

void PSliceSavingLightConeData(pmparticletype *pmparticles,int np,
		int nx,int ny, int nz, int mstep, float amax, float a, 
		float astep,float nextastep, float boxsize,float omei,
		float omep,float omepb,float omeplam,float wlam0, float wlam1, float hubble,int myid,int nid){
	pmparticletype *bp;
	float redshift;
	double r,dr1,dr2;
	double xp,yp,zp,distsq,dist;
	float rminsq,rmaxsq;
	int nbox;
	int i,j,k;
	float ax,ay,az;
	int *nnp;
	slcparticletype *tp;
	lightconetype *lightcones;
	float ainv,app,apsq,afact,bfact,vfact1h,vfact2h,rng;
	float asteph;
	float vfact1,vfact2;
	float gv1,gv2;
	MPI_Status status;



	OBSTEMREAD(slcparticles,snlcp,myid,nid,10);
	MPI_Reduce(&snlcp,&tsnlcp,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&tsnlcp,1,MPI_INT,0,MPI_COMM_WORLD);
	if(tsnlcp==0) return;


	asteph = astep*0.5;
	rng = nx;
	ainv = 1./a;
	/*
    app  = -4.*M_PI/3.*(ainv*ainv+(1.+3*wlam0)*omeplam/omep/pow(a,2.+3*wlam0)*pow(amax, 3*wlam0));
    apsq = 8.*M_PI/3.*(ainv+1./omei-1.+omeplam/omep*(pow(a,-1.-3*wlam0)-1.)*pow(amax,3*wlam0));
	*/
	void appapsq_(float *, float *, float *, float *, float *, float *, float *, float *);
	appapsq_(&omep, &omeplam, &wlam0, &wlam1, &amax, &a, &app, &apsq);
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


	nnp = (int *)Malloc(sizeof(int)*3,PPTR(nnp));
	nnp[0] = nnp[1] = nnp[2] =0;
	tp = slcparticles;
	for(i=0;i<snlcp;i++){
		bp = tp->bp;
		/* write down the half-step advanced velocity data for identifying particles
		 * in the next step */
		tp->nvx = bp->vx;
		tp->nvy = bp->vy;
		tp->nvz = bp->vz;
		/*
		ax = (bp->vx-vfact1*tp->ovx)/vfact2;
		ay = (bp->vy-vfact1*tp->ovy)/vfact2;
		az = (bp->vz-vfact1*tp->ovz)/vfact2;
		tp->vx = vfact1h*tp->vx + vfact2h*ax;
		tp->vy = vfact1h*tp->vy + vfact2h*ay;
		tp->vz = vfact1h*tp->vz + vfact2h*az;
		*/
		tp->vx = gv1*tp->ovx + gv2*(tp->nvx-vfact1*tp->ovx);
		tp->vy = gv1*tp->ovy + gv2*(tp->nvy-vfact1*tp->ovy);
		tp->vz = gv1*tp->ovz + gv2*(tp->nvz-vfact1*tp->ovz);

		tp++;
	}
	for(i=0;i<snlcp;i++){
		nnp[slcparticles[i].type]++;
	}
	if(snlcp>0) qsort(slcparticles,snlcp,sizeof(slcparticletype),ssorttype);
	if(snlcp>0) lightcones = (lightconetype *)Malloc(sizeof(lightconetype)*snlcp,
			PPTR(lightcones));
	for(i=0;i<snlcp;i++){
		lightcones[i].x = slcparticles[i].x;
		lightcones[i].y = slcparticles[i].y;
		lightcones[i].z = slcparticles[i].z;
		lightcones[i].vx = slcparticles[i].vx;
		lightcones[i].vy = slcparticles[i].vy;
		lightcones[i].vz = slcparticles[i].vz;
		lightcones[i].ovx = slcparticles[i].ovx;
		lightcones[i].ovy = slcparticles[i].ovy;
		lightcones[i].ovz = slcparticles[i].ovz;
		lightcones[i].nvx = slcparticles[i].nvx;
		lightcones[i].nvy = slcparticles[i].nvy;
		lightcones[i].nvz = slcparticles[i].nvz;
#ifdef INDEX
		lightcones[i].indx = slcparticles[i].indx;
#endif
	}
	{
		lightconetype *pp;

		sprintf(surveyfilename,"Pslightcone.main.%.5d",mstep);
		pp = lightcones;
		writedownlightconedata(surveyfilename,pp,nnp[0]);

		sprintf(surveyfilename,"Pslightcone.outer.%.5d",mstep);
		pp += nnp[0];
		writedownlightconedata(surveyfilename,pp,nnp[1]);

		sprintf(surveyfilename,"Pslightcone.inner.%.5d",mstep);
		pp += nnp[1];
		writedownlightconedata(surveyfilename,pp,nnp[2]);
	}
	/*
	for(i=0;i<nid;i++){
		if(myid==i){
			lightconetype *pp;
			sprintf(surveyfilename,"lightcone.main.%.5d",mstep);
			pp = lightcones;
			writedownlightconedata(surveyfilename,pp,nnp[0]);
	
			sprintf(surveyfilename,"lightcone.outer.%.5d",mstep);
			pp += nnp[0];
			writedownlightconedata(surveyfilename,pp,nnp[1]);
	
			sprintf(surveyfilename,"lightcone.inner.%.5d",mstep);
			pp += nnp[1];
			writedownlightconedata(surveyfilename,pp,nnp[2]);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	*/
	if(snlcp>0) Free(lightcones);
	Free(nnp);
	Free(slcparticles);
}

#undef npstep
