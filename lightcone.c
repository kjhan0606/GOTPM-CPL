/*
 * 이 코드는 임의의 X0,Y0,Z0 에 위치한 관측자가 시뮬레이션 박스를 기반으로
 * 시뮬레이션 입자들을 관측했을때  입자의 분포를 나타낸다.
 * 최총 결과는 comoving space 에서의 시뮬레이션에서 사용한 거리단위로 나온다.
 * 현재는 theta, phi  의 값이 1사분면에서만 계산이 가능하다.
 * 만약 2,3,4 사분면에서도 가능하게 할려면, min,max 의 값들을 적절히
 * 바꾸어야 한다. 그리고 마지막 부분에 해당되 입자를 찾을 때 phi 값을 
 * 계산하는데,
 * 여기에서는 그냥 1,4사분면이라고 가정을 했다. 
 * wflag 의 값은 2진법으로 표현해야 한다.
 * 09/08/2002 김주한 
 */
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <mpi.h>
#include "Memory.h"
#include "pmheader.h"
#include "lightcone.h"
#define MIN(A,B) ((A) <(B) ? (A):(B))
#define PI 3.14159265354979L
double Theta,Phi,DTheta,DPhi;
char surveyfilename[100];
int *freework;
enum boolean { NO = 00,YES = 01};
typedef struct Spherical{
	double r,theta,phi;
} Spherical;
unsigned int setbits(unsigned int x,int p) {
    unsigned int  b;
	b = x|(~((~0)<<1) << (p-1));
	return b;
}
InBox *box;
double r1, r2;
double theta1,phi1,theta2,phi2;
float omega0,H,lambda0;
float func1(float red){
    return 2997.92458L/sqrt((1.+red)*(1.+red)*(1.+red)*omega0
             +(1.+red)*(1.+red)*(1.-omega0-lambda0) +lambda0);
}
float qsimp(float (*)(float),float,float);
float trapzd(float (*)(float), float , float , int );
void FindcomovingpixelWidth(double *R, double *dR1, double *dR2,int nx, 
		float amax,float a, float astep1,float astep2,
		float omep,float omeplam,float size, float h){
	float initz = 0;
	float shellwidth1,shellwidth2;
	float redshift = amax/a - 1.L;
	shellwidth1 = amax/(a+astep2)-1.L;
	shellwidth2 = amax/(a-astep1)-1.L;
	if(omep == 1.0) {
		*R = 2997.92458L*2.L * (1.L-1.L/sqrt(1.L+redshift));
		*dR1 = 2997.92458L*2.L*(1.L-1.L/sqrt(1.L+shellwidth1));
		*dR2 = 2997.92458L*2.L*(1.L-1.L/sqrt(1.L+shellwidth2));
	}
	else {
		omega0 = omep;
		lambda0 = omeplam;
		H = h;
		*R = qsimp(func1,initz,redshift);
		*dR1 = qsimp(func1,initz,shellwidth1);
		*dR2 = qsimp(func1,initz,shellwidth2);
	}
	*R = *R/size*nx;
	*dR1 = *dR1/size*nx;
	*dR2 = *dR2/size*nx;
	*dR1 = *R - *dR1;
	*dR2 = *dR2 - *R;
}

void a2comovingpixel(double *R, double *dR1, double *dR2,int nx, 
		float amax,float a, float astep1,float astep2,
		float omep,float omeplam,float size, float h){
	float initz = 0;
	float shellwidth1,shellwidth2;
	float redshift = amax/a - 1.L;
	shellwidth1 = amax/(a+astep2*0.5L)-1.L;
	shellwidth2 = amax/(a-astep1*0.5L)-1.L;
	if(omep == 1.0) {
		*R = 2997.92458L*2.L * (1.L-1.L/sqrt(1.L+redshift));
		*dR1 = 2997.92458L*2.L*(1.L-1.L/sqrt(1.L+shellwidth1));
		*dR2 = 2997.92458L*2.L*(1.L-1.L/sqrt(1.L+shellwidth2));
	}
	else {
		omega0 = omep;
		lambda0 = omeplam;
		H = h;
		*R = qsimp(func1,initz,redshift);
		*dR1 = qsimp(func1,initz,shellwidth1);
		*dR2 = qsimp(func1,initz,shellwidth2);
	}
	*R = *R/size*nx;
	*dR1 = *dR1/size*nx;
	*dR2 = *dR2/size*nx;
	*dR1 = *R - *dR1;
	*dR2 = *dR2 - *R;
}
typedef struct lcparticletype{
#ifdef XYZDBL
	double x,y,z;
#else
	float x,y,z;
#endif
	float vx,vy,vz;
	pmparticletype *bp;
	float dist;
#ifdef INDEX
	indextype indx;
#endif
	int type;
} lcparticletype;
typedef struct lightconetype{
#ifdef XYZDBL
	double x,y,z;
#else
	float x,y,z;
#endif
	float vx,vy,vz;
#ifdef INDEX
	indextype indx;
#endif
} lightconetype;

int sorttype(const void *a,const void *b){
	lcparticletype *aa,*bb;
	aa = (lcparticletype *)a;
	bb = (lcparticletype *)b;
	if(aa->type < bb->type) return -1;
	else if(aa->type > bb->type) return +1;
	else return 0;
}
#define NPSTEP 10000

static float rmin2,rmax2;
static lcparticletype *lcparticles;
static int nlcp,maxlcp,tnlcp;
/* if dist >= r1-maxd && dist < r1 --> located in inner buffer zone 
   if dist >= r1 && dist < r2 --> located in the main zone
   if dist >= r2 && dist < r2+maxd --> located in the outter buffer zone
   */
void ExtractingLightConeParticles(pmparticletype *pmparticles,
		int np,int nx,int ny,int nz,int nspace,
		int mstep, float amax, float a, float astep,float nextastep,
		float boxsize,
		float omep,float omeplam,float hubble,int myid,int nid,float maxd){
	pmparticletype *bp;
	float redshift;
	double r,dr1,dr2;
	double xp,yp,zp;
	double distsq,dist;
	float rminsq,rmaxsq;
	int nbox;
	int i,j,k;


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
	if(myid==0) printf("SURVEY %g %g %g at observer=%g %g %g\n",r,dr1,dr2,
			X0,Y0,Z0);
	nlcp = 0;
	maxlcp = NPSTEP;
	lcparticles = (lcparticletype *)Malloc(sizeof(lcparticletype)*maxlcp,
			PPTR(lcparticles));
	bp = pmparticles;
	for(i=0;i<np;i++){
		xp = XofP(bp)-X0;
		yp = YofP(bp)-Y0;
		zp = ZofP(bp)-Z0;
		distsq = xp*xp+yp*yp+zp*zp;
		if(distsq >= rminsq && distsq < rmaxsq){
			dist = sqrt(distsq);
			lcparticles[nlcp].x = xp;
			lcparticles[nlcp].y = yp;
			lcparticles[nlcp].z = zp;
			lcparticles[nlcp].vx = bp->vx;
			lcparticles[nlcp].vy = bp->vy;
			lcparticles[nlcp].vz = bp->vz;
			lcparticles[nlcp].bp = bp;
			lcparticles[nlcp].dist = dist;
#ifdef INDEX
			lcparticles[nlcp].indx = bp->indx;
#endif
			nlcp ++;
			if(nlcp >= maxlcp){
				maxlcp += NPSTEP;
				lcparticles = (lcparticletype*)Realloc(lcparticles,
						sizeof(lcparticletype)*maxlcp);

			}
		}
		bp++;
	}
	lcparticles = (lcparticletype*)Realloc(lcparticles,
		sizeof(lcparticletype)*nlcp);
	for(i=0;i<nlcp;i++){
		dist = lcparticles[i].dist;
		if(dist >=r1 && dist < r2) lcparticles[i].type = 0;
		else if(dist >=r2 && dist < r2+maxd) lcparticles[i].type = 1;
		else lcparticles[i].type = 2;
	}
	if(myid==0) printf("SURVEY %d detected\n",nlcp);
}
#define writedownlightconedata(A,B,C) {\
	FILE *wp;\
	float x0,y0,z0;\
	int nmlcp;\
	lightconetype *ttmp;\
	x0 = X0;\
	y0 = Y0;\
	z0 = Z0;\
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
			fwrite(&rmin2,sizeof(float),1,wp);\
			fwrite(&rmax2,sizeof(float),1,wp);\
			fwrite(&C,sizeof(int),1,wp);\
			fwrite(B,sizeof(lightconetype),C,wp);\
			for(i=1;i<nid;i++){	\
				MPI_Recv(&nmlcp,1,MPI_INT,i,i,MPI_COMM_WORLD,&status);\
				ttmp=(lightconetype*)Malloc(sizeof(lightconetype)*nmlcp,PPTR(ttmp));\
				MPI_Recv(ttmp,sizeof(lightconetype)*nmlcp,MPI_BYTE,i,i,\
						MPI_COMM_WORLD,&status);\
				fwrite(&nmlcp,sizeof(int),1,wp);\
				fwrite(ttmp,sizeof(lightconetype),nmlcp,wp);\
				Free(ttmp);\
			}\
			fclose(wp);\
		}\
		else  {\
			MPI_Send(&C,1,MPI_INT,0,myid,MPI_COMM_WORLD);\
			MPI_Send(B,sizeof(lightconetype)*C,MPI_BYTE,0,myid,MPI_COMM_WORLD);\
		}\
	}\
}
void SavingLightConeData(pmparticletype *pmparticles,int np,int nx,int ny,
		int nz, int mstep, float amax, float a, float astep,float nextastep,
		float boxsize,float omei,
		float omep,float omepb,float omeplam,float wlam0, float wlam1,float hubble,int myid,int nid){
	pmparticletype *bp;
	float redshift;
	double r,dr1,dr2;
	double xp,yp,zp,distsq,dist;
	float rminsq,rmaxsq;
	int nbox;
	int i,j,k;
	float ax,ay,az;
	int *nnp;
	lcparticletype *tp;
	lightconetype *lightcones;
	float ainv,app,apsq,afact,bfact,vfact1h,vfact2h,rng;
	float asteph;
	float vfact1,vfact2;
	MPI_Status status;

	MPI_Reduce(&nlcp,&tnlcp,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&tnlcp,1,MPI_INT,0,MPI_COMM_WORLD);
	if(tnlcp==0) return;


	asteph = astep*0.5;
	rng = nx;
	ainv = 1./a;
	/*
	app  = -4.*M_PI/3.*(ainv*ainv+(1.+3*wlam0)*omeplam/omep/pow(a,2.+3*wlam0)*pow(amax, 3*wlam0));
	apsq = 8.*M_PI/3.*(ainv+1./omei-1.+omeplam/omep*(pow(a,-1.-3*wlam0)-1.)*pow(amax,3*wlam0));
	*/
	appapsq_(&omep, &omeplam, &wlam0, &wlam1, &amax, &a, &app, &apsq);

	afact = (2.+app*a/apsq)/2.*ainv;
	bfact = ainv*ainv*ainv/apsq;
	vfact1 = (1.-afact*astep)/(1.+afact*astep);
	vfact2 = bfact*astep/(1.+afact*astep)*rng;
	vfact1h = (1.-afact*asteph)/(1.+afact*asteph);
	vfact2h = bfact*asteph/(1.+afact*asteph)*rng;

	vfact1h = (1.-afact*astep);
	vfact2h = bfact*astep*0.5*rng;

	nnp = (int *)Malloc(sizeof(int)*3,PPTR(nnp));
	nnp[0] = nnp[1] = nnp[2] =0;
	tp = lcparticles;
	for(i=0;i<nlcp;i++){
		bp = tp->bp;
		ax = (bp->vx-vfact1*tp->vx)/vfact2;
		ay = (bp->vy-vfact1*tp->vy)/vfact2;
		az = (bp->vz-vfact1*tp->vz)/vfact2;
		tp->vx = vfact1h*tp->vx + vfact2h*ax;
		tp->vy = vfact1h*tp->vy + vfact2h*ay;
		tp->vz = vfact1h*tp->vz + vfact2h*az;
		tp++;
	}
	for(i=0;i<nlcp;i++){
		nnp[lcparticles[i].type]++;
	}
	if(nlcp>0) qsort(lcparticles,nlcp,sizeof(lcparticletype),sorttype);
	if(nlcp>0) lightcones = (lightconetype *)Malloc(sizeof(lightconetype)*nlcp,
			PPTR(lightcones));
	for(i=0;i<nlcp;i++){
		lightcones[i].x = lcparticles[i].x;
		lightcones[i].y = lcparticles[i].y;
		lightcones[i].z = lcparticles[i].z;
		lightcones[i].vx = lcparticles[i].vx;
		lightcones[i].vy = lcparticles[i].vy;
		lightcones[i].vz = lcparticles[i].vz;
#ifdef INDEX
		lightcones[i].indx = lcparticles[i].indx;
#endif
	}
	{
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
	if(nlcp>0) Free(lightcones);
	Free(nnp);
	Free(lcparticles);
}
#undef npstep
#include <math.h>
#define EPS 1.0e-6
#define JMAX 20

float qsimp(float (*func)(float), float a, float b)
{
	float trapzd(float (*func)(float), float a, float b, int n);
	int j;
	float s,st,ost,os;

	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=trapzd(func,a,b,j);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		os=s;
		ost=st;
	}
	return 0.0;
}
#undef EPS
#undef JMAX
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
#define FUNC(x) ((*func)(x))

float trapzd(float (*func)(float), float a, float b, int n)
{
	float x,tnm,sum,del;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
