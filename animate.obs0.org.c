#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include<sys/times.h>
#include "Memory.h"
#include "pmheader.h"
#include "lightcone.h"
#include "animate.h"

double ODScaleHeight=50; /* Scale height of the optical depth */

ImageSize Image = Planetarium;

#define MIN(A,B) ((A) <(B) ? (A):(B))
#define PI 3.1415926535L
#define hPI 1.57079632677L

double Theta,Phi,DTheta,DPhi;
char surveyfilename[100];
int *freework;
enum boolean { NO = 00,YES = 01};
typedef struct Spherical{
	double r,theta,phi;
} Spherical;
InBox *box;
float omega0,H,lambda0;
float qsimp(float (*)(float),float,float);
float trapzd(float (*)(float), float , float , int );
int ssorttype(const void *, const void *);

#define NPSTEP 10000
#define twoPI (6.283185306L)
#define deg2rad (0.017453292516L)


#define ObsX(view)   ((view).pos.x)
#define ObsY(view)   ((view).pos.y)
#define ObsZ(view)   ((view).pos.z)

#define ObsVWX(view)   ((view).E3.x)
#define ObsVWY(view)   ((view).E3.y)
#define ObsVWZ(view)   ((view).E3.z)

#define ObsVWHX(view)   ((view).E1.x)
#define ObsVWHY(view)   ((view).E1.y)
#define ObsVWHZ(view)   ((view).E1.z)

void ObsAnimate(pmparticletype *pmparticles,int np, Viewer *view, int stepcount){
	float HOpticalDepth;
	int nx,ny,nz, nspace;
	int mstep, amax, a, astep, boxsize;
	float omep,omeplam,hubble;
	int myid,nid;
	pmparticletype *bp;
	float redshift;
	double xp,yp,zp;
	double rminsq,rmaxsq;
	int nbox,ibox;
	int i,j,k;
	long ix,iy,iz,nxy;
	double ox,oy,oz;
	float *img,*omg;
	double dtheta, dphi;
	double vwx,vwy,vwz;
	double vwhx,vwhy,vwhz;
	double vwvx,vwvy,vwvz;
	char outfile[190];
	int nowframe;

	nowframe = stepcount-1;
	TIMER_START(44)

	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;
	nspace = simpar.nspace;
	mstep = simpar.nstep;
	a = simpar.anow;
	astep = simpar.astep;
	boxsize = simpar.boxsize;
	omep = simpar.omep;
	omeplam = simpar.omeplam;
	hubble = simpar.hubble;
	myid = simpar.myid;
	nid = simpar.nid;

	HOpticalDepth = view[nowframe].HOpticalDepth;


	nbox = 27;
	box = (InBox *)Malloc(sizeof(InBox)*nbox,PPTR(box));
	img = (float *) Malloc(sizeof(float)*Image.nx*Image.ny,PPTR(img));
	omg = (float *) Malloc(sizeof(float)*Image.nx*Image.ny,PPTR(omg));

	for(i=0;i<Image.nx*Image.ny;i++) img[i] = 0.;
	for(k=0;k<3;k++) for(j=0;j<3;j++) for(i=0;i<3;i++) {
		box[i+3*(j+3*k)].nx = nx*(i-1);
		box[i+3*(j+3*k)].ny = ny*(j-1);
		box[i+3*(j+3*k)].nz = nz*(k-1);
	}
	/* (ox,oy,oz) position of the observer */
	ox = ObsX(view[nowframe]);
	oy = ObsY(view[nowframe]);
	oz = ObsZ(view[nowframe]);
	/* (vwx,vwy,vwz) directional cosine of view point */
	vwx = ObsVWX(view[nowframe]);
	vwy = ObsVWY(view[nowframe]);
	vwz = ObsVWZ(view[nowframe]);
	
	/* (vwhx,vwhy,vwhz) directional cosine of horizontal view frame */
	vwhx = ObsVWHX(view[nowframe]);
	vwhy = ObsVWHY(view[nowframe]);
	vwhz = ObsVWHZ(view[nowframe]);

	/* (vwvx,vwvy,vwvz) directional cosine of vertical view frame */
	vwvx = vwy*vwhz - vwz*vwhy;
	vwvy = vwz*vwhx - vwx*vwhz;
	vwvz = vwx*vwhy - vwy*vwhx;

	if(myid==0){
		printf("nowframe %d \n",nowframe);
		printf("WKWW %g %g %g  \n",ox,oy,oz);
		printf("AAAA %g %g %g  \n",vwx,vwy,vwz);
		printf("BBBB %g %g %g  \n",vwhx,vwhy,vwhz);
		printf("CCCC %g %g %g  \n",vwvx,vwvy,vwvz);
	}

	dtheta = (Image.nx/2)/(PI/2.);
	bp = pmparticles;
	{
		double theta, phi,r,r2;
		double a, b, c;
		double Xp,Yp,Zp;
		int mx,my;
		for(i=0;i<np;i++){
			Xp = XofP(bp);
			Yp = YofP(bp);
			Zp = ZofP(bp);
			for(j=0;j<nbox;j++){
				xp = Xp - (ox + box[j].nx);
				yp = Yp - (oy + box[j].ny);
				zp = Zp - (oz + box[j].nz);
				r2 = xp*xp+yp*yp+zp*zp;
				r = sqrt(r2);
				if(r > HOpticalDepth*5) continue;
				c = (xp*vwx + yp*vwy + zp*vwz);
				if(r ==0) theta = 0;
				else theta = acos(c/r);
				if(theta > hPI) continue;
				a = (xp*vwhx + yp*vwhy + zp*vwhz);
				b = (xp*vwvx + yp*vwvy + zp*vwvz);
				if(b ==0) {
					if(a > 0) phi = hPI;
					else if(a<0) phi = PI + hPI;
					else phi = 0;
				}
				else {
					phi = atan(a/b);
					if(b < 0) phi = phi + PI;
				}
				my = -(dtheta * theta)*cos(phi) + (Image.ny/2);
				mx = (dtheta * theta)*sin(phi) + (Image.nx/2);
				img[mx + Image.nx*my] += exp(-r/HOpticalDepth)/r2;
			}
			bp++;
		}
	}
	if(myid == 0){
		FILE *wp;
		MPI_Status status;
		sprintf(outfile,"%s.%.5d",view[nowframe].outfile,stepcount);
		wp = fopen(outfile,"w");
		fwrite(&Image.nx,sizeof(int),1,wp);
		fwrite(&Image.ny,sizeof(int),1,wp);
		for(i=1;i<nid;i++){
			MPI_Recv(omg,Image.nx*Image.ny,MPI_FLOAT,i,i,MPI_COMM_WORLD,&status);
			for(j=1;j<Image.nx*Image.ny;j++){
				img[j] += omg[j];
			}
		}
		fwrite(img,sizeof(float),Image.nx*Image.ny,wp);
		fclose(wp);
	}
	else {
		MPI_Send(img,Image.nx*Image.ny, MPI_FLOAT,0,myid,MPI_COMM_WORLD);
	}
	Free(omg); Free(img); Free(box);
	TIMER_STOP(44)
	if(simpar.myid==0) printf("CPU(making animated image) = %g second\n",ELAPSED_TIME(44));  
	if(myid==0) printf("Single frame is made at file=  %s\n",outfile);
	/*
	MPI_Finalize();
	exit(99);
	*/
}
