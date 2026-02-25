#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include<omp.h>

#include "pmheader.h"
#include "Memory.h"
#define pi (3.14159265354979L)

#define den(i,j,k) (den[(i)+mx*((j)+ny*(k))])
#define PTR(den,i,j,k) (den+(i)+mx*((j)+ny*(k)))
/*
#define icopysign(a) (signbit(a)? -1:1)
*/
#define MAX(a,b) ((a)>(b)?(a):(b))
#define icopysign(a) copysignf(1.,a)

#ifdef XYZDBL
#define FABS(a) fabs(a)
#define RINT(a) rint(a)
#else
#define FABS(a) fabsf(a)
#define RINT(a) rintf(a)
#endif


void fda4(pmparticletype *pmp, int np, float *den, int local_nz,int nstart_z, int ihflag1, int ihflag2){
	float anow,astep,amax;
	float omep,omeplam;
	long long i,j,k,ii,mx,slicesize;
	float rngx,rngy,rngz,rng;
	long long nx,ny,nz;
	float start_z;
	float fact1,fact2;
	float a,ainv,app,apsq,afact,bfact,fact3;
	float fact1h, fact2h,ffact1,ffact2;
	int ng,ngx,ngy,ngz;
	int *xabp0,*xabp1,*xabm1;
	int *yabp0,*yabp1,*yabm1;
	int *zabp0,*zabp1,*zabm1;

	int nz_per,dest,src,sendsize,recvsize,stag,rtag;
	MPI_Status status;
	int myid,nid,ierror;
	int *nxperiodic,*nyperiodic,*nzperiodic;
	int *nxp,*nyp,*nzp;
	int nspace;
	float omei,pfact;


	/* Retrieve values from the simpar */
	nx = simpar.nx; ny = simpar.ny; nz = simpar.nz;
	anow = simpar.anow; astep = simpar.astep;
	amax = simpar.amax; 
	omep = simpar.omep;
	omeplam = simpar.omeplam;
	omei = simpar.omei;
	nspace = simpar.nspace;

	mx = 2*(nx/2+1);
	slicesize = mx*ny;




	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);

	TIMER_START(55);
	nz_per = local_nz+7;


	nxperiodic = (int *)Malloc(sizeof(int)*(nx+32),PPTR(nxperiodic)); nxp = nxperiodic+16;
	nyperiodic = (int *)Malloc(sizeof(int)*(ny+32),PPTR(nyperiodic)); nyp = nyperiodic+16;
	nzperiodic = (int *)Malloc(sizeof(int)*(nz_per+32),PPTR(nzperiodic)); nzp = nzperiodic+16;
	for(i=-16;i<nx+16;i++) nxp[i] = (i+nx)%nx;
	for(i=-16;i<ny+16;i++) nyp[i] = (i+ny)%ny;
	for(i=-16;i<nz_per+16;i++) nzp[i] = (i+nz_per)%nz_per;

	xabm1 = nxp-1; xabp1 = nxp+1; 
	yabm1 = nyp-1; yabp1 = nyp+1; 
	zabm1 = nzp-1; zabp1 = nzp+1;

#ifndef VarPM
	sendsize = recvsize = slicesize*3;
	dest = (myid+1+nid)%nid;
	src = (myid-1+nid)%nid;
	stag = rtag = 1;
	MPI_Sendrecv(den+slicesize*(local_nz-3),sendsize,MPI_FLOAT,dest,
			stag,den+slicesize*(nz_per-3),recvsize,MPI_FLOAT,src,tag,
			MPI_COMM_WORLD,&status);
	sendsize = recvsize = slicesize*4;
	dest = (myid-1+nid)%nid;
	src = (myid+1+nid)%nid;
	MPI_Sendrecv(den, sendsize,MPI_FLOAT,dest,stag,
			den+slicesize*local_nz,recvsize,MPI_FLOAT,src,rtag,
			MPI_COMM_WORLD,&staus);
#endif
	ngx = ng = nx;
	ngy = ny;
	ngz = nz;
	rng = ng;
	rngx = ngx;
	rngy = ngy;
	rngz = ngz;
	float wlam0 = simpar.wlam0;
	float wlam1 = simpar.wlam1;

	if(anow == 1.) astep =  0.5*astep;
	else if(fabs(anow-1.) < 1.e-5) astep =  0.5*astep;

	ainv = 1./anow;

	/*
	app = -4.*pi/3.*(ainv*ainv-2.*omeplam/omep*anow/amax/amax/amax);
	apsq = 8.*pi/3.*(ainv+1./omei-1.+omeplam/omep*(anow*anow-1.)/amax/amax/amax);
	*/

	/*
	app = -4.*M_PI/3.*(ainv*ainv+(1.+3*wlam0)*omeplam/omep/pow(anow,2.+3*wlam0)*pow(amax,3*wlam0));
	apsq = 8.*M_PI/3.*(ainv+1./omei-1.+omeplam/omep*(pow(anow,-1.-3*wlam0)-1.)*pow(amax,3*wlam0));
	*/
	
	appapsq_(&omep, &omeplam, &wlam0, &wlam1, &amax, &anow, &app, &apsq);
	afact = (2.+app*anow/apsq)/2.*ainv;
	bfact = ainv*ainv*ainv/apsq;
	fact1 = (1.-afact*astep)/(1.+afact*astep);
	fact1h = (1.-afact*astep*0.5)/(1.+afact*astep*0.5);
	fact2 = bfact*astep/(1.+afact*astep)*rng;
	fact2h = bfact*astep*0.5/(1.+afact*astep*0.5)*rng;
	fact3 = fact2/fact1;
	if(ihflag1 == 1 && ihflag2 == 0) {
		/*
		fact1 = fact1h;
		fact2 = fact2h;
		*/
		fact1 = (1.-afact*astep);
		fact2 = bfact*astep*0.5*rng;
	}
	else if(ihflag1 ==1 && ihflag2 == 1){
		/*
		fact1 = fact1/fact1h;
		fact2 = (fact2-fact1*fact2h);
		*/
		fact1 = 1./(1.+afact*astep);
		fact2 = bfact*astep*0.5/(1.+afact*astep)*rng;
	}
	if(anow == 1.) astep = 2.*astep;
	else if(fabs(anow-1.) < 1.e-5) astep =  2.*astep;

	pfact = rng*astep;

	/* Update the simulation parameters */
	simpar.pfact = pfact;
	simpar.fact1 = fact1;
	simpar.fact2 = fact2;

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<np;i++){
		pmp[i].vx = pmp[i].vx*fact1;
		pmp[i].vy = pmp[i].vy*fact1;
		pmp[i].vz = pmp[i].vz*fact1;
	}
	start_z = nstart_z;
	ffact1 = (4./3.)/2.;
	ffact2 = (-1./3.)/4.; 

#ifdef NOT_USE_GPU
	{
		void call_cuda_fda4(int,int,int ,int , int , float *,
				pmparticletype *, int , float , float );
		call_cuda_fda4(myid,nx,mx,ny, local_nz, den,pmp, np, start_z, fact2);
	}

#else

#	ifdef _OPENMP
#	pragma omp parallel for 
#	endif
	for(i=0;i<np;i++){
#		ifdef XYZDBL
		double xp,yp,zp;
#		else
		float xp,yp,zp;
#		endif
		float xmin,ymin,zmin,xd1,yd1,zd1;
		float wx1,wy1,wz1,wx2,wy2,wz2,wx3,wy3,wz3;
		float wt1,wt2,wt3,wt4,wt5,wt6,wt7,wt8,wt9,wt10;
		float wt11,wt12,wt13,wt14,wt15,wt16,wt17,wt18,wt19,wt20;
		float wt21,wt22,wt23,wt24,wt25,wt26,wt27;
		float fx,fy,fz;
		float fx1,fx2,fx3,fx4,fx5,fx6,fx7,fx8,fx9,fx10;
		float fx11,fx12,fx13,fx14,fx15,fx16,fx17,fx18,fx19,fx20;
		float fx21,fx22,fx23,fx24,fx25,fx26,fx27;
		float fy1,fy2,fy3,fy4,fy5,fy6,fy7,fy8,fy9,fy10;
		float fy11,fy12,fy13,fy14,fy15,fy16,fy17,fy18,fy19,fy20;
		float fy21,fy22,fy23,fy24,fy25,fy26,fy27;
		float fz1,fz2,fz3,fz4,fz5,fz6,fz7,fz8,fz9,fz10;
		float fz11,fz12,fz13,fz14,fz15,fz16,fz17,fz18,fz19,fz20;
		float fz21,fz22,fz23,fz24,fz25,fz26,fz27;
		int nearx,neary,nearz,xsign,ysign,zsign;
		int i1,j1,k1,i2,j2,k2,i3,j3,k3;
		float wx1wy1, wx2wy1, wx3wy1, wx1wy2, wx2wy2, wx3wy2, wx1wy3, wx2wy3, wx3wy3;
		int i1m1,i1m2,i1p1,i1p2,i2m1,i2m2,i2p1,i2p2,i3m1,i3m2,i3p1;
		int i3p2,j1m1,j1m2,j1p1,j1p2,j2m1,j2m2,j2p1,j2p2,j3m1,j3m2;
		int j3p1,j3p2,k1m1,k1m2,k1p1,k1p2,k2m1,k2m2,k2p1,k2p2,k3m1;
		int k3m2,k3p1,k3p2;
		int jk;

		xp = XofP(pmp+i);
		yp = YofP(pmp+i);
		zp = ZofP(pmp+i)-start_z;
		nearx = RINT(xp);
		neary = RINT(yp);
		nearz = RINT(zp);
		xmin = xp - nearx;
		ymin = yp - neary;
		zmin = zp - nearz;
		xsign = icopysign(xmin);
		ysign = icopysign(ymin);
		zsign = icopysign(zmin);

		i1 = nxp[nearx];
		i2 = nxp[nearx+xsign];
		i3 = nxp[nearx-xsign];

		j1 = nyp[neary];
		j2 = nyp[neary+ysign];
		j3 = nyp[neary-ysign];

		k1 = nzp[nearz];
		k2 = nzp[nearz+zsign];
		k3 = nzp[nearz-zsign];

		xd1 = FABS(xmin);
		yd1 = FABS(ymin);
		zd1 = FABS(zmin);


		wx1 =  0.75-xd1*xd1;
		wy1 =  0.75-yd1*yd1;
		wz1 =  0.75-zd1*zd1;

		wx3 = 0.5*(0.25+xd1*(xd1-1.));
		wx2 = wx3 + xd1;

		wy3 = 0.5*(0.25+yd1*(yd1-1.));
		wy2 = wy3 + yd1;

		wz3 = 0.5*(0.25+zd1*(zd1-1.));
		wz2 = wz3 + zd1;

		wx1wy1 = wx1*wy1;
		wx2wy1 = wx2*wy1;
		wx3wy1 = wx3*wy1;
		wx1wy2 = wx1*wy2;
		wx2wy2 = wx2*wy2;
		wx3wy2 = wx3*wy2;
		wx1wy3 = wx1*wy3;
		wx2wy3 = wx2*wy3;
		wx3wy3 = wx3*wy3;

		wt1 = wx1wy1*wz1;
		wt2 = wx2wy1*wz1;
		wt3 = wx3wy1*wz1;
		wt4 = wx1wy2*wz1;
		wt5 = wx2wy2*wz1;
		wt6 = wx3wy2*wz1;
		wt7 = wx1wy3*wz1;
		wt8 = wx2wy3*wz1;
		wt9 = wx3wy3*wz1;

		wt10 = wx1wy1*wz2;
		wt11 = wx2wy1*wz2;
		wt12 = wx3wy1*wz2;
		wt13 = wx1wy2*wz2;
		wt14 = wx2wy2*wz2;
		wt15 = wx3wy2*wz2;
		wt16 = wx1wy3*wz2;
		wt17 = wx2wy3*wz2;
		wt18 = wx3wy3*wz2;

		wt19 = wx1wy1*wz3;
		wt20 = wx2wy1*wz3;
		wt21 = wx3wy1*wz3;
		wt22 = wx1wy2*wz3;
		wt23 = wx2wy2*wz3;
		wt24 = wx3wy2*wz3;
		wt25 = wx1wy3*wz3;
		wt26 = wx2wy3*wz3;
		wt27 = wx3wy3*wz3;

		i1m1 = xabm1[i1];
		i1m2 = xabm1[i1m1];
		i1p1 = xabp1[i1];
		i1p2 = xabp1[i1p1];
		i2m1 = xabm1[i2];
		i2m2 = xabm1[i2m1];
		i2p1 = xabp1[i2];
		i2p2 = xabp1[i2p1];
		i3m1 = xabm1[i3];
		i3m2 = xabm1[i3m1];
		i3p1 = xabp1[i3];
		i3p2 = xabp1[i3p1];

		j1m1 = yabm1[j1];
		j1m2 = yabm1[j1m1];
		j1p1 = yabp1[j1];
		j1p2 = yabp1[j1p1];
		j2m1 = yabm1[j2];
		j2m2 = yabm1[j2m1];
		j2p1 = yabp1[j2];
		j2p2 = yabp1[j2p1];
		j3m1 = yabm1[j3];
		j3m2 = yabm1[j3m1];
		j3p1 = yabp1[j3];
		j3p2 = yabp1[j3p1];

		k1m1 = zabm1[k1];
		k1m2 = zabm1[k1m1];
		k1p1 = zabp1[k1];
		k1p2 = zabp1[k1p1];
		k2m1 = zabm1[k2];
		k2m2 = zabm1[k2m1];
		k2p1 = zabp1[k2];
		k2p2 = zabp1[k2p1];
		k3m1 = zabm1[k3];
		k3m2 = zabm1[k3m1];
		k3p1 = zabp1[k3];
		k3p2 = zabp1[k3p1];

#ifndef FAST
		{
			long j1k1,j1k2,j1k3,j2k1,j2k2,j2k3,j3k1,j3k2,j3k3;
			j1k1 = mx*(j1+ny*k1);
			fx1 = ffact1*(den[i1m1+j1k1]-den[i1p1+j1k1])+ffact2*(den[i1m2+j1k1]-den[i1p2+j1k1]);
			fx2 = ffact1*(den[i2m1+j1k1]-den[i2p1+j1k1])+ffact2*(den[i2m2+j1k1]-den[i2p2+j1k1]);
			fx3 = ffact1*(den[i3m1+j1k1]-den[i3p1+j1k1])+ffact2*(den[i3m2+j1k1]-den[i3p2+j1k1]);
			j2k1 = mx*(j2+ny*k1);
			fx4 = ffact1*(den[i1m1+j2k1]-den[i1p1+j2k1])+ffact2*(den[i1m2+j2k1]-den[i1p2+j2k1]);
			fx5 = ffact1*(den[i2m1+j2k1]-den[i2p1+j2k1])+ffact2*(den[i2m2+j2k1]-den[i2p2+j2k1]);
			fx6 = ffact1*(den[i3m1+j2k1]-den[i3p1+j2k1])+ffact2*(den[i3m2+j2k1]-den[i3p2+j2k1]);
			j3k1 = mx*(j3+ny*k1);
			fx7 = ffact1*(den[i1m1+j3k1]-den[i1p1+j3k1])+ffact2*(den[i1m2+j3k1]-den[i1p2+j3k1]);
			fx8 = ffact1*(den[i2m1+j3k1]-den[i2p1+j3k1])+ffact2*(den[i2m2+j3k1]-den[i2p2+j3k1]);
			fx9 = ffact1*(den[i3m1+j3k1]-den[i3p1+j3k1])+ffact2*(den[i3m2+j3k1]-den[i3p2+j3k1]);
			j1k2 = mx*(j1+ny*k2);
			fx10 = ffact1*(den[i1m1+j1k2]-den[i1p1+j1k2])+ffact2*(den[i1m2+j1k2]-den[i1p2+j1k2]);
			fx11 = ffact1*(den[i2m1+j1k2]-den[i2p1+j1k2])+ffact2*(den[i2m2+j1k2]-den[i2p2+j1k2]);
			fx12 = ffact1*(den[i3m1+j1k2]-den[i3p1+j1k2])+ffact2*(den[i3m2+j1k2]-den[i3p2+j1k2]);
			j2k2 = mx*(j2+ny*k2);
			fx13 = ffact1*(den[i1m1+j2k2]-den[i1p1+j2k2])+ffact2*(den[i1m2+j2k2]-den[i1p2+j2k2]);
			fx14 = ffact1*(den[i2m1+j2k2]-den[i2p1+j2k2])+ffact2*(den[i2m2+j2k2]-den[i2p2+j2k2]);
			fx15 = ffact1*(den[i3m1+j2k2]-den[i3p1+j2k2])+ffact2*(den[i3m2+j2k2]-den[i3p2+j2k2]);
			j3k2 = mx*(j3+ny*k2);
			fx16 = ffact1*(den[i1m1+j3k2]-den[i1p1+j3k2])+ffact2*(den[i1m2+j3k2]-den[i1p2+j3k2]);
			fx17 = ffact1*(den[i2m1+j3k2]-den[i2p1+j3k2])+ffact2*(den[i2m2+j3k2]-den[i2p2+j3k2]);
			fx18 = ffact1*(den[i3m1+j3k2]-den[i3p1+j3k2])+ffact2*(den[i3m2+j3k2]-den[i3p2+j3k2]);
			j1k3 = mx*(j1+ny*k3);
			fx19 = ffact1*(den[i1m1+j1k3]-den[i1p1+j1k3])+ffact2*(den[i1m2+j1k3]-den[i1p2+j1k3]);
			fx20 = ffact1*(den[i2m1+j1k3]-den[i2p1+j1k3])+ffact2*(den[i2m2+j1k3]-den[i2p2+j1k3]);
			fx21 = ffact1*(den[i3m1+j1k3]-den[i3p1+j1k3])+ffact2*(den[i3m2+j1k3]-den[i3p2+j1k3]);
			j2k3 = mx*(j2+ny*k3);
			fx22 = ffact1*(den[i1m1+j2k3]-den[i1p1+j2k3])+ffact2*(den[i1m2+j2k3]-den[i1p2+j2k3]);
			fx23 = ffact1*(den[i2m1+j2k3]-den[i2p1+j2k3])+ffact2*(den[i2m2+j2k3]-den[i2p2+j2k3]);
			fx24 = ffact1*(den[i3m1+j2k3]-den[i3p1+j2k3])+ffact2*(den[i3m2+j2k3]-den[i3p2+j2k3]);
			j3k3 = mx*(j3+ny*k3);
			fx25 = ffact1*(den[i1m1+j3k3]-den[i1p1+j3k3])+ffact2*(den[i1m2+j3k3]-den[i1p2+j3k3]);
			fx26 = ffact1*(den[i2m1+j3k3]-den[i2p1+j3k3])+ffact2*(den[i2m2+j3k3]-den[i2p2+j3k3]);
			fx27 = ffact1*(den[i3m1+j3k3]-den[i3p1+j3k3])+ffact2*(den[i3m2+j3k3]-den[i3p2+j3k3]);
		}
		{
			long jk1,jk2,jk3,jk4;
			jk1 = mx*(j1m1+ny*k1);
			jk2 = mx*(j1p1+ny*k1);
			jk3 = mx*(j1m2+ny*k1);
			jk4 = mx*(j1p2+ny*k1);
			fy1 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fy2 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fy3 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j2m1+ny*k1);
			jk2 = mx*(j2p1+ny*k1);
			jk3 = mx*(j2m2+ny*k1);
			jk4 = mx*(j2p2+ny*k1);
			fy4 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fy5 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fy6 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j3m1+ny*k1);
			jk2 = mx*(j3p1+ny*k1);
			jk3 = mx*(j3m2+ny*k1);
			jk4 = mx*(j3p2+ny*k1);
			fy7 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fy8 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fy9 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j1m1+ny*k2);
			jk2 = mx*(j1p1+ny*k2);
			jk3 = mx*(j1m2+ny*k2);
			jk4 = mx*(j1p2+ny*k2);
			fy10 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fy11 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fy12 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j2m1+ny*k2);
			jk2 = mx*(j2p1+ny*k2);
			jk3 = mx*(j2m2+ny*k2);
			jk4 = mx*(j2p2+ny*k2);
			fy13 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fy14 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fy15 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j3m1+ny*k2);
			jk2 = mx*(j3p1+ny*k2);
			jk3 = mx*(j3m2+ny*k2);
			jk4 = mx*(j3p2+ny*k2);
			fy16 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fy17 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fy18 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j1m1+ny*k3);
			jk2 = mx*(j1p1+ny*k3);
			jk3 = mx*(j1m2+ny*k3);
			jk4 = mx*(j1p2+ny*k3);
			fy19 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fy20 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fy21 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j2m1+ny*k3);
			jk2 = mx*(j2p1+ny*k3);
			jk3 = mx*(j2m2+ny*k3);
			jk4 = mx*(j2p2+ny*k3);
			fy22 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fy23 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fy24 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j3m1+ny*k3);
			jk2 = mx*(j3p1+ny*k3);
			jk3 = mx*(j3m2+ny*k3);
			jk4 = mx*(j3p2+ny*k3);
			fy25 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fy26 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fy27 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
		}
		{
			long jk1,jk2,jk3,jk4;
			jk1 = mx*(j1+ny*k1m1);
			jk2 = mx*(j1+ny*k1p1);
			jk3 = mx*(j1+ny*k1m2);
			jk4 = mx*(j1+ny*k1p2);
			fz1 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fz2 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fz3 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j2+ny*k1m1);
			jk2 = mx*(j2+ny*k1p1);
			jk3 = mx*(j2+ny*k1m2);
			jk4 = mx*(j2+ny*k1p2);
			fz4 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fz5 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fz6 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j3+ny*k1m1);
			jk2 = mx*(j3+ny*k1p1);
			jk3 = mx*(j3+ny*k1m2);
			jk4 = mx*(j3+ny*k1p2);
			fz7 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fz8 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fz9 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j1+ny*k2m1);
			jk2 = mx*(j1+ny*k2p1);
			jk3 = mx*(j1+ny*k2m2);
			jk4 = mx*(j1+ny*k2p2);
			fz10 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fz11 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fz12 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j2+ny*k2m1);
			jk2 = mx*(j2+ny*k2p1);
			jk3 = mx*(j2+ny*k2m2);
			jk4 = mx*(j2+ny*k2p2);
			fz13 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fz14 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fz15 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j3+ny*k2m1);
			jk2 = mx*(j3+ny*k2p1);
			jk3 = mx*(j3+ny*k2m2);
			jk4 = mx*(j3+ny*k2p2);
			fz16 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fz17 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fz18 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j1+ny*k3m1);
			jk2 = mx*(j1+ny*k3p1);
			jk3 = mx*(j1+ny*k3m2);
			jk4 = mx*(j1+ny*k3p2);
			fz19 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fz20 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fz21 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j2+ny*k3m1);
			jk2 = mx*(j2+ny*k3p1);
			jk3 = mx*(j2+ny*k3m2);
			jk4 = mx*(j2+ny*k3p2);
			fz22 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fz23 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fz24 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
			jk1 = mx*(j3+ny*k3m1);
			jk2 = mx*(j3+ny*k3p1);
			jk3 = mx*(j3+ny*k3m2);
			jk4 = mx*(j3+ny*k3p2);
			fz25 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
			fz26 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
			fz27 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
		}
#else
		fx1 = ffact1*(den(i1m1,j1,k1)-den(i1p1,j1,k1))+ffact2*(den(i1m2,j1,k1)-den(i1p2,j1,k1));
		fx2 = ffact1*(den(i2m1,j1,k1)-den(i2p1,j1,k1))+ffact2*(den(i2m2,j1,k1)-den(i2p2,j1,k1));
		fx3 = ffact1*(den(i3m1,j1,k1)-den(i3p1,j1,k1))+ffact2*(den(i3m2,j1,k1)-den(i3p2,j1,k1));
		fx4 = ffact1*(den(i1m1,j2,k1)-den(i1p1,j2,k1))+ffact2*(den(i1m2,j2,k1)-den(i1p2,j2,k1));
		fx5 = ffact1*(den(i2m1,j2,k1)-den(i2p1,j2,k1))+ffact2*(den(i2m2,j2,k1)-den(i2p2,j2,k1));
		fx6 = ffact1*(den(i3m1,j2,k1)-den(i3p1,j2,k1))+ffact2*(den(i3m2,j2,k1)-den(i3p2,j2,k1));
		fx7 = ffact1*(den(i1m1,j3,k1)-den(i1p1,j3,k1))+ffact2*(den(i1m2,j3,k1)-den(i1p2,j3,k1));
		fx8 = ffact1*(den(i2m1,j3,k1)-den(i2p1,j3,k1))+ffact2*(den(i2m2,j3,k1)-den(i2p2,j3,k1));
		fx9 = ffact1*(den(i3m1,j3,k1)-den(i3p1,j3,k1))+ffact2*(den(i3m2,j3,k1)-den(i3p2,j3,k1));
		fx10 = ffact1*(den(i1m1,j1,k2)-den(i1p1,j1,k2))+ffact2*(den(i1m2,j1,k2)-den(i1p2,j1,k2));
		fx11 = ffact1*(den(i2m1,j1,k2)-den(i2p1,j1,k2))+ffact2*(den(i2m2,j1,k2)-den(i2p2,j1,k2));
		fx12 = ffact1*(den(i3m1,j1,k2)-den(i3p1,j1,k2))+ffact2*(den(i3m2,j1,k2)-den(i3p2,j1,k2));
		fx13 = ffact1*(den(i1m1,j2,k2)-den(i1p1,j2,k2))+ffact2*(den(i1m2,j2,k2)-den(i1p2,j2,k2));
		fx14 = ffact1*(den(i2m1,j2,k2)-den(i2p1,j2,k2))+ffact2*(den(i2m2,j2,k2)-den(i2p2,j2,k2));
		fx15 = ffact1*(den(i3m1,j2,k2)-den(i3p1,j2,k2))+ffact2*(den(i3m2,j2,k2)-den(i3p2,j2,k2));
		fx16 = ffact1*(den(i1m1,j3,k2)-den(i1p1,j3,k2))+ffact2*(den(i1m2,j3,k2)-den(i1p2,j3,k2));
		fx17 = ffact1*(den(i2m1,j3,k2)-den(i2p1,j3,k2))+ffact2*(den(i2m2,j3,k2)-den(i2p2,j3,k2));
		fx18 = ffact1*(den(i3m1,j3,k2)-den(i3p1,j3,k2))+ffact2*(den(i3m2,j3,k2)-den(i3p2,j3,k2));
		fx19 = ffact1*(den(i1m1,j1,k3)-den(i1p1,j1,k3))+ffact2*(den(i1m2,j1,k3)-den(i1p2,j1,k3));
		fx20 = ffact1*(den(i2m1,j1,k3)-den(i2p1,j1,k3))+ffact2*(den(i2m2,j1,k3)-den(i2p2,j1,k3));
		fx21 = ffact1*(den(i3m1,j1,k3)-den(i3p1,j1,k3))+ffact2*(den(i3m2,j1,k3)-den(i3p2,j1,k3));
		fx22 = ffact1*(den(i1m1,j2,k3)-den(i1p1,j2,k3))+ffact2*(den(i1m2,j2,k3)-den(i1p2,j2,k3));
		fx23 = ffact1*(den(i2m1,j2,k3)-den(i2p1,j2,k3))+ffact2*(den(i2m2,j2,k3)-den(i2p2,j2,k3));
		fx24 = ffact1*(den(i3m1,j2,k3)-den(i3p1,j2,k3))+ffact2*(den(i3m2,j2,k3)-den(i3p2,j2,k3));
		fx25 = ffact1*(den(i1m1,j3,k3)-den(i1p1,j3,k3))+ffact2*(den(i1m2,j3,k3)-den(i1p2,j3,k3));
		fx26 = ffact1*(den(i2m1,j3,k3)-den(i2p1,j3,k3))+ffact2*(den(i2m2,j3,k3)-den(i2p2,j3,k3));
		fx27 = ffact1*(den(i3m1,j3,k3)-den(i3p1,j3,k3))+ffact2*(den(i3m2,j3,k3)-den(i3p2,j3,k3));

		fy1 = ffact1*(den(i1,j1m1,k1)-den(i1,j1p1,k1))+ffact2*(den(i1,j1m2,k1)-den(i1,j1p2,k1));
		fy2 = ffact1*(den(i2,j1m1,k1)-den(i2,j1p1,k1))+ffact2*(den(i2,j1m2,k1)-den(i2,j1p2,k1));
		fy3 = ffact1*(den(i3,j1m1,k1)-den(i3,j1p1,k1))+ffact2*(den(i3,j1m2,k1)-den(i3,j1p2,k1));
		fy4 = ffact1*(den(i1,j2m1,k1)-den(i1,j2p1,k1))+ffact2*(den(i1,j2m2,k1)-den(i1,j2p2,k1));
		fy5 = ffact1*(den(i2,j2m1,k1)-den(i2,j2p1,k1))+ffact2*(den(i2,j2m2,k1)-den(i2,j2p2,k1));
		fy6 = ffact1*(den(i3,j2m1,k1)-den(i3,j2p1,k1))+ffact2*(den(i3,j2m2,k1)-den(i3,j2p2,k1));
		fy7 = ffact1*(den(i1,j3m1,k1)-den(i1,j3p1,k1))+ffact2*(den(i1,j3m2,k1)-den(i1,j3p2,k1));
		fy8 = ffact1*(den(i2,j3m1,k1)-den(i2,j3p1,k1))+ffact2*(den(i2,j3m2,k1)-den(i2,j3p2,k1));
		fy9 = ffact1*(den(i3,j3m1,k1)-den(i3,j3p1,k1))+ffact2*(den(i3,j3m2,k1)-den(i3,j3p2,k1));
		fy10 = ffact1*(den(i1,j1m1,k2)-den(i1,j1p1,k2))+ffact2*(den(i1,j1m2,k2)-den(i1,j1p2,k2));
		fy11 = ffact1*(den(i2,j1m1,k2)-den(i2,j1p1,k2))+ffact2*(den(i2,j1m2,k2)-den(i2,j1p2,k2));
		fy12 = ffact1*(den(i3,j1m1,k2)-den(i3,j1p1,k2))+ffact2*(den(i3,j1m2,k2)-den(i3,j1p2,k2));
		fy13 = ffact1*(den(i1,j2m1,k2)-den(i1,j2p1,k2))+ffact2*(den(i1,j2m2,k2)-den(i1,j2p2,k2));
		fy14 = ffact1*(den(i2,j2m1,k2)-den(i2,j2p1,k2))+ffact2*(den(i2,j2m2,k2)-den(i2,j2p2,k2));
		fy15 = ffact1*(den(i3,j2m1,k2)-den(i3,j2p1,k2))+ffact2*(den(i3,j2m2,k2)-den(i3,j2p2,k2));
		fy16 = ffact1*(den(i1,j3m1,k2)-den(i1,j3p1,k2))+ffact2*(den(i1,j3m2,k2)-den(i1,j3p2,k2));
		fy17 = ffact1*(den(i2,j3m1,k2)-den(i2,j3p1,k2))+ffact2*(den(i2,j3m2,k2)-den(i2,j3p2,k2));
		fy18 = ffact1*(den(i3,j3m1,k2)-den(i3,j3p1,k2))+ffact2*(den(i3,j3m2,k2)-den(i3,j3p2,k2));
		fy19 = ffact1*(den(i1,j1m1,k3)-den(i1,j1p1,k3))+ffact2*(den(i1,j1m2,k3)-den(i1,j1p2,k3));
		fy20 = ffact1*(den(i2,j1m1,k3)-den(i2,j1p1,k3))+ffact2*(den(i2,j1m2,k3)-den(i2,j1p2,k3));
		fy21 = ffact1*(den(i3,j1m1,k3)-den(i3,j1p1,k3))+ffact2*(den(i3,j1m2,k3)-den(i3,j1p2,k3));
		fy22 = ffact1*(den(i1,j2m1,k3)-den(i1,j2p1,k3))+ffact2*(den(i1,j2m2,k3)-den(i1,j2p2,k3));
		fy23 = ffact1*(den(i2,j2m1,k3)-den(i2,j2p1,k3))+ffact2*(den(i2,j2m2,k3)-den(i2,j2p2,k3));
		fy24 = ffact1*(den(i3,j2m1,k3)-den(i3,j2p1,k3))+ffact2*(den(i3,j2m2,k3)-den(i3,j2p2,k3));
		fy25 = ffact1*(den(i1,j3m1,k3)-den(i1,j3p1,k3))+ffact2*(den(i1,j3m2,k3)-den(i1,j3p2,k3));
		fy26 = ffact1*(den(i2,j3m1,k3)-den(i2,j3p1,k3))+ffact2*(den(i2,j3m2,k3)-den(i2,j3p2,k3));
		fy27 = ffact1*(den(i3,j3m1,k3)-den(i3,j3p1,k3))+ffact2*(den(i3,j3m2,k3)-den(i3,j3p2,k3));

		fz1 = ffact1*(den(i1,j1,k1m1)-den(i1,j1,k1p1))+ffact2*(den(i1,j1,k1m2)-den(i1,j1,k1p2));
		fz2 = ffact1*(den(i2,j1,k1m1)-den(i2,j1,k1p1))+ffact2*(den(i2,j1,k1m2)-den(i2,j1,k1p2));
		fz3 = ffact1*(den(i3,j1,k1m1)-den(i3,j1,k1p1))+ffact2*(den(i3,j1,k1m2)-den(i3,j1,k1p2));
		fz4 = ffact1*(den(i1,j2,k1m1)-den(i1,j2,k1p1))+ffact2*(den(i1,j2,k1m2)-den(i1,j2,k1p2));
		fz5 = ffact1*(den(i2,j2,k1m1)-den(i2,j2,k1p1))+ffact2*(den(i2,j2,k1m2)-den(i2,j2,k1p2));
		fz6 = ffact1*(den(i3,j2,k1m1)-den(i3,j2,k1p1))+ffact2*(den(i3,j2,k1m2)-den(i3,j2,k1p2));
		fz7 = ffact1*(den(i1,j3,k1m1)-den(i1,j3,k1p1))+ffact2*(den(i1,j3,k1m2)-den(i1,j3,k1p2));
		fz8 = ffact1*(den(i2,j3,k1m1)-den(i2,j3,k1p1))+ffact2*(den(i2,j3,k1m2)-den(i2,j3,k1p2));
		fz9 = ffact1*(den(i3,j3,k1m1)-den(i3,j3,k1p1))+ffact2*(den(i3,j3,k1m2)-den(i3,j3,k1p2));
		fz10 = ffact1*(den(i1,j1,k2m1)-den(i1,j1,k2p1))+ffact2*(den(i1,j1,k2m2)-den(i1,j1,k2p2));
		fz11 = ffact1*(den(i2,j1,k2m1)-den(i2,j1,k2p1))+ffact2*(den(i2,j1,k2m2)-den(i2,j1,k2p2));
		fz12 = ffact1*(den(i3,j1,k2m1)-den(i3,j1,k2p1))+ffact2*(den(i3,j1,k2m2)-den(i3,j1,k2p2));
		fz13 = ffact1*(den(i1,j2,k2m1)-den(i1,j2,k2p1))+ffact2*(den(i1,j2,k2m2)-den(i1,j2,k2p2));
		fz14 = ffact1*(den(i2,j2,k2m1)-den(i2,j2,k2p1))+ffact2*(den(i2,j2,k2m2)-den(i2,j2,k2p2));
		fz15 = ffact1*(den(i3,j2,k2m1)-den(i3,j2,k2p1))+ffact2*(den(i3,j2,k2m2)-den(i3,j2,k2p2));
		fz16 = ffact1*(den(i1,j3,k2m1)-den(i1,j3,k2p1))+ffact2*(den(i1,j3,k2m2)-den(i1,j3,k2p2));
		fz17 = ffact1*(den(i2,j3,k2m1)-den(i2,j3,k2p1))+ffact2*(den(i2,j3,k2m2)-den(i2,j3,k2p2));
		fz18 = ffact1*(den(i3,j3,k2m1)-den(i3,j3,k2p1))+ffact2*(den(i3,j3,k2m2)-den(i3,j3,k2p2));
		fz19 = ffact1*(den(i1,j1,k3m1)-den(i1,j1,k3p1))+ffact2*(den(i1,j1,k3m2)-den(i1,j1,k3p2));
		fz20 = ffact1*(den(i2,j1,k3m1)-den(i2,j1,k3p1))+ffact2*(den(i2,j1,k3m2)-den(i2,j1,k3p2));
		fz21 = ffact1*(den(i3,j1,k3m1)-den(i3,j1,k3p1))+ffact2*(den(i3,j1,k3m2)-den(i3,j1,k3p2));
		fz22 = ffact1*(den(i1,j2,k3m1)-den(i1,j2,k3p1))+ffact2*(den(i1,j2,k3m2)-den(i1,j2,k3p2));
		fz23 = ffact1*(den(i2,j2,k3m1)-den(i2,j2,k3p1))+ffact2*(den(i2,j2,k3m2)-den(i2,j2,k3p2));
		fz24 = ffact1*(den(i3,j2,k3m1)-den(i3,j2,k3p1))+ffact2*(den(i3,j2,k3m2)-den(i3,j2,k3p2));
		fz25 = ffact1*(den(i1,j3,k3m1)-den(i1,j3,k3p1))+ffact2*(den(i1,j3,k3m2)-den(i1,j3,k3p2));
   		fz26 = ffact1*(den(i2,j3,k3m1)-den(i2,j3,k3p1))+ffact2*(den(i2,j3,k3m2)-den(i2,j3,k3p2));
  		fz27 = ffact1*(den(i3,j3,k3m1)-den(i3,j3,k3p1))+ffact2*(den(i3,j3,k3m2)-den(i3,j3,k3p2));
#endif

		fx=(wt1*fx1+wt2*fx2+wt3*fx3+wt4*fx4  +wt5 *fx5 +wt6 *fx6+
				wt7 *fx7 +wt8 *fx8 +wt9 *fx9 +wt10*fx10+wt11*fx11+wt12*fx12+
				wt13*fx13+wt14*fx14+wt15*fx15+wt16*fx16+wt17*fx17+wt18*fx18+
				wt19*fx19+wt20*fx20+wt21*fx21+wt22*fx22+wt23*fx23+wt24*fx24+
				wt25*fx25+wt26*fx26+wt27*fx27);
		fy=(wt1*fy1+wt2*fy2+wt3*fy3+wt4*fy4  +wt5 *fy5 +wt6 *fy6+
				wt7 *fy7 +wt8 *fy8 +wt9 *fy9 +wt10*fy10+wt11*fy11+wt12*fy12+
				wt13*fy13+wt14*fy14+wt15*fy15+wt16*fy16+wt17*fy17+wt18*fy18+
				wt19*fy19+wt20*fy20+wt21*fy21+wt22*fy22+wt23*fy23+wt24*fy24+
				wt25*fy25+wt26*fy26+wt27*fy27);
		fz=(wt1*fz1+wt2*fz2+wt3*fz3+wt4*fz4  +wt5 *fz5 +wt6 *fz6+
				wt7 *fz7 +wt8 *fz8 +wt9 *fz9 +wt10*fz10+wt11*fz11+wt12*fz12+
				wt13*fz13+wt14*fz14+wt15*fz15+wt16*fz16+wt17*fz17+wt18*fz18+
				wt19*fz19+wt20*fz20+wt21*fz21+wt22*fz22+wt23*fz23+wt24*fz24+
				wt25*fz25+wt26*fz26+wt27*fz27);
//		/* The correction factor is moved before correlation estimation. - kjhan 
		pmp[i].vx += fact2*fx*simpar.pcorr;
		pmp[i].vy += fact2*fy*simpar.pcorr;
		pmp[i].vz += fact2*fz*simpar.pcorr;
//		*/ minseong_2024_6_15 pmmain also changed
//		pmp[i].vx += fact2*fx;
//		pmp[i].vy += fact2*fy;
//		pmp[i].vz += fact2*fz;
	}
#endif
	Free(nzperiodic); Free(nyperiodic); Free(nxperiodic);
	TIMER_STOP(55);
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid == 0) {
		fprintf(stdout,"CPU(FDA4) %g\n",ELAPSED_TIME(55)); fflush(stdout);
	}
}
void onestepforwardposition(pmparticletype *pmp,int np){
	float rngx,rngy,rngz,tdmax=0;
	int myid;
	long long i,j,k;
	float dmax;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	rngx = simpar.nx;
	rngy = simpar.ny;
	rngz = simpar.nz;
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		float ddmax=-9.e29;
#ifdef _OPENMP
#pragma omp  for  private(i)
#endif
		for(i=0;i<np;i++){
			float xd,yd,zd,rd;
			xd = pmp[i].vx*simpar.pfact;
			yd = pmp[i].vy*simpar.pfact;
			zd = pmp[i].vz*simpar.pfact;
			rd = xd*xd + yd*yd + zd*zd;
			ddmax = MAX(ddmax,rd);
#ifdef XYZDBL
			pmp[i].x += xd;
			pmp[i].y += yd;
			pmp[i].z += zd;
#else
			pmp[i].x = fmod(pmp[i].x+xd+rngx,rngx);
			pmp[i].y = fmod(pmp[i].y+yd+rngy,rngy);
			pmp[i].z = fmod(pmp[i].z+zd+rngz,rngz);
#endif
		}
#pragma omp critical
		{
			if(ddmax > tdmax) tdmax = ddmax;
		}
	}
	dmax = sqrt(tdmax);
	MPI_Reduce(&dmax,&tdmax,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
	if(myid ==0) printf("Maximum displacement = %d %g\n ",simpar.stepcount,tdmax);
}
#undef den
#undef PTR
