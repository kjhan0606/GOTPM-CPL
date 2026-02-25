#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include<omp.h>

#include "pmheader.h"
#include "Memory.h"

/*
#define icopysign(a) (signbit(a)? -1:1)
*/
#define icopysign(a) copysignf(1.,a)

#define denvel(i,j,k) (denvel[(i)+mx*((j)+ny*(k))])
#ifdef XYZDBL
#define FABS(a) fabs(a)
#define RINT(a) rint(a)
#else
#define FABS(a) fabsf(a)
#define RINT(a) rintf(a)
#endif

static int nid,myid;
#define NZWIDTH 4


void veltscX(float *denvel, 
		int local_nz, int nstart_z, pmparticletype *pmp, int np, float pmas){
	float  start_z;
	float pmas0,p05;
	int xsign,ysign,zsign,nz_per;
	long long i,j,k,ii;
	int nearx,neary,nearz,i1,j1,k1,i2,j2,k2,i3,j3,k3;
#ifdef XYZDBL
	double xp,yp,zp;
#else
	float xp,yp,zp;
#endif
	float xmin,ymin,zmin;

	float wx1,wy1,wz1,wx2,wy2,wz2,wx3,wy3,wz3;
	float xd1,yd1,zd1;
	double rngc;
	int usndsize,dgetsize,dsndsize,ugetsize,dest,src;
	MPI_Status status;
	int ierror,stag,rtag,tag;
	int local_z_start,local_ny_after_transpose,local_y_start_after_transpose,
		total_local_size;
	long long nx,ny,nz,mx,mgrid,ngrid,slicesize;
	int *nxperiodic,*nyperiodic,*nzperiodic;
	int *nxp,*nyp,*nzp,nspace;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	TIMER_START(59);
	nx = simpar.nx; ny = simpar.ny; nz = simpar.nz;
	nspace = simpar.nspace;
	start_z = nstart_z;
	rngc = (double) nx *(double)ny *(double)nz;
	nz_per = local_nz+NZWIDTH; /* local cyclic width in z direction */


	nxperiodic = (int *)Malloc(sizeof(int)*(nx+10),PPTR(nxperiodic)); nxp = nxperiodic+5;
	nyperiodic = (int *)Malloc(sizeof(int)*(ny+10),PPTR(nyperiodic)); nyp = nyperiodic+5;
	nzperiodic = (int *)Malloc(sizeof(int)*(nz_per+10),PPTR(nzperiodic)); nzp = nzperiodic+5;
	for(i=-5;i<nx+5;i++) nxp[i] = (i+nx)%nx;
	for(i=-5;i<ny+5;i++) nyp[i] = (i+ny)%ny;
	for(i=-5;i<nz_per+5;i++) nzp[i] = (i+nz_per)%nz_per;


	mx = 2*(nx/2+1);
	slicesize = mx * ny;
	mgrid = slicesize * nz_per;
	ngrid = slicesize * local_nz;

#ifdef VarPM
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<mgrid;i++) denvel[i] = 0;
#else
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<ngrid;i++) denvel[i] = 0.;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=ngrid;i<mgrid;i++) denvel[i] = 0;
#endif

#ifndef MULTIMASS
	pmas0 = pmas;
	p05 = pmas0 * 0.5;
#endif
	for(i=0;i<np;i++){
		float wx1wy1,wx2wy1,wx3wy1;
		float wx1wy2,wx2wy2,wx3wy2;
		float wx1wy3,wx2wy3,wx3wy3;
#if defined(MULTIMASS)
		pmas0= pmp[i].mass * rngc;
		p05 = pmas0 * 0.5;
#endif
		pmas0 = pmas*pmp[i].vx;
		p05 = pmas0*0.5;

		xp = XofP(pmp+i);
		yp = YofP(pmp+i);
		zp = ZofP(pmp+i) - start_z;
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


		wx1 = (0.75-xd1*xd1)*pmas0;
		wy1 =  0.75-yd1*yd1;
		wz1 =  0.75-zd1*zd1;

		wx3 = p05*(0.25+xd1*(xd1-1.));
		wx2 = wx3 + pmas0*xd1;

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

		denvel(i1,j1,k1) += wx1wy1*wz1;
		denvel(i2,j1,k1) += wx2wy1*wz1;
		denvel(i3,j1,k1) += wx3wy1*wz1;
		denvel(i1,j2,k1) += wx1wy2*wz1;
		denvel(i2,j2,k1) += wx2wy2*wz1;
		denvel(i3,j2,k1) += wx3wy2*wz1;
		denvel(i1,j3,k1) += wx1wy3*wz1;
		denvel(i2,j3,k1) += wx2wy3*wz1;
		denvel(i3,j3,k1) += wx3wy3*wz1;

		denvel(i1,j1,k2) += wx1wy1*wz2;
		denvel(i2,j1,k2) += wx2wy1*wz2;
		denvel(i3,j1,k2) += wx3wy1*wz2;
		denvel(i1,j2,k2) += wx1wy2*wz2;
		denvel(i2,j2,k2) += wx2wy2*wz2;
		denvel(i3,j2,k2) += wx3wy2*wz2;
		denvel(i1,j3,k2) += wx1wy3*wz2;
		denvel(i2,j3,k2) += wx2wy3*wz2;
		denvel(i3,j3,k2) += wx3wy3*wz2;

		denvel(i1,j1,k3) += wx1wy1*wz3;
		denvel(i2,j1,k3) += wx2wy1*wz3;
		denvel(i3,j1,k3) += wx3wy1*wz3;
		denvel(i1,j2,k3) += wx1wy2*wz3;
		denvel(i2,j2,k3) += wx2wy2*wz3;
		denvel(i3,j2,k3) += wx3wy2*wz3;
		denvel(i1,j3,k3) += wx1wy3*wz3;
		denvel(i2,j3,k3) += wx2wy3*wz3;
		denvel(i3,j3,k3) += wx3wy3*wz3;
	}

#ifdef VarPM
#ifdef INTEL /* intel compiler may sometimes experience problems 
				in pushing memory copy backward.
			This is an optimization issue. */
	i=slicesize*(local_nz+2);
	while(i >=0){
		denvel[i+slicesize] = denvel[i]; i--;
	}
	i=0;
	while(i<slicesize){
		denvel[i] = denvel[i + (slicesize*(local_nz+NZWIDTH-1))]; i++;
	}
#else
	for(i=slicesize*(local_nz+2);i>=0;i--){
		denvel[i+slicesize] = denvel[i];
	}
	for(i=0;i<slicesize;i++){
		denvel[i] = denvel[i + (slicesize*(local_nz+NZWIDTH-1))];
	}
#endif
#else
	dest = (myid+1+nid)%nid;
	src = (myid-1+nid)%nid;
	stag = rtag = 0;
	usndsize = slicesize*2;
	dgetsize = usndsize;
	MPI_Sendrecv_replace((denvel+slicesize*local_nz),usndsize,
			MPI_REAL,dest,stag,src,rtag,MPI_COMM_WORLD,&status);
	for(i=0;i<slicesize*2;i++){
		denvel[i] += denvel[i+slicesize*local_nz];
	}
	ugetsize = (dsndsize = slicesize);
	dest = (myid-1+nid)%nid;
	src = (myid+1+nid)%nid;
	stag = rtag = 0;
	MPI_Sendrecv_replace((denvel+slicesize*(local_nz+NZWIDTH-1)),dsndsize,
			MPI_REAL,dest,stag,src,rtag,MPI_COMM_WORLD,&status);
	for(i=0;i<slicesize;i++){
		denvel[i+slicesize*(local_nz-1)] += denvel[i+slicesize*(local_nz+NZWIDTH-1)];
	}
#endif
	Free(nzperiodic); Free(nyperiodic); Free(nxperiodic);
	TIMER_STOP(59);
	if(myid==0) fprintf(stdout,"CPU(VELX TSC)   = %f\n",ELAPSED_TIME(59));
}

void veltscY(float *denvel, 
		int local_nz, int nstart_z, pmparticletype *pmp, int np, float pmas){
	float  start_z;
	float pmas0,p05;
	int xsign,ysign,zsign,nz_per;
	long long i,j,k,ii;
	int nearx,neary,nearz,i1,j1,k1,i2,j2,k2,i3,j3,k3;
#ifdef XYZDBL
	double xp,yp,zp;
#else
	float xp,yp,zp;
#endif
	float xmin,ymin,zmin;

	float wx1,wy1,wz1,wx2,wy2,wz2,wx3,wy3,wz3;
	float xd1,yd1,zd1;
	double rngc;
	int usndsize,dgetsize,dsndsize,ugetsize,dest,src;
	MPI_Status status;
	int ierror,stag,rtag,tag;
	int local_z_start,local_ny_after_transpose,local_y_start_after_transpose,
		total_local_size;
	long long nx,ny,nz,mx,mgrid,ngrid,slicesize;
	int *nxperiodic,*nyperiodic,*nzperiodic;
	int *nxp,*nyp,*nzp,nspace;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	TIMER_START(59);
	nx = simpar.nx; ny = simpar.ny; nz = simpar.nz;
	nspace = simpar.nspace;
	start_z = nstart_z;
	rngc = (double) nx *(double)ny *(double)nz;
	nz_per = local_nz+NZWIDTH; /* local cyclic width in z direction */


	nxperiodic = (int *)Malloc(sizeof(int)*(nx+10),PPTR(nxperiodic)); nxp = nxperiodic+5;
	nyperiodic = (int *)Malloc(sizeof(int)*(ny+10),PPTR(nyperiodic)); nyp = nyperiodic+5;
	nzperiodic = (int *)Malloc(sizeof(int)*(nz_per+10),PPTR(nzperiodic)); nzp = nzperiodic+5;
	for(i=-5;i<nx+5;i++) nxp[i] = (i+nx)%nx;
	for(i=-5;i<ny+5;i++) nyp[i] = (i+ny)%ny;
	for(i=-5;i<nz_per+5;i++) nzp[i] = (i+nz_per)%nz_per;


	mx = 2*(nx/2+1);
	slicesize = mx * ny;
	mgrid = slicesize * nz_per;
	ngrid = slicesize * local_nz;

#ifdef VarPM
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<mgrid;i++) denvel[i] = 0;
#else
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<ngrid;i++) denvel[i] = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=ngrid;i<mgrid;i++) denvel[i] = 0;
#endif
#ifndef MULTIMASS
	pmas0 = pmas;
	p05 = pmas0 * 0.5;
#endif
	for(i=0;i<np;i++){
		float wx1wy1,wx2wy1,wx3wy1;
		float wx1wy2,wx2wy2,wx3wy2;
		float wx1wy3,wx2wy3,wx3wy3;
#if defined(MULTIMASS)
		pmas0= pmp[i].mass * rngc;
		p05 = pmas0 * 0.5;
#endif
		pmas0 = pmas*pmp[i].vy;
		p05 = pmas0*0.5;

		xp = XofP(pmp+i);
		yp = YofP(pmp+i);
		zp = ZofP(pmp+i) - start_z;
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


		wx1 = (0.75-xd1*xd1)*pmas0;
		wy1 =  0.75-yd1*yd1;
		wz1 =  0.75-zd1*zd1;

		wx3 = p05*(0.25+xd1*(xd1-1.));
		wx2 = wx3 + pmas0*xd1;

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

		denvel(i1,j1,k1) += wx1wy1*wz1;
		denvel(i2,j1,k1) += wx2wy1*wz1;
		denvel(i3,j1,k1) += wx3wy1*wz1;
		denvel(i1,j2,k1) += wx1wy2*wz1;
		denvel(i2,j2,k1) += wx2wy2*wz1;
		denvel(i3,j2,k1) += wx3wy2*wz1;
		denvel(i1,j3,k1) += wx1wy3*wz1;
		denvel(i2,j3,k1) += wx2wy3*wz1;
		denvel(i3,j3,k1) += wx3wy3*wz1;

		denvel(i1,j1,k2) += wx1wy1*wz2;
		denvel(i2,j1,k2) += wx2wy1*wz2;
		denvel(i3,j1,k2) += wx3wy1*wz2;
		denvel(i1,j2,k2) += wx1wy2*wz2;
		denvel(i2,j2,k2) += wx2wy2*wz2;
		denvel(i3,j2,k2) += wx3wy2*wz2;
		denvel(i1,j3,k2) += wx1wy3*wz2;
		denvel(i2,j3,k2) += wx2wy3*wz2;
		denvel(i3,j3,k2) += wx3wy3*wz2;

		denvel(i1,j1,k3) += wx1wy1*wz3;
		denvel(i2,j1,k3) += wx2wy1*wz3;
		denvel(i3,j1,k3) += wx3wy1*wz3;
		denvel(i1,j2,k3) += wx1wy2*wz3;
		denvel(i2,j2,k3) += wx2wy2*wz3;
		denvel(i3,j2,k3) += wx3wy2*wz3;
		denvel(i1,j3,k3) += wx1wy3*wz3;
		denvel(i2,j3,k3) += wx2wy3*wz3;
		denvel(i3,j3,k3) += wx3wy3*wz3;
	}
#ifdef VarPM
#ifdef INTEL /* intel compiler may experience problems in pushing memory copy toward address. 
			This is a optimization issue. */
	i=slicesize*(local_nz+2);
	while(i >=0){
		denvel[i+slicesize] = denvel[i]; i--;
	}
	i=0;
	while(i<slicesize){
		denvel[i] = denvel[i + (slicesize*(local_nz+NZWIDTH-1))]; i++;
	}
#else
	for(i=slicesize*(local_nz+2);i>=0;i--){
		denvel[i+slicesize] = denvel[i];
	}
	for(i=0;i<slicesize;i++){
		denvel[i] = denvel[i + (slicesize*(local_nz+NZWIDTH-1))];
	}
#endif
#else
	dest = (myid+1+nid)%nid;
	src = (myid-1+nid)%nid;
	stag = rtag = 0;
	usndsize = slicesize*2;
	dgetsize = usndsize;
	MPI_Sendrecv_replace((denvel+slicesize*local_nz),usndsize,
			MPI_REAL,dest,stag,src,rtag,MPI_COMM_WORLD,&status);
	for(i=0;i<slicesize*2;i++){
		denvel[i] += denvel[i+slicesize*local_nz];
	}
	ugetsize = (dsndsize = slicesize);
	dest = (myid-1+nid)%nid;
	src = (myid+1+nid)%nid;
	stag = rtag = 0;
	MPI_Sendrecv_replace((denvel+slicesize*(local_nz+NZWIDTH-1)),dsndsize,
			MPI_REAL,dest,stag,src,rtag,MPI_COMM_WORLD,&status);
	for(i=0;i<slicesize;i++){
		denvel[i+slicesize*(local_nz-1)] += denvel[i+slicesize*(local_nz+NZWIDTH-1)];
	}
#endif
	Free(nzperiodic); Free(nyperiodic); Free(nxperiodic);
	TIMER_STOP(59);
	if(myid==0) fprintf(stdout,"CPU(VELY TSC)   = %f\n",ELAPSED_TIME(59));
}

void veltscZ(float *denvel, 
		int local_nz, int nstart_z, pmparticletype *pmp, int np, float pmas){
	float  start_z;
	float pmas0,p05;
	int xsign,ysign,zsign,nz_per;
	long long i,j,k,ii;
	int nearx,neary,nearz,i1,j1,k1,i2,j2,k2,i3,j3,k3;
#ifdef XYZDBL
	double xp,yp,zp;
#else
	float xp,yp,zp;
#endif
	float xmin,ymin,zmin;

	float wx1,wy1,wz1,wx2,wy2,wz2,wx3,wy3,wz3;
	float xd1,yd1,zd1;
	double rngc;
	int usndsize,dgetsize,dsndsize,ugetsize,dest,src;
	MPI_Status status;
	int ierror,stag,rtag,tag;
	int local_z_start,local_ny_after_transpose,local_y_start_after_transpose,
		total_local_size;
	long long nx,ny,nz,mx,mgrid,ngrid,slicesize;
	int *nxperiodic,*nyperiodic,*nzperiodic;
	int *nxp,*nyp,*nzp,nspace;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	TIMER_START(59);
	nx = simpar.nx; ny = simpar.ny; nz = simpar.nz;
	nspace = simpar.nspace;
	start_z = nstart_z;
	rngc = (double) nx *(double)ny *(double)nz;
	nz_per = local_nz+NZWIDTH; /* local cyclic width in z direction */


	nxperiodic = (int *)Malloc(sizeof(int)*(nx+10),PPTR(nxperiodic)); nxp = nxperiodic+5;
	nyperiodic = (int *)Malloc(sizeof(int)*(ny+10),PPTR(nyperiodic)); nyp = nyperiodic+5;
	nzperiodic = (int *)Malloc(sizeof(int)*(nz_per+10),PPTR(nzperiodic)); nzp = nzperiodic+5;
	for(i=-5;i<nx+5;i++) nxp[i] = (i+nx)%nx;
	for(i=-5;i<ny+5;i++) nyp[i] = (i+ny)%ny;
	for(i=-5;i<nz_per+5;i++) nzp[i] = (i+nz_per)%nz_per;


	mx = 2*(nx/2+1);
	slicesize = mx * ny;
	mgrid = slicesize * nz_per;
	ngrid = slicesize * local_nz;

#ifdef VarPM
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<mgrid;i++) denvel[i] = 0;
#else
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<ngrid;i++) denvel[i] = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=ngrid;i<mgrid;i++) denvel[i] = 0;
#endif
#ifndef MULTIMASS
	pmas0 = pmas;
	p05 = pmas0 * 0.5;
#endif
	for(i=0;i<np;i++){
		float wx1wy1,wx2wy1,wx3wy1;
		float wx1wy2,wx2wy2,wx3wy2;
		float wx1wy3,wx2wy3,wx3wy3;
#if defined(MULTIMASS)
		pmas0= pmp[i].mass * rngc;
		p05 = pmas0 * 0.5;
#endif
		pmas0 = pmas*pmp[i].vz;
		p05 = pmas0*0.5;

		xp = XofP(pmp+i);
		yp = YofP(pmp+i);
		zp = ZofP(pmp+i) - start_z;
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


		wx1 = (0.75-xd1*xd1)*pmas0;
		wy1 =  0.75-yd1*yd1;
		wz1 =  0.75-zd1*zd1;

		wx3 = p05*(0.25+xd1*(xd1-1.));
		wx2 = wx3 + pmas0*xd1;

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

		denvel(i1,j1,k1) += wx1wy1*wz1;
		denvel(i2,j1,k1) += wx2wy1*wz1;
		denvel(i3,j1,k1) += wx3wy1*wz1;
		denvel(i1,j2,k1) += wx1wy2*wz1;
		denvel(i2,j2,k1) += wx2wy2*wz1;
		denvel(i3,j2,k1) += wx3wy2*wz1;
		denvel(i1,j3,k1) += wx1wy3*wz1;
		denvel(i2,j3,k1) += wx2wy3*wz1;
		denvel(i3,j3,k1) += wx3wy3*wz1;

		denvel(i1,j1,k2) += wx1wy1*wz2;
		denvel(i2,j1,k2) += wx2wy1*wz2;
		denvel(i3,j1,k2) += wx3wy1*wz2;
		denvel(i1,j2,k2) += wx1wy2*wz2;
		denvel(i2,j2,k2) += wx2wy2*wz2;
		denvel(i3,j2,k2) += wx3wy2*wz2;
		denvel(i1,j3,k2) += wx1wy3*wz2;
		denvel(i2,j3,k2) += wx2wy3*wz2;
		denvel(i3,j3,k2) += wx3wy3*wz2;

		denvel(i1,j1,k3) += wx1wy1*wz3;
		denvel(i2,j1,k3) += wx2wy1*wz3;
		denvel(i3,j1,k3) += wx3wy1*wz3;
		denvel(i1,j2,k3) += wx1wy2*wz3;
		denvel(i2,j2,k3) += wx2wy2*wz3;
		denvel(i3,j2,k3) += wx3wy2*wz3;
		denvel(i1,j3,k3) += wx1wy3*wz3;
		denvel(i2,j3,k3) += wx2wy3*wz3;
		denvel(i3,j3,k3) += wx3wy3*wz3;
	}
#ifdef VarPM
#ifdef INTEL /* intel compiler may experience problems in pushing memory copy toward address. 
			This is a optimization issue. */
	i=slicesize*(local_nz+2);
	while(i >=0){
		denvel[i+slicesize] = denvel[i]; i--;
	}
	i=0;
	while(i<slicesize){
		denvel[i] = denvel[i + (slicesize*(local_nz+NZWIDTH-1))]; i++;
	}
#else
	for(i=slicesize*(local_nz+2);i>=0;i--){
		denvel[i+slicesize] = denvel[i];
	}
	for(i=0;i<slicesize;i++){
		denvel[i] = denvel[i + (slicesize*(local_nz+NZWIDTH-1))];
	}
#endif
#else
	dest = (myid+1+nid)%nid;
	src = (myid-1+nid)%nid;
	stag = rtag = 0;
	usndsize = slicesize*2;
	dgetsize = usndsize;
	MPI_Sendrecv_replace((denvel+slicesize*local_nz),usndsize,
			MPI_REAL,dest,stag,src,rtag,MPI_COMM_WORLD,&status);
	for(i=0;i<slicesize*2;i++){
		denvel[i] += denvel[i+slicesize*local_nz];
	}
	ugetsize = (dsndsize = slicesize);
	dest = (myid-1+nid)%nid;
	src = (myid+1+nid)%nid;
	stag = rtag = 0;
	MPI_Sendrecv_replace((denvel+slicesize*(local_nz+NZWIDTH-1)),dsndsize,
			MPI_REAL,dest,stag,src,rtag,MPI_COMM_WORLD,&status);
	for(i=0;i<slicesize;i++){
		denvel[i+slicesize*(local_nz-1)] += denvel[i+slicesize*(local_nz+NZWIDTH-1)];
	}
#endif
	Free(nzperiodic); Free(nyperiodic); Free(nxperiodic);
	TIMER_STOP(59);
	if(myid==0) fprintf(stdout,"CPU(VELZ TSC)   = %f\n",ELAPSED_TIME(59));
}

