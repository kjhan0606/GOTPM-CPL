#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stddef.h>
#include<mpi.h>
#include<omp.h>
#include "pmheader.h"
#include "Memory.h"
#include "varpm.h"

void diffphi_(float *,float *,int *);
#define max(a,b) (a)>(b)? (a):(b)
#define ffx1(i,j,k) ffx1[(i)+mx*((j)+ny*(k))]
#define ffx2(i,j,k) ffx2[(i)+mx*((j)+ny*(k))]
#define f1(i,j,k) f1[(i)+mx*((j)+ny*(k))]
#define f2(i,j,k) f2[(i)+mx*((j)+ny*(k))]
void TwoLPT_getinnerparticlesforce(long *mp,pmvrparticletype *pvr,
				int nzparticlewidth, int nzprefixel, float zstart, float zheight,
				float *phi1, float *phi2, float vamp1, float vamp2){
	long mx;
	int nspace;
	long nx,ny,nz;
	long np,i,j,k;
	long ii,jj,kk;
	int ixyz;
	long long densize, fftwsize;
	float zfinal,vfact;

	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;
	mx = 2*(nx/2+1);
	nspace = simpar.nspace;

	densize = 2*(simpar.nx/2+1) * (long)simpar.ny * (long)simpar.local_nz;
	if(densize ==0) densize = 10;
	fftwsize = 2*(simpar.nx/2+1)*(long)simpar.ny*(long)(ceil(zstart+zheight)-(int)zstart+1);

	zfinal = zstart + zheight;

	vfact = vamp1;

	for(ixyz=1;ixyz<=3;ixyz++){
		float *f1, *ffx1;
		f1 = (float *)Malloc(sizeof(float)*densize,PPTR(f1));
		ffx1 = (float*)Malloc(sizeof(float)*fftwsize, PPTR(ffx1));
		diffphi_(phi1,f1,&ixyz);
		adVarPM_fda4_mesh(f1,simpar.local_nz,simpar.local_z_start,ffx1,zstart,zfinal,simpar.nx,simpar.ny,simpar.nz);
		Free(f1);
		np = 0;
		for(k=1;k<(nzparticlewidth)-1;k++){
			kk = nspace*k + (nzprefixel);
			for(j=0;j<ny/nspace;j++){
				jj = nspace*j;
				for(i=0;i<nx/nspace;i++){
					ii = nspace*i;
					if(ixyz==1)      pvr[np].vx = vfact*ffx1(ii,jj,kk);
					else if(ixyz==2) pvr[np].vy = vfact*ffx1(ii,jj,kk);
					else             pvr[np].vz = vfact*ffx1(ii,jj,kk);
					np++;
				}
			}
		}
		Free(ffx1);
	}
	vfact = vamp2;
	for(ixyz=1;ixyz<=3;ixyz++){
		float *f2, *ffx2;
		f2 = (float *)Malloc(sizeof(float)*densize,PPTR(f2));
		ffx2 = (float*)Malloc(sizeof(float)*fftwsize, PPTR(ffx2));
		diffphi_(phi2,f2,&ixyz);
		adVarPM_fda4_mesh(f2,simpar.local_nz,simpar.local_z_start,ffx2,zstart,zfinal,simpar.nx,simpar.ny,simpar.nz);
		Free(f2);
		np = 0;
		for(k=1;k<(nzparticlewidth)-1;k++){
			kk = nspace*k + (nzprefixel);
			for(j=0;j<ny/nspace;j++){
				jj = nspace*j;
				for(i=0;i<nx/nspace;i++){
					ii = nspace*i;
					if(ixyz==1)      pvr[np].vx +=                      vfact *ffx2(ii,jj,kk);
					else if(ixyz==2) pvr[np].vy +=                      vfact *ffx2(ii,jj,kk);
					else             pvr[np].vz +=                      vfact *ffx2(ii,jj,kk);
					np++;
				}
			}
		}
		Free(ffx2);
	}


	MPI_Barrier(MPI_COMM_WORLD);
	*mp = np;
}
void TwoLPT_setinnerparticlesforce(long *mp,pmparticletype *particles,
				     int nzparticlewidth,int nzprefixel,float zstartpos, float *ken,
					 float zstart,float zheight, float omei, float *phi1, float *phi2, 
					 float damp1, float damp2){
	long i,j,k,ii,jj,kk;
	float zfinal;
	int myid,nid;
	float rnx,rny,rnz;
	float pfact,rng;
	float ai = 1;
	int nx,ny,nz;
	long np,tnp;
	double delx,displ,mdispl,tdelx,tmdispl;
	long long densize,fftwsize;
	int nspace = simpar.nspace;
	long mx;

	zfinal = zstart + zheight;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	rnx = nx = simpar.nx;
	rny = ny = simpar.ny;
	rnz = nz = simpar.nz;
	rng = rnx;
	mx = 2*(nx/2+1);

	np = (*mp);
	pfact = rng*ai*damp1;
	for(k=1;k<(nzparticlewidth)-1;k++){
		for(j=0;j<ny/(nspace);j++){
			for(i=0;i<nx/(nspace);i++){
#ifdef XYZDBL
				particles[np].x = 0;
				particles[np].y = 0;
				particles[np].z = 0;
#else
				particles[np].x = (nspace)*i;
				particles[np].y = (nspace)*j;
				particles[np].z = (nspace)*k+(zstartpos);
#endif
				np ++;
			}
		}
	}
#ifdef INDEX
	void indexingip_(pmparticletype *,long *,long *,int *,int *, int *,int *,int *);
	indexingip_(particles,&np,mp,&nzparticlewidth,&nspace,&nx,&ny,&nz);
#endif
	densize = 2*(simpar.nx/2+1) * (long) simpar.ny * (long)simpar.local_nz;
	if(densize==0) densize = 10;
	fftwsize = 2*(simpar.nx/2+1)*(long)simpar.ny*(long)(ceil(zstart+zheight)-(int)zstart+1);

	int ixyz;
	for(ixyz=1;ixyz<=3;ixyz++){
		float *f1,*ffx1;
		f1 = (float *)Malloc(sizeof(float)*densize,PPTR(f1));
		ffx1 = (float*)Malloc(sizeof(float)*fftwsize, PPTR(ffx1));
		diffphi_(phi1,f1,&ixyz);
		adVarPM_fda4_mesh(f1,simpar.local_nz,simpar.local_z_start,ffx1,zstart,zfinal,simpar.nx,simpar.ny,simpar.nz);
		Free(f1);

		np = *mp;
		for(k=1;k<(nzparticlewidth)-1;k++){
			kk = (nspace)*k + (nzprefixel);
			for(j=0;j<ny/(nspace);j++){
				jj = (nspace)*j;
				for(i=0;i<nx/(nspace);i++){
					ii = (nspace)*i;
					if(ixyz==1)      particles[np].x += pfact*ffx1(ii,jj,kk);
					else if(ixyz==2) particles[np].y += pfact*ffx1(ii,jj,kk);
					else             particles[np].z += pfact*ffx1(ii,jj,kk);
					np ++;
				}
			}
		}
		Free(ffx1);
	}
	pfact = rng*ai*damp2;
	for(ixyz=1;ixyz<=3;ixyz++){
		float *f2,*ffx2;
		f2 = (float *)Malloc(sizeof(float)*densize,PPTR(f2));
		ffx2 = (float*)Malloc(sizeof(float)*fftwsize, PPTR(ffx2));
		diffphi_(phi2,f2,&ixyz);
		adVarPM_fda4_mesh(f2,simpar.local_nz,simpar.local_z_start,ffx2,zstart,zfinal,simpar.nx,simpar.ny,simpar.nz);
		Free(f2);

		np = *mp;
		for(k=1;k<(nzparticlewidth)-1;k++){
			kk = (nspace)*k + (nzprefixel);
			for(j=0;j<ny/(nspace);j++){
				jj = (nspace)*j;
				for(i=0;i<nx/(nspace);i++){
					ii = (nspace)*i;
					if(ixyz==1){
						particles[np].x += pfact*ffx2(ii,jj,kk);
#ifndef XYZDBL
						particles[np].x = fmod(particles[np].x+rnx,rnx);
#endif
					}
					else if(ixyz==2){
						particles[np].y += pfact*ffx2(ii,jj,kk);
#ifndef XYZDBL
						particles[np].y = fmod(particles[np].y+rny,rny);
#endif
					}
					else{
						particles[np].z += pfact*ffx2(ii,jj,kk);
#ifndef XYZDBL
						particles[np].z = fmod(particles[np].z+rnz,rnz);
#endif
					}
					np ++;
				}
			}
		}
		Free(ffx2);
	}


	pfact = rng*ai*damp1;
	delx = displ = mdispl = 0;
	for(i=0;i<np;i++){
		double vamp;
		vamp  = particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy 
				+particles[i].vz*particles[i].vz;
		delx+= vamp;
		displ = sqrt(vamp)*pfact;
		mdispl = max(displ,mdispl);
	}
	MPI_Reduce(&delx,&tdelx,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	if(myid==0) delx = tdelx;
	MPI_Bcast(&delx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Reduce(&np,&tnp,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&tnp,1,MPI_LONG,0,MPI_COMM_WORLD);
	MPI_Reduce(&mdispl,&tmdispl,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	if(myid==0) {
		double pi=3.141592653L,constant,poten,ktot;
		poten = ktot = 0;
		printf("P%d %g %g %ld\n",myid,delx,pfact,tnp);
		printf("rms displacement = %g\n",sqrt(delx)*pfact/sqrt(tnp)),
		printf("Maximum displacement = %g\n",tmdispl);
		printf("sum of dx/da**2=%g\n",tdelx);
		tdelx = 0.5*8.*pi/3.L*(1.L/ai+1.L/(omei)-1.L)*tdelx;
		constant = ai*ai*ai * tdelx + poten+ktot;
		printf("poten,ken,const = %g %g %g\n",poten,tdelx,constant);
	}
	*ken = tdelx;
	*mp  = np;
}

void TwoLPT_setboundaryparticlesforce(long *mp,pmparticletype *pmparticles,
				int nzparticlewidth, float zstartpos, int nzprefixel, float zstart, float zwidth,
				float *phi1,float *phi2,float damp1, float damp2, float vamp1, float vamp2){
	float ai;
	int myid,nid;
	float rng,pfact;
	float rnx,rny,rnz;
	long np,i,j,k,mx;
	long ii,jj,kk;
	long jidx,jkidx,kidx;
	int nzstride;
	int nx,ny,nz;
	long long densize,fftwsize;
	float zheight, zfinal;
	int nspace = simpar.nspace;
	ai = 1;



	/*
	printf("P%d has  damp1/damp2 = %g %g  : damp2/vamp2 %g %g\n",simpar.myid,damp1,damp2,vamp1,vamp2);
	*/

	zheight = zwidth;
	zfinal = zstart + zwidth;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	rng = simpar.nx;
	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;
	rnx = nx; rny =ny; rnz = nz;
	mx = 2*(nx/2+1);
	nzstride  = max((nzparticlewidth)-1,1);
	densize = 2*(simpar.nx/2+1) * (long)simpar.ny * (long) simpar.local_nz;
	if(densize ==0) densize = 10;
	fftwsize = 2*(simpar.nx/2+1)*(long)simpar.ny*(long)(ceil(zstart+zheight)-(int)zstart+1);
	int ixyz;
	float vfact;
	pfact = ai*simpar.nx*damp1;
	vfact = vamp1;

	for(ixyz=1;ixyz<=3;ixyz++){
		float *f1, *ffx1;
		f1 = (float *)Malloc(sizeof(float)*densize,PPTR(f1));
		ffx1 = (float*)Malloc(sizeof(float)*fftwsize, PPTR(ffx1));
		diffphi_(phi1,f1,&ixyz);
		adVarPM_fda4_mesh(f1,simpar.local_nz,simpar.local_z_start,ffx1,zstart,zfinal,simpar.nx,simpar.ny,simpar.nz);
		if(simpar.myid==0){
			printf("P%d is now after differencing the phi1 and phi2 in the %d direction\n",simpar.myid,ixyz);
			printf("P%d has %d %d %d %d : vamp1/damp1 = %g %g\n",simpar.myid,nzparticlewidth,nzstride,
					nspace,nzprefixel,vamp1,damp1);
		}
		Free(f1);
		np = 0;
		for(k=0;k<(nzparticlewidth);k+=nzstride){
			kk = (nspace)*k + (nzprefixel);
			for(j=0;j<ny/(nspace);j++){
				jj = (nspace)*j;
				for(i=0;i<nx/(nspace);i++){
					ii = (nspace)*i;
					if(ixyz ==1 ) {
#ifdef XYZDBL
						pmparticles[np].x = 0;
#else
						pmparticles[np].x = (nspace)*i;
#endif
						pmparticles[np].vx = vfact*ffx1(ii,jj,kk);
						pmparticles[np].x += pfact*ffx1(ii,jj,kk);
					}
					else if(ixyz ==2){
#ifdef XYZDBL
						pmparticles[np].y = 0;
#else
						pmparticles[np].y = (nspace)*j;
#endif
						pmparticles[np].vy = vfact*ffx1(ii,jj,kk);
						pmparticles[np].y += pfact*ffx1(ii,jj,kk);
					}
					else {
#ifdef XYZDBL
						pmparticles[np].z = 0;
#else
						pmparticles[np].z = (nspace)*k+(zstartpos);
#endif
						pmparticles[np].vz = vfact*ffx1(ii,jj,kk);
						pmparticles[np].z += pfact*ffx1(ii,jj,kk);
					}
					np ++;
				}
			}
		}
		Free(ffx1);
	}
	pfact = ai*simpar.nx*damp2;
	vfact = vamp2;
	for(ixyz=1;ixyz<=3;ixyz++){
		float *f2,*ffx2;
		f2 = (float *)Malloc(sizeof(float)*densize,PPTR(f2));
		ffx2 = (float*)Malloc(sizeof(float)*fftwsize, PPTR(ffx2));
		diffphi_(phi2,f2,&ixyz);
		adVarPM_fda4_mesh(f2,simpar.local_nz,simpar.local_z_start,ffx2,zstart,zfinal,simpar.nx,simpar.ny,simpar.nz);
		Free(f2);
		np = 0;
		for(k=0;k<(nzparticlewidth);k+=nzstride){
			kk = (nspace)*k + (nzprefixel);
			for(j=0;j<ny/(nspace);j++){
				jj = (nspace)*j;
				for(i=0;i<nx/(nspace);i++){
					ii = (nspace)*i;
					if(ixyz ==1 ) {
						pmparticles[np].vx +=                      vfact*ffx2(ii,jj,kk);
						pmparticles[np].x +=                       pfact*ffx2(ii,jj,kk);
					}
					else if(ixyz ==2){
						pmparticles[np].vy +=                       vfact*ffx2(ii,jj,kk);
						pmparticles[np].y +=                       pfact*ffx2(ii,jj,kk);
					}
					else {
						pmparticles[np].vz +=                       vfact*ffx2(ii,jj,kk);
						pmparticles[np].z +=                       pfact*ffx2(ii,jj,kk);
					}
					np ++;
				}
			}
		}
		Free(ffx2);
	}

#ifdef INDEX
	void indexingbp_(pmparticletype *,long *,int *,int *, int *,int *,int *);
	indexingbp_(pmparticles,&np,&nzparticlewidth,&nspace,&nx,&ny,&nz);
#endif
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	for(i=0;i<np;i++){
#ifdef XYZDBL 
#else 
		pmparticles[i].x = fmod(pmparticles[i].x+rnx,rnx); 
		pmparticles[i].y = fmod(pmparticles[i].y+rny,rny); 
		pmparticles[i].z = fmod(pmparticles[i].z+rnz,rnz);
#endif

	}
	*mp = np;
}

