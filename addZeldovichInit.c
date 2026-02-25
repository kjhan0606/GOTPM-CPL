#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stddef.h>
#include<mpi.h>
#include "pmheader.h"
#include "Memory.h"

#define max(a,b) (a)>(b)? (a):(b)
#define fx(i,j,k) fx[(i)+mx*((j)+ny*(k))]
#define fy(i,j,k) fy[(i)+mx*((j)+ny*(k))]
#define fz(i,j,k) fz[(i)+mx*((j)+ny*(k))]
void getinnerparticlesforce_(long *mp,pmvrparticletype *pvr,int *nn,int *nxpos,
				int *nzparticlewidth, float *zstartpos,
				int *nzprefixel,float *omep,float *omeplam,
				float *ken, float *zstart, float *zwidth, float *vamp,
				int *mspace,float *omei,float *amax, float *fx,float *fy,float *fz){
		long mx;
		int nspace;
		long nx,ny,nz;
		float vfact;
		long np,i,j,k;
		long ii,jj,kk;

		nx = simpar.nx;
		ny = simpar.ny;
		nz = simpar.nz;
		mx = 2*(nx/2+1);
		nspace = simpar.nspace;

		vfact = *vamp;
		np = 0;
		MPI_Barrier(MPI_COMM_WORLD);
		/*
		printf("P%d has nzparticlewidth and nx, ny nzprefixel %d %d %d, %d, %g: %g %g %g\n",simpar.myid,
						(*nzparticlewidth),nx,nspace,(*nzprefixel),*amax,
						fx(nx-1,ny-1,(*nzparticlewidth)-1), 
						fy(nx-1,ny-1,(*nzparticlewidth)-1), 
						fz(nx-1,ny-1,(*nzparticlewidth)-1));
						*/
		for(k=1;k<(*nzparticlewidth)-1;k++){
			kk = nspace*k + (*nzprefixel);
			for(j=0;j<ny/nspace;j++){
				jj = nspace*j;
				/*
				if(simpar.myid == 290) printf("now %ld %ld   : %ld\n",jj,kk,np);
				*/
				for(i=0;i<nx/nspace;i++){
					ii = nspace*i;
					pvr[np].vx = vfact*fx(ii,jj,kk);
					pvr[np].vy = vfact*fy(ii,jj,kk);
					pvr[np].vz = vfact*fz(ii,jj,kk);
					np++;
				}
			}
		}
		/*
		printf("P%.4d is now exiting with np=%ld\n",simpar.myid,np);
		*/
		MPI_Barrier(MPI_COMM_WORLD);
		*mp = np;
}
void setinnerparticlesforce_(long *mp,pmparticletype *particles,int *nn,int *nxpos,
				     int *nzparticlewidth,float *zstartpos,
					 int *nzprefixel,float *omep,float *omeplam,float *ken,
					 float *zstart,float *zwidth,float *vamp,int *nspace,
					float *omei,float *amax){
		long i,j,k,ii,jj,kk;
		float *zfinal;
		int myid,nid;
		float rnx,rny,rnz;
		float pfact,rng;
		float ai = 1;
		int nx,ny,nz;
		long np,tnp;
		double delx,displ,mdispl,tdelx,tmdispl;

		MPI_Comm_rank(MPI_COMM_WORLD,&myid);
		MPI_Comm_size(MPI_COMM_WORLD,&nid);
		rnx = nx = simpar.nx;
		rny = ny = simpar.ny;
		rnz = nz = simpar.nz;
		rng = rnx;
		np = (*mp);
		pfact = rng*ai;
		for(k=1;k<(*nzparticlewidth)-1;k++){
			for(j=0;j<ny/(*nspace);j++){
				for(i=0;i<nx/(*nspace);i++){
#ifdef XYZDBL
					particles[np].x = 0;
					particles[np].y = 0;
					particles[np].z = 0;
#else
					particles[np].x = (*nspace)*i;
					particles[np].y = (*nspace)*j;
					particles[np].z = (*nspace)*k+(*zstartpos);
#endif
					np ++;
				}
			}
		}
#ifdef INDEX
		void indexingip_(pmparticletype *,long *,long *,int *,int *, int *,int *,int *);
		indexingip_(particles,&np,mp,nzparticlewidth,nspace,&nx,&ny,&nz);
#endif
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
		for(i=(*mp);i<np;i++){
#ifdef XYZDBL
			particles[i].x += particles[i].vx*pfact; 
			particles[i].y += particles[i].vy*pfact; 
			particles[i].z += particles[i].vz*pfact; 
#else
			particles[i].x = fmod(particles[i].x +particles[i].vx*pfact+rnx,rnx); 
			particles[i].y = fmod(particles[i].y +particles[i].vy*pfact+rny,rny); 
			particles[i].z = fmod(particles[i].z +particles[i].vz*pfact+rnz,rnz);
#endif
		}
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
			tdelx = 0.5*8.*pi/3.L*(1.L/ai+1.L/(*omei)-1.L)*tdelx;
			constant = ai*ai*ai * tdelx + poten+ktot;
			printf("poten,ken,const = %g %g %g\n",poten,tdelx,constant);
		}
		*ken = tdelx;

		*mp  = np;
}

void setboundaryparticlesforce_(long *mp,pmparticletype *pmparticles,
				int *nn, int *nxpos,int *nzparticlewidth,
				float *zstartpos, int *nzprefixel,
				int *omep,int *omeplam,
				float *ken, float *zstart, float *zwidth,
				float *vamp,int *nspace,
				float *omei,float *amax,
				float *fx,float *fy,float *fz){
		float zfinal,ai;
		int myid,nid;
		float vfact,rng,pfact;
		float rnx,rny,rnz;
		long np,i,j,k,mx;
		long ii,jj,kk;
		long jidx,jkidx,kidx;
		int nzstride;
		int nx,ny,nz;
		ai = 1;

		MPI_Comm_rank(MPI_COMM_WORLD,&myid);
		MPI_Comm_size(MPI_COMM_WORLD,&nid);
		vfact = *vamp;
		rng = simpar.nx;
		nx = simpar.nx;
		ny = simpar.ny;
		nz = simpar.nz;
		rnx = nx; rny =ny; rnz = nz;
		mx = 2*(nx/2+1);
		np = 0;
		nzstride  = max((*nzparticlewidth)-1,1);
		for(k=0;k<(*nzparticlewidth);k+=nzstride){
			kk = (*nspace)*k + (*nzprefixel);
			for(j=0;j<ny/(*nspace);j++){
				jj = (*nspace)*j;
				for(i=0;i<nx/(*nspace);i++){
#ifdef XYZDBL
					pmparticles[np].x = 0;
					pmparticles[np].y = 0;
					pmparticles[np].z = 0;
#else
					pmparticles[np].x = (*nspace)*i;
					pmparticles[np].y = (*nspace)*j;
					pmparticles[np].z = (*nspace)*k+(*zstartpos);
#endif
					ii = (*nspace)*i;
					pmparticles[np].vx = vfact*fx(ii,jj,kk);
					pmparticles[np].vy = vfact*fy(ii,jj,kk);
					pmparticles[np].vz = vfact*fz(ii,jj,kk);
					np ++;
				}
			}
		}
#ifdef INDEX
		void indexingbp_(pmparticletype *,long *,int *,int *, int *,int *,int *);
		indexingbp_(pmparticles,&np,nzparticlewidth,nspace,&nx,&ny,&nz);
#endif
		pfact = (rng)*ai;
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
		for(i=0;i<np;i++){
#ifdef XYZDBL
			pmparticles[i].x += pmparticles[i].vx*pfact; 
			pmparticles[i].y += pmparticles[i].vy*pfact; 
			pmparticles[i].z += pmparticles[i].vz*pfact; 
#else
			pmparticles[i].x += pmparticles[i].vx*pfact; 
			pmparticles[i].x = fmod(pmparticles[i].x+rnx,rnx);
			pmparticles[i].y += pmparticles[i].vy*pfact; 
			pmparticles[i].y = fmod(pmparticles[i].y+rny,rny);
			pmparticles[i].z += pmparticles[i].vz*pfact; 
			pmparticles[i].z = fmod(pmparticles[i].z+rnz,rnz);
#endif
		}
		*mp = np;
}

