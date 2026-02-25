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

#define den(i,j,k) (den[(i)+mx*((j)+ny*(k))])
#ifdef XYZDBL
#define FABS(a) fabs(a)
#define RINT(a) rint(a)
#else
#define FABS(a) fabsf(a)
#define RINT(a) rintf(a)
#endif

#ifndef _OPENMP
#define omp_get_thread_num() nullfct0()
#define omp_get_num_threads() nullfct1()
#endif


void tsc_omp(float *den, int local_nz, int nstart_z, pmparticletype *pmp, int np, float pmas){
	int nid,myid;
	int NZWIDTH=4;
	float  start_z;
	float pmas0,p05;
//	long long i,j,k;

//	double rngc;
	int usndsize,dgetsize,dsndsize,ugetsize,dest,src;
	MPI_Status status;
	int ierror,stag,rtag,tag;
	int local_z_start,local_ny_after_transpose,local_y_start_after_transpose,
		total_local_size;
	long long nx,ny,nz,mx,mgrid,ngrid,slicesize;
	long long *nxperiodic,*nyperiodic,*nzperiodic;
	long long *nxp,*nyp,*nzp,nspace;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	TIMER_START(59);
	nx = simpar.nx; ny = simpar.ny; nz = simpar.nz;
	nspace = simpar.nspace;
	start_z = nstart_z;
//	rngc = (double) nx *(double)ny *(double)nz;
	long long nz_per = local_nz+NZWIDTH; /* local cyclic width in z direction */

//	if(myid==0) printf("Now in the OMP tsc: %d %d %d %g\n", local_nz,nstart_z,np, pmas);

	nxperiodic = (long long *)Malloc(sizeof(long long)*(nx+10),PPTR(nxperiodic)); nxp = nxperiodic+5;
	nyperiodic = (long long *)Malloc(sizeof(long long)*(ny+10),PPTR(nyperiodic)); nyp = nyperiodic+5;
	nzperiodic = (long long *)Malloc(sizeof(long long)*(nz_per+10),PPTR(nzperiodic)); nzp = nzperiodic+5;
	{
		long long i;
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
		for(i=0;i<mgrid;i++) den[i] = 0;
//		if(myid==0) printf("Now in the OMP tsc: 1\n");
#else
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(i=0;i<ngrid;i++) den[i] = -1.;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(i=ngrid;i<mgrid;i++) den[i] = 0;
#endif
	}
#ifndef MULTIMASS
	pmas0 = pmas;
	p05 = pmas0 * 0.5;
#endif


	int nthreads, mthreads;
	float hwidth,xwidth;
#ifdef _OPENMP
#pragma omp parallel shared(nthreads)
#endif
	{
#ifdef _OPENMP
#pragma omp master
#endif
		{
			nthreads = omp_get_num_threads();
		}
	}
	xwidth = (float)nx/(float)nthreads;
	hwidth = 0.5*xwidth;
	if(hwidth < 5) {
		mthreads = (float)nx/(2*5.);
	}
	else {
		mthreads = nthreads;
	}

//	if(myid==0) printf("Now in the OMP tsc with thread= %d\n", mthreads);
	signed char *target = (signed char*)Malloc(sizeof(short)*np, PPTR(target));
#ifdef _OPENMP
	{
#ifdef XYZDBL
		double Xwidth = (double)nx/(double)mthreads;
#else
		float Xwidth = (float)nx/(float)mthreads;
#endif
		long long i;
#pragma omp parallel for num_threads(mthreads)
		for(i=0;i<np;i++){
#ifdef XYZDBL
			double xp;
#else
			float xp;
#endif
			xp = XofP(pmp+i);
			int j = (int)(xp/Xwidth) + 1;
			int k = (int)(xp/Xwidth+0.5L) + 1;
			if(j >= mthreads+1){
				target[i] = mthreads;
			}
			else if(k ==j) {
				target[i] = -j;
			}
			else {
				target[i] = j;
			}
		}
	}
#endif

#if defined(MULTIMASS)
#ifdef _OPENMP
#pragma omp parallel private(pmas0,p05) num_threads(mthreads)
#endif
#else
#ifdef _OPENMP
#pragma omp parallel num_threads(mthreads)
#endif
#endif
	{
		int idthread = omp_get_thread_num() + 1;
		long long i;
		for(i=0;i<np;i++){
			if(target[i] != -idthread) continue;
#ifdef XYZDBL
			double xp,yp,zp;
#else
			float xp,yp,zp;
#endif
			float wx1wy1,wx2wy1,wx3wy1;
			float wx1wy2,wx2wy2,wx3wy2;
			float wx1wy3,wx2wy3,wx3wy3;
			float xmin,ymin,zmin;
			long long xsign,ysign,zsign;
			long long nearx,neary,nearz,i1,j1,k1,i2,j2,k2,i3,j3,k3;
			float wx1,wy1,wz1,wx2,wy2,wz2,wx3,wy3,wz3;
			float xd1,yd1,zd1;
#if defined(MULTIMASS)
			pmas0= pmp[i].mass;
			p05 = pmas0 * 0.5;
#endif
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
			den(i1,j1,k1) += wx1wy1*wz1; 
			den(i2,j1,k1) += wx2wy1*wz1; 
			den(i3,j1,k1) += wx3wy1*wz1; 
			den(i1,j2,k1) += wx1wy2*wz1; 
			den(i2,j2,k1) += wx2wy2*wz1; 
			den(i3,j2,k1) += wx3wy2*wz1; 
			den(i1,j3,k1) += wx1wy3*wz1; 
			den(i2,j3,k1) += wx2wy3*wz1; 
			den(i3,j3,k1) += wx3wy3*wz1;
	
		    den(i1,j1,k2) += wx1wy1*wz2; 
			den(i2,j1,k2) += wx2wy1*wz2; 
			den(i3,j1,k2) += wx3wy1*wz2; 
			den(i1,j2,k2) += wx1wy2*wz2; 
			den(i2,j2,k2) += wx2wy2*wz2; 
			den(i3,j2,k2) += wx3wy2*wz2; 
			den(i1,j3,k2) += wx1wy3*wz2; 
			den(i2,j3,k2) += wx2wy3*wz2; 
			den(i3,j3,k2) += wx3wy3*wz2;
	
		    den(i1,j1,k3) += wx1wy1*wz3; 
			den(i2,j1,k3) += wx2wy1*wz3; 
			den(i3,j1,k3) += wx3wy1*wz3; 
			den(i1,j2,k3) += wx1wy2*wz3; 
			den(i2,j2,k3) += wx2wy2*wz3; 
			den(i3,j2,k3) += wx3wy2*wz3; 
			den(i1,j3,k3) += wx1wy3*wz3; 
			den(i2,j3,k3) += wx2wy3*wz3; 
			den(i3,j3,k3) += wx3wy3*wz3;
	
		}
	}
#if defined(MULTIMASS)
#ifdef _OPENMP
#pragma omp parallel private(pmas0,p05) num_threads(mthreads)
#endif
#else
#ifdef _OPENMP
#pragma omp parallel num_threads(mthreads)
#endif
#endif
	{
		int idthread = omp_get_thread_num()+1;
		long long i;
		for(i=0;i<np;i++){
			if(target[i] != idthread) continue;
#ifdef XYZDBL
			double xp,yp,zp;
#else
			float xp,yp,zp;
#endif
			float wx1wy1,wx2wy1,wx3wy1;
			float wx1wy2,wx2wy2,wx3wy2;
			float wx1wy3,wx2wy3,wx3wy3;
			float xmin,ymin,zmin;
			long long xsign,ysign,zsign;
			long long nearx,neary,nearz,i1,j1,k1,i2,j2,k2,i3,j3,k3;
			float wx1,wy1,wz1,wx2,wy2,wz2,wx3,wy3,wz3;
			float xd1,yd1,zd1;
#if defined(MULTIMASS)
			pmas0= pmp[i].mass;
			p05 = pmas0 * 0.5;
#endif
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
			den(i1,j1,k1) += wx1wy1*wz1; 
			den(i2,j1,k1) += wx2wy1*wz1; 
			den(i3,j1,k1) += wx3wy1*wz1; 
			den(i1,j2,k1) += wx1wy2*wz1; 
			den(i2,j2,k1) += wx2wy2*wz1; 
			den(i3,j2,k1) += wx3wy2*wz1; 
			den(i1,j3,k1) += wx1wy3*wz1; 
			den(i2,j3,k1) += wx2wy3*wz1; 
			den(i3,j3,k1) += wx3wy3*wz1;
	
		    den(i1,j1,k2) += wx1wy1*wz2; 
			den(i2,j1,k2) += wx2wy1*wz2; 
			den(i3,j1,k2) += wx3wy1*wz2; 
			den(i1,j2,k2) += wx1wy2*wz2; 
			den(i2,j2,k2) += wx2wy2*wz2; 
			den(i3,j2,k2) += wx3wy2*wz2; 
			den(i1,j3,k2) += wx1wy3*wz2; 
			den(i2,j3,k2) += wx2wy3*wz2; 
			den(i3,j3,k2) += wx3wy3*wz2;
	
		    den(i1,j1,k3) += wx1wy1*wz3; 
			den(i2,j1,k3) += wx2wy1*wz3; 
			den(i3,j1,k3) += wx3wy1*wz3; 
			den(i1,j2,k3) += wx1wy2*wz3; 
			den(i2,j2,k3) += wx2wy2*wz3; 
			den(i3,j2,k3) += wx3wy2*wz3; 
			den(i1,j3,k3) += wx1wy3*wz3; 
			den(i2,j3,k3) += wx2wy3*wz3; 
			den(i3,j3,k3) += wx3wy3*wz3;
	
		}
	}
	if(myid==0) printf("Now exit the main tsc\n");


	{
		long long i;

#ifdef VarPM
#ifdef INTEL /* intel compiler may experience problems in pushing memory copy toward address. 
			This is a optimization issue. */
		i=slicesize*(local_nz+2);
		while(i >=0){
			den[i+slicesize] = den[i]; i--;
		}
		i=0;
		while(i<slicesize){
			den[i] = den[i + (slicesize*(local_nz+NZWIDTH-1))]; i++;
		}
#else
		for(i=slicesize*(local_nz+2);i>=0;i--){
			den[i+slicesize] = den[i];
		}
		for(i=0;i<slicesize;i++){
			den[i] = den[i + (slicesize*(local_nz+NZWIDTH-1))];
		}
#endif
#else
		dest = (myid+1+nid)%nid;
		src = (myid-1+nid)%nid;
		stag = rtag = 0;
		usndsize = slicesize*2;
		dgetsize = usndsize;
		MPI_Sendrecv_replace((den+slicesize*local_nz),usndsize,
				MPI_REAL,dest,stag,src,rtag,MPI_COMM_WORLD,&status);
		for(i=0;i<slicesize*2;i++){
			den[i] += den[i+slicesize*local_nz];
		}
		ugetsize = (dsndsize = slicesize);
		dest = (myid-1+nid)%nid;
		src = (myid+1+nid)%nid;
		stag = rtag = 0;
		MPI_Sendrecv_replace((den+slicesize*(local_nz+NZWIDTH-1)),dsndsize,
				MPI_REAL,dest,stag,src,rtag,MPI_COMM_WORLD,&status);
		for(i=0;i<slicesize;i++){
			den[i+slicesize*(local_nz-1)] += den[i+slicesize*(local_nz+NZWIDTH-1)];
		}
#endif
	}
	Free(target);
	Free(nzperiodic); Free(nyperiodic); Free(nxperiodic);
	MPI_Barrier(MPI_COMM_WORLD);
	TIMER_STOP(59);
	if(myid==0) fprintf(stdout,"CPU(TSC)   = %f\n",ELAPSED_TIME(59));
}

#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))

void preBinDen_nurion(float *fftwmesh, int mode){
	long long i,j,k;
	int mx,my,mz;
	char filename[100];
	MPI_Status status;
	float *tmpden;
	int nbin = 4;
	float dfact = 1./(nbin*nbin*nbin);
	int myid, nid;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);

	if(mode==0) sprintf(filename,"BinnedDen.%.5d.%.5d",simpar.stepcount,myid);
	else if(mode==1) sprintf(filename,"BinnedVx.%.5d.%.5d",simpar.stepcount,myid);
	else if(mode==2) sprintf(filename,"BinnedVy.%.5d.%.5d",simpar.stepcount,myid);
	else if(mode==3) sprintf(filename,"BinnedVz.%.5d.%.5d",simpar.stepcount,myid);
	else {
		fprintf(stderr,"Error in the mode of the binned output \n");
		exit(999);
	}

	long nnx = 2*(simpar.nx/2+1);

	mx = simpar.nx/nbin;
	my = simpar.ny/nbin;
	mz = simpar.nz/nbin;

	int itag = 0;
	int isend,iget,WGroupSize;
	int src, tgt;
	isend = iget = 1;
	WGroupSize = WGROUPSIZE;
	src = myid-1;
	tgt = myid+1;
	if(RANKINGROUP(myid,WGroupSize) != 0 ) MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);
	if(simpar.local_nz >0){
		FILE *fp = fopen(filename,"w");
		fwrite(&mx, sizeof(int),1,fp);
		fwrite(&my, sizeof(int),1,fp);
		fwrite(&mz, sizeof(int),1,fp);
		long nPixSlice = mx*my;
		int islice,jslice;
		float *binden = (float*)Malloc(sizeof(float)*nPixSlice,PPTR(binden));
		int zstart, zfinal;
		zstart = simpar.local_z_start;
		zfinal = zstart + simpar.local_nz;
		for(islice=zstart;islice<zfinal;islice++){
			jslice = islice-simpar.local_z_start;
			if(islice==zstart || islice%nbin==0){
				for(j=0;j<nPixSlice;j++) binden[j] = 0;
			}
			int j2;
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for(j2=0;j2<my;j2++){
				long j1 = j2 * nbin;
				long i2;
				for(i2=0;i2<mx;i2++){
					long i1 = i2*nbin;
					long j,i;
					long k2 = i1 + nnx*(j1+simpar.ny*jslice);
					float sum=0;
					for(j=0;j<nbin;j++) for(i=0;i<nbin;i++){
						sum += (1.+fftwmesh[k2+i+nnx*j])*dfact;
					}
					binden[i2+mx*j2] += sum;
				}
			}
			if( (islice+1)%nbin ==0 || islice == zfinal-1){
				int ix = islice/nbin;
				fwrite(&ix, sizeof(int), 1,fp);
				fwrite(binden, sizeof(float), mx*my,fp);
			}
		}
		fclose(fp);
		Free(binden);
	}

	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) 
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);
}

