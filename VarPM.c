/* This modifies the paddingpotentialmesh for the better 
 * performance in speed. 08 May 2008 */
#include <stdio.h>
#include <stdlib.h>
//#include<srfftw_mpi.h>
#include<math.h>
#include<omp.h>
#include<mpi.h>
#include "pmheader.h"
#include "fortran.h"
#include "varpm.h"
#include "Memory.h"
#include "slicemigrate.h"
#define NOBOUND -999999
static int myid,nid;
static int nx,ny,nz;
static int *nslicerank,*tmpnslicerank;
static NZBound *nzbound;
#define MIN(a,b) ((a)<(b)? (a):(b))
#define MAX(a,b) ((a)>(b)? (a):(b))

void VarPM_tsc2(float *den,int Nx,int Ny,int Nz, float lzfloor, float lzceil, 
		pmparticletype *particles, int nn,int nxpos,int np, float pmas, 
		float *fftwmesh, int after_local_nz, int after_z_start){
	int nzfloor,nzceil,nzwidth,mode;
	int before_local_nz;
	MPI_Status status;
	int i,j,k;
	int now_local_nz;
	int nslicesize,now_z_start,nslicefinal;
	double sum,tsum;

	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	nx = Nx;ny=Ny;nz=Nz;
	nslicesize = 2*(nx/2+1)*ny;

	nzfloor = (int) lzfloor;
	nzceil = (int)(lzceil);
	if(((float)nzceil) == lzceil) nzceil = nzceil - 1;
	before_local_nz = (nzceil-nzfloor)+1;


#ifdef TSC_OLD
	tsc_(den,&before_local_nz,&nzfloor,
			particles,&nn,&nxpos,&np,&pmas,&myid,&nid);
#else
	{
		void tsc(float *, int, int, pmparticletype *,int,float);
		tsc(den,before_local_nz,nzfloor,particles,np,pmas);
	}
#endif
	now_z_start = nzfloor-1; /* Due to the TSC, lowering one slice is needed */
	now_local_nz = before_local_nz+3; /* Three padding slices are needed */

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<after_local_nz*(2*(nx/2+1))*ny;i++) fftwmesh[i] = 0.;

	mode = 0; /* This is the density mode */
	tsc_fftw_mesh(den,now_z_start,now_local_nz,fftwmesh,after_z_start,
			after_local_nz,nx,ny,nz,mode);
}
void VarPM_veltsc2(float *vel, int Nx,int Ny,int Nz, float lzfloor, float lzceil, 
		pmparticletype *particles, int nn,int nxpos,int np, float pmas, 
		float *fftwmeshvel, int after_local_nz, int after_z_start,int mode){
	int nzfloor,nzceil,nzwidth;
	int before_local_nz;
	MPI_Status status;
	int i,j,k;
	int now_local_nz;
	int nslicesize,now_z_start,nslicefinal;
	double sum,tsum;

	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	nx = Nx;ny=Ny;nz=Nz;
	nslicesize = 2*(nx/2+1)*ny;

	nzfloor = (int) lzfloor;
	nzceil = (int)(lzceil);
	if(((float)nzceil) == lzceil) nzceil = nzceil - 1;
	before_local_nz = (nzceil-nzfloor)+1;


#ifdef TSC_OLD
	fprintf(stderr,"Not yet implemented for veltsc\n");
	exit(99);
	
#else
	{
		void veltscX(float *, int, int, pmparticletype *,int,float);
		void veltscY(float *, int, int, pmparticletype *,int,float);
		void veltscZ(float *, int, int, pmparticletype *,int,float);
		if(mode==1)      veltscX(vel,before_local_nz,nzfloor,particles,np,pmas);
		else if(mode==2) veltscY(vel,before_local_nz,nzfloor,particles,np,pmas);
		else if(mode==3) veltscZ(vel,before_local_nz,nzfloor,particles,np,pmas);
	}
#endif
	now_z_start = nzfloor-1; /* Due to the TSC, lowering one slice is needed */
	now_local_nz = before_local_nz+3; /* Three padding slices are needed */

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<after_local_nz*(2*(nx/2+1))*ny;i++) fftwmeshvel[i] = 0.;
	tsc_fftw_mesh(vel,now_z_start,now_local_nz,fftwmeshvel,after_z_start,
			after_local_nz,nx,ny,nz,mode);
}

void fftw_fda4_mesh(float *,int, int,float *,int,int,int,int,int, int);
void VarPM_fda4_mesh( float *fftwmesh, int now_local_nz, int now_z_start,
		float *den, float lzfloor, float lzceil, int Nx,int Ny,int Nz){
	int after_z_start,nzwidth,after_z_final;
	int after_local_nz;
	MPI_Status status;
	int i,j,k;
	int nslicesize,nslicefinal;
	double sum;

	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	nx = Nx;ny=Ny;nz=Nz;
	nslicesize = 2*(nx/2+1)*ny;

	after_z_start = (int) lzfloor-3;
	after_z_final = (int) lzceil+4;

	after_local_nz = (after_z_final-after_z_start)+1;
	/*
    printf("P%d has zslice from %d + %d\n",myid,now_z_start,now_local_nz);
    printf("P%d will have zslice from %d + %d\n",myid,after_z_start,
		after_local_nz);
	*/
	fftw_fda4_mesh(fftwmesh,now_z_start,now_local_nz,den,after_z_start,
			after_local_nz,nx,ny,nz,1);
}

/* this routine outputs den[*] from input fftwmesh[*] */
void adVarPM_fda4_mesh( float *fftwmesh, int now_local_nz, int now_z_start,
		float *den, float lzfloor, float lzceil, int Nx,int Ny,int Nz){
	int after_z_start,nzwidth,after_z_final;
	int after_local_nz;
	MPI_Status status;
	int i,j,k;
	int nslicesize,nslicefinal;
	double sum;

	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0) printf("P%d is in adVarPM_fda4_mesh\n",myid);
	nx = Nx;ny=Ny;nz=Nz;
	nslicesize = 2*(nx/2+1)*ny;

	after_z_start = (int) lzfloor;
	after_z_final = (int) lzceil;
	if((float)after_z_final == lzceil) after_z_final --;

	after_local_nz = (after_z_final-after_z_start)+1;
	/*
    printf("P%d has zslice from %d to %d (fftw->fda4)\n",
			myid,now_z_start,now_z_start+now_local_nz-1);
    printf("P%d will have zslice from %d to %d (fftw->fda4)\n",
			myid,after_z_start,after_z_start+after_local_nz-1);
			*/
	fftw_fda4_mesh(fftwmesh,now_z_start,now_local_nz,den,after_z_start,
			after_local_nz,nx,ny,nz,0);
	/* free fftwmesh */
}
typedef struct IZDomain{
	int izstart,izfinal;
} IZDomain;
/* Beforehand, den array should be expanded -3, ..., +4 .
 * This routine is rather slow than that in the slicemigrate.c  which
 * uses the ring interconnecting network structure.
 * This routine uses nid blocked communications but the "stride" sending
 * only uses nid-1 in the worst case. Generally, "striding" can 
 * transfer data in 3 blocked communications.
 * */
#ifdef SLOW_PADDING
void paddingpotentialmesh(float *den,int izstart, int izwidth,
		int local_nz_varpm, int Nx,int Ny, int Nz, int m3z,int p4z){
	int izfinal;
	int i,j,k,nowz,slicepixel,ii;
	int mpitag=1;
	int src,dest;
	int myXizstart,myXizfinal;
	MPI_Status mpistatus;
	IZDomain *izdomain,*orgizdomain;
	
	nx = Nx; ny = Ny; nz = Nz;
	slicepixel = 2*(nx/2+1)*ny;
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	izdomain = (IZDomain *) Malloc(sizeof(IZDomain)*nid,PPTR(izdomain));
	orgizdomain = (IZDomain *) Malloc(sizeof(IZDomain)*nid,PPTR(orgizdomain));
	izfinal = izstart + izwidth - 1;
	izdomain[myid].izstart = izstart;
	izdomain[myid].izfinal = izfinal;
	if(myid==0){
		for(i=1;i<nid;i++)
			MPI_Recv(izdomain+i,sizeof(IZDomain),MPI_BYTE,i,mpitag,
					MPI_COMM_WORLD,&mpistatus);
		for(i=0;i<nid;i++)
			orgizdomain[i] = izdomain[i];
	}

	else
		MPI_Send(izdomain+myid,sizeof(IZDomain),MPI_BYTE,0,mpitag,
				MPI_COMM_WORLD);
	/* This is to avoid of overlapping of a slice between adjacent domains */
	if(myid==0){
		for(i=0;i<nid-1;i++){
			izdomain[i+1].izstart = izdomain[i].izfinal+1;
		}
	}
	MPI_Bcast(izdomain,sizeof(IZDomain)*nid,MPI_BYTE,0,MPI_COMM_WORLD);
	MPI_Bcast(orgizdomain,sizeof(IZDomain)*nid,MPI_BYTE,0,MPI_COMM_WORLD);
	myXizstart = izdomain[myid].izstart;
	myXizfinal = izdomain[myid].izfinal;
	if(myid==0) printf("P%d has zslices %d %d from %d %d (padding 7 slices)\n",
			myid,myXizstart,myXizfinal,izstart, izfinal);

	for(i=0;i<nid;i++){
		int targetizstart,targetizfinal;
		targetizstart = orgizdomain[i].izstart;
		targetizfinal = orgizdomain[i].izfinal;
		if(myid==i){
			ii = 0;
			for(j=1;j<p4z+1;j++){
				nowz = (j+izfinal+nz)%nz;
				for(k=0;k<nid;k++){
					if(nowz >= izdomain[k].izstart && nowz<=izdomain[k].izfinal)
						break;
				}
				MPI_Recv(den+(local_nz_varpm+ii)*slicepixel,slicepixel,MPI_FLOAT,k,
						nowz,MPI_COMM_WORLD,&mpistatus);
				ii++;
			}
			for(j=m3z;j<0;j++){
				nowz = (j+izstart+nz)%nz;
				for(k=0;k<nid;k++){
					if(nowz >= izdomain[k].izstart && nowz<=izdomain[k].izfinal)
						break;
				}
				MPI_Recv(den+(local_nz_varpm+ii)*slicepixel,slicepixel,MPI_FLOAT,k,
						nowz,MPI_COMM_WORLD,&mpistatus);
				ii++;
			}
		}
		else {
			for(j=targetizfinal+1;j<targetizfinal+p4z+1;j++){
				nowz = (j+nz)%nz;
				if(nowz >= myXizstart && nowz <= myXizfinal){
					MPI_Send(den+(nowz-izstart)*slicepixel,slicepixel,
							MPI_FLOAT,i,nowz,MPI_COMM_WORLD);
				}
			}
			for(j=targetizstart+m3z;j<targetizstart;j++){
				nowz = (j+nz)%nz;
				if(nowz >= myXizstart && nowz <= myXizfinal){
					MPI_Send(den+(nowz-izstart)*slicepixel,slicepixel,
							MPI_FLOAT,i,nowz,MPI_COMM_WORLD);
				}
			}
		}
		/*
		if(i%10 == 1) MPI_Barrier(MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		*/
	}
	if(myid==0) printf("P%d has passed\n",myid);
	Free(orgizdomain);
	Free(izdomain);
}
#else
//#error The function of Paddingpotentialmesh is Not Yet Fully Tested! \
 Please define SLOW_PADDING as -DSLOW_PADDING
#define PaddingPair(istep,step){\
	for(i=istep;i<nid;i+=step){\
		int targetizstart,targetizfinal;\
		targetizstart = orgizdomain[i].izstart;\
		targetizfinal = orgizdomain[i].izfinal;\
		if(myid==i){\
			ii = 0;\
			for(j=1;j<p4z+1;j++){\
				nowz = (j+izfinal+nz)%nz;\
				for(k=0;k<nid;k++)\
					if(nowz >= izdomain[k].izstart && nowz<=izdomain[k].izfinal) break;\
				MPI_Recv(den+(local_nz_varpm+ii)*slicepixel,slicepixel,MPI_FLOAT,k,\
						nowz,MPI_COMM_WORLD,&mpistatus); ii++;\
			}\
			for(j=m3z;j<0;j++){\
				nowz = (j+izstart+nz)%nz;\
				for(k=0;k<nid;k++)\
					if(nowz >= izdomain[k].izstart && nowz<=izdomain[k].izfinal) break;\
				MPI_Recv(den+(local_nz_varpm+ii)*slicepixel,slicepixel,MPI_FLOAT,k,\
						nowz,MPI_COMM_WORLD,&mpistatus); ii++;\
			}\
		}\
		else {\
			for(j=targetizfinal+1;j<targetizfinal+p4z+1;j++){\
				nowz = (j+nz)%nz;\
				if(nowz >= myXizstart && nowz <= myXizfinal)\
					MPI_Send(den+(nowz-izstart)*slicepixel,slicepixel,\
							MPI_FLOAT,i,nowz,MPI_COMM_WORLD);\
			}\
			for(j=targetizstart+m3z;j<targetizstart;j++){\
				nowz = (j+nz)%nz;\
				if(nowz >= myXizstart && nowz <= myXizfinal)\
					MPI_Send(den+(nowz-izstart)*slicepixel,slicepixel,\
							MPI_FLOAT,i,nowz,MPI_COMM_WORLD);\
			}\
		}\
	}\
}
void paddingpotentialmesh(float *den,int izstart, int izwidth,
		int local_nz_varpm, int Nx,int Ny, int Nz, int m3z,int p4z){
	int izfinal;
	int i,j,k,nowz,slicepixel,ii,jj;
	int mpitag=1;
	int src,dest;
	int myXizstart,myXizfinal;
	MPI_Status mpistatus;
	IZDomain *izdomain,*orgizdomain;
	
	nx = Nx; ny = Ny; nz = Nz;
	slicepixel = 2*(nx/2+1)*ny;
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	izdomain = (IZDomain *) Malloc(sizeof(IZDomain)*nid,PPTR(izdomain));
	orgizdomain = (IZDomain *) Malloc(sizeof(IZDomain)*nid,PPTR(orgizdomain));
	izfinal = izstart + izwidth - 1;
	izdomain[myid].izstart = izstart;
	izdomain[myid].izfinal = izfinal;
	if(myid==0){
		for(i=1;i<nid;i++)
			MPI_Recv(izdomain+i,sizeof(IZDomain),MPI_BYTE,i,mpitag,
					MPI_COMM_WORLD,&mpistatus);
		for(i=0;i<nid;i++)
			orgizdomain[i] = izdomain[i];
	}

	else
		MPI_Send(izdomain+myid,sizeof(IZDomain),MPI_BYTE,0,mpitag,
				MPI_COMM_WORLD);
	/* This is to avoid of overlapping of a slice between adjacent domains */
	if(myid==0){
		for(i=0;i<nid-1;i++){
			izdomain[i+1].izstart = izdomain[i].izfinal+1;
		}
	}
	MPI_Bcast(izdomain,sizeof(IZDomain)*nid,MPI_BYTE,0,MPI_COMM_WORLD);
	MPI_Bcast(orgizdomain,sizeof(IZDomain)*nid,MPI_BYTE,0,MPI_COMM_WORLD);
	myXizstart = izdomain[myid].izstart;
	myXizfinal = izdomain[myid].izfinal;
	if(myid==0) printf("P%d has zslices %d %d from %d %d (padding 7 slices)\n",
			myid,myXizstart,myXizfinal,izstart, izfinal);
	for(jj=0;jj<4;jj++) PaddingPair(jj,4);

	if(myid==0) printf("P%d has passed\n",myid);
	Free(orgizdomain);
	Free(izdomain);
}

#endif
/* izslice should be -3,-2,-1, and local_nz_varpm+1,...,local_nz_varpm+4 */
/* upperORbelow should be 1 (upper sending) or -1 (lower sending) */
void paddingonepotentialmesh(float *den,int izstart, int izwidth,
		int local_nz_varpm, int Nx,int Ny, int Nz, int izslice,int upperORbelow){
	int izfinal;
	int i,j,k,nowz,slicepixel,ii;
	int mpitag=1;
	int src,dest;
	int myXizstart,myXizfinal;
	MPI_Status mpistatus;
	IZDomain *izdomain,*orgizdomain;
	
	nx = Nx; ny = Ny; nz = Nz;
	slicepixel = 2*(nx/2+1)*ny;
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	izdomain = (IZDomain *) Malloc(sizeof(IZDomain)*nid,PPTR(izdomain));
	orgizdomain = (IZDomain *) Malloc(sizeof(IZDomain)*nid,PPTR(orgizdomain));
	izfinal = izstart + izwidth - 1;
	izdomain[myid].izstart = izstart;
	izdomain[myid].izfinal = izfinal;
	if(myid==0){
		for(i=1;i<nid;i++)
			MPI_Recv(izdomain+i,sizeof(IZDomain),MPI_BYTE,i,mpitag,
					MPI_COMM_WORLD,&mpistatus);
		for(i=0;i<nid;i++)
			orgizdomain[i] = izdomain[i];
	}
	else
		MPI_Send(izdomain+myid,sizeof(IZDomain),MPI_BYTE,0,mpitag,
				MPI_COMM_WORLD);
	/* This is to avoid of overlapping of a slice between adjacent domains */
	if(myid==0){
		for(i=0;i<nid-1;i++){
			izdomain[i+1].izstart = izdomain[i].izfinal+1;
		}
	}
	MPI_Bcast(izdomain,sizeof(IZDomain)*nid,MPI_BYTE,0,MPI_COMM_WORLD);
	MPI_Bcast(orgizdomain,sizeof(IZDomain)*nid,MPI_BYTE,0,MPI_COMM_WORLD);
	myXizstart = izdomain[myid].izstart;
	myXizfinal = izdomain[myid].izfinal;
	/*
	printf("P%d has zslices %d %d from %d %d (for one slice padding)\n",
			myid,myXizstart,myXizfinal,izstart, izfinal);
			*/

	for(i=0;i<nid;i++){
		int targetizstart,targetizfinal;
		/*
		if(myid==i){
			targetizstart = izstart;
			targetizfinal = izfinal;
		}
		MPI_Bcast(&targetizstart,1,MPI_INT,i,MPI_COMM_WORLD);
		MPI_Bcast(&targetizfinal,1,MPI_INT,i,MPI_COMM_WORLD);
		*/
		targetizstart = orgizdomain[i].izstart;
		targetizfinal = orgizdomain[i].izfinal;
		if(myid==i){
			j=izslice;
			{
				if(upperORbelow == UpperSending) nowz = (j+izfinal+nz)%nz;
				else if(upperORbelow == BelowSending) nowz = (j+izstart+nz)%nz;
				for(k=0;k<nid;k++){
					if(nowz >= izdomain[k].izstart && nowz<=izdomain[k].izfinal)
						break;
				}
				MPI_Recv(den+(local_nz_varpm)*slicepixel,slicepixel,MPI_FLOAT,k,
						nowz,MPI_COMM_WORLD,&mpistatus);
				/*
				printf("P%d has a padding slice at %d from %d proc.\n",
						myid,local_nz_varpm+ii,k);
						*/
			}
		}
		else {
			j=izslice;
			{
				if(upperORbelow==UpperSending) nowz = (j+targetizfinal+nz)%nz;
				else if(upperORbelow==BelowSending) nowz = (j+targetizstart+nz)%nz;
				if(nowz >= myXizstart && nowz <= myXizfinal){
					/*
					printf("P%d is sending nowz=%d to i=%d\n",myid,nowz,i);
					*/
					MPI_Send(den+(nowz-izstart)*slicepixel,slicepixel,
							MPI_FLOAT,i,nowz,MPI_COMM_WORLD);
				}
			}
		}
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	Free(orgizdomain);
	Free(izdomain);
}
