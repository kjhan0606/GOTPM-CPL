#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include "pmheader.h"
#include "Memory.h"
#define max(a,b) ((a)>(b)? (a):(b))

#ifdef INDEX
void indexingbp_(pmparticletype *pmparticles,long *mp,int *mzp,int *mspace,
		int *mx,int *my,int *mz){
	long np,ii,i,j;
	long nzp;
	int nspace,nx,ny,nz;
	int myid,nid;
	long nzpcum;
	long *nzps,npinslice;
	indextype iindx,indxoffset;
	MPI_Status mpi_status;

	np = *mp; nzp = *mzp; nspace = *mspace;
	nx = *mx;ny = *my; nz= *mz;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	npinslice = (nx/nspace)*(ny/nspace);
	nzps = (long *)malloc(sizeof(long)*nid);
	MPI_Gather(&nzp,1,MPI_LONG,nzps,1,MPI_LONG,0,MPI_COMM_WORLD);
	if(myid==0){
		nzpcum = 0;
		for(i=1;i<nid;i++){
			nzpcum += nzps[i-1];
			MPI_Send(&nzpcum,1,MPI_LONG,i,i,MPI_COMM_WORLD);
		}
		nzpcum = 0; /* reinitialize for myid==0 */
	}
	else 
		MPI_Recv(&nzpcum,1,MPI_LONG,0,myid,MPI_COMM_WORLD,&mpi_status);
	iindx = 0; ii = 0;
	for(i=0;i<nzp;i+=max(nzp-1,1)){
		indxoffset = (indextype) npinslice *(indextype)(nzpcum+i);
		/*
		printf("S%d has z= %d & indxoffset =%ld\n",myid,nzpcum+i,indxoffset);
		*/
		for(j=0;j<npinslice;j++){
			/*
			pmparticles[ii].indx = indxoffset + iindx;
			iindx ++; ii++;
			*/
			pmparticles[ii].indx = indxoffset + j;
			ii++;
		}
	}
	free(nzps);
}
void indexingip_(pmparticletype *pmparticles,long *mp,long *premp,int *mzp,
		int *mspace, int *mx,int *my,int *mz){
	long np,nzp;
	int nspace,nx,ny,nz;
	long prenp,i,j;
	int myid,nid;
	long nzpcum;
	long *nzps,npinslice;
	long ii;
	indextype iindx,indxoffset;
	MPI_Status mpi_status;

	np = *mp; nzp = *mzp; nspace = *mspace;
	nx = *mx;ny = *my; nz= *mz;prenp = *premp;


	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	npinslice = (nx/nspace)*(ny/nspace);
	nzps = (long *)malloc(sizeof(long)*nid);
	MPI_Gather(&nzp,1,MPI_LONG,nzps,1,MPI_LONG,0,MPI_COMM_WORLD);
	if(myid==0){
		nzpcum = 0;
		for(i=1;i<nid;i++){
			nzpcum += nzps[i-1];
			MPI_Send(&nzpcum,1,MPI_LONG,i,i,MPI_COMM_WORLD);
		}
		nzpcum = 0; /* reinitialize for myid==0 */
	}
	else 
		MPI_Recv(&nzpcum,1,MPI_LONG,0,myid,MPI_COMM_WORLD,&mpi_status);
	iindx = 0; ii = 0;
	for(i=1;i<nzp-1;i++){
		indxoffset = (indextype) npinslice *(indextype)(nzpcum+i);
		for(j=0;j<npinslice;j++){
			/*
			pmparticles[prenp+ii].indx = indxoffset + iindx;
			iindx ++; ii++;
			*/
			pmparticles[prenp+ii].indx = indxoffset + j;
			ii++;
		}
	}
	free(nzps);
}
#endif
