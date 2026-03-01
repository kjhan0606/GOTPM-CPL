/* Version 0.5 of parallel sorting */
/* Copyright to Juhan Kim (kjhan@kias.re.kr) */
/* onewaypqsort(), twowaypqsort(), pasort() */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<limits.h>
#include<mpi.h>
#define __MAIN_PQSORT__
#include "pqsort.SMG.h"
#undef __MAIN_PQSORT__
/*
#include "tree.h"
#include "peak.h"
*/
#include "merger.h"
#include "mkmerger.h"
#include "SMG.h"




/*
GenerateParallelComp(MBP,IDTYPE,nowhid);
*/


extern int iredshift;
GenerateParallelComp(GID,IDTYPE,gid);
GenerateParallelComp(GID,IDTYPE,newgid);
GenerateParallelComp(MajorGid,IDTYPE,gid);
GenerateParallelComp(MajorGid,IDTYPE,mgid);
GenerateParallelComp(GidMergeHistory,IDTYPE,gid);
/*
GenerateParallelComp(HID,IDTYPE,mbp);
GenerateParallelComp(HID,IDTYPE,gid);
GenerateParallelComp(MBP,IDTYPE,mbp);
*/


#define SWAP(a,b,stmp) do{\
	memcpy(stmp,a,n_size);\
	memcpy(a,b,n_size);\
	memcpy(b,stmp,n_size);\
}while(0)



void my_MPI_Sendrecv(void *sbase, long smem, MPI_Datatype datatype1, int dest, int stag,
	void *rbase, long rmem, MPI_Datatype datatype2, int src, int rtag, MPI_Comm comm, MPI_Status *status){
	long nmax, nchunk = INT_MAX/2, max;
	long i;
	MPI_Reduce(&smem, &nmax, 1,MPI_LONG, MPI_MAX, 0, comm);
	max = nmax;
	MPI_Reduce(&rmem, &nmax, 1,MPI_LONG, MPI_MAX, 0, comm);
	if(nmax > max) max = nmax;
	MPI_Bcast(&max, 1, MPI_LONG,0,comm);
	for(i=0;i<max;i+=nchunk){
		int ismem, irmem;
		long lsmem, lrmem;
		lsmem = smem -i;
		lrmem = rmem -i;
		if(lsmem > nchunk) ismem = nchunk;
		else if(lsmem <=0 ) ismem = 0;
		else ismem = lsmem;
		if(lrmem > nchunk) irmem = nchunk;
		else if(lrmem <=0 ) irmem = 0;
		else irmem = lrmem;
		MPI_Sendrecv((char*)sbase+i, ismem, datatype1,dest,stag, 
			(char*)rbase+i,irmem, datatype2, src, rtag, comm, status);
	}
}
void my_MPI_Send(void *sbase, long smem, MPI_Datatype datatype, int dest, int stag, MPI_Comm comm){
	long nchunk = INT_MAX/2;
	long i;
	for(i=0;i<smem;i+=nchunk){
		int ismem;
		long lsmem;
		lsmem = smem - i;
		if(lsmem > nchunk) ismem = nchunk;
		else if (lsmem <=0) ismem = 0;
		else ismem = lsmem;
		MPI_Send((char*) sbase+i, ismem, datatype, dest, stag, comm);
	}
}
void my_MPI_Recv(void *rbase, long smem, MPI_Datatype datatype, int src, int rtag, MPI_Comm comm, MPI_Status *status){
	long nchunk = INT_MAX/2;
	long i;
	for(i=0;i<smem;i+=nchunk){
		int ismem;
		long lsmem;
		lsmem = smem - i;
		if(lsmem > nchunk) ismem = nchunk;
		else if (lsmem <=0) ismem = 0;
		else ismem = lsmem;
		MPI_Recv((char*) rbase+i, ismem, datatype, src, rtag, comm, status);
	}
}





void *pqsort(void *base,size_t *mmem, size_t n_size, 
	int (*CompR)(const void *, const void *, const void *), 
	int (*CompL)(const void *, const void *, const void *),
	int (*Comp)(const void *, const void *),
	void *localmin, void *localmax, MPI_Comm Comm){
	MPI_Status status;
	int src,dest;
	long nrecv,nsend;
	long i,j,k;
	int myid,nid;
	size_t nmem = *mmem;
	char swaptmp[n_size];
	void *tlocalmin, *tlocalmax;
	int ipair;
	size_t nowmem;
	char *rbase = NULL;




	MPI_Comm_rank(Comm,&myid);
	MPI_Comm_size(Comm,&nid);


	tlocalmin = (void*)malloc(n_size);
	tlocalmax = (void*)malloc(n_size);

	nowmem = 0;
	{
		char *left,*right;
		left = (char*)base;
		right = (char*)base + nmem*n_size;
		memcpy(tlocalmin, localmin, n_size);
		memcpy(tlocalmax, localmax, n_size);
		for(;left<right;){
			if(CompR(left,tlocalmin,tlocalmax)){
				right -= n_size;
				SWAP(left,right,swaptmp);
			}
			else left += n_size;
		}
		long nsend = ( (char*)right - (char*)base)/n_size;
		if(myid==0) {
			nrecv = nsend;
			printf("+: %d ::: %ld :: %ld\n",0, nrecv, nsend); fflush(stdout);
		}
		rbase = (char*) base;
		base = (char*)malloc((nmem-nsend)*n_size);
		memcpy(base, right, (nmem-nsend)*n_size);
		nmem = nmem -(nsend);
		nowmem += nsend;
		rbase = (char*)realloc(rbase, nowmem*n_size);
	}


	for(ipair=1;ipair<nid;ipair++){
		dest = (myid+nid+ipair)%nid;
		src = (myid+nid-ipair)%nid;
		char *left,*right;
		left = (char*)base;
		right = (char*)base + nmem*n_size;
		if(ipair ==0){
			memcpy(tlocalmin, localmin, n_size);
			memcpy(tlocalmax, localmax, n_size);
		}
		else {
			MPI_Sendrecv(localmin,n_size, MPI_BYTE,src, 0, tlocalmin, n_size, MPI_BYTE,dest,0,Comm,&status);
			MPI_Sendrecv(localmax,n_size, MPI_BYTE,src, 1, tlocalmax, n_size, MPI_BYTE,dest,1,Comm,&status);
		}
		if(ipair !=0){
			for(;left<right;){
				if(!CompL(left,tlocalmin,tlocalmax)){
					right -= n_size;
					SWAP(left,right,swaptmp);
				}
				else left += n_size;
			}
		}
		else {
			for(;left<right;){
				if(!CompR(left,tlocalmin,tlocalmax)){
					right -= n_size;
					SWAP(left,right,swaptmp);
				}
				else left += n_size;
			}
		}
		long  nsend = nmem-((char*)right-(char *)base)/n_size;
		if(ipair ==0){
			nrecv = nsend;
		}
		else {
			MPI_Sendrecv(&nsend,1,MPI_LONG,dest,0,&nrecv,1,MPI_LONG,src,0,Comm,&status);
		}
		if(myid==0) {
			printf("+: %d ::: %ld :: %ld\n",ipair, nrecv, nsend); fflush(stdout);
		}
		if(rbase == NULL) rbase = (char*)malloc(sizeof(char)*nrecv*n_size);
		else rbase = (char*)realloc(rbase,sizeof(char)*(nrecv + nowmem) *n_size);
		if(ipair ==0){
			memcpy(rbase + nowmem*n_size, right, nrecv*n_size);
		}
		else {
			my_MPI_Sendrecv(right,nsend*sizeof(char)*n_size,MPI_BYTE,dest,0,
					rbase+nowmem*n_size,nrecv*sizeof(char)*n_size,MPI_BYTE,src,0,Comm,&status);
		}
		nmem = nmem-nsend;
		base = (void*)realloc(base, nmem * n_size);
		nowmem += nrecv;
	}

	free(base);
	free(tlocalmin);free(tlocalmax);
	qsort(rbase,nowmem,n_size,Comp);
	*mmem = nowmem;
	return rbase;
} 

void *pqsort_old(void *base,size_t *mmem, size_t n_size, 
	int (*CompR)(const void *, const void *, const void *), 
	int (*CompL)(const void *, const void *, const void *),
	int (*Comp)(const void *, const void *),
	void *localmin, void *localmax, MPI_Comm Comm){
	MPI_Status status;
	int src,dest;
	long nrecv,nsend;
	long i,j,k;
	int myid,nid;
	size_t nmem = *mmem;
	char swaptmp[n_size];
	void *tlocalmin, *tlocalmax;
	int ipair;
	size_t nowmem;
	char *rbase = NULL;




	MPI_Comm_rank(Comm,&myid);
	MPI_Comm_size(Comm,&nid);


	tlocalmin = (void*)malloc(n_size);
	tlocalmax = (void*)malloc(n_size);

	nowmem = 0;


	for(ipair=0;ipair<nid;ipair++){
		dest = (myid+nid+ipair)%nid;
		src = (myid+nid-ipair)%nid;
		char *left,*right;
		left = (char*)base;
		right = (char*)base + nmem*n_size;
		if(ipair ==0){
			memcpy(tlocalmin, localmin, n_size);
			memcpy(tlocalmax, localmax, n_size);
		}
		else {
			MPI_Sendrecv(localmin,n_size, MPI_BYTE,src, 0, tlocalmin, n_size, MPI_BYTE,dest,0,Comm,&status);
			MPI_Sendrecv(localmax,n_size, MPI_BYTE,src, 1, tlocalmax, n_size, MPI_BYTE,dest,1,Comm,&status);
		}
		if(ipair !=0){
			for(;left<right;){
				if(!CompL(left,tlocalmin,tlocalmax)){
					right -= n_size;
					SWAP(left,right,swaptmp);
				}
				else left += n_size;
			}
		}
		else {
			for(;left<right;){
				if(!CompR(left,tlocalmin,tlocalmax)){
					right -= n_size;
					SWAP(left,right,swaptmp);
				}
				else left += n_size;
			}
		}
		long  nsend = nmem-((char*)right-(char *)base)/n_size;
		if(ipair ==0){
			nrecv = nsend;
		}
		else {
			MPI_Sendrecv(&nsend,1,MPI_LONG,dest,0,&nrecv,1,MPI_LONG,src,0,Comm,&status);
		}
		if(myid==0) {
			printf("+: %d ::: %ld :: %ld\n",ipair, nrecv, nsend); fflush(stdout);
		}
		if(rbase == NULL) rbase = (char*)malloc(sizeof(char)*nrecv*n_size);
		else rbase = (char*)realloc(rbase,sizeof(char)*(nrecv + nowmem) *n_size);
		if(ipair ==0){
			memcpy(rbase + nowmem*n_size, right, nrecv*n_size);
		}
		else {
			my_MPI_Sendrecv(right,nsend*sizeof(char)*n_size,MPI_BYTE,dest,0,
					rbase+nowmem*n_size,nrecv*sizeof(char)*n_size,MPI_BYTE,src,0,Comm,&status);
		}
		nmem = nmem-nsend;
		base = (void*)realloc(base, nmem * n_size);
		nowmem += nrecv;
	}

	free(base);
	free(tlocalmin);free(tlocalmax);
	qsort(rbase,nowmem,n_size,Comp);
	*mmem = nowmem;
	return rbase;
} 

