#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<limits.h>
#include<mpi.h>
#define __MAIN_PQSORT__
#include "pqsort.h"
#undef __MAIN_PQSORT__
#include "merger.h"

/*
int PsortTrHaloIDTYPE(const void *a, const void *b){
	TrHalo *aa = (TrHalo*)a;
	TrHalo *bb = (TrHalo*)b;
	if(aa->aid < bb->aid) return -1;
	else if(aa->aid > bb->aid) return 1;
	else return 0;
}

int PsortTrHaloIDTYPE_L(const void *a, const void *lb,const void *lc){
	if(PsortTrHaloIDTYPE(a,lc) >= 0 || PsortTrHaloIDTYPE(a,lb) == -1) return 1;
	else return 0;
}
int PsortTrHaloIDTYPE_R(const void *a, const void *lb,const void *lc,  void *blc){
	if(PsortTrHaloIDTYPE(a,lc) >=0 || PsortTrHaloIDTYPE(a,lb) == -1) {
		TrHalo *aa = (TrHalo*) a;
		TrHalo *lbb = (TrHalo*) lb;
		TrHalo *lcc = (TrHalo*) lc;
		TrHalo *blcc = (TrHalo*) blc;
		IDTYPE rdist,ldist;
		rdist = fabs(lcc->aid - aa->aid);
		rdist = min(rdist, blcc->aid-rdist);
		ldist = fabs(lbb->aid - aa->aid);
		ldist = min(ldist, blcc->aid-ldist);
		if(rdist <= ldist ) return 1;
		else return 0;
	}
	else return 0;
}
*/

GenerateParallelComp(TrHalo,INDXTYPE,mbp);
GenerateParallelComp(UpLink,INDXTYPE,mbp);
GenerateParallelComp(TrHalo,IDTYPE,nowhid);
GenerateParallelComp(TrHalo,IDTYPE,nexthid);
GenerateParallelComp(UpLink,IDTYPE,nexthid);
GenerateParallelComp(UpLink,IDTYPE,nowhid);
GenerateParallelComp(TrHalo,IDTYPE,aid);


#define SWAP(a,b,stmp) do{\
	memcpy(stmp,a,n_size);\
	memcpy(a,b,n_size);\
	memcpy(b,stmp,n_size);\
}while(0)



void my_MPI_Sendrecv(void *sbase, long long smem, MPI_Datatype datatype1, int dest, int stag,
	void *rbase, long long rmem, MPI_Datatype datatype2, int src, int rtag, MPI_Comm comm, MPI_Status *status){
	long nmax, nchunk = INT_MAX/2, max;
	long i;
	MPI_Reduce(&smem, &nmax, 1,MPI_LONG, MPI_MAX, 0, comm);
	max = nmax;
	MPI_Reduce(&rmem, &nmax, 1,MPI_LONG, MPI_MAX, 0, comm);
	if(nmax > max) max = nmax;
	MPI_Bcast(&max, 1, MPI_LONG,0,comm);
	for(i=0;i<max;i+=nchunk){
		int ismem, irmem;
		long long lsmem, lrmem;
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

void my_MPI_Send(void *sbase, long long smem, MPI_Datatype datatype, int dest, int stag, MPI_Comm comm){
	long long nchunk = INT_MAX/2;
	long long i;
	for(i=0;i<smem;i+=nchunk){
		int ismem;
		long long lsmem;
		lsmem = smem - i;
		if(lsmem > nchunk) ismem = nchunk;
		else if (lsmem <=0) ismem = 0;
		else ismem = lsmem;
		MPI_Send((char*) sbase+i, ismem, datatype, dest, stag, comm);
	}
}
void my_MPI_Recv(void *rbase, long long smem, MPI_Datatype datatype, int src, int rtag, MPI_Comm comm, MPI_Status *status){
	int nchunk = INT_MAX/2;
	long long i;
	for(i=0;i<smem;i+=nchunk){
		int ismem;
		long long lsmem;
		lsmem = smem - i;
		if(lsmem > nchunk) ismem = nchunk;
		else if (lsmem <=0) ismem = 0;
		else ismem = lsmem;
		MPI_Recv((char*) rbase+i, ismem, datatype, src, rtag, comm, status);
	}
}





void *pqsort(void *base,size_t *mmem, size_t n_size, 
	int (*CompR)(const void *, const void *, const void *, int), 
	int (*CompL)(const void *, const void *, const void *),
	int (*Comp)(const void *, const void *),
	void *localmin, void *localmax, MPI_Comm Comm){
	MPI_Status status;
	int src,dest;
	long nrecv,nsend;
	long commsize,i,j,k;
	int myid,nid;
	size_t nmem = *mmem;
	char swaptmp[n_size];
	MPI_Comm_rank(Comm,&myid);
	MPI_Comm_size(Comm,&nid);
	dest = (myid+nid+1)%nid;
	src = (myid+nid-1)%nid;
	char *left,*right;
	commsize =10;
	left = (char*)base;
	while(commsize){
		right = (char*)base + nmem*n_size;
		for(;left<right;){
			if(CompR(left,localmin,localmax,nid)){
				right -= n_size;
				SWAP(left,right,swaptmp);
			}
			else left += n_size;
		}
		char *bp2send = right;
		long  nsend = nmem-(bp2send-(char *)base)/n_size;
		MPI_Sendrecv(&nsend,1,MPI_LONG,dest,0,&nrecv,1,MPI_LONG,src,0,Comm,&status);
		MPI_Reduce(&nsend,&commsize,1,MPI_LONG,MPI_SUM,0,Comm);
		MPI_Bcast(&commsize,1,MPI_LONG,0,Comm);
		if(myid==0) {
			printf("+: Total commsize = %ld\n",commsize);
		}
		long net = nrecv-nsend;
		char *tmp = (char*)malloc(sizeof(char)*nsend *n_size);
/*
		i = 0;
		for(right = bp2send;right<(char*)base+nmem*n_size;right++) tmp[i++]=*right;
*/
		memcpy(tmp,right,nsend*n_size);
		if(nmem==0) base = (void*) malloc(sizeof(char)*n_size*net);
		else base = (void*)realloc(base,sizeof(char)*n_size*(nmem+net));
		left = (char *)base+(nmem-nsend)*n_size;
		my_MPI_Sendrecv(tmp,nsend*sizeof(char)*n_size,MPI_BYTE,dest,0,
				left,nrecv*sizeof(char)*n_size,MPI_BYTE,src,0,Comm,&status);
		free(tmp);
		nmem = nmem+net;
	}
	commsize = 100;
	dest = (myid+nid-1)%nid;
	src = (myid+nid+1)%nid;
	left=(char*)base;
	while(commsize){
		right = (char*)base + nmem*n_size;
		for(;left<right;){
			if(CompL(left,localmin,localmax)){
				right-= n_size;
				SWAP(left,right,swaptmp);
			}
			else left += n_size;
		}
		char *bp2send = right;
		long nsend = nmem - (bp2send-(char *)base)/n_size;
		MPI_Sendrecv(&nsend,1,MPI_LONG,dest,0,&nrecv,1,MPI_LONG,src,0,Comm,&status);
		MPI_Reduce(&nsend,&commsize,1,MPI_LONG,MPI_SUM,0,Comm);
		MPI_Bcast(&commsize,1,MPI_LONG,0,Comm);
		if(myid==0) {
			printf("-: Total commsize = %ld\n",commsize);
		}
		long net = nrecv-nsend;
		char *tmp = (char*)malloc(sizeof(char)*nsend*n_size);
/*
		i = 0;
		for(right = bp2send;right<(char*)base+nmem*n_size;right++) tmp[i++]=*right;
*/
		memcpy(tmp,right,nsend*n_size);
		if(nmem==0) base = (void*) malloc(sizeof(char)*n_size*net);
		else base = (void*)realloc(base,sizeof(char)*n_size*(nmem+net));
		left = (char *)base+(nmem-nsend)*n_size;
		my_MPI_Sendrecv(tmp,nsend*sizeof(char)*n_size,MPI_BYTE,dest,0,
				left,nrecv*sizeof(char)*n_size,MPI_BYTE,src,0,Comm,&status);
		free(tmp);
		nmem = nmem+net;
	}
	qsort(base,nmem,n_size,Comp);
	*mmem = nmem;
	return base;
} 

/*
size_t MyPqsort(void *base, size_t nmem, size_t n_size, int (*CompR)(const void *, const void *, const void *, const void *), 
	int (*CompL)(const void *, const void *, const void *), int (*Comp)(const void *, const void *),MPI_Comm Comm){
	
	nmem= pqsort(base,nmem, n_size, CompR, CompL, Comp,&localmin,&localmax,&basemax,&Comm);
	return nmem;
}
*/
