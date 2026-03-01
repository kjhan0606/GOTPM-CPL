#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>





int my_MPI_Bcast(void *base, int count, MPI_Datatype datatype, int root, MPI_Comm comm){

	int myid, nid;

	MPI_Comm_size(comm, &nid);
	MPI_Comm_rank(comm, &myid);
	int *map;
	int i,j,k,nstep;
	MPI_Status status;


	map = (int*)malloc(sizeof(int)*nid);

	for(i=0;i<nid;i++) map[i] = -1;
	nstep = 0;
	while((nid >> nstep) >0 ) nstep ++;
	if((nid &(1<<nstep)-1) == 0) nstep --;
	for(i=0;i<nstep;i++){
		int istart = (1<<i);
		int iwidth = (1<<(i+1));
		for(j=0;j<iwidth;j++){
			if(j > nid) break;
			map[j] = (istart + j)&(iwidth-1);
		}
		int ipos = (myid - root +nid)%nid;
		if (ipos <=iwidth && map[ipos] >=0){
			if(ipos >= istart ) {
				int src = (map[ipos]+root+nid)%nid;
				if(map[ipos] < nid) MPI_Recv(base, count, datatype, src,src,comm, &status);
			}
			else if(ipos < istart) {
				int tgt = (map[ipos]+root+nid)%nid;
				if(map[ipos] < nid) MPI_Send(base, count, datatype, tgt,myid,comm); 
			}
		}
		MPI_Barrier(comm);
	}
	free(map);
	return 1;
}


#define NSIZE 100000000
int main(int argc, char **argv){
	int myid, nid;
	float a[NSIZE];
	MPI_Init(&argc, &argv);
	int i;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	if(myid==1) for(i=0;i<NSIZE;i++) a[i] = i*3.45/2.34;


	for(i=0;i<NSIZE;i+= 1000000){
		if(myid==0) printf("################# %d ############\n",i);
		double s1,s2;
		s1 = MPI_Wtime();
		my_MPI_Bcast(a, i, MPI_FLOAT, 1, MPI_COMM_WORLD);
		s2 = MPI_Wtime();
		if(myid==0) printf("my wtime = %g\n",(s2-s1));
		s1 = MPI_Wtime();
		MPI_Bcast(a, i, MPI_FLOAT, 0, MPI_COMM_WORLD);
		s2 = MPI_Wtime();
		if(myid==0) printf("Default wtime = %g\n",(s2-s1));
	}
	MPI_Finalize();
	return 0;
}
