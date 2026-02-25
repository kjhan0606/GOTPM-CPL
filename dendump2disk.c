#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<stddef.h>
#include<mpi.h>



int dendump2dist_(float *den,int *nx,int *ny,int *local_nz, float *amax){
	int myid,nid;
	int i,j,k;
	int mx,my,mz;
	size_t nsize;
	FILE *wp;
	char outfile[100];
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	mx = *nx; my = *ny; mz = *local_nz;
	nsize = (size_t) mx*(size_t)my*(size_t)mz;

	for(i=0;i<nid;i++){
		sprintf(outfile,"INITIAL_ZEL_POT.%.5d",myid);
		if(myid==i){
			wp = fopen(outfile,"w");
			fwrite(&mx,sizeof(int),1,wp);
			fwrite(&my,sizeof(int),1,wp);
			fwrite(&mz,sizeof(int),1,wp);
			fwrite(amax,sizeof(float),1,wp);
			fwrite(den,sizeof(float),nsize,wp);
			fclose(wp);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	return 1;
}
