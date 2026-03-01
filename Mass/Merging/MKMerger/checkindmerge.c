/* icc -o check check.c -DINDEX -DXYZDBL -g */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include "merger.h"
#include "mkmerger.h"



float *reds;
float boxsize, hubble, npower, omep, omeplam, bias,astep,amax,anow,omepb;
int nx, nspace;
int msteps,*nsteps;


int main(int argc, char *argv[]){
	IDTYPE i,j,k;
	SavedHaloType *bp;
	IDTYPE np;


	/*
	FILE *wp = fopen("MergingTree.dat","r");
	FILE *wp = fopen("sub2_MergingTree.dat","r");
	*/
	FILE *fp = fopen(argv[1],"r");

	bp = (SavedHaloType*)malloc(sizeof(SavedHaloType)*100000);
	size_t tcount,icount = 0;
	tcount = 0;

	while((np = fread(bp,sizeof(SavedHaloType),100000,fp))>0){
		for(i=0;i<np;i++){
			tcount ++;
			if(bp[i].mbp <= 0) {
				icount ++;
			}
			if(tcount%1000000 == 0) {
				float frac = (float)icount/(float) tcount;
				printf("statistics %ld in %ld :: %g\n",icount, tcount, frac);
			}
		}
	}
	{
				float frac = (float)icount/(float) tcount;
				printf("statistics %ld in %ld :: %g\n",icount, tcount, frac);
	}
	fclose(fp);

}
