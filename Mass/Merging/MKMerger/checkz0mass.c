/* icc -o check check.c -DINDEX -DXYZDBL -g */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include "merger.h"
#include "mkmerger.h"



typedef struct MergeHistory{
	SavedHaloType a[75];
} MergeHistory;

int *nsteps;
float *reds;
float boxsize, hubble, npower, omep, omeplam, bias,astep,amax,anow,omepb;
int nx, nspace;
int msteps,*nsteps;


int main(int argc, char *argv[]){
	IDTYPE i,j,k;
	MergeHistory *bp;
	IDTYPE np;


	/*
	FILE *wp = fopen("MergingTree.dat","r");
	FILE *wp = fopen("sub2_MergingTree.dat","r");
	*/
	FILE *wp = fopen(argv[1],"r");
	fread(&boxsize, sizeof(float), 1,wp);/* simulation box size in h^-1 Mpc */
	fread(&hubble, sizeof(float), 1,wp); /* hubble parameter, for example, h=0.72 */
	fread(&npower, sizeof(float), 1,wp); /* power spectral index */
	fread(&omep, sizeof(float), 1,wp); /* omega_matter at present */
	fread(&omepb, sizeof(float), 1,wp); /* Omega_baryon0 at present */
	fread(&omeplam, sizeof(float), 1,wp); /* Omega_lambda at present */
	fread(&bias, sizeof(float), 1,wp); /* bias factor or inverse of sigma_8 */
	fread(&astep, sizeof(float), 1,wp); /* simulation astep */
	fread(&msteps, sizeof(float), 1,wp); /* number of steps in the merging tree */
	nsteps = (int*)malloc(sizeof(int)*msteps);
	reds = (float*)malloc(sizeof(float)*msteps);
	fread(nsteps, sizeof(int), msteps,wp); /* step number save to merging tree */
	fread(reds, sizeof(float), msteps,wp); /* redshift info. save to merging tree */


	bp = (MergeHistory*)malloc(sizeof(MergeHistory)*100000);
	size_t tcount,icount = 0;
	tcount = 0;

	while((np = fread(bp,sizeof(MergeHistory),100000,wp))>0){
		for(i=0;i<np;i++){
			tcount ++;
			if(bp[i].a[msteps-1].mass < 1.e2) {
				icount ++;
			}
			if(tcount%100000 == 0) {
				float frac = (float)icount/(float) tcount;
				printf("statistics %ld in %ld :: %g\n",icount, tcount, frac);
			}
		}
	}
	fclose(wp);

}
