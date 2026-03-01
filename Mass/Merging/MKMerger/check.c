/* icc -o check check.c -DINDEX -DXYZDBL -g */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include "merger.h"
#include "mkmerger.h"



typedef struct MergeHistory{
	/*
	SavedHaloType a[133];
	*/
	SavedHaloType a[130];
} MergeHistory;

int *nsteps;
float *reds;
float boxsize, hubble, npower, omep, omeplam, bias,astep,amax,anow,omepb;
float wlam1, wlam0;
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
	FILE *fp = fopen(argv[1],"r");
	FILE *wp = fopen(argv[2],"w");
	fread(&boxsize, sizeof(float), 1,fp);/* simulation box size in h^-1 Mpc */
	fread(&hubble, sizeof(float), 1,fp); /* hubble parameter, for example, h=0.72 */
	fread(&npower, sizeof(float), 1,fp); /* power spectral index */
	fread(&omep, sizeof(float), 1,fp); /* omega_matter at present */
	fread(&omepb, sizeof(float), 1,fp); /* Omega_baryon0 at present */
	fread(&omeplam, sizeof(float), 1,fp); /* Omega_lambda at present */
	fread(&wlam0, sizeof(float), 1,fp); /* Omega_lambda at present */
	fread(&wlam1, sizeof(float), 1,fp); /* Omega_lambda at present */
	fread(&bias, sizeof(float), 1,fp); /* bias factor or inverse of sigma_8 */
	fread(&astep, sizeof(float), 1,fp); /* simulation astep */
	fread(&msteps, sizeof(float), 1,fp); /* number of steps in the merging tree */
	nsteps = (int*)malloc(sizeof(int)*msteps);
	reds = (float*)malloc(sizeof(float)*msteps);
	fread(nsteps, sizeof(int), msteps,fp); /* step number save to merging tree */
	fread(reds, sizeof(float), msteps,fp); /* redshift info. save to merging tree */


	bp = (MergeHistory*)malloc(sizeof(MergeHistory)*1000000);
	while((np = fread(bp,sizeof(MergeHistory),1000000,fp))>0){
		printf(".\n");
		for(i=0;i<np;i++){
		}
	}
}
