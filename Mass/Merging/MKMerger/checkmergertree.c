/* icc -o check check.c -DINDEX -DXYZDBL -g */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include "merger.h"
#include "mkmerger.h"


#define NSTEPS 131

typedef struct MergeHistory{
	SavedHaloType a[NSTEPS];
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
	FILE *fp = fopen(argv[1],"r");
	fread(&boxsize, sizeof(float), 1,fp);/* simulation box size in h^-1 Mpc */
	fread(&hubble, sizeof(float), 1,fp); /* hubble parameter, for example, h=0.72 */
	fread(&npower, sizeof(float), 1,fp); /* power spectral index */
	fread(&omep, sizeof(float), 1,fp); /* omega_matter at present */
	fread(&omepb, sizeof(float), 1,fp); /* Omega_baryon0 at present */
	fread(&omeplam, sizeof(float), 1,fp); /* Omega_lambda at present */
	fread(&bias, sizeof(float), 1,fp); /* bias factor or inverse of sigma_8 */
	fread(&astep, sizeof(float), 1,fp); /* simulation astep */
	fread(&msteps, sizeof(float), 1,fp); /* number of steps in the merging tree */
	nsteps = (int*)malloc(sizeof(int)*msteps);
	reds = (float*)malloc(sizeof(float)*msteps);
	fread(nsteps, sizeof(int), msteps,fp); /* step number save to merging tree */
	fread(reds, sizeof(float), msteps,fp); /* redshift info. save to merging tree */


	bp = (MergeHistory*)malloc(sizeof(MergeHistory)*1000000);
	size_t tcount,icount = 0;
	tcount = 0;
	int nblock = 1000000,iblock=0;
	float block[nblock];
	float rminmass,rmaxmass;
	int nbin = 30;
	double bin[nbin+1],massbin[nbin+1],err[nbin+1];
	rminmass = log10(2.7E11);
	rmaxmass = log10(5E15);
	float binstep = (rmaxmass-rminmass)/(float)nbin;
	for(i=0;i<=nbin;i++) {
		bin[i] = 0;
		massbin[i] = rminmass + (rmaxmass-rminmass)/(double) nbin * i;
	}
	long long gid= 0;
	long long ncount = 0;
	double boxsize = 3150.;

	while((np = fread(bp,sizeof(MergeHistory),1000000,fp))>0){
		for(i=0;i<np;i++){
			int ii=0;
			for(j=0;j<NSTEPS;j++){
				if(bp[i].a[j].majorglobalid ==0 && bp[i].a[j+1].majorglobalid !=0) {
					ii ++;
				}
				else if(bp[i].a[j].majorglobalid !=0 && bp[i].a[j+1].majorglobalid ==0) {
					ii++;
				}
			}
			if(ii>=3){
				fprintf(stderr,"Error strange \n");
			}
			gid++;
		}
		printf("ruuning cumulative number is %ld ::: missing= %ld\n",gid,ncount);fflush(stdout);
	}
}
