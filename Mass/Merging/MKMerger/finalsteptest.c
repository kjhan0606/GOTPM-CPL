/* icc -o makeone makeone.c -g -lm -DINDEX -DXYZDBL
 * */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<string.h>
#include"merger.h"
#include"mkmerger.h"


#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

int *nsteps;
float *reds;
float boxsize, hubble, npower, omep, omeplam, bias,astep,amax,anow,omepb;
int nx, nspace;
int msteps,*nsteps;


int main (int argc, char **argv){
	FILE *fp,*wp;
	IDTYPE i,j,k,npmax;
	IDTYPE gidblock = 2000000,gidshift;
	SavedHaloType *tmp = (SavedHaloType*)malloc(sizeof(SavedHaloType)*gidblock);
	int np,mp,rflag;



	gidshift = 0;
	npmax = gidblock;

	wp = fopen("MergingTree.dat","r");
	fread(&boxsize, sizeof(float), 1,wp);
	fread(&hubble, sizeof(float), 1,wp);
	fread(&npower, sizeof(float), 1,wp);
	fread(&omep, sizeof(float), 1,wp);
	fread(&omepb, sizeof(float), 1,wp);
	fread(&omeplam, sizeof(float), 1,wp);
	fread(&bias, sizeof(float), 1,wp);
	fread(&astep, sizeof(float), 1,wp);
	fread(&msteps, sizeof(float), 1,wp);
	nsteps = (int*)malloc(sizeof(int)*msteps);
	reds = (float*)malloc(sizeof(float)*msteps);
	fread(nsteps, sizeof(int), msteps,wp);
	fread(reds, sizeof(float), msteps,wp);

	fflush(wp);
	IDTYPE ncount = 0,gid=0;

	for(k=0; ; k+= npmax){
		npmax = 0;
		i=74;
		{
			char infile[190]; sprintf(infile,"MergingTree.dat.%.5d",(int)i);
			fp = fopen(infile,"r");
			fseek(fp, k*sizeof(SavedHaloType),SEEK_SET);
			np = fread(tmp,sizeof(SavedHaloType),gidblock,fp);
			fflush(stdout);
			fclose(fp);
			if(np>0){
				for(j=0;j<np;j++){
					if(tmp[j].mbp <=0) {
						ncount ++;
					};
				}
			}
			npmax = MAX(npmax,np);
		}
		if(npmax ==0) break;
		gid += gidblock;
		float frac = ncount/(float)gid;
		printf("zero mbp : %ld : %ld ::: %g\n",ncount,gid,frac);
	}
	fclose(wp);
}

