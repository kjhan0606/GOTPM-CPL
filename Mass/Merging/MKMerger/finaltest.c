/* icc -o finaltest finaltest.c -g -lm -DINDEX -DXYZDBL
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
float boxsize, wlam, hubble, npower, omep, omeplam, bias,astep,amax,anow,omepb;
int nx, nspace;
int msteps,*nsteps;


int main (int argc, char **argv){
	FILE *fp,*wp;
	IDTYPE i,j,k,npmax;
	IDTYPE gidblock = 500000,gidshift;
	SavedHaloType *tmp = (SavedHaloType*)malloc(sizeof(SavedHaloType)*gidblock);
	int np,mp,rflag;



	gidshift = 0;
	npmax = gidblock;

	wp = fopen("wMergingTree.dat","r");
	fread(&boxsize, sizeof(float), 1,wp);
	fread(&hubble, sizeof(float), 1,wp);
	fread(&npower, sizeof(float), 1,wp);
	fread(&omep, sizeof(float), 1,wp);
	fread(&omepb, sizeof(float), 1,wp);
	fread(&omeplam, sizeof(float), 1,wp);
	fread(&wlam, sizeof(float), 1,wp);
	fread(&bias, sizeof(float), 1,wp);
	fread(&astep, sizeof(float), 1,wp);
	fread(&msteps, sizeof(float), 1,wp);
	nsteps = (int*)malloc(sizeof(int)*msteps);
	reds = (float*)malloc(sizeof(float)*msteps);
	fread(nsteps, sizeof(int), msteps,wp);
	fread(reds, sizeof(float), msteps,wp);

	IDTYPE ncount = 0,gid=0;

	for(k=0; ; k+= npmax){
		npmax = 0;
		SavedHaloType *bp = (SavedHaloType*)calloc(gidblock*msteps,sizeof(SavedHaloType));
		for(i=0;i<msteps*gidblock;i++){
			bp[i].mbp = 0;
		}
		{
			char infile[190]; sprintf(infile,"MergingTree.dat.%.5d",(int)i);
			np = fread(bp,sizeof(SavedHaloType),gidblock*msteps,wp)/msteps;
			npmax = MAX(npmax,np);
		}
		if(npmax ==0) break;
		for(i=0;i<npmax;i++){
			int ii=0;
			for(j=0;j<msteps-1;j++){
				if(bp[j + i*msteps].mbp <=0 && bp[j+1 + i*msteps].mbp >0) {
					ii ++;
				}
				else if(bp[j + i*msteps].mbp >0 && bp[j+1 + i*msteps].mbp <=0) {
					ii++;
				}
			}
			if(ii>=3){
				fprintf(stderr,"Error strange \n");
			}
		}
		{
			for(j=0;j<npmax;j++){
				IDTYPE mbp;
				IDTYPE ii = 0;
				while((mbp=bp[j*msteps + ii].mbp )== 0) {
					if(ii == msteps) break;
					ii++; 
				}
				if(ii== msteps) continue;
				for(i=ii+1;i<msteps;i++){
					if(bp[i+j*msteps].mbp !=0  && mbp != bp[i+j*msteps].mbp){
						/*
						fprintf(stderr,"Error here %ld %ld  : %ld %ld\n",i,j,mbp,bp[i+j*msteps].mbp);
						exit(999);
						*/
					}
					else if(bp[i+j*msteps].mbp ==0){
						ncount ++;
						break;
					}
				}
			}
			gid += gidblock;
			float frac = (float)ncount/(float) gid;
			printf("Missing halos %ld : %ld :::: %g\n",ncount,gid,frac);
		}
		free(bp);
	}
	fclose(wp);
}

