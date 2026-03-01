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
float boxsize, hubble, npower, omep, omeplam, bias,astep,amax,anow,omepb, wlam0, wlam1;
float wlam0, wlam1;
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
	fread(&wlam0, sizeof(float), 1,wp);
	fread(&wlam1, sizeof(float), 1,wp);
	fread(&bias, sizeof(float), 1,wp);
	fread(&astep, sizeof(float), 1,wp);
	fread(&msteps, sizeof(float), 1,wp);
	nsteps = (int*)malloc(sizeof(int)*msteps);
	reds = (float*)malloc(sizeof(float)*msteps);
	fread(nsteps, sizeof(int), msteps,wp);
	fread(reds, sizeof(float), msteps,wp);
	fclose(wp);
	wp = fopen("MergingTree.dat","w");
	fwrite(&boxsize, sizeof(float), 1,wp);
	fwrite(&hubble, sizeof(float), 1,wp);
	fwrite(&npower, sizeof(float), 1,wp);
	fwrite(&omep, sizeof(float), 1,wp);
	fwrite(&omepb, sizeof(float), 1,wp);
	fwrite(&omeplam, sizeof(float), 1,wp);
	fwrite(&wlam0, sizeof(float), 1,wp);
	fwrite(&wlam1, sizeof(float), 1,wp);
	fwrite(&bias, sizeof(float), 1,wp);
	fwrite(&astep, sizeof(float), 1,wp);
	fwrite(&msteps, sizeof(float), 1,wp);
	fwrite(nsteps, sizeof(int), msteps,wp);
	fwrite(reds, sizeof(float), msteps,wp);

	fflush(wp);

	for(k=0; ; k+= npmax){
		npmax = 0;
		SavedHaloType *bp = (SavedHaloType*)calloc(gidblock*msteps,sizeof(SavedHaloType));
		for(i=0;i<msteps*gidblock;i++){
			bp[i].mbp = 0;
		}
		for(i=0;i<msteps;i++){
			char infile[190]; sprintf(infile,"MergingTree.dat.%.5d",(int)i);
			fp = fopen(infile,"r");
			fseek(fp, k*sizeof(SavedHaloType),SEEK_SET);
			np = fread(tmp,sizeof(SavedHaloType),gidblock,fp);
			printf("%d : %s is being read with np = %d\n",(int)k, infile, (int)np);
			fflush(stdout);
			fclose(fp);
			if(np>0){
				for(j=0;j<np;j++){
					bp[i + j*msteps] = tmp[j];
				}
			}
			npmax = MAX(npmax,np);
		}
		if(npmax ==0) break;
		fwrite(bp,sizeof(SavedHaloType),npmax*msteps,wp);
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
					if(bp[i+j*msteps].mbp > 0  && mbp != bp[i+j*msteps].mbp){
						fprintf(stderr,"Error here %ld %ld  : %ld %ld\n",i,j,mbp,bp[i+j*msteps].mbp);
						exit(999);
					}
				}
			}
		}
		free(bp);
	}
	fclose(wp);
}

