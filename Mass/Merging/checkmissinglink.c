#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include "merger.h"

float minmass,maxmass;
#define nbin 32


long bin[nbin],mbin[nbin];



int main(int argc, char **argv){
	IDTYPE nowmbp,i,j;
	float redshift;
	TrHalo *nowmbplist;


	for(i=0;i<nbin;i++) {
		bin[i] =0;
		mbin[i] =0;
	}
	minmass = log10(2.7e11L);
	maxmass = log10(5.e15L);
	float mstep = (maxmass-minmass)/nbin;



	FILE *fp;

	fp = fopen(argv[1],"r");

	fread(&redshift, sizeof(float),1,fp);
	fread(&nowmbp, sizeof(IDTYPE),1,fp);

	nowmbplist = (TrHalo*)malloc(sizeof(TrHalo)*nowmbp);
	fread(nowmbplist,sizeof(TrHalo),nowmbp,fp);
	fclose(fp);
	IDTYPE count = 0;
	for(i=0;i<nowmbp;i++){
		{
			int k = (log10(nowmbplist[i].mass)-minmass)/mstep;
			if(k>=0 && k < nbin) mbin[k] ++;
		}
		if(nowmbplist[i].upaid >0) count ++;
		else {
			int k = (log10(nowmbplist[i].mass)-minmass)/mstep;
			if(k>=0 && k < nbin) bin[k] ++;
		}
	}
	printf("%s is being checked\n",argv[1]);
	for(i=0;i<nbin;i++){
		float mass = pow(10., minmass + (i+0.5) * mstep);
		float frac = (float)bin[i]/(float)mbin[i];
		printf("mass = %g count= %ld in tcount= %ld : %g\n",mass, bin[i],mbin[i], frac);
	}
	float frac = (float) count / (float) nowmbp;
	
	printf("%s  %ld : %ld   ::::: %g\n",argv[1],nowmbp,count,frac);
}
