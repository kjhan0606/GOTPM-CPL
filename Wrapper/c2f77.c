#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include "pmheader.h"
#include "params.h"

SimParameters read_header(FILE *);

SimParameters simpar;


void intchecknp_(int *np, char *string, int string_len){
	char infile[19];
	int i;

	for(i=0;i<string_len;i++){
		infile[i] = string[i];
	}
	infile[i] = '\0';
	FILE *fp = fopen(infile,"r");
	simpar = read_header(fp);

	fclose(fp);
	*np = simpar.np;
}

void simparameters_(char *string, int string_len,float *omep,float *omepb, float *omeplam,
		float *hubble,float *red, float *nps, float *sigma8, float *bxsize){
	char infile[19];
	int i;

	for(i=0;i<string_len;i++){
		infile[i] = string[i];
	}
	infile[i] = '\0';
	FILE *fp = fopen(infile,"r");
	simpar = read_header(fp);
	fclose(fp);
	*omep = simpar.omep;
	*omepb = simpar.omepb;
	*omeplam = simpar.omeplam;
	*hubble = simpar.hubble;
	*red = simpar.amax /simpar.anow;
	*nps = simpar.npow;
	*sigma8 = 1./simpar.bias;
	*bxsize = simpar.boxsize;

}
void readparticles_(char *string, int string_len,float *x,float *y, float *z,
		float *vx,float *vy, float *vz, long long *index){
	char infile[19];
	int i;
	pmparticletype *bp;
	double r2kineticfact;

	for(i=0;i<string_len;i++){
		infile[i] = string[i];
	}
	infile[i] = '\0';
	FILE *fp = fopen(infile,"r");
	simpar = read_header(fp);

	double Hsub;
	float omepk =  1. - simpar.omep - simpar.omeplam;
	Hsub = sqrt(simpar.omep*pow(simpar.amax/simpar.anow,3.L)+simpar.omeplam+omepk*pow(simpar.amax/simpar.anow,2.L));
	r2kineticfact = simpar.boxsize/simpar.amax*simpar.anow*simpar.anow*100.L*Hsub;

	double rscale = simpar.boxsize /simpar.nx;

	bp = (pmparticletype*)malloc(sizeof(pmparticletype)*simpar.np);
	fread(bp,sizeof(pmparticletype),simpar.np,fp);
	for(i=0;i<simpar.np;i++){
		x[i] = XofP(bp+i) * rscale;
		y[i] = YofP(bp+i) * rscale;
		z[i] = ZofP(bp+i) * rscale;
		vx[i] = bp[i].vx * r2kineticfact;
		vy[i] = bp[i].vy * r2kineticfact;
		vz[i] = bp[i].vz * r2kineticfact;
		index[i] = bp[i].indx;
	}

	fclose(fp);

	free(bp);
}
