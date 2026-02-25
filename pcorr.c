#if DEMODEL == 1
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>

#include "pmheader.h"


#define mkh 256

static float pcorr[mkh], pcork[mkh], rkh[mkh];
static float rscale, kscale;
static float kmax;
static float rmax;


float interppcorr(float dist){
	float r = dist*rscale;
	int ir = floor(r);
	if(ir < 255){
		return ((pcorr[ir+1]-pcorr[ir]) * (r-ir) + pcorr[ir]);
	}
	else {
		return pcorr[mkh-1];
	}
}

float interppcork(float k){
	float rk = k*kscale;
	int irk = floor(rk);
	if(irk < 255){
		return ((pcork[irk+1]-pcork[irk]) * (rk-irk) + pcork[irk]);
	}
	else {
		return pcork[mkh-1];
	}
}


extern void Fgetpcorr(int *, float*, float *, float *, float *, float *, float *, char *, float *,
        float *, float *, float *, float *, float *, float *, int *, float *, float *);


void CalPcorK(int nkh, float *kh, float *pfactor, float boxsize, float nx){
	int i,j,k;

	kmax = 2*M_PI * sqrt(3.)*nx/boxsize;

	kscale = mkh/(sqrt(3.L)*nx); /* := mkh/(kmax) * 2PI /(boxsize) */


	rkh[0] = 0;
	pcork[0] = 1;


	for(i=1;i<mkh;i++){
		rkh[i] = i * kmax/mkh;
		for(j=0;j<nkh-1;j++){
			if(rkh[i] >= kh[j]  && rkh[i] < kh[j+1]){
				pcork[i] = pfactor[j] + (pfactor[j+1]-pfactor[j])/(kh[j+1]-kh[j])*(rkh[i]-kh[j]);
				break;
			}
		}
		if(j==nkh-1) pcork[i] = pfactor[nkh-1];
		printf("k-space %d %g %g\n", j, rkh[i], pcork[i]);
	}
}

void CalPcorR(int nkh, float *kh, float *pfactor, float boxsize, float nx){
	int i,j,k;
	float rr;

	rmax = 10;
	rscale = mkh/rmax;
	for(i=1;i<mkh;i++){
		rr = i / rscale *(boxsize/nx);
		float cr = 0;
		for(j=0;j<nkh;j++){
		    float kR = rr*(rkh[j]+rkh[j+1])*0.5L;
			cr += 1./(2*M_PI*M_PI)*sin(kR)*rkh[j]/rr *(rkh[j+1]-rkh[j]) * (pfactor[j]-1);
		}
		pcorr[i] = cr + 1;
		printf("r-space %g %g\n", rr, pcorr[i]);
	}
}


void InitializePoissonCorrection(){
	float rng;
    int nkh;
    char DarkEnergyModel[20];

    float *factpoisson, *kh;

	nkh = 1024;
	factpoisson = (float*)malloc(sizeof(float)*nkh);
    kh = (float*)malloc(sizeof(float)*nkh);

	rng = simpar.nx;
	Fgetpcorr(&simpar.nx, &simpar.boxsize, &simpar.hubble, &simpar.npow, &simpar.omep, &simpar.omepb, 
			&simpar.omeplam, simpar.DarkEnergyModel, &simpar.wlam0, &simpar.wlam1, &simpar.CsDE,
			 &simpar.As, &rng, &simpar.amax, &simpar.anow, &nkh, kh, factpoisson);

	int i;
	for(i=0;i<nkh;i++){
		printf("%d : %g %g\n", i, kh[i], factpoisson[i]);
	}


	CalPcorK(nkh, kh, factpoisson, simpar.boxsize, simpar.nx);
	CalPcorR(nkh, kh, factpoisson, simpar.boxsize, simpar.nx);

	free(kh);
	free(factpoisson);


}
#endif
