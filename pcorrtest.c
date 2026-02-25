#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>


extern void Fgetpcorr(int *, float*, float *, float *, float *, float *, float *, char *, float *,
		float *, float *, float *, float *, float *, float *, float *, int *, float *, float *);


int main(int argc, char**argv){
	float boxsize, hubble, npower, omep, omeplam, omepb,w0,wa,cs2_lam, bias, As, rng, amax, anow;
	int nkh;
	char DarkEnergyModel[20];

	float *factpoisson, *kh;


	nkh = 1024;
	factpoisson = (float*)malloc(sizeof(float)*nkh);
	kh = (float*)malloc(sizeof(float)*nkh);

	boxsize = 1024;
	hubble = 0.72;
	npower = 0.96;
	omep = 0.26;
	omepb = 0.044;
	omeplam = 0.74;
	sprintf(DarkEnergyModel,"PPF");
	w0 = -1;
	wa = 0.1;
	bias = 1.26;
	As = 2.154E-9;
	rng = 1024;
	amax = 48;
	anow = 1;
	cs2_lam = 1;
	int nx = 1024;

	Fgetpcorr(&nx,&boxsize, &hubble, &npower, &omep,&omepb, &omeplam, DarkEnergyModel, &w0,
			&wa, &cs2_lam, &bias, &As,&rng, &amax, &anow, &nkh, kh, factpoisson);
	int i;
	for(i=0;i<nkh;i++){
		printf("%d %g %g\n", i, kh[i], factpoisson[i]);
	}

	return 0;

}
