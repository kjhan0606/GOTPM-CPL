#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>

#include "pmheader.h"
#include "params.h"

int read_viewers(char *paramfile,SimParameters *simpar){
	int isize,ncnt=0;
	int nviewer, mview = 0;
	char line[MAX_LINE_LENGTH];
	FILE *fp;
	fp  = fopen(paramfile,"r");
   	while(strcmp(P_Closing,fgets(line, MAX_LINE_LENGTH,fp)) != 0){ /* Loop until it reaches P_Closing */ 
		int ncnt = 0;
		if(line[0] != '#'){
			if(strstr(line,"define") != NULL ) {
				VIEWER2(sscanf,line, &, simpar,ncnt);
				if(ncnt == 1) {
					nviewer = simpar->nviewer;
					break;
				}
			}
		}
	}
	fclose(fp);
	if(mview != nviewer){
		fprintf(stderr,"Strange in the number of viewers\n");
	}
	return nviewer;
}
SimParameters  read_sim_parameter_file(FILE *);

void getsimparm_(char *paramfile,float *boxsize,float *hubble,float *npower, 
		float *omep, float *omepb, float*omeplam, float *wlam0,
		float *bias8,  float *rng, 
		float *amax, float *da, float *anow, int *iflag, char *outfile){ 
	SimParameters simpar; 
	char cpar[190];
	FILE *fp; 
	sprintf(cpar,"%s",paramfile);
	fp = fopen(cpar,"r");
   	simpar=read_sim_parameter_file(fp);
   	*boxsize = simpar.boxsize;
   	*hubble = simpar.hubble;
   	*npower = simpar.npow;
   	*omep = simpar.omep;
   	*omepb = simpar.omepb;
   	*omeplam = simpar.omeplam;
   	*bias8 = simpar.bias8;
   	*rng = simpar.nx;
   	*amax = simpar.amax;
   	*da = simpar.astep;
   	*anow = simpar.anow;
   	*iflag = simpar.powreadflag;
	*wlam0 = simpar.wlam0;
	if(simpar.powreadflag ==1) sprintf(outfile,"%s", simpar.powfilename);
	else if(simpar.powreadflag ==2) sprintf(outfile,"%s", simpar.inpapkfilename);
   	fclose(fp);
}

void getwsimparm_(char *paramfile,float *boxsize,float *hubble,float *npower, 
		float *omep, float *omepb, float*omeplam, char *DarkEnergyModel,
		float *wlam0,float *wlama, float *cs2_lam,
		float *bias8, float *As, float *rng, 
		float *amax, float *da, float *anow, int *iflag, char *outfile){ 
	SimParameters simpar; 
	char cpar[190];
	FILE *fp; 
	sprintf(cpar,"%s",paramfile);
	fp = fopen(cpar,"r");
   	simpar=read_sim_parameter_file(fp);
   	*boxsize = simpar.boxsize;
   	*hubble = simpar.hubble;
   	*npower = simpar.npow;
   	*omep = simpar.omep;
   	*omepb = simpar.omepb;
   	*omeplam = simpar.omeplam;
   	*bias8 = simpar.bias8;
   	*As = simpar.As;
   	*rng = simpar.nx;
   	*amax = simpar.amax;
   	*da = simpar.astep;
   	*anow = simpar.anow;
   	*iflag = simpar.powreadflag;
	*wlam0 = simpar.wlam0;
	*wlama = simpar.wlam1;
	*cs2_lam = simpar.CsDE;
//   	sprintf(outfile,"%s", simpar.powfilename);
   	sprintf(outfile,"%s", simpar.inpapkfilename);
	sprintf(DarkEnergyModel,"%s", simpar.DarkEnergyModel);
   	fclose(fp);
}

