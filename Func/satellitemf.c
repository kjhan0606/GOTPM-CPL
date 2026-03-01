#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "sed.h"
#define MAX_HALO_NUM 20000000
int sorttype(const void *a,const void *b){
	Halo *aa,*bb;
	aa = (Halo *)a;
	bb = (Halo *)b;
	if(aa->indx < bb->indx) return -1;
	else if(aa->indx > bb->indx) return +1;
	else {
		if(aa->mass < bb->mass) return +1;
		else if(aa->mass > bb->mass)return -1;
		else return 0;
	}
}

/*
#define MIN(A,B) ((A) > (B) ? (B):(A))
#define MAX(A,B) ((A) < (B) ? (B):(A))
*/

int main(int argc, char **argv){
	int i,j,k;
	char filename[80];
	float mass[MAX_HALO_NUM];
	int np;
	float minmass=1.E23,maxmass=-1.E23;
	float *massbin,*numbin, *err;
	float binstep,boxsize;
	Halo *halo;
	FILE *fp,*wp;
	char *outfile;
	char onechar;
	int nbin;
	float hubble;
	float x,y,z,index,mpeak,peakindex;

	if(argc != 6) {
		fprintf(stderr,
				"input : massfunc2 [filename] [# of bin] [boxsize] [h] [outfile]\n");
		exit(0);
	}
	fp = fopen(argv[1],"r");
	nbin = atoi(argv[2]);
	boxsize = atof(argv[3]);
	hubble = atof(argv[4]);
	outfile = argv[5];
	halo = (Halo *)malloc(sizeof(Halo)*MAX_HALO_NUM);
	i = 0;
	/*
	while(fscanf(fp,"%f %f %f %f %d %d %d\n",mass+i,&x,&y,&z,
				&index,&mpeak,&peakindex) != EOF){
				*/
	while(fscanf(fp,"%f %f %f %f %d %d %d ",&(halo[i].mass),&(halo[i].x),&(halo[i].y),&(halo[i].z),
		    &(halo[i].indx),&(halo[i].mpeak),&(halo[i].peakindx))!= EOF){
		onechar = fgetc(fp);
		while(onechar != '\n' && onechar != EOF){
			onechar = fgetc(fp);
		}

		/* multiply by hubble  because massfunc1.c outputs 
		 * mass including hubble constant */
		halo[i].mass *= hubble;
		minmass = MIN(minmass,halo[i].mass);
		maxmass = MAX(maxmass,halo[i].mass);
		/*
		*(mass+i) = *(mass+i)*hubble;
		minmass = MIN(minmass,*(mass+i));
		maxmass = MAX(maxmass,*(mass+i));
		*/
		i++;
	}
	np = i;
	printf("Total halo np= %d\n",np);
	fclose(fp);
	{
		int mp,iflag;
		qsort(halo,np,sizeof(Halo),sorttype);
		printf("End of qsort\n");
		mp = 0;
		for(i=0;i<np;){
			for(j=i+1;j<np;j++){
				if(halo[j].indx != halo[i].indx){
					break;
				}
			}
			for(k=i+1;k<j;k++){
				halo[mp] = halo[k];
				mp++;
			}
			i = j;
		}
		np = mp;
	}
	massbin = (float *) malloc(sizeof(float)*(nbin+1));
	numbin = (float *) malloc(sizeof(float)*(nbin+1));
	err = (float *) malloc(sizeof(float)*(nbin+1));
	minmass = log10(minmass);
	maxmass = log10(maxmass);
	for(i=0;i<=nbin;i++){
		massbin[i] = minmass + (maxmass-minmass)/(float) nbin *i;
		numbin[i] = 0;
	}
	binstep = (maxmass-minmass)/(float) nbin;
	for(i=0;i<np;i++){
		float tmpmass;
		int ntmp;
		tmpmass = log10(halo[i].mass);
		ntmp = (tmpmass-minmass)/binstep;
		numbin[ntmp] += 1;
	}
	/* centering of bin's to take account of the [dex] */
	wp = fopen(outfile,"w");
	for(i=0;i<nbin;i++) {
		err[i] = sqrt(numbin[i])/pow(boxsize,3);
		numbin[i] = numbin[i]/pow(boxsize,3);
	}
	for(i=0;i<nbin;i++){
		float binsize;
		float cmass;
		float sigma,jenkins;
		binsize = pow(10,massbin[i+1]) - pow(10,massbin[i]);
		massbin[i] = (pow(10,massbin[i]) + pow(10,massbin[i+1]))*0.5;
		cmass = 0;
		for(j=i;j<nbin;j++){
			cmass += numbin[j];
		}
		printf("%g %g %g %g %g\n",massbin[i],
				numbin[i]/binsize*massbin[i]*log(10.),
				(numbin[i]+err[i])/binsize*massbin[i]*log(10.),
				(numbin[i]-err[i])/binsize*massbin[i]*log(10.),
				cmass);
		fprintf(wp,"%g %g %g %g %g\n",massbin[i],
				numbin[i]/binsize*massbin[i]*log(10.),
				(numbin[i]+err[i])/binsize*massbin[i]*log(10.),
				(numbin[i]-err[i])/binsize*massbin[i]*log(10.),
				cmass);
	}
	fclose(wp);
	return 0;
}
