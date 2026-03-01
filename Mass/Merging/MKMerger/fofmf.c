#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
typedef struct Halo{
	size_t nmem;
	double x,y,z;
	float vx,vy,vz;
} Halo;
#define MAX_HALO_NUM 100000000

#define MIN(A,B) ((A) > (B) ? (B):(A))
#define MAX(A,B) ((A) < (B) ? (B):(A))

int main(int argc, char **argv){
	int j,k;
	size_t i;
	char filename[80];
	size_t nummass;
	float minmass=1.E23,maxmass=-1.E23;
	double *massbin,*numbin, *err;
	float binstep,boxsize;
	float hmass;
	FILE *fp,*wp;
	char *outfile;
	int nbin;
	float x,y,z,index,mpeak,peakindex;
	float size,hubble,npower,omep,omepl,bias,amax,astep,anow;
	float omepb,unitmass;
	int nx,nspace;

	if(argc != 4) {
		fprintf(stderr,
				"input : massfunc2 [filename] [# of bin] [outfile]\n");
		exit(0);
	}
	fp = fopen(argv[1],"r");
	nbin = atoi(argv[2]);
	wp = fopen(argv[3],"w");
	printf("nbin=%d\n",nbin);
	massbin = (double *) malloc(sizeof(double)*(nbin+1));
	numbin = (double *) malloc(sizeof(double)*(nbin+1));
	err = (double *) malloc(sizeof(double)*(nbin+1));


	int iflag = 0;
	long long ntot = 0;

	float rminmass,rmaxmass;

	minmass = 2.7E11;
	maxmass = 5.E15;
	rminmass = log10(minmass);
	rmaxmass = log10(maxmass);
	binstep = (rmaxmass-rminmass)/(float) nbin;
	for(i=0;i<=nbin;i++){
		massbin[i] = rminmass + (rmaxmass-rminmass)/(double) nbin *i;
		numbin[i] = 0;
	}



	while((nummass = fread(&hmass,sizeof(Halo),MAX_HALO_NUM,fp))>0){
		float tmpmass;
		int ntmp;
		tmpmass = log10(hmass);
		ntmp = (tmpmass-rminmass)/binstep;
		if(ntmp >=0 && ntmp < nbin) numbin[ntmp] += 1;
		ntot += nummass;
	}


	printf("minum %g and maximum %g\n",minmass,maxmass);
	printf("total number of halo =%ld\n",ntot);
	fclose(fp);


	/* centering of bin's to take account of the [dex] */
	for(i=0;i<nbin;i++) {
		err[i] = sqrt(numbin[i])/pow(size,3);
		numbin[i] = numbin[i]/pow(size,3);
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
		printf("%lg %lg %lg %lg %lg\n",massbin[i],
				numbin[i]/binsize*massbin[i]*log(10.),
				(numbin[i]+err[i])/binsize*massbin[i]*log(10.),
				(numbin[i]-err[i])/binsize*massbin[i]*log(10.),
				cmass);
		fprintf(wp,"%lg %lg %lg %lg %lg\n",massbin[i],
				numbin[i]/binsize*massbin[i]*log(10.),
				(numbin[i]+err[i])/binsize*massbin[i]*log(10.),
				(numbin[i]-err[i])/binsize*massbin[i]*log(10.),
				cmass);
	}
	fclose(wp);
	return 0;
}
