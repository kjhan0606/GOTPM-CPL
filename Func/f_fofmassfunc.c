#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
typedef struct Halo{
	size_t nmem;
	float x,y,z;
	float vx,vy,vz;
} Halo;
#define MAX_HALO_NUM 90000000

#define MIN(A,B) ((A) > (B) ? (B):(A))
#define MAX(A,B) ((A) < (B) ? (B):(A))

int main(int argc, char **argv){
	int j,k;
	size_t i;
	char filename[80];
	size_t nummass;
	float minmass=1.E23,maxmass=-1.E23;
	float *massbin,*numbin, *err;
	float binstep,boxsize;
	Halo *halo,*ptr;
	FILE *fp,*wp;
	char *outfile;
	int nbin;
	float x,y,z,index,mpeak,peakindex;
	float size,hubble,npower,omep,omepl,bias,amax,astep,anow;
	float omepb,unitmass;
	int nx,nspace;

	if(argc != 6) {
		fprintf(stderr,
				"input : massfunc2 [filename] [# of bin] [min mass] [max mass] [outfile]\n");
		exit(0);
	}
	fp = fopen(argv[1],"r");
	nbin = atoi(argv[2]);
	float rminmass,rmaxmass;
	rminmass = atof(argv[3]);
	rmaxmass = atof(argv[4]);
	wp = fopen(argv[5],"w");
	printf("nbin=%d\n",nbin);
	fread(&size,sizeof(float),1,fp);
	fread(&hubble,sizeof(float),1,fp);
	fread(&npower,sizeof(float),1,fp);
	fread(&omep,sizeof(float),1,fp);
	fread(&omepb,sizeof(float),1,fp);
	fread(&omepl,sizeof(float),1,fp);
	fread(&bias,sizeof(float),1,fp);
	fread(&nx,sizeof(int),1,fp);
	fread(&nspace,sizeof(int),1,fp);
	fread(&amax,sizeof(int),1,fp);
	fread(&astep,sizeof(int),1,fp);
	fread(&anow,sizeof(int),1,fp);
	unitmass = 2.7755E11*omep*pow(size/nx*nspace,3);
	printf("size=%g hubble=%g omep=%g nx=%d nspace=%d\n",size,hubble,omep,nx,nspace);
	printf("particle mass=%g\n",unitmass);
	halo = (Halo *)malloc(sizeof(Halo)*MAX_HALO_NUM);
	nummass = fread(halo,sizeof(Halo),MAX_HALO_NUM,fp);
	printf("total number of halo =%ld\n",nummass);
	fclose(fp);
	ptr = halo;
	for(i=0;i<nummass;i++){
		if(ptr->nmem < 10 ) continue;
		minmass = MIN(minmass,ptr->nmem*unitmass);
		maxmass = MAX(maxmass,ptr->nmem*unitmass);
		ptr++;
	}
	printf("minum %g and maximum %g\n",minmass,maxmass);
	minmass = rminmass;
	maxmass = rmaxmass;
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
	for(i=0;i<nummass;i++){
		float tmpmass;
		int ntmp;
		if(halo[i].nmem < 10 ) continue;
		tmpmass = log10(halo[i].nmem*unitmass);
		ntmp = (tmpmass-minmass)/binstep;
		if(ntmp >=0 && ntmp < nbin) numbin[ntmp] += 1;
	}
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
