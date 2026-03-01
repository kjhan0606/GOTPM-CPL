#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>

#define indxtype long
typedef struct Halo{
	int nmem;
	/*
	float x,y,z;
	float vx,vy,vz;
	*/
} Halo;
#define MAX_HALO_NUM 50000000

#define MIN(A,B) ((A) > (B) ? (B):(A))
#define MAX(A,B) ((A) < (B) ? (B):(A))
float mscale;

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

	if(argc != 9) {
		fprintf(stderr,
				"input : massfunc2 [size] [omep] [nx]  [nspace] [hubble] [infilename] [# of bin] [outfile]\n");
		exit(0);
	}
	size = atof(argv[1]);
	omep = atof(argv[2]);
	nx = atoi(argv[3]);
	nspace = atoi(argv[4]);
	hubble = atof(argv[5]);
	fp = fopen(argv[6],"r");
	nbin = atoi(argv[7]);
	wp = fopen(argv[8],"w");
	printf("nbin=%d\n",nbin);
	halo = (Halo *)malloc(sizeof(Halo)*MAX_HALO_NUM);
	unitmass = 2.7755E11*omep*pow(size/nx*nspace,3);
	printf("size=%g hubble=%g omep=%g nx=%d nspace=%d\n",size,hubble,omep,nx,nspace);
	printf("particle mass=%g\n",unitmass);
	{
		int mpeak,np,nm;
		float *r,*vr;
		indxtype *indx;
		nm = 0;
		while(fread(&mpeak,sizeof(int),1,fp)==1){
			for(i=0;i<mpeak;i++){
				fread(&np,sizeof(int),1,fp);
				if(np > 0){
					r = (float*)malloc(sizeof(float)*np*3);
					vr = (float*)malloc(sizeof(float)*np*3);
					indx = (indxtype*)malloc(sizeof(indxtype)*np);
					fread(r,sizeof(float),3*np,fp);
					fread(vr,sizeof(float),3*np,fp);
					fread(indx,sizeof(indxtype),np,fp);
					halo[nm].nmem = np;
					free(indx);
					free(r);
					free(vr);
					nm ++;
				}
			}
		}
		nummass = nm;
	}
	printf("total number of halo =%ld\n",nummass);
	fclose(fp);
	ptr = halo;
	minmass = 1.e30;
	maxmass = -1.e30;
	for(i=0;i<nummass;i++){
		minmass = MIN(minmass,ptr->nmem*unitmass);
		maxmass = MAX(maxmass,ptr->nmem*unitmass);
		ptr++;
	}
	printf("minum %g and maximum %g\n",minmass,maxmass);
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
		tmpmass = log10(halo[i].nmem*unitmass);
		ntmp = (tmpmass-minmass)/binstep;
		numbin[ntmp] += 1;
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
