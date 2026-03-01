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

	if(argc != 2) {
		fprintf(stderr,
				"input : massfunc2 [filename] \n");
		exit(0);
	}
	fp = fopen(argv[1],"r");
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
	printf("%g %g %g \n",halo[0].vx,halo[0].vy,halo[0].vz);
	return 0;
}
