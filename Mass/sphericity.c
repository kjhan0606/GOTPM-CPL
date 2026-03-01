#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include "merger.h"





HaloProp halop[1000000];
int nx=128;
int ny=128;
double *img;
float *fimg;



int main(int argc, char **argv)
{
	size_t np;
	FILE *fp;
	FILE *wp;


	char infile[190],outfile[100];
	if(argc !=6) {
		printf("bin2ascii.exe infile outfile minmass maxmass hubble\n");
		exit(99);
	}
	sprintf(infile,"%s",argv[1]);
	sprintf(outfile,"%s",argv[2]);
	float minmass = atof(argv[3]);
	float maxmass = atof(argv[4]);
	float hubble = atof(argv[5]);
	int iaxis = atoi(argv[6]);

	img = (double*)malloc(sizeof(double)*nx*ny);
	fimg = (float*)malloc(sizeof(float)*nx*ny);
	int i,j;
	for(i=0;i<nx*ny;i++){
		img[i] = fimg[i] = 0;
	}
	long long ntot = 0;
	fp = fopen(infile,"r");

	while((np = fread(halop,sizeof(HaloProp),1000000,fp))>0){
		for(i=0;i<np;i++){
			halop[i].mass *= hubble;
			if(halop[i].mass >=minmass && halop[i].mass < maxmass){
				int ix,iy;
				ix = (log10(halop[i].mass)-log10(minmass))*nx / (log10(maxmass)-log10(minmass));
				if(ix<0 || ix >=nx)continue;
				double aa = sqrt(halop[i].q * halop[i].s);
				iy = (aa)*ny;
				if(iy<0 || iy >=ny)continue;
				img[ix+nx*iy] += 1.L;
				ntot ++;
			}
		}
	}
	for(i=0;i<nx*ny;i++) fimg[i] = img[i];
	wp = fopen(outfile,"w");
	fwrite(&nx,sizeof(int),1,wp);
	fwrite(&ny,sizeof(int),1,wp);
	fwrite(fimg,sizeof(float),nx*ny,wp);
	fclose(wp);
}
