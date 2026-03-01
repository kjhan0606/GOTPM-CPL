#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>


enum boolean {YES=1,NO=0};





float searchimg(float *img, int *count, int nx,int ny, int icx,int icy, float thres){
	int i,j;
	float icount=0;
	for(j=icy-1;j<=icy+1;j++) {
		if(j<0 || j >= ny ) continue;
		if(img[icx+nx*j] >= thres && count[j] == NO) {
			count[j] = YES;
			icount += img[icx+nx*j];
			icount += searchimg(img, count, nx,ny,icx,j,thres);
		}
	}
	return icount;
}



int main(int argc, char **argv)
{
	int i,j,k;
	FILE *fp;
	FILE *wp;
	int nx,ny;
	float *img,*wimg;
	int *count;

	fp = fopen(argv[1],"r");
	wp = fopen(argv[2],"w");

	fread(&nx,sizeof(int),1,fp);
	fread(&ny,sizeof(int),1,fp);
	img = (float*)malloc(sizeof(float)*nx*ny);
	wimg = (float*)malloc(sizeof(float)*nx*ny);
	count = (int*)malloc(sizeof(int)*ny);
	fread(img,sizeof(float),nx*ny,fp);

	for(i=0;i<nx*ny;i++) wimg[i] = 0.;



	float maxval,minval,totval=0;
	maxval = -1e20;
	minval =  1e20;
	int icy[nx];

	for(i=0;i<nx;i++){
		maxval = -1e20;
		totval = 0;
		for(j=0;j<ny;j++){ 
			if(img[i+nx*j]> maxval ) {
				icy[i] = j; 
				maxval = img[i+nx*j];
			}
			totval += img[i+nx*j];
		}
		for(k=0;k<256;k++){
			float thres = maxval * k/256.;
			for(j=0;j<ny;j++) count[j] = NO;
			float icount = searchimg(img,count, nx,ny,i,icy[i],thres);
			float percentage = (1.-icount/totval);
			for(j=0;j<ny;j++){
				if(count[j] == YES) wimg[i+nx*j] = percentage;
			}
		}
		/*
		for(j=0;j<ny;j++){
			wimg[i+nx*j] = img[i+nx*j]/totval;
		}
		*/
	}
	fwrite(&nx,sizeof(int),1,wp);
	fwrite(&ny,sizeof(int),1,wp);
	fwrite(wimg,sizeof(float),nx*ny,wp);
	fclose(wp);
}
