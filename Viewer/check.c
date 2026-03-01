#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>


typedef struct HaloQ{
	size_t np;
	float x,y,z;
	float vx,vy,vz;
} HaloQ;

HaloQ halo;


int main( int argc, char ** argv){
	int i,j,k;
	int np;


	FILE *fp;

	fp = fopen("FoF_member_particle.01602","r");

	for(i=0;i<100;i++){
		fread(&halo,sizeof(HaloQ),1,fp);
		np = halo.np;
		printf("%d kdfkjdkfdk %ld %d\n",i,halo.np,np);
	}
	fclose(fp);
}
