#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>


#define NMAX (1<<14)


typedef struct Particle{
	float x,y,z,vx,vy,vz;
	long long indx;
} Particle;

Particle *a,*bp;


int main(int argc, char **argv){
	int tnp,np,size;
	FILE *fp;

	fp = fopen("tmp.dat","r");

	bp = (Particle*)malloc(sizeof(Particle)*NMAX);
	size = NMAX;
	a = bp;
	tnp = 0;
	while((np=fread(a,sizeof(Particle),NMAX,fp))>0){
		tnp += np;
		a += np;
		if(size-tnp < NMAX){
			size += 2*NMAX;
			bp = (Particle*)realloc(bp,sizeof(Particle)*size);
			a = bp + tnp;
		}
	}
	printf("Total np=%d are read\n",tnp);
	free(bp);
}
