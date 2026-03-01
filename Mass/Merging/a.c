#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<string.h>
#include<mpi.h>
#include "pqsort.h"


typedef struct halo{
	float x,y,z;
	float vx,vy,vz;
	int hid;
} Halo;

Halo halo[100];


PQSORT(base,Halo,int,hid,nmem,Comp,min,max,Comm)

int comp(const void *a, const void *b){
	Halo *aa = (Halo*)a;
	Halo *bb = (Halo*)b;
	if(aa->hid > bb->hid) return 1;
	else if(aa->hid < bb->hid) return -1;
	return 0;
}




int main(int argc, char *argv[]){
	return 1;
}
