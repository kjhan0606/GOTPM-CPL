#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include "Memory.h"
#include "pmheader.h"
#include "animate.h"


int mkviewer(char *, int ,Viewer **);

Viewer *viewer;
int nstep;
char infile[190];


int MAIN_(int argc, char *argv[]){
	int i,j,k;
	if(Make_Total_Memory() == 0 ){
		fprintf(stderr,"Error opening memory\n");
		exit(99);
	}


	nstep = 400;
	sprintf(infile,"Viewer1.dat");

	mkviewer(infile,nstep,&viewer);

	for(i=0;i<nstep;i++){
		printf("%d %d |    %g %g %g : %g %g %g ::: %g\n",
				viewer[i].frame,viewer[i].nstep,
				viewer[i].pos.x,viewer[i].pos.y,viewer[i].pos.z,
				viewer[i].view.x,viewer[i].view.y,viewer[i].view.z,viewer[i].rot);
	}

}
