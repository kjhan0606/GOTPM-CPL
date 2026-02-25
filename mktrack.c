#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include "Memory.h"
#include "pmheader.h"
#include "animate.h"


int mkviewer(char *, Viewer **);

Viewer *viewer;
int nstep,nframe;
char infile[190],outfile[190];
FILE *wp;

#define rad2deg (180.L/3.14159265354979L)


int MAIN_(int argc, char *argv[]){
	int i,j,k;
	if(Make_Total_Memory() == 0 ){
		fprintf(stderr,"Error opening memory\n");
		exit(99);
	}

	if(argc != 3){
		fprintf(stderr,"Error in the number of arg's. It needs 2 arg's\n");
		fprintf(stderr,"mktrack infile outfile\n");
		exit(99);
	}

	sprintf(infile,"%s",argv[1]);
	sprintf(outfile,"%s",argv[2]);

	viewer = (Viewer *) Malloc(sizeof(Viewer)*MAXFRAME,PPTR(viewer));

	nframe = mkviewer(infile,&viewer);

	wp = fopen(outfile,"w");

	for(i=0;i<nframe;i++){
		printf("%d %d |    %g %g %g : %g %g %g ::: %g %g %g || %g %g %g\n",
				viewer[i].frame,viewer[i].nstep,
				viewer[i].pos.x,viewer[i].pos.y,viewer[i].pos.z,
				viewer[i].E1.x,viewer[i].E1.y,viewer[i].E1.z,
				viewer[i].E2.x,viewer[i].E2.y, viewer[i].E2.z,
				viewer[i].E3.x,viewer[i].E3.y, viewer[i].E3.z
				);
		fprintf(wp,"%d %d %g %g %g %g %g %g %g %g %g %g %g %g\n",
				viewer[i].frame,viewer[i].nstep,
				viewer[i].pos.x,viewer[i].pos.y,viewer[i].pos.z,
				viewer[i].E3.x,viewer[i].E3.y,viewer[i].E3.z,
				viewer[i].E2.x,viewer[i].E2.y, viewer[i].E2.z,
				viewer[i].alpha,viewer[i].beta,viewer[i].gamma);
		/*
		{
			double a,b,c;
			a = viewer[i].E1.x*viewer[i].E3.x+viewer[i].E1.y*viewer[i].E3.y +viewer[i].E1.z*viewer[i].E3.z;
			b = viewer[i].E2.x*viewer[i].E3.x+viewer[i].E2.y*viewer[i].E3.y +viewer[i].E2.z*viewer[i].E3.z;
			c = viewer[i].E1.x*viewer[i].E2.x+viewer[i].E1.y*viewer[i].E2.y +viewer[i].E1.z*viewer[i].E2.z;
			printf(" %g %g %g\n",a,b,c);
		}
		*/
	}
	Free(viewer);

}
