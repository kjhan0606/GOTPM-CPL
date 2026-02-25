#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>

#include "Memory.h"
#include "pmheader.h"
#include "animate.h"
#define rad2deg (90./3.14159265)


Viewer *viewer;
int mode;
int mkviewer(char *, Viewer **);

int mkanimate(pmparticletype *pmparticles, int np){
	char viewfile[190];
	char infile[190];
	int nframe=0,i,j,k;
	FILE *fp,*fp1;
	MPI_Status status;


	viewer = (Viewer *)Malloc(sizeof(Viewer)*MAXFRAME,PPTR(viewer));

	if(simpar.myid ==0) { 
		if((fp = fopen(simpar.Viewerfile,"r")) != NULL) {
			printf("making animation from viewer file= %s\n",simpar.Viewerfile);
			while((fscanf(fp,"%s %d",viewfile,&mode))!= EOF){
				float HOpticalDepth;
				nframe = mkviewer(viewfile, &viewer);
				printf("making animation for viewfile= %s,  nframe= %d out= %s\n",
						viewfile,nframe,viewer[0].outfile);
				for(i=1;i<simpar.nid;i++){
					MPI_Send(&nframe,1,MPI_INT,i,i,MPI_COMM_WORLD);
					MPI_Send(viewer,nframe*sizeof(Viewer),MPI_BYTE,i,i,MPI_COMM_WORLD);
					MPI_Send(&mode,1,MPI_INT,i,i,MPI_COMM_WORLD);
				}
				ObsAnimate(pmparticles,np, viewer, nframe,simpar.stepcount,mode);
			}
			fclose(fp);
		}
		nframe = 0;
		for(i=1;i<simpar.nid;i++) MPI_Send(&nframe,1,MPI_INT,i,i,MPI_COMM_WORLD);
	}
	else {
		MPI_Recv(&nframe,1,MPI_INT,0,simpar.myid,MPI_COMM_WORLD,&status);
		while(nframe > 0) {
			MPI_Recv(viewer,nframe*sizeof(Viewer),MPI_BYTE,0,simpar.myid,MPI_COMM_WORLD,&status);
			MPI_Recv(&mode,1,MPI_INT,0,simpar.myid,MPI_COMM_WORLD,&status);
			ObsAnimate(pmparticles,np, viewer, nframe,simpar.stepcount,mode);
			MPI_Recv(&nframe,1,MPI_INT,0,simpar.myid,MPI_COMM_WORLD,&status);
		}
	}

	Free(viewer);
	return 0;
}
