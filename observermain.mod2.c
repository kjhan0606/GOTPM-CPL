#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "pmheader.h"
#include "pm_common.h"
#include "lightcone.boss.h"


#define NOBS 64
double xobs[NOBS],yobs[NOBS],zobs[NOBS];



void GottPMeshDump(float *,int ,int , double *, double *, double *,int , float );
void GottGrid(float *fftwmesh){
	int i,j,k;
	float xshift,yshift,zshift,shift;
	float redi;
	int iobs;


	shift = simpar.nx/4./2.;
	iobs = 0;
	for(k=0;k<4;k++) {
		for(j=0;j<4;j++) {
			for(i=0;i<4;i++){ 
				xobs[iobs] = (float)simpar.nx * (float)i/4. + shift;
				yobs[iobs] = (float)simpar.ny * (float)j/4. + shift;
				zobs[iobs] = (float)simpar.nz * (float)k/4. + shift;
				iobs ++;
			}
		}
	}
	redi = 0.8;

	TIMER_START(62);
	GottPMeshDump(fftwmesh,simpar.stepcount,0,xobs,yobs, zobs,iobs,redi);
	TIMER_STOP(62);
	if(simpar.myid==0) fprintf(stdout,"BOSS Grid Survey CPU= %f \n", ELAPSED_TIME(62));
}




float maxd;

void observerextmain(){
	long i,j,k;
	int iobs;
	float xshift,yshift,zshift,shift;
	float redi;
	int flagsync = 0;

	maxd= simpar.nspace*0.1;
	
	/* This is for the d = 1.2 Mpc/h
	iobs = 0;
	shift = simpar.nx/4.;
	for(k=0;k<4;k++) for(j=0;j<4;j++){
		xshift = ((j+k)%2)*shift;
		for(i=0;i<2;i++){
			xobs[iobs] = i *simpar.nx/2. + xshift;
			yobs[iobs] = j *simpar.ny/4.;
			zobs[iobs] = k *simpar.nz/4.;
			iobs ++;
		}
	}
	redi = 0.8;
	*/
	shift = simpar.nx/4./2.;
	iobs = 0;
	for(k=0;k<4;k++) {
		for(j=0;j<4;j++) {
			for(i=0;i<4;i++){ 
				xobs[iobs] = (float)simpar.nx * (float)i/4. + shift;
				yobs[iobs] = (float)simpar.ny * (float)j/4. + shift;
				zobs[iobs] = (float)simpar.nz * (float)k/4. + shift;
				iobs ++;
			}
		}
	}
	redi = 0.8;


	TIMER_START(62);

	BossESLightConeData(pmparticles,simpar.np,simpar.stepcount,simpar.anow,
			simpar.astep,simpar.astep,maxd,flagsync,0,xobs,yobs, zobs,iobs,redi);
	TIMER_STOP(62);
	if(simpar.myid==0) fprintf(stdout,"BOSS ES Survey CPU= %f \n", ELAPSED_TIME(62));
	TIMER_START(62);
	DiskESLightConeData(pmparticles,simpar.np,simpar.stepcount,simpar.anow,
			simpar.astep,simpar.astep,maxd,flagsync);
	TIMER_STOP(62);
	if(simpar.myid==0) fprintf(stdout,"DISK ES Survey CPU= %f \n", ELAPSED_TIME(62));
}
void observersavemain(){
	float redi;
	int iobs = NOBS;
	redi = 0.8;
	TIMER_START(62);
	BossSavingLightConeData(pmparticles,simpar.np, simpar.stepcount,
			simpar.amax,simpar.anow,
			simpar.astep,simpar.astep,0,xobs,yobs, zobs,iobs,redi);
	TIMER_STOP(62);
	if(simpar.myid==0) fprintf(stdout,"BOSS S Survey CPU= %f \n", ELAPSED_TIME(62));
	TIMER_START(62);
	DiskSavingLightConeData(pmparticles,simpar.np,simpar.stepcount,simpar.amax,simpar.anow,
			simpar.astep,simpar.astep);
	TIMER_STOP(62);
	if(simpar.myid==0) fprintf(stdout,"DISK S Survey CPU= %f \n", ELAPSED_TIME(62));

}
