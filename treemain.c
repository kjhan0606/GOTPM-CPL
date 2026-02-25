#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<math.h>

#include "pmheader.h"
#include "pm_common.h"
#include "flow.h"


void tree2pmtype(long);
void pm2treetype(long);
static float treetime;

void TreeAndPreFoFmain(int iflagPreFoF){
	long nplong;
	int myid;
	float treecorrection(treeparticletype *, long);

	TIMER_START(59);
	nplong = np;
	myid = simpar.myid;
	if(simpar.stepcount < 10) simpar.theta = 0.3;
	else if(simpar.stepcount < 20) simpar.theta = 0.35;
	else if(simpar.stepcount < 60) simpar.theta = 0.4;
	else simpar.theta = 0.45;
	treeparticles = (treeparticletype *) pmparticles;
	pm2treetype(nplong);
	treetime=treecorrection(treeparticles,nplong);
	TIMER_STOP(59);
	if(myid==0) fprintf(stdout,"Treecorrectingtime CPU= %f \n", ELAPSED_TIME(59));
	if(halfstep.first==1 && halfstep.second==0 && iflagPreFoF){
		float PreFoF(treeparticletype *, long);
		if(myid==0) printf("P%d entering the PreFoF\n",myid);
		TIMER_START(58);
		PreFoF(treeparticles,nplong);
		TIMER_STOP(58);
		if(myid==0) fprintf(stdout,"PreFoF CPU= %f \n", ELAPSED_TIME(58));
	}
	tree2pmtype(nplong);
}

