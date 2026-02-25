#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>

#include "pmheader.h"
#include "pm_common.h"



void CObs0ExtractingLightConeParticles(pmparticletype *, int ,float , float);
void Obs0ExtractingLightConeParticles(pmparticletype *, int ,float , float);
void Obs1ExtractingLightConeParticles(pmparticletype *, int ,float , float);
void Obs2ExtractingLightConeParticles(pmparticletype *, int ,float , float);
void Obs3ExtractingLightConeParticles(pmparticletype *, int ,float , float);
void Obs4ExtractingLightConeParticles(pmparticletype *, int ,float , float);
void Obs5ExtractingLightConeParticles(pmparticletype *, int ,float , float);
void Obs6ExtractingLightConeParticles(pmparticletype *, int ,float , float);
void Obs7ExtractingLightConeParticles(pmparticletype *, int ,float , float);
void PSliceExtractingLightConeParticles(pmparticletype *, int ,float , float);
void PBeamExtractingLightConeParticles(pmparticletype *, int ,float , float);


void observerextmain(){
	float nextastep,maxd;
	int myid;
	{/* This is to extract particles in the lightcone space */
		float maxd;
		nextastep = simpar.astep;
		maxd = 0.1*simpar.nspace;
		myid = simpar.myid;
		TIMER_START(82);
		CObs0ExtractingLightConeParticles(pmparticles, np, nextastep, maxd);
		Obs0ExtractingLightConeParticles(pmparticles, np, nextastep, maxd);
		Obs1ExtractingLightConeParticles(pmparticles, np, nextastep, maxd);
		Obs2ExtractingLightConeParticles(pmparticles, np, nextastep, maxd);
		Obs3ExtractingLightConeParticles(pmparticles, np, nextastep, maxd);
		Obs4ExtractingLightConeParticles(pmparticles, np, nextastep, maxd);
		Obs5ExtractingLightConeParticles(pmparticles, np, nextastep, maxd);
		Obs6ExtractingLightConeParticles(pmparticles, np, nextastep, maxd);
		Obs7ExtractingLightConeParticles(pmparticles, np, nextastep, maxd);
		PSliceExtractingLightConeParticles(pmparticles, np, nextastep, maxd);
		PBeamExtractingLightConeParticles(pmparticles, np, nextastep, maxd);
		TIMER_STOP(82);
		if(myid==0) fprintf(stdout,"Setting Observation CPU= %f \n", ELAPSED_TIME(82));
	}
}
void CObs0SavingLightConeData(pmparticletype *,int ,float); 
void Obs0SavingLightConeData(pmparticletype *,int ,float); 
void Obs1SavingLightConeData(pmparticletype *,int ,float); 
void Obs2SavingLightConeData(pmparticletype *,int ,float); 
void Obs3SavingLightConeData(pmparticletype *,int ,float); 
void Obs4SavingLightConeData(pmparticletype *,int ,float); 
void Obs5SavingLightConeData(pmparticletype *,int ,float); 
void Obs6SavingLightConeData(pmparticletype *,int ,float); 
void Obs7SavingLightConeData(pmparticletype *,int ,float); 
void PBeamSavingLightConeData(pmparticletype *,int ,float); 
void PSliceSavingLightConeData(pmparticletype *,int ,float); 
void observersavemain(){
	float nextastep,maxd;
	int myid;
	{
		TIMER_START(82);
		myid = simpar.myid;
		nextastep = simpar.astep;
		CObs0SavingLightConeData(pmparticles,np,nextastep);
		Obs0SavingLightConeData(pmparticles,np,nextastep);
		Obs1SavingLightConeData(pmparticles,np,nextastep);
		Obs2SavingLightConeData(pmparticles,np,nextastep);
		Obs3SavingLightConeData(pmparticles,np,nextastep);
		Obs4SavingLightConeData(pmparticles,np,nextastep);
		Obs5SavingLightConeData(pmparticles,np,nextastep);
		Obs6SavingLightConeData(pmparticles,np,nextastep);
		Obs7SavingLightConeData(pmparticles,np,nextastep);
		PBeamSavingLightConeData(pmparticles,np,nextastep);
		PSliceSavingLightConeData(pmparticles,np,nextastep);
		TIMER_STOP(82);
		if(myid==0) fprintf(stdout,"Saving Observation CPU= %f \n", ELAPSED_TIME(82));
	}
}

