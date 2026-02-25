#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<stdint.h>
#include<math.h>

#include "pmheader.h"
#include "pm_common.h"
#include "flow.h"
#include "kjhtree.h"

#define MAX(a,b) ((a)>(b) ? (a):(b))


EvolFact GetEvolFactor(float anow,float astep){

	float ainv,app,apsq;
	double afact,bfact,fact1,fact2;
	EvolFact evolfactor;
	ainv = 1./anow; 
	/*
	app = -4.*M_PI/3.*(ainv*ainv+(1.+3*simpar.wlam)*simpar.omeplam/simpar.omep/pow(anow,2.+3*simpar.wlam)*pow(simpar.amax, 3*simpar.wlam)); 
	apsq = 8.*M_PI/3.*(ainv+1./simpar.omei-1.+simpar.omeplam/simpar.omep*(pow(anow,-1.-3*simpar.wlam)-1.)*pow(simpar.amax,3*simpar.wlam));
	*/
	appapsq_(&simpar.omep, &simpar.omeplam, &simpar.wlam0, &simpar.wlam1, &simpar.amax, &anow, &app, &apsq);

	afact = (2+app*anow/apsq)/2*ainv;
	bfact = ainv*ainv*ainv/apsq;
	fact1 = (1-afact*astep)/(1+afact*astep);
	fact2 = bfact*astep/(1+afact*astep);
	evolfactor.afact = afact;
	evolfactor.bfact = bfact;
	evolfactor.fact1 = fact1;
	evolfactor.fact2 = fact2;
	evolfactor.pfact = astep;
	evolfactor.fact1_push = 1/(1+afact*astep*0.5);
	evolfactor.fact2_push = bfact*astep*0.5/(1+afact*astep*0.5);
	evolfactor.fact1_pull = (1-afact*astep*0.5);
	evolfactor.fact2_pull = bfact*astep*0.5;
	return evolfactor;
}

#ifdef INDTIME
void DetermineIndT(){
	int i,nsub,isub;
	float astep;
	double vr,dist;
	pmparticletype *pp;
	astep = simpar.astep;
	pp = pmparticles;
	maxsubT = 0;
	for(i=0;i<np;i++){
		vr = pp->vx*pp->vx + pp->vy*pp->vy + pp->vz*pp->vz;
		vr = sqrt(vr);
		dist = vr*astep;
		nsub = GetNSub(dist);
		isub = ceil(log(nsub)/log(2.));
		pp->tflag = isub;
		maxsubT = MAX(maxsubT,isub);
	}
}
#endif
