#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>



int main(int argc, char **argv){
	float app,apsq;
	int i,j,k;



	float astep = 0.05;
	float anow = 1;
	float amax =  100;

	int nstep = 1981;


	float omeplam = 0.74;
	float omep = 0.26;


	float afact,bfact;

	float wlam0 = -1.5;
	float zi = amax-1.;
	float rng = 2048;



	for(i=0;i<nstep;i++){
		float ainv = 1./anow;
		float fact1,fact2;
		float red = amax/anow-1;
		float omei=omep*pow(1+zi,3)/(omep*pow(1+zi,3) + omeplam+(1-omep-omeplam)*pow(1+zi,2));


		app = -4.*M_PI/3.*(ainv*ainv-2.*omeplam/omep*anow/amax/amax/amax);
		apsq = 8.*M_PI/3.*(ainv+1./omei-1.+omeplam/omep*(anow*anow-1.)/amax/amax/amax);
		afact = (2.+app*anow/apsq)/2.*ainv;
		bfact = ainv*ainv*ainv/apsq;
		fact1 = (1.-afact*astep)/(1.+afact*astep);
		fact2 = bfact*astep/(1.+afact*astep)*rng;
		printf("%g : %g %g %g %g %g %g",red,app,apsq,afact,bfact,fact1,fact2);


		omei=omep*pow(1+zi,3)/(omep*pow(1+zi,3) + omeplam*pow(1+zi,3*(1+wlam0))+(1-omep-omeplam)*pow(1+zi,2));
		app = -4.*M_PI/3.*(ainv*ainv+(1.+3*wlam0)*omeplam/omep/pow(anow,2.+3*wlam0)*pow(amax,3*wlam0));
		apsq = 8.*M_PI/3.*(ainv+1./omei-1.+omeplam/omep*(pow(anow,-1.-3*wlam0)-1.)*pow(amax,3*wlam0));
		afact = (2.+app*anow/apsq)/2.*ainv;
		bfact = ainv*ainv*ainv/apsq;
		fact1 = (1.-afact*astep)/(1.+afact*astep);
		fact2 = bfact*astep/(1.+afact*astep)*rng;
		printf(" :: %g %g %g %g %g %g\n",app,apsq,afact,bfact,fact1,fact2);
		anow += astep;
	}
}
