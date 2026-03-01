#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<omp.h>
#include "pmheader.h"
#include "kjhtree.h"

#define DEFINE_TREE_COMMON
#include "force_spline.h"
#undef DEFINE_TREE_COMMON

#define RANGE 6
#define MIN(a,b) (a)>(b) ? (b): (a)
double LL;
double pmcf(double);
double AA[8] = {
	1.L/256.L/256.L,
	0.78224E+00L,
	0.37971E-06L,
	-0.60338E+00L,
	-0.34419E-07L,
	-0.39741E+01L,
	-0.10607E+01L,
	-0.38145E+00L
};
void i_force_spline(int nx, float rsphere){
	double r;
	int i,j,k;
	double e,dr;
	double rp2,rp1,rm1,rm2;
	double  fp, fpp,fr;

	LL = nx*nx;
	ran2nran=(double)RANGE/(double)NSPLINE;
	invran2nran = 1./ran2nran;

	for(i=0;i<3;i++){
		forcecorrectslope(0,i) = forcecorrectdiff(0,i) = 0.;
	}
	dr = ran2nran*0.1;
#ifdef _OPENMP
#pragma omp parallel for private(i,r,rm1,rm2,rp1,rp2,fp,fpp,fr)
#endif
	for(i=1;i<NSPLINE;i++){
		r = (double)(i) * ran2nran;
		rm1 = r - dr; rm2 = r - 2*dr; rp1 = r + dr; rp2 = r + 2*dr;
		fp = (pmcf(rp1) - pmcf(rm1))/(dr+dr);
		fpp = (pmcf(rp2) - 2*pmcf(r) + pmcf(rm2))/(4*dr*dr);
		fr = pmcf(r);
		forcecorrectdiff(i,0) = -(double)(fr/r);
		forcecorrectdiff(i,1) = -(double)(0.5L*(fp*r-fr)/r/r/r);
		forcecorrectdiff(i,2) = -(double)(0.5L*(fpp*r*r - 3*fp*r + 3*fr)/r/r/r/r/r);
		/*
		diff[i][0] = -fr;
		diff[i][1] = -0.5*(fp*r-fr);
		diff[i][2] = -0.5L*(fpp*r*r - 3.L*fp*r + 3.L*fr);
		*/
	}
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	for(i=0;i<NSPLINE;i++){
		forcecorrectslope(i,0) = (forcecorrectdiff(i+1,0)-forcecorrectdiff(i,0))/ran2nran;
		forcecorrectslope(i,1) = (forcecorrectdiff(i+1,1)-forcecorrectdiff(i,1))/ran2nran;
		forcecorrectslope(i,2) = (forcecorrectdiff(i+1,2)-forcecorrectdiff(i,2))/ran2nran;
	}
	for(i=0;i<3;i++){
		forcecorrectdiff(NSPLINE-1,i)  =0 ;
		forcecorrectslope(NSPLINE-1,i) =0 ;
	}

	return;
}
double pmcf(double r){
	double fr;
	double tmp1;
	double e;
	/*
	if(r == 0.L) {
		fr = 0.;
	}
	else 
	*/
	{
		tmp1 = cosh(AA[1]*r);
		fr = r/pow(r*r+EPSILON*EPSILON,1.5L)
			- (
					(1L/(r*r)*tanh(AA[1]*r)-AA[1]/(tmp1*tmp1*r)
					-2L*AA[2]/AA[0]*r*exp(AA[3]*r*r)
					-2L*AA[2]*AA[3]/AA[0]*r*r*r*exp(AA[3]*r*r)
					+AA[4]/AA[0]*(1L+AA[5]*r*r+AA[6]*r*r*r*r)*exp(AA[7]*r*r))
				);
		fr = fr/LL; /* scaling for the PM force \propto 1/L**2 */
	}
	return fr;
}
void i_potent_spline(){
	double x;
	double xstep;
	int i,j,k;
	double e,dx;
	double g(double);
	double xp2,xp1,xm1,xm2;
	double  gp, gpp;
	ran2nran=(double)RANGE/(double)NSPLINE;
	invran2nran=1./ran2nran;

	for(i=0;i<3;i++){
		forcecorrectslope(0,i) = forcecorrectdiff(0,i) = 0.;
	}
	dx = ran2nran;

#ifdef _OPENMP
#pragma omp parallel for private(i,x,xm1,xm2,xp1,xp2,gp,gpp)
#endif
	for(i=1;i<NSPLINE;i++){
		x = (double)(i) * ran2nran;
		xm1 = x - dx; xm2 = x - 2*dx; xp1 = x + dx; xp2 = x + 2*dx;
		gp = (g(xp1) - g(xm1))/(2*dx);
		gpp = (g(xp2) - 2*g(x) + g(xm2))/(4.*dx*dx);
		forcecorrectdiff(i,0) = g(x);
		forcecorrectdiff(i,1) = 0.5L*gp/x;
		forcecorrectdiff(i,2) = 0.5L*(gpp/x/x - gp/x/x/x);
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=1;i<NSPLINE;i++){
		forcecorrectslope(i,0) = (forcecorrectdiff(i+1,0)-forcecorrectdiff(i,0))/dx;
		forcecorrectslope(i,1) = (forcecorrectdiff(i+1,1)-forcecorrectdiff(i,1))/dx;
		forcecorrectslope(i,2) = (forcecorrectdiff(i+1,2)-forcecorrectdiff(i,2))/dx;
	}
	return;
}
double g(double x){
	double e,g;
	g = -1.L/sqrt(x*x+EPSILON*EPSILON);
	return g;
}
