#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include "Memory.h"
#include "nr.h"
#include "nrutil.h"
#define pi 3.1415926535
#define NDIM 3
#define NAN -9.E20
#define N 200
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

typedef struct vector3d{
	float x,y,z;
} vector3d;
double a,b,c,A,B,C;
double R[4][4];
double invR[4][4];
double xprime,yprime,zprime;
void determinABC(){
	A = (1.L/(a*a)+c*c*R[1][3]*R[1][3]/R[3][3]/R[3][3]/(a*a*a*a));
	C = (1.L/(b*b)+c*c*R[2][3]*R[2][3]/R[3][3]/R[3][3]/(b*b*b*b));
	B = c*c/(a*a*b*b)*R[1][3]*R[2][3]/R[3][3]/R[3][3]/A;

}
double determinexprime(int n){
	double tmp;
	tmp = 1.L+yprime*yprime*(A*B*B-C);
	if(tmp < 0.L) return NAN;
	xprime = -B*yprime + n*sqrt(tmp)/sqrt(A);
	if(xprime >= -a && xprime <= a){
		return xprime;
	}
	else {
		return NAN;
	}
}
double determinezprime(){
	return -c*c/R[3][3]*(R[1][3]*xprime/(a*a) + R[2][3]*yprime/(b*b));
}
void inverseR(){
	int i,j,k;
	for(i=1;i<=3;i++){
		for(j=1;j<=3;j++){
			invR[i][j] = R[j][i];
		}
	}
}
#define REPS 1.E-7L
double findyprime(int n){
	int i,j,k;
	double res;
	double ymin,ymax,localxprime; 
	if(n == 1) {
		ymax = b+REPS;ymin = 0.L;
	}
	else if(n == -1) {
		ymin = -(b+REPS);ymax = 0.L;
	}
	else {
		fprintf(stderr,"Error input findyprime\n");exit(999);
	}
	yprime = (ymax+ymin)*0.5L;
	if(n == 1){
		while(fabs(yprime-ymin) > REPS){
			if((localxprime = determinexprime(-n)) != NAN){
				ymin = yprime;
			}
			else {
				ymax = yprime;
			}
			yprime = (ymax+ymin)*0.5L;
		}
		return yprime;
	}
	else {
		while(fabs(yprime-ymin) > REPS){
			if((localxprime = determinexprime(-n)) != NAN){
				ymax = yprime;
			}
			else {
				ymin = yprime;
			}
			yprime = (ymax+ymin)*0.5L;
		}
		return yprime;
	}
}
#undef REPS
int shape(float r,float q, float s,float *r1,
		float *r2,float *r3, float *x,float *y, float *z){
	int i,j,k;
	FILE *wp;
	double yprimestep;
	double minyprime,maxyprime;
	double oxprime;
	int np;
	a = r; b = a*q; c = a*s;
	R[1][1] = r1[0]; R[1][2] = r1[1]; R[1][3] = r1[2];
	R[2][1] = r2[0]; R[2][2] = r2[1]; R[2][3] = r2[2];
	R[3][1] = r3[0]; R[3][2] = r3[1]; R[3][3] = r3[2];
	inverseR();
	determinABC();
	/*
	yprimestep =  2.L * b /(double)N;
	*/
	yprimestep = pi/(double)N;
	yprime = -b; findyprime(-1);
	np = 0;
	/*
	while(yprime <= b){
	*/
	for(i=0;i<N;i++){
		xprime = determinexprime(1);
		if(xprime == NAN) {
			/*
			yprime += yprimestep;
			*/
			yprime = b*cos(acos(yprime/b)-yprimestep);
			continue;
		}
		zprime = determinezprime();
		x[np] = invR[1][1]*xprime + invR[1][2]*yprime + invR[1][3]*zprime;
		y[np] = invR[2][1]*xprime + invR[2][2]*yprime + invR[2][3]*zprime;
		z[np] = invR[3][1]*xprime + invR[3][2]*yprime + invR[3][3]*zprime;
		np ++;
		/*
		yprime += yprimestep;
		*/
		yprime = b*cos(acos(yprime/b)-yprimestep);
	}

	yprime = b; findyprime(1);
	/*
	while(yprime >= -b){
	*/
	for(i=0;i<N;i++){
		xprime = determinexprime(-1);
		if(xprime == NAN) {
			yprime = b*cos(acos(yprime/b)+yprimestep);
			/*
			yprime -= yprimestep;
			*/
			continue;
		}
		zprime = determinezprime();
		x[np] = invR[1][1]*xprime + invR[1][2]*yprime + invR[1][3]*zprime;
		y[np] = invR[2][1]*xprime + invR[2][2]*yprime + invR[2][3]*zprime;
		z[np] = invR[3][1]*xprime + invR[3][2]*yprime + invR[3][3]*zprime;
		np++;
		yprime = b*cos(acos(yprime/b)+yprimestep);
		/*
		yprime -= yprimestep;
		*/
	}
	x[np] = x[0]; y[np] = y[0]; z[np] = z[0];
	np++;
	return np;
}
#define REPS 1.E-4
#define NEPS 0.98
extern int myid;
float find_r(vector3d *r,int np, float q, float s, float *r1,float *r2,float *r3,
		double cx,double cy,double cz){
	int i,j,k;
	float rmin,rmax,radius;
	float tmpx,tmpy,tmpz,tmp;
	float xx,yy,zz;
	float *dist;
	int mp,ncount;
	mp = MAX(np*NEPS,1);

	/*
	cx = cy = cz = 0.;
	for(i=0;i<np;i++){
		cx += r[i].x; cy += r[i].y; cz += r[i].z;
	}
	cx = cx / (float)np; cy = cy / (float)np; cz = cz / (float)np;
	*/

	R[1][1] = r1[0]; R[1][2] = r1[1]; R[1][3] = r1[2];
	R[2][1] = r2[0]; R[2][2] = r2[1]; R[2][3] = r2[2];
	R[3][1] = r3[0]; R[3][2] = r3[1]; R[3][3] = r3[2];
	dist = (float *)Malloc(sizeof(float)*np,PPTR(dist));
	rmin = rmax = 0.;
	for(i=0;i<np;i++){
		tmpx = r[i].x-cx; tmpy = r[i].y-cy; tmpz = r[i].z-cz;
		xx = R[1][1]*tmpx + R[1][2]*tmpy + R[1][3]*tmpz;
		yy = R[2][1]*tmpx + R[2][2]*tmpy + R[2][3]*tmpz;
		zz = R[3][1]*tmpx + R[3][2]*tmpy + R[3][3]*tmpz;
		tmp = sqrt(xx*xx + yy*yy/(q*q) + zz*zz/(s*s));
		rmax = MAX(rmax,tmp);
		dist[i] = tmp;
	}
	radius = (rmax+rmin)*0.5;
	while(fabs(radius-rmin)/radius > REPS){
		ncount = 0;
		for(i=0;i<np;i++){
			if(dist[i] < radius) ncount ++;
		}
		if(ncount > mp) rmax = radius;
		else rmin = radius;
		radius = (rmax+rmin)*0.5;
	}
	Free(dist);
	return radius;
}
#undef NEPS
#undef REPS
