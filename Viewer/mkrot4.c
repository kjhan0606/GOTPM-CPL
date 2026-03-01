/* These are additional frames inserted between the Viewer1 and Viewer2.dat 
 *  for 1.5 + 1.5 seconds.
 * */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>


#define x0  40
#define y0  500
#define z0  420
#define PI (3.1415926535L)

#define dist 14.



int main ( int argc, char **argv){
	int i,j,k;
	int nstep1,stepcount,nstep0;
	int nstep2,nstep3,nstep4,nstep5;
	int nstep6,nstep7,nstep8,nstep9;
	float rx,ry,rz;
	float px,py,pz,p0;
	float vx,vy,vz,vv;
	float dy,dx,dz;
	double yrot,xrot,zrot,rot,drot;
	double dvz,amp,nsmooth;
	int ii,jj;
	float rij,rfact;


	stepcount = 1601;

	nstep0 = 0;
	nstep1 = 45;
	printf("MO4 150\n");

	/* This is to go further */
	nsmooth = 20.;
	rot = 0;
	dz = 0.075;
	pz = z0;
	for(i=nstep0;i<nstep1;i++){
		double ddrot;
		px = x0;
		py = y0;
		pz += dz * erf((nstep1-i)/nsmooth);
		printf("%d %d %g %g %g 1 0 1 1 0 -1 %g 0 0\n",i+1,stepcount,px,py,pz,rot);
	}

	nstep2 = 90;
	for(i=nstep1;i<nstep2;i++){
		printf("%d %d %g %g %g 1 0 1 1 0 -1 %g 0 0\n",i+1,stepcount,px,py,pz,rot);
	}
}
