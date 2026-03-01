#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>


#define x0  40
#define y0  500
#define z0  432
#define PI (3.1415926535L)

#define dist 12.



int main ( int argc, char **argv){
	int i,j,k;
	int nstep1,stepcount;
	int nstep2,nstep3,nstep4,nstep5;
	int nstep6,nstep7,nstep8,nstep9;
	float rx,ry,rz;
	float px,py,pz;
	float vx,vy,vz,vv;
	float dy,dx,dz;
	double rot,drot;
	double dvz;


	stepcount = 1601;

	nstep1 = 900;
	drot = 2*PI/nstep1;
	printf("MO3 150\n");

	/* This is to rotate around the cluster halo */
	for(i=0;i<nstep1;i++){
		rot = drot*i;
		px = x0;
		py = y0 + dist*sin(rot);
		pz = z0 - dist*cos(rot);
		vx = x0-px;
		vy = y0-py;
		vz = z0-pz;
		vv = sqrt(vx*vx+vy*vy+vz*vz);
		vx = vx/vv;
		vy = vy/vv;
		vz = vz/vv;
		vx = vx + 1;
		printf("%d %d %g %g %g %g %g %g 0\n",i+1,stepcount,px,py,pz,vx,vy,vz);
	}
	rot = 0;
	/* This is to enter into the cluster and pass through it */
	nstep2 = 1300;
	dz = dist*2./(float)(nstep2-nstep1);
	for(i=nstep1;i<nstep2;i++){
		pz = pz + dz;
		printf("%d %d %g %g %g %g %g %g %g\n",i+1,stepcount,px,py,pz,vx,vy,vz,rot);
	}

	/* This is to turn around the observer's view point */
	nstep3 = 1400;
	drot = PI/(nstep3-nstep2);
	for(i=nstep2;i<nstep3;i++){
		vx = 1.;
		vy = sin(drot*(i-nstep2));
		vz = cos(drot*(i-nstep2));
		printf("%d %d %g %g %g %g %g %g %g\n",i+1,stepcount,px,py,pz,vx,vy,vz,rot);
	}

	/* This is to go back */
	nstep4 = 1800;
	dz = dist*2./(float)(nstep4-nstep3);
	for(i=nstep3;i<nstep4;i++){
		pz = pz - dz;
		vx = 1.;
		vy = 0.;
		vz = -1.;
		printf("%d %d %g %g %g %g %g %g %g\n",i+1,stepcount,px,py,pz,vx,vy,vz,rot);
	}

	/* This is to turn around the observer's view point */
	nstep5 = 1900;
	drot = PI/(nstep5-nstep4);
	for(i=nstep4;i<nstep5;i++){
		vx = 1.;
		vy = -sin(drot*(i-nstep4));
		vz = -cos(drot*(i-nstep4));
		printf("%d %d %g %g %g %g %g %g %g\n",i+1,stepcount,px,py,pz,vx,vy,vz,rot);
	}
	printf("1901 1601 40 500 420 1 0 1 0 \n");
	px = 40;
	py = 500;
	pz = 420;
	vx = 1.;
	vy = 0;
	rot = 0;

	/* This is to enter into the cluster and pass through it */
	nstep5 = 1901;
	nstep6 = 2000;

	dvz = (sqrt(3.)-1.)/(nstep6-nstep5);
	dz = dist*0.5/(float)(nstep6-nstep5);
	for(i=nstep5;i<nstep6;i++){
		pz = pz + dz;
		vz = dvz*(i-nstep5) + 1.;
		printf("%d %d %g %g %g %g %g %g %g\n",i+1,stepcount,px,py,pz,vx,vy,vz,rot);
	}
}
