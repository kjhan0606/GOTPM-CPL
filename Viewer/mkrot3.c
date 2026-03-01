#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>


#define x0  40
#define y0  500
#define z0  434
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
	nstep1 = 900;
	drot = 2*PI/(nstep1-nstep0);
	printf("MO3 150\n");

	/* This is to rotate around the cluster halo */
	amp = 0;
	nsmooth = 80.;
	for(i=nstep0;i<nstep1;i++){
		amp += erf((i+1.-nstep0)/nsmooth)*erf((nstep1-i)/nsmooth);
	}
	drot = 2*PI/amp;
	rot = 0;
	for(i=nstep0;i<nstep1;i++){
		double ddrot;
		rot += drot*erf((i+1.-nstep0)/nsmooth)*erf((nstep1-i)/nsmooth);
		px = x0;
		py = y0 + dist*sin(rot);
		pz = z0 - dist*cos(rot);
		printf("%d %d %g %g %g 1 0 1 1 0 -1 %g 0 0\n",i+1,stepcount,px,py,pz,rot);
	}

	rot = 0;
	/* This is to enter into the cluster and pass through it */
	nstep2 = 1300;
	dz = dist*2./(float)(nstep2-nstep1);

	amp = 0;
	nsmooth = 80.;
	for(i=nstep1;i<nstep2;i++){
		amp += erf((i+1.-nstep1)/nsmooth)*erf((nstep2-i)/nsmooth);
	}
	dz = dist*2./amp;
	for(i=nstep1;i<nstep2;i++){
		pz += dz*erf((i+1.-nstep1)/nsmooth)*erf((nstep2-i)/nsmooth);
		printf("%d %d %g %g %g 1 0 1 1 0 -1 %g 0 0\n",i+1,stepcount,px,py,pz,rot);
	}

	/* This is to turn around the observer's view point */
	nstep3 = 1400;

	amp = 0;
	nsmooth = 80.;
	for(i=nstep2;i<nstep3;i++){
		amp += erf((i+1.-nstep2)/nsmooth)*erf((nstep3-i)/nsmooth);
	}
	drot = PI/amp;

	for(i=nstep2;i<nstep3;i++){
		vx = 1.;
		rot += drot*erf((i+1.-nstep2)/nsmooth)*erf((nstep3-i)/nsmooth);
		vy = sin(rot);
		vz = cos(rot);
		rot = drot *(i-nstep2);
		printf("%d %d %g %g %g 1 0 1 1 0 -1 %g 0 0\n",i+1,stepcount,px,py,pz,rot);
	}

	/* This is to go back */
	nstep4 = 1800;

	amp = 0;
	nsmooth = 80.;
	for(i=nstep3;i<nstep4;i++){
		amp += erf((i+1.-nstep3)/nsmooth)*erf((nstep4-i)/nsmooth);
	}
	dz = -dist*2./amp;
	for(i=nstep3;i<nstep4;i++){
		pz += dz*erf((i+1.-nstep3)/nsmooth)*erf((nstep4-i)/nsmooth);
		vx = 1.;
		vy = 0.;
		vz = -1.;
		printf("%d %d %g %g %g 1 0 1 1 0 -1 %g 0 0\n",i+1,stepcount,px,py,pz,rot);
	}

	/* This is to turn around the observer's view point */
	nstep5 = 1900;

	amp = 0;
	nsmooth = 15.;
	for(i=nstep4;i<nstep5;i++){
		amp += erf((i+1.-nstep4)/nsmooth)*erf((nstep5-i)/nsmooth);
	}


	drot = PI/amp;
	rot = 0;
	for(i=nstep4;i<nstep5;i++){
		vx = 1.;
		vy = -sin(rot);
		vz = -cos(rot);
		rot += drot*erf((i+1.-nstep4)/nsmooth)*erf((nstep5-i)/nsmooth);;
		printf("%d %d %g %g %g 1 0 1 1 0 -1 %g 0 0\n",i+1,stepcount,px,py,pz,rot);
	}
	printf("1901 1601 40 500 420 1 0 1 1 0 -1 0 0 0 \n");
	px = 40;
	py = 500;
	pz = 420;
	vx = 1.;
	vy = 0;
	rot = 0;

	/* This is to approach to the cluster and lookup up slightly. */
	nstep5 = 1901;
	nstep6 = 2000;

	amp = 0;
	nsmooth = 15.;
	for(i=nstep5;i<nstep6;i++){
		amp += erf((i+1.-nstep5)/nsmooth)*erf((nstep6-i)/nsmooth);
	}
	dz = dist*0.5/amp;
	drot = -15.*PI/180./amp;
	rot = 0;


	for(i=nstep5;i<nstep6;i++){
		pz += dz*erf((i+1.-nstep5)/nsmooth)*erf((nstep6-i)/nsmooth);
		rot += drot*erf((i+1.-nstep5)/nsmooth)*erf((nstep6-i)/nsmooth);
		printf("%d %d %g %g %g 1 0 1 1 0 -1 0 %g 0\n",i+1,stepcount,px,py,pz,rot);
	}


	/* This is to shift the observer's position */

	nstep7 = 2100;
	amp = 0;
	nsmooth = 15.;
	for(i=nstep6;i<nstep7;i++){
		amp += erf((i+1.-nstep6)/nsmooth)*erf((nstep7-i)/nsmooth);
	}

	dx = (40.-x0)/amp;
	dy = (502.-y0)/amp;
	for(i=nstep6;i<nstep7;i++){
		px += dx*erf((i+1.-nstep6)/nsmooth)*erf((nstep7-i)/nsmooth);
		py += dy*erf((i+1.-nstep6)/nsmooth)*erf((nstep7-i)/nsmooth);
		printf("%d %d %g %g %g 1 0 1 1 0 -1 0 %g 0\n",i+1,stepcount,px,py,pz,rot);
	}

	nstep8 = 2796;
	amp = 0;
	nsmooth = 100.;
	for(i=nstep7;i<nstep8;i++){
		amp += erf((i+1.-nstep7)/nsmooth)*erf((nstep8-i)/nsmooth);
	}


	drot = 2*PI/amp;
	yrot = rot;
	xrot = zrot = 0;
	dz = dist*0.15/amp;
	for(i=nstep7;i<nstep8;i++){
		xrot += drot*erf((i+1.-nstep7)/nsmooth)*erf((nstep8-i)/nsmooth);
		zrot -= drot*erf((i+1.-nstep7)/nsmooth)*erf((nstep8-i)/nsmooth);
		yrot += drot*erf((i+1.-nstep7)/nsmooth)*erf((nstep8-i)/nsmooth);
		pz += dz*erf((i+1.-nstep7)/nsmooth)*erf((nstep8-i)/nsmooth);
		printf("%d %d %g %g %g 1 0 1 1 0 -1 %g %g %g\n",i+1,stepcount,px,py,pz,
				xrot,yrot,zrot);
	}
}
