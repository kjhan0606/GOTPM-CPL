#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "Memory.h"
#include "pmheader.h"
#include "animate.h"


#define LMAX 1000
void splintnr_(float *, float *, float *, int *, float *, float *, float *);
void polintnr_(float *, float *, int *, float *, float *, float *);
#define deg2rad (3.14159265354979L/180.L)
/*
#define INTPOL(a1,a2,a3,a4,a5,a6,a7) polintnr_(a1,a2,a4,a5,a6,a7)
*/
#define INTPOL(a1,a2,a3,a4,a5,a6,a7) splintnr_(a1,a2,a3,a4,a5,a6,a7)

int intint(int *x, int *y, int np, int xx){
	int i,j,k;
	int lo,hi,ires;
	float res;

	for(i=0;i<np-1;i++){
		if((xx-x[i])*(xx-x[i+1]) <= 0){
			res = (float)(y[i+1]-y[i])/(float)(x[i+1]-x[i])*(float)(xx-x[i]) + y[i];
			ires = rint(res);
			return ires;
		}
	}
}

int mkviewer(char *infile, Viewer **viewer){
	FILE *fp;
	int nframe;
	int i,j,k,np,mp;
	float ri;
	float *jframe;
	int *iframe;
	int *istep;
	float *x,*y,*z,*vx,*vy,*vz,*rot;
	float *x2,*y2,*z2,*vx2,*vy2,*vz2,*rot2;
	float ans,ans1;
	float odepth;
	double alpha,beta;
	Viewer *view;
	char outfile[190];


	jframe = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(jframe));
	iframe = (int *) Malloc(sizeof(int)*MAXFRAME,PPTR(iframe));
	istep = (int *) Malloc(sizeof(int)*MAXFRAME,PPTR(istep));
	x = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(x));
	y = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(y));
	z = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(z));
	vx = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(vx));
	vy = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(vy));
	vz = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(vz));
	rot = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(rot));
	x2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(x2));
	y2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(y2));
	z2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(z2));
	vx2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(vx2));
	vy2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(vy2));
	vz2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(vz2));
	rot2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(rot2));

	fp = fopen(infile,"r");
	fscanf(fp,"%s %g\n",outfile,&odepth);
	np = 0;
	while((mp = fscanf(fp,"%d %d %g %g %g %g %g %g %g\n",
					iframe+np,istep+np,x+np,y+np,z+np,vx+np,vy+np,vz+np,rot+np))!= EOF){
		/*
	    printf("Now iframe= %d with np= %d %g %g %g %g %g %g\n", iframe[np],np,
				x[np],y[np],z[np],vx[np],vy[np],vz[np]);
				*/
		np ++;
	}
	fclose(fp);
	nframe = iframe[np-1];

	view = *viewer;

	for(i=0;i<np;i++){
		jframe[i] = iframe[i];
	}
	{
		float yp1,ypn;
		void splinenr_(float *,float*, int *,float *, float *, float *);
		yp1 = ypn = 1.E31;
		splinenr_(jframe,x,&np,&yp1,&ypn,x2);
		splinenr_(jframe,y,&np,&yp1,&ypn,y2);
		splinenr_(jframe,z,&np,&yp1,&ypn,z2);
		splinenr_(jframe,vx,&np,&yp1,&ypn,vx2);
		splinenr_(jframe,vy,&np,&yp1,&ypn,vy2);
		splinenr_(jframe,vz,&np,&yp1,&ypn,vz2);
		splinenr_(jframe,rot,&np,&yp1,&ypn,rot2);
	}

	for(i=0;i<nframe;i++){
		view[i].frame = i+1;
		ri = i+1;
		view[i].HOpticalDepth = odepth;
		sprintf(view[i].outfile,"%s",outfile);

		view[i].nstep = intint(iframe,istep,np,ri);

		INTPOL(jframe,x,x2,&np,&ri,&ans,&ans1);
		view[i].pos.x = ans;
		INTPOL(jframe,y,y2,&np,&ri,&ans,&ans1);
		view[i].pos.y = ans;
		INTPOL(jframe,z,z2,&np,&ri,&ans,&ans1);
		view[i].pos.z = ans;
		INTPOL(jframe,vx,vx2,&np,&ri,&ans,&ans1);
		view[i].view.x = ans;
		INTPOL(jframe,vy,vy2,&np,&ri,&ans,&ans1);
		view[i].view.y = ans;
		INTPOL(jframe,vz,vz2,&np,&ri,&ans,&ans1);
		view[i].view.z = ans;
		INTPOL(jframe,rot,rot2,&np,&ri,&ans,&ans1);
		ans = deg2rad*ans;
		view[i].rot = ans;
	}
	for(i=0;i<nframe;i++){
		double r2,xx,yy,zz;
		r2 = view[i].view.x*view[i].view.x + 
			view[i].view.y*view[i].view.y+
			view[i].view.z*view[i].view.z;
		r2 = sqrt(r2);
		view[i].view.x = view[i].view.x/r2;
		view[i].view.y = view[i].view.y/r2;
		view[i].view.z = view[i].view.z/r2;
		view[i].E3.x = view[i].view.x;
		view[i].E3.y = view[i].view.y;
		view[i].E3.z = view[i].view.z;

		/*
		xx = sqrt(view[i].E3.y*view[i].E3.y/
				(view[i].E3.x*view[i].E3.x + view[i].E3.y*view[i].E3.y));
		yy = sqrt(view[i].E3.x*view[i].E3.x/
				(view[i].E3.x*view[i].E3.x + view[i].E3.y*view[i].E3.y));
				*/
		alpha = acos();
		beta = acos();
		xx = view[i].E3.x;
		zz = view[i].E3.z;
		if(zz != 0 ) {
			double r2;
			view[i].E1.x = 1;
			view[i].E1.z = -xx/zz;
			r2 = 1 + xx*xx/zz/zz;
			r2 = sqrt(r2);
			view[i].E1.x = view[i].E1.x/r2;
			view[i].E1.y = 0;
			view[i].E1.z = view[i].E1.z/r2;
		}
		else {
			view[i].E1.x = 0;
			view[i].E1.y = 0;
			view[i].E1.z = 1;

		}




		view[i].E2.x =  view[i].E3.y*view[i].E1.z - view[i].E3.z*view[i].E1.y;
		view[i].E2.y =  view[i].E3.z*view[i].E1.x - view[i].E3.x*view[i].E1.z;
		view[i].E2.z =  view[i].E3.x*view[i].E1.y - view[i].E3.y*view[i].E1.x;
		/* rotating around E3 */
		{
			double a11,a12,a13,a21,a22,a23,a31,a32,a33;
			double crot,srot;
			a11 = view[i].E2.y*view[i].E3.z - view[i].E3.y*view[i].E2.z;
			a12 = view[i].E1.z*view[i].E3.y - view[i].E1.y*view[i].E3.z;
			a13 = view[i].E1.y*view[i].E2.z - view[i].E2.y*view[i].E1.z;
			a21 = view[i].E3.x*view[i].E2.z - view[i].E2.x*view[i].E3.z;
			a22 = view[i].E1.x*view[i].E3.z - view[i].E3.x*view[i].E1.z;
			a23 = view[i].E2.x*view[i].E1.z - view[i].E1.x*view[i].E2.z;
			a31 = view[i].E2.x*view[i].E3.y - view[i].E3.x*view[i].E2.y;
			a32 = view[i].E3.x*view[i].E1.y - view[i].E1.x*view[i].E3.y;
			a33 = view[i].E1.x*view[i].E2.y - view[i].E2.x*view[i].E1.y;
			crot = cos(view[i].rot); srot = sin(view[i].rot);
			view[i].E1.x =  a11*crot + a12*srot;
			view[i].E1.y =  a21*crot + a22*srot;
			view[i].E1.z =  a31*crot + a32*srot;
			view[i].E2.x = -a11*srot + a12*crot;
			view[i].E2.y = -a21*srot + a22*crot;
			view[i].E2.z = -a31*srot + a32*crot;
		}
	}
	Free(rot2); Free(vz2); Free(vy2); Free(vx2);
	Free(z2); Free(y2); Free(x2);
	Free(rot); Free(vz); Free(vy); Free(vx);
	Free(z); Free(y); Free(x);
	Free(istep);Free(iframe);Free(jframe);
	return nframe;
}
