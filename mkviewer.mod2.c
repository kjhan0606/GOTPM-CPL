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
double Rxyz[3][3];

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
	return 0;
}

int mkviewer(char *infile, Viewer **viewer){
	FILE *fp;
	int nframe;
	int i,j,k,np,mp;
	float ri;
	float *jframe;
	int *iframe;
	int *istep;
	float *x,*y,*z,*rx,*ry,*rz,*rot;
	float *vx,*vy,*vz;
	float *x2,*y2,*z2,*rx2,*ry2,*rz2,*rot2;
	float *vx2,*vy2,*vz2;
	float *alpha,*beta,*gamma;
	float *alpha2,*beta2,*gamma2;
	float ans,ans1;
	float odepth;
	double al,be,ga;
	Viewer *view;
	char outfile[190];


	jframe = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(jframe));
	iframe = (int *) Malloc(sizeof(int)*MAXFRAME,PPTR(iframe));
	istep = (int *) Malloc(sizeof(int)*MAXFRAME,PPTR(istep));
	x = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(x));
	y = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(y));
	z = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(z));
	rx = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(rx));
	ry = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(ry));
	rz = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(rz));
	vx = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(rx));
	vy = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(ry));
	vz = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(rz));
	alpha = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(alpha));
	beta = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(beta));
	gamma = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(gamma));

	x2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(x2));
	y2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(y2));
	z2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(z2));
	rx2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(rx2));
	ry2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(ry2));
	rz2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(rz2));
	vx2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(rx2));
	vy2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(ry2));
	vz2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(rz2));
	alpha2= (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(alpha2));
	beta2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(beta2));
	gamma2 = (float *) Malloc(sizeof(float)*MAXFRAME,PPTR(gamma2));





	fp = fopen(infile,"r");
	fscanf(fp,"%s %g\n",outfile,&odepth);
	np = 0;
	while((mp = fscanf(fp,"%d %d %g %g %g %g %g %g %g %g %g %g %g %g\n",
					iframe+np,istep+np,x+np,y+np,z+np,rx+np,ry+np,rz+np,vx+np,
					vy+np,vz+np,alpha+np,beta+np,gamma+np))!= EOF){
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

		splinenr_(jframe,rx,&np,&yp1,&ypn,rx2);
		splinenr_(jframe,ry,&np,&yp1,&ypn,ry2);
		splinenr_(jframe,rz,&np,&yp1,&ypn,rz2);

		splinenr_(jframe,vx,&np,&yp1,&ypn,vx2);
		splinenr_(jframe,vy,&np,&yp1,&ypn,vy2);
		splinenr_(jframe,vz,&np,&yp1,&ypn,vz2);

		splinenr_(jframe,alpha,&np,&yp1,&ypn,alpha2);
		splinenr_(jframe,beta,&np,&yp1,&ypn,beta2);
		splinenr_(jframe,gamma,&np,&yp1,&ypn,gamma2);
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

		INTPOL(jframe,rx,rx2,&np,&ri,&ans,&ans1);
		view[i].E3.x = view[i].view.x = ans;
		INTPOL(jframe,ry,ry2,&np,&ri,&ans,&ans1);
		view[i].E3.y = view[i].view.y = ans;
		INTPOL(jframe,rz,rz2,&np,&ri,&ans,&ans1);
		view[i].E3.z = view[i].view.z = ans;

		INTPOL(jframe,vx,vx2,&np,&ri,&ans,&ans1);
		view[i].E1.x = ans;
		INTPOL(jframe,vy,vy2,&np,&ri,&ans,&ans1);
		view[i].E1.y = ans;
		INTPOL(jframe,vz,vz2,&np,&ri,&ans,&ans1);
		view[i].E1.z = ans;
		{
			double r2;
			/* Normalize them */
			r2 = sqrt(view[i].E3.x*view[i].E3.x+view[i].E3.y*view[i].E3.y+
					view[i].E3.z*view[i].E3.z);
			view[i].E3.x = view[i].view.x = view[i].E3.x/r2;
			view[i].E3.y = view[i].view.y = view[i].E3.y/r2;
			view[i].E3.z = view[i].view.z = view[i].E3.z/r2;
			r2 = sqrt(view[i].E1.x*view[i].E1.x+view[i].E1.y*view[i].E1.y+
					view[i].E1.z*view[i].E1.z);
			view[i].E1.x = view[i].E1.x/r2;
			view[i].E1.y = view[i].E1.y/r2;
			view[i].E1.z = view[i].E1.z/r2;
		}



		view[i].E2.x =  view[i].E3.y*view[i].E1.z - view[i].E3.z*view[i].E1.y;
		view[i].E2.y =  view[i].E3.z*view[i].E1.x - view[i].E3.x*view[i].E1.z;
		view[i].E2.z =  view[i].E3.x*view[i].E1.y - view[i].E3.y*view[i].E1.x;


		INTPOL(jframe,alpha,alpha2,&np,&ri,&ans,&ans1);
		al = view[i].alpha = ans;
		INTPOL(jframe,beta,beta2,&np,&ri,&ans,&ans1);
		be = view[i].beta = ans;
		INTPOL(jframe,gamma,gamma2,&np,&ri,&ans,&ans1);
		ga = view[i].gamma = ans;


		/*
		Rxyz[0][0] =  cos(ga)*cos(be);
		Rxyz[0][1] =  sin(al)*sin(be)*cos(ga) + cos(al)*sin(ga);
		Rxyz[0][2] = -cos(al)*sin(be)*cos(ga) + sin(al)*sin(ga);
		Rxyz[1][0] = -cos(be)*sin(ga);
		Rxyz[1][1] = -sin(al)*sin(be)*sin(ga) + cos(al)*cos(ga);
		Rxyz[1][2] =  cos(al)*sin(be)*sin(ga) + sin(al)*cos(ga);
		Rxyz[2][0] =  sin(be);
		Rxyz[2][1] = -sin(al)*cos(be);
		Rxyz[2][2] =  cos(al)*cos(be);
		*/
		Rxyz[0][0] =  cos(be)*cos(ga);
		Rxyz[0][1] =  sin(al)*sin(be)*cos(ga) - cos(al)*sin(ga);
		Rxyz[0][2] =  cos(al)*sin(be)*cos(ga) + sin(al)*sin(ga);
		Rxyz[1][0] =  cos(be)*sin(ga);
		Rxyz[1][1] =  sin(al)*sin(be)*sin(ga) + cos(al)*cos(ga);
		Rxyz[1][2] =  cos(al)*sin(be)*sin(ga) - sin(al)*cos(ga);
		Rxyz[2][0] = -sin(be);
		Rxyz[2][1] =  sin(al)*cos(be);
		Rxyz[2][2] =  cos(al)*cos(be);
		{
			double a00,a01,a02,a10,a11,a12,a20,a21,a22;
			a00 = Rxyz[0][0]*view[i].E1.x+Rxyz[0][1]*view[i].E1.y+Rxyz[0][2]*view[i].E1.z;
			a10 = Rxyz[1][0]*view[i].E1.x+Rxyz[1][1]*view[i].E1.y+Rxyz[1][2]*view[i].E1.z;
			a20 = Rxyz[2][0]*view[i].E1.x+Rxyz[2][1]*view[i].E1.y+Rxyz[2][2]*view[i].E1.z;
			a01 = Rxyz[0][0]*view[i].E2.x+Rxyz[0][1]*view[i].E2.y+Rxyz[0][2]*view[i].E2.z;
			a11 = Rxyz[1][0]*view[i].E2.x+Rxyz[1][1]*view[i].E2.y+Rxyz[1][2]*view[i].E2.z;
			a21 = Rxyz[2][0]*view[i].E2.x+Rxyz[2][1]*view[i].E2.y+Rxyz[2][2]*view[i].E2.z;
			a02 = Rxyz[0][0]*view[i].E3.x+Rxyz[0][1]*view[i].E3.y+Rxyz[0][2]*view[i].E3.z;
			a12 = Rxyz[1][0]*view[i].E3.x+Rxyz[1][1]*view[i].E3.y+Rxyz[1][2]*view[i].E3.z;
			a22 = Rxyz[2][0]*view[i].E3.x+Rxyz[2][1]*view[i].E3.y+Rxyz[2][2]*view[i].E3.z;
			view[i].E1.x = a00;
			view[i].E1.y = a10;
			view[i].E1.z = a20;
			view[i].E2.x = a01;
			view[i].E2.y = a11;
			view[i].E2.z = a21;
			view[i].E3.x = a02;
			view[i].E3.y = a12;
			view[i].E3.z = a22;
		}
	}




	Free(gamma2);Free(beta2);Free(alpha2); 
	Free(vz2);Free(vy2);Free(vx2); 
	Free(rz2); Free(ry2); Free(rx2);
	Free(z2); Free(y2); Free(x2);

	Free(gamma);Free(beta);Free(alpha);
	Free(vz);Free(vy);Free(vx);
	Free(rz); Free(ry); Free(rx);

	Free(z); Free(y); Free(x);
	Free(istep);Free(iframe);Free(jframe);
	return nframe;
}
