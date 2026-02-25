#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "pmheader.h"

#define GetCIC(mtype,den,zstart,nz_per) do {\
	for(i=0;i<lvsimpar.##mtype##.np;i++) { \
		float xp,yp,zp,pmas;\
		pmas = MofP(bp+i,lvsimpar);\
		xp = XofP(bp+i)- 0.5L;\
		yp = YofP(bp+i)- 0.5L;\
		zp = ZofP(bp+i)- 0.5L - zstart;\
		int nearx = rint(xp);\
		int neary = rint(yp);\
		int nearz = rint(zp);\
		int ic,jc,kc;\
		ic = (nearx + ng)%ng;\
		jc = (neary + ng)%ng;\
		kc = (nearz + nz_per)%nz_per;\
		xmin = xp - nearx;\
		ymin = yp - neary;\
		zmin = zp - nearz;\
		int icc = (ng+nearx + (int)(copysign(1.,xmin))) %ng;\
		int jcc = (ng+neary + (int)(copysign(1.,ymin))) %ng;\
		int kcc = (ng+nearz + (int)(copysign(1.,zmin))) %ng;\
		xmin = fabs(xmin);\
		ymin = fabs(ymin);\
		zmin = fabs(zmin);\
		xmax = (1.-xmin)*pmas;\
		xmin = xmin*pmas;\
		ymax = 1. -ymin;\
		zmax = 1. -zmin;\
		den[ic+nx*(jc+ny*kc)] += xmax*ymax*zmax;\
		den[icc+nx*(jc+ny*kc)] += xmin*ymax*zmax;\
		den[ic+nx*(jcc+ny*kc)] += xmax*ymin*zmax;\
		den[icc+nx*(jcc+ny*kc)] += xmin*ymin*zmax;\
		den[ic+nx*(jc+ny*kcc)] += xmax*ymax*zmin;\
		den[icc+nx*(jc+ny*kcc)] += xmin*ymax*zmin;\
		den[ic+nx*(jcc+ny*kcc)] += xmax*ymin*zmin;\
		den[icc+nx*(jcc+ny*kcc)] += xmin*ymin*zmin;\
	}\
}while(0)

void cic(lvsimpar){
	float xmin,ymin,zmin;
	int ng;
	float *den;

	ng = lvsimpar.nx;

	lvsimpar.mass = lvsimpar.nspace*lvsimpar.nspace*lvsimpar.nspace;

}

