#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "near.h"
#include "Memory.h"
#include "defs.h"
#define MIN(a,b) (a)<(b)? (a): (b);
#define MAX(a,b) (a)>(b)? (a): (b);
particle *p;
float x,y,z;
double f(double);
float mindist;

/* inverse of pi */
#define AMPW41 (-0.2387324)
#define AMPW42 (0.07957747)
inline float W4(float dist,float h){
	float twomr,r;
	r = dist/h;
	if(r<1.){
		/*
		return 1/PI*(1.-1.5*dist*dist+0.75*dist*dist*dist);
		*/
		return AMPW41*((2.-r)*r*r-1.333333)/(h*h*h);
	}
	else if(r < 2.){
		twomr = (2.-r)/h;
		/*
		return 1/PI*0.25*twomr*twomr*twomr;
		*/
		return AMPW42*twomr*twomr*twomr;
	}
	else {
		return 0.;
	}
}
#undef AMPW41
#undef AMPW42
void findsphdensity(vector3d *bp,int np,int *nearindex, int Numnear, 
		float *densph){
	float *h;
	double std,mean;
	int ntmp;
	float tmpx,tmpy,tmpz;
	float *dist2;
	float fplmf,ptlmass;
	int i,j,k;
	int N,M;
	TPtlStruct *ptl;
	BeginEndTree beginend;
	Box box;
	TStruct *TREE;
	float theta = 1.;
	float wtime;
	long iseed=-9;
	float x0,y0,z0,pscale;
	float Rmin,Rmax,size;
	int index, *tmpindx;

	tmpindx = (int *)Malloc(sizeof(int)*Numnear,PPTR(tmpindx));

	p = (particle *) Malloc(sizeof(particle)*np,PPTR(p));
	for(i=0;i<np;i++){
		p[i].x = bp[i].x;
		p[i].y = bp[i].y;
		p[i].z = bp[i].z;
	}
	ptl = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
	for(i=0;i<np;i++){
		ptl[i].type = TYPE_PTL;
		ptl[i].r[0] = bp[i].x;
		ptl[i].r[1] = bp[i].y;
		ptl[i].r[2] = bp[i].z;
		ptl[i].sibling = ptl+(i+1);
		ptl[i].mass = 1;
	}
	ptl[np-1].sibling = NULL;
	TREE = (TStruct *) Malloc(sizeof(TStruct)*np,PPTR(TREE));
	h = (float*) Malloc(sizeof(float)*np,PPTR(h));
	dist2 = (float*) Malloc(sizeof(float)*np*Numnear,PPTR(dist2));

	/*
	(void) WALLCLOCK();
	*/
	box = findsphbox(ptl,np);
#ifdef DEBUG
	printf("after findphsbox rmin =%g %g %g with width=%g\n",box.x,box.y,box.z,box.width);
#endif
	Make_Tree_Near(TREE,ptl,np,box);
#ifdef DEBUG
	printf("after Make_Tree_near\n");
#endif
	/*
	printf("CPU(Make_Tree_test)=%g\n",WALLCLOCK());
	*/
	mindist = 2.e25;
	for(i=0;i<np;i++){
		int res;
		float neardist;
		res = Find_Near(p+i,Numnear,TREE,ptl,&neardist,tmpindx,dist2+i*Numnear);
		h[i] = neardist;
		for(j=0;j<Numnear;j++){
			nearindex[i*Numnear+j] = tmpindx[j];
		}
	}
	for(i=0;i<np*Numnear;i++) dist2[i] = sqrt(dist2[i]);
#ifdef DEBUG
	printf("after Find_near\n");
#endif

	k = 0;
	for(i=0;i<np;i++){
		densph[i] = 0.;
		for(j=0;j<Numnear;j++){
			densph[i] += W4(dist2[k],h[i]);
			k++;
		}
	}
	Free(dist2); Free(h);
	Free(TREE);Free(ptl);Free(p);Free(tmpindx);
}
void DenPeakCenter(vector3d *bp,int np,double *cx,double *cy,double *cz,int Numnear){
	int i,j,k;
	int *nearindex;
	float *densph;
	float pden;
	Numnear = MIN(Numnear,np);
	nearindex = (int *)Malloc(sizeof(int)*np*Numnear,PPTR(nearindex));
	densph = (float*)Malloc(sizeof(float)*np,PPTR(densph));
	findsphdensity(bp,np,nearindex,Numnear,densph);
	*cx = *cy = *cz = 0.L;
	pden = 0.L;
	for(i=0;i<np;i++){
		if(densph[i]> pden){
			j  =i;
			pden = densph[i];
		}
	}
	*cx = bp[j].x; *cy = bp[j].y; *cz = bp[j].z;
	Free(densph);
	Free(nearindex);
}
