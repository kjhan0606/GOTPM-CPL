#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>

#include "pmheader.h"
#include "Memory.h"
#include "cutil.h"
#include "cuda.h"
#include "nvgpu.h"
#define pi (3.14159265354979L)

#define den(i,j,k) (den[(i)+mx*((j)+ny*(k))])
#define PTR(den,i,j,k) (den+(i)+mx*((j)+ny*(k)))

#define MAX(a,b) ((a)>(b)?(a):(b))
#define icopysign(a) copysignf(1.,a);

#ifdef XYZDBL
#define FABS(a) fabs(a)
#define RINT(a) rint(a)
#else
#define FABS(a) fabsf(a)
#define RINT(a) rintf(a)
#endif

#define nxp(a) ((a+nx)%nx)
#define nyp(a) ((a+ny)%ny)
#define nzp(a) ((a+nz_per)%nz_per)

#define xabm1(a) ((a-1+nx)%nx)
#define xabp1(a) ((a+1+nx)%nx)
#define yabm1(a) ((a-1+ny)%ny)
#define yabp1(a) ((a+1+ny)%ny)
#define zabm1(a) ((a-1+nz_per)%nz_per)
#define zabp1(a) ((a+1+nz_per)%nz_per)



__global__ void Get_Fda4(float *den,int mx,int nx,int ny, int local_nz,pmparticletype *pmp, int np,
		float start_z, float fact2){
	int idx,nz_per;
	float ffact1,ffact2;


	ffact1 = (4./3.)/2.;
	ffact2 = (-1./3.)/4.; 
	nz_per = local_nz + 7;

/*
	idx = (gridDim.x*blockIdx.y + blockIdx.x)*blockDim.x + threadIdx.x;
	*/
	idx = blockIdx.x*blockDim.x + threadIdx.x;
	if(idx<np){
#ifdef XYZDBL
		double xp,yp,zp;
#else
		float xp,yp,zp;
#endif
		float3 min;
		float fx1,fx2,fx3,fx4,fx5,fx6,fx7,fx8,fx9,fx10;
		float fx11,fx12,fx13,fx14,fx15,fx16,fx17,fx18,fx19,fx20;
		float fx21,fx22,fx23,fx24,fx25,fx26,fx27;
		float fy1,fy2,fy3,fy4,fy5,fy6,fy7,fy8,fy9,fy10;
		float fy11,fy12,fy13,fy14,fy15,fy16,fy17,fy18,fy19,fy20;
		float fy21,fy22,fy23,fy24,fy25,fy26,fy27;
		float fz1,fz2,fz3,fz4,fz5,fz6,fz7,fz8,fz9,fz10;
		float fz11,fz12,fz13,fz14,fz15,fz16,fz17,fz18,fz19,fz20;
		float fz21,fz22,fz23,fz24,fz25,fz26,fz27;
		int3 near,sign;
		int i1,j1,k1,i2,j2,k2,i3,j3,k3;
		int i1m1,i1m2,i1p1,i1p2,i2m1,i2m2,i2p1,i2p2,i3m1,i3m2,i3p1;
		int i3p2,j1m1,j1m2,j1p1,j1p2,j2m1,j2m2,j2p1,j2p2,j3m1,j3m2;
		int j3p1,j3p2,k1m1,k1m2,k1p1,k1p2,k2m1,k2m2,k2p1,k2p2,k3m1;
		int k3m2,k3p1,k3p2;

		xp = XofP(pmp+idx);
		yp = YofP(pmp+idx);
		zp = ZofP(pmp+idx)-start_z;


		near.x = RINT(xp);
		near.y = RINT(yp);
		near.z = RINT(zp);
		min.x = xp - near.x;
		min.y = yp - near.y;
		min.z = zp - near.z;

		sign.x = icopysign(min.x);
		sign.y = icopysign(min.y);
		sign.z = icopysign(min.z);

		i1 = nxp(near.x);
		i2 = nxp(near.x+sign.x);
		i3 = nxp(near.x-sign.x);

		j1 = nyp(near.y);
		j2 = nyp(near.y+sign.y);
		j3 = nyp(near.y-sign.y);

		k1 = nzp(near.z);
		k2 = nzp(near.z+sign.z);
		k3 = nzp(near.z-sign.z);




		i1m1 = xabm1(i1);
		i1m2 = xabm1(i1m1);
		i1p1 = xabp1(i1);
		i1p2 = xabp1(i1p1);
		i2m1 = xabm1(i2);
		i2m2 = xabm1(i2m1);
		i2p1 = xabp1(i2);
		i2p2 = xabp1(i2p1);
		i3m1 = xabm1(i3);
		i3m2 = xabm1(i3m1);
		i3p1 = xabp1(i3);
		i3p2 = xabp1(i3p1);

		j1m1 = yabm1(j1);
		j1m2 = yabm1(j1m1);
		j1p1 = yabp1(j1);
		j1p2 = yabp1(j1p1);
		j2m1 = yabm1(j2);
		j2m2 = yabm1(j2m1);
		j2p1 = yabp1(j2);
		j2p2 = yabp1(j2p1);
		j3m1 = yabm1(j3);
		j3m2 = yabm1(j3m1);
		j3p1 = yabp1(j3);
		j3p2 = yabp1(j3p1);

		k1m1 = zabm1(k1);
		k1m2 = zabm1(k1m1);
		k1p1 = zabp1(k1);
		k1p2 = zabp1(k1p1);
		k2m1 = zabm1(k2);
		k2m2 = zabm1(k2m1);
		k2p1 = zabp1(k2);
		k2p2 = zabp1(k2p1);
		k3m1 = zabm1(k3);
		k3m2 = zabm1(k3m1);
		k3p1 = zabp1(k3);
		k3p2 = zabp1(k3p1);
        {
            long jk1,jk2,jk3;
            jk1 = mx*(j1+ny*k1);
            fx1 = ffact1*(den[i1m1+jk1]-den[i1p1+jk1])+ffact2*(den[i1m2+jk1]-den[i1p2+jk1]);
            fx2 = ffact1*(den[i2m1+jk1]-den[i2p1+jk1])+ffact2*(den[i2m2+jk1]-den[i2p2+jk1]);
            fx3 = ffact1*(den[i3m1+jk1]-den[i3p1+jk1])+ffact2*(den[i3m2+jk1]-den[i3p2+jk1]);
            jk1 = mx*(j2+ny*k1);
            fx4 = ffact1*(den[i1m1+jk1]-den[i1p1+jk1])+ffact2*(den[i1m2+jk1]-den[i1p2+jk1]);
            fx5 = ffact1*(den[i2m1+jk1]-den[i2p1+jk1])+ffact2*(den[i2m2+jk1]-den[i2p2+jk1]);
            fx6 = ffact1*(den[i3m1+jk1]-den[i3p1+jk1])+ffact2*(den[i3m2+jk1]-den[i3p2+jk1]);
            jk1 = mx*(j3+ny*k1);
            fx7 = ffact1*(den[i1m1+jk1]-den[i1p1+jk1])+ffact2*(den[i1m2+jk1]-den[i1p2+jk1]);
            fx8 = ffact1*(den[i2m1+jk1]-den[i2p1+jk1])+ffact2*(den[i2m2+jk1]-den[i2p2+jk1]);
            fx9 = ffact1*(den[i3m1+jk1]-den[i3p1+jk1])+ffact2*(den[i3m2+jk1]-den[i3p2+jk1]);
            jk2 = mx*(j1+ny*k2);
            fx10 = ffact1*(den[i1m1+jk2]-den[i1p1+jk2])+ffact2*(den[i1m2+jk2]-den[i1p2+jk2]);
            fx11 = ffact1*(den[i2m1+jk2]-den[i2p1+jk2])+ffact2*(den[i2m2+jk2]-den[i2p2+jk2]);
            fx12 = ffact1*(den[i3m1+jk2]-den[i3p1+jk2])+ffact2*(den[i3m2+jk2]-den[i3p2+jk2]);
            jk2 = mx*(j2+ny*k2);
            fx13 = ffact1*(den[i1m1+jk2]-den[i1p1+jk2])+ffact2*(den[i1m2+jk2]-den[i1p2+jk2]);
            fx14 = ffact1*(den[i2m1+jk2]-den[i2p1+jk2])+ffact2*(den[i2m2+jk2]-den[i2p2+jk2]);
            fx15 = ffact1*(den[i3m1+jk2]-den[i3p1+jk2])+ffact2*(den[i3m2+jk2]-den[i3p2+jk2]);
            jk2 = mx*(j3+ny*k2);
            fx16 = ffact1*(den[i1m1+jk2]-den[i1p1+jk2])+ffact2*(den[i1m2+jk2]-den[i1p2+jk2]);
            fx17 = ffact1*(den[i2m1+jk2]-den[i2p1+jk2])+ffact2*(den[i2m2+jk2]-den[i2p2+jk2]);
            fx18 = ffact1*(den[i3m1+jk2]-den[i3p1+jk2])+ffact2*(den[i3m2+jk2]-den[i3p2+jk2]);
            jk3 = mx*(j1+ny*k3);
            fx19 = ffact1*(den[i1m1+jk3]-den[i1p1+jk3])+ffact2*(den[i1m2+jk3]-den[i1p2+jk3]);
            fx20 = ffact1*(den[i2m1+jk3]-den[i2p1+jk3])+ffact2*(den[i2m2+jk3]-den[i2p2+jk3]);
            fx21 = ffact1*(den[i3m1+jk3]-den[i3p1+jk3])+ffact2*(den[i3m2+jk3]-den[i3p2+jk3]);
            jk3 = mx*(j2+ny*k3);
            fx22 = ffact1*(den[i1m1+jk3]-den[i1p1+jk3])+ffact2*(den[i1m2+jk3]-den[i1p2+jk3]);
            fx23 = ffact1*(den[i2m1+jk3]-den[i2p1+jk3])+ffact2*(den[i2m2+jk3]-den[i2p2+jk3]);
            fx24 = ffact1*(den[i3m1+jk3]-den[i3p1+jk3])+ffact2*(den[i3m2+jk3]-den[i3p2+jk3]);
            jk3 = mx*(j3+ny*k3);
            fx25 = ffact1*(den[i1m1+jk3]-den[i1p1+jk3])+ffact2*(den[i1m2+jk3]-den[i1p2+jk3]);
            fx26 = ffact1*(den[i2m1+jk3]-den[i2p1+jk3])+ffact2*(den[i2m2+jk3]-den[i2p2+jk3]);
            fx27 = ffact1*(den[i3m1+jk3]-den[i3p1+jk3])+ffact2*(den[i3m2+jk3]-den[i3p2+jk3]);
        }
        {
            long jk1,jk2,jk3,jk4;
            jk1 = mx*(j1m1+ny*k1);
            jk2 = mx*(j1p1+ny*k1);
            jk3 = mx*(j1m2+ny*k1);
            jk4 = mx*(j1p2+ny*k1);
            fy1 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fy2 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fy3 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j2m1+ny*k1);
            jk2 = mx*(j2p1+ny*k1);
            jk3 = mx*(j2m2+ny*k1);
            jk4 = mx*(j2p2+ny*k1);
            fy4 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fy5 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fy6 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j3m1+ny*k1);
            jk2 = mx*(j3p1+ny*k1);
            jk3 = mx*(j3m2+ny*k1);
            jk4 = mx*(j3p2+ny*k1);
            fy7 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fy8 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fy9 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j1m1+ny*k2);
            jk2 = mx*(j1p1+ny*k2);
            jk3 = mx*(j1m2+ny*k2);
            jk4 = mx*(j1p2+ny*k2);
            fy10 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fy11 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fy12 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j2m1+ny*k2);
            jk2 = mx*(j2p1+ny*k2);
            jk3 = mx*(j2m2+ny*k2);
            jk4 = mx*(j2p2+ny*k2);
            fy13 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fy14 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fy15 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j3m1+ny*k2);
            jk2 = mx*(j3p1+ny*k2);
            jk3 = mx*(j3m2+ny*k2);
            jk4 = mx*(j3p2+ny*k2);
            fy16 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fy17 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fy18 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j1m1+ny*k3);
            jk2 = mx*(j1p1+ny*k3);
            jk3 = mx*(j1m2+ny*k3);
            jk4 = mx*(j1p2+ny*k3);
            fy19 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fy20 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fy21 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j2m1+ny*k3);
            jk2 = mx*(j2p1+ny*k3);
            jk3 = mx*(j2m2+ny*k3);
            jk4 = mx*(j2p2+ny*k3);
            fy22 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fy23 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fy24 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j3m1+ny*k3);
            jk2 = mx*(j3p1+ny*k3);
            jk3 = mx*(j3m2+ny*k3);
            jk4 = mx*(j3p2+ny*k3);
            fy25 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fy26 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fy27 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
        }
        {
            long jk1,jk2,jk3,jk4;
            jk1 = mx*(j1+ny*k1m1);
            jk2 = mx*(j1+ny*k1p1);
            jk3 = mx*(j1+ny*k1m2);
            jk4 = mx*(j1+ny*k1p2);
            fz1 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fz2 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fz3 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j2+ny*k1m1);
            jk2 = mx*(j2+ny*k1p1);
            jk3 = mx*(j2+ny*k1m2);
            jk4 = mx*(j2+ny*k1p2);
            fz4 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fz5 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fz6 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j3+ny*k1m1);
            jk2 = mx*(j3+ny*k1p1);
            jk3 = mx*(j3+ny*k1m2);
            jk4 = mx*(j3+ny*k1p2);
            fz7 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fz8 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fz9 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j1+ny*k2m1);
            jk2 = mx*(j1+ny*k2p1);
            jk3 = mx*(j1+ny*k2m2);
            jk4 = mx*(j1+ny*k2p2);
            fz10 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fz11 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fz12 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j2+ny*k2m1);
            jk2 = mx*(j2+ny*k2p1);
            jk3 = mx*(j2+ny*k2m2);
            jk4 = mx*(j2+ny*k2p2);
            fz13 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fz14 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fz15 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j3+ny*k2m1);
            jk2 = mx*(j3+ny*k2p1);
            jk3 = mx*(j3+ny*k2m2);
            jk4 = mx*(j3+ny*k2p2);
            fz16 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fz17 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fz18 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j1+ny*k3m1);
            jk2 = mx*(j1+ny*k3p1);
            jk3 = mx*(j1+ny*k3m2);
            jk4 = mx*(j1+ny*k3p2);
            fz19 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fz20 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fz21 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j2+ny*k3m1);
            jk2 = mx*(j2+ny*k3p1);
            jk3 = mx*(j2+ny*k3m2);
            jk4 = mx*(j2+ny*k3p2);
            fz22 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fz23 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fz24 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
            jk1 = mx*(j3+ny*k3m1);
            jk2 = mx*(j3+ny*k3p1);
            jk3 = mx*(j3+ny*k3m2);
            jk4 = mx*(j3+ny*k3p2);
            fz25 = ffact1*(den[i1+jk1]-den[i1+jk2])+ffact2*(den[i1+jk3]-den[i1+jk4]);
            fz26 = ffact1*(den[i2+jk1]-den[i2+jk2])+ffact2*(den[i2+jk3]-den[i2+jk4]);
            fz27 = ffact1*(den[i3+jk1]-den[i3+jk2])+ffact2*(den[i3+jk3]-den[i3+jk4]);
        }
		{
			float wx1wy1, wx2wy1, wx3wy1, wx1wy2, wx2wy2, wx3wy2, wx1wy3, wx2wy3, wx3wy3;
			float3 rd1,wr1,wr2,wr3;
			float wt1,wt2,wt3,wt4,wt5,wt6,wt7,wt8,wt9,wt10;
			float wt11,wt12,wt13,wt14,wt15,wt16,wt17,wt18,wt19,wt20;
			float wt21,wt22,wt23,wt24,wt25,wt26,wt27;

			rd1.x = FABS(min.x);
			rd1.y = FABS(min.y);
			rd1.z = FABS(min.z);
	
	
			wr1.x =  0.75-rd1.x*rd1.x;
			wr1.y =  0.75-rd1.y*rd1.y;
			wr1.z =  0.75-rd1.z*rd1.z;
	
			wr3.x = 0.5*(0.25+rd1.x*(rd1.x-1.));
			wr2.x = wr3.x + rd1.x;
	
			wr3.y = 0.5*(0.25+rd1.y*(rd1.y-1.));
			wr2.y = wr3.y + rd1.y;
	
			wr3.z = 0.5*(0.25+rd1.z*(rd1.z-1.));
			wr2.z = wr3.z + rd1.z;

			
			wx1wy1 = wr1.x*wr1.y;
			wx2wy1 = wr2.x*wr1.y;
			wx3wy1 = wr3.x*wr1.y;
			wx1wy2 = wr1.x*wr2.y;
			wx2wy2 = wr2.x*wr2.y;
			wx3wy2 = wr3.x*wr2.y;
			wx1wy3 = wr1.x*wr3.y;
			wx2wy3 = wr2.x*wr3.y;
			wx3wy3 = wr3.x*wr3.y;
	
			wt1 = wx1wy1*wr1.z;
			wt2 = wx2wy1*wr1.z;
			wt3 = wx3wy1*wr1.z;
			wt4 = wx1wy2*wr1.z;
			wt5 = wx2wy2*wr1.z;
			wt6 = wx3wy2*wr1.z;
			wt7 = wx1wy3*wr1.z;
			wt8 = wx2wy3*wr1.z;
			wt9 = wx3wy3*wr1.z;
	
			wt10 = wx1wy1*wr2.z;
			wt11 = wx2wy1*wr2.z;
			wt12 = wx3wy1*wr2.z;
			wt13 = wx1wy2*wr2.z;
			wt14 = wx2wy2*wr2.z;
			wt15 = wx3wy2*wr2.z;
			wt16 = wx1wy3*wr2.z;
			wt17 = wx2wy3*wr2.z;
			wt18 = wx3wy3*wr2.z;
	
			wt19 = wx1wy1*wr3.z;
			wt20 = wx2wy1*wr3.z;
			wt21 = wx3wy1*wr3.z;
			wt22 = wx1wy2*wr3.z;
			wt23 = wx2wy2*wr3.z;
			wt24 = wx3wy2*wr3.z;
			wt25 = wx1wy3*wr3.z;
			wt26 = wx2wy3*wr3.z;
			wt27 = wx3wy3*wr3.z;
	
			float fx,fy,fz;

			fx=(wt1*fx1+wt2*fx2+wt3*fx3+wt4*fx4  +wt5 *fx5 +wt6 *fx6+
					wt7 *fx7 +wt8 *fx8 +wt9 *fx9 +wt10*fx10+wt11*fx11+wt12*fx12+
					wt13*fx13+wt14*fx14+wt15*fx15+wt16*fx16+wt17*fx17+wt18*fx18+
					wt19*fx19+wt20*fx20+wt21*fx21+wt22*fx22+wt23*fx23+wt24*fx24+
					wt25*fx25+wt26*fx26+wt27*fx27);
			fy=(wt1*fy1+wt2*fy2+wt3*fy3+wt4*fy4  +wt5 *fy5 +wt6 *fy6+
					wt7 *fy7 +wt8 *fy8 +wt9 *fy9 +wt10*fy10+wt11*fy11+wt12*fy12+
					wt13*fy13+wt14*fy14+wt15*fy15+wt16*fy16+wt17*fy17+wt18*fy18+
					wt19*fy19+wt20*fy20+wt21*fy21+wt22*fy22+wt23*fy23+wt24*fy24+
					wt25*fy25+wt26*fy26+wt27*fy27);
			fz=(wt1*fz1+wt2*fz2+wt3*fz3+wt4*fz4  +wt5 *fz5 +wt6 *fz6+
					wt7 *fz7 +wt8 *fz8 +wt9 *fz9 +wt10*fz10+wt11*fz11+wt12*fz12+
					wt13*fz13+wt14*fz14+wt15*fz15+wt16*fz16+wt17*fz17+wt18*fz18+
					wt19*fz19+wt20*fz20+wt21*fz21+wt22*fz22+wt23*fz23+wt24*fz24+
					wt25*fz25+wt26*fz26+wt27*fz27);
			pmp[idx].vx += fact2*fx;
			pmp[idx].vy += fact2*fy;
			pmp[idx].vz += fact2*fz;
		}
	}
#ifdef __DEVICE_EMULATION__
    __syncthreads();
#endif
}

#define MIN(a,b) (a)<(b)? (a):(b)
void cuda_setting(int, int);

extern "C" void call_cuda_fda4(int myid, int nx, int mx, int ny, int local_nz, 
		float *den, pmparticletype *pmp, int np, float start_z, float fact2){
	int i;
	float *devden;
	extern cudaDeviceProp mydev;
	pmparticletype *devpmp;
	cuda_setting(myid,GPUSPERNODE);
	CUDA_SAFE_CALL(cudaMalloc(PPTR(devden),sizeof(float)*mx*ny*(local_nz+7)));
	CUDA_SAFE_CALL(cudaMemcpy(devden, den, mx*ny*(local_nz+7)*sizeof(float), cudaMemcpyHostToDevice));



	int blocksize = BlockSize;
    int numthreads = (mydev.maxGridSize[0]/2)*blocksize;
	if(numthreads < 1024) {
		fprintf(stderr,"Error in the maxGridSize[0] = %d\n",mydev.maxGridSize[0]);
		exit(99);
	}
	int nsplit = (np+numthreads-1)/numthreads;
	CUDA_SAFE_CALL(cudaMalloc(PPTR(devpmp),sizeof(pmparticletype)*numthreads));
	for(i=0;i<nsplit;i++){
		int mp = MIN(numthreads,(np-i*numthreads));
		int nblock = (mp+blocksize-1)/blocksize;
		CUDA_SAFE_CALL(cudaMemcpy(devpmp, &(pmp[i*numthreads]), 
				mp*sizeof(pmparticletype), cudaMemcpyHostToDevice));
		Get_Fda4<<<nblock,blocksize>>>(devden,mx,nx,ny,local_nz,devpmp,mp,start_z,fact2);
		CUDA_SAFE_CALL(cudaMemcpy(&(pmp[i*numthreads]),devpmp, mp*sizeof(pmparticletype), 
				cudaMemcpyDeviceToHost));
	}
	cudaThreadSynchronize();
	if(cudaGetLastError() != cudaSuccess){
		printf("Error in the cudaFDA4: %s\n",cudaGetErrorString(cudaGetLastError()));
	}

	CUDA_SAFE_CALL(cudaFree(devpmp));
	CUDA_SAFE_CALL(cudaFree(devden));

	/*
	for(i=0;i<np;i++){
		float xd,yd,zd,rd;
		xd = pmp[i].vx*simpar.pfact;
		yd = pmp[i].vy*simpar.pfact;
		zd = pmp[i].vz*simpar.pfact;
		rd = xd*xd+yd*yd+zd*zd;
		rd = sqrt(rd);
		if(rd > 1.) printf("PPP%d had %g %g %g : %ld\n",i,pmp[i].vx,pmp[i].vy,pmp[i].vz,pmp[i].indx);
	}
	*/


}
#undef MIN

