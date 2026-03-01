#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#define DEFINE_SIM_PARA
#include "pmheader.h"
#undef DEFINE_SIM_PARA
#include "kjhtree.h"
#include "force_spline.h"
#include "Memory.h"
#include<sys/times.h>

void Direct_Nbody(particle *, int, TPtlStruct *, int, int, int, int);
void DirectSummation(long ,TPtlStruct *,long ,particle *);
#define SQRT3(a,b,c) (sqrt(a*a+b*b+c*c))

#define N 40

float ran3(long *);
int npos,ngrv;
int nx=1024;
int GpusPerNode=1;
int myid=0;

int initflag=0;
#define rsphere 6.
#define rspheresq (rsphere*rsphere)


int main(int argc, char **argv){
	int i,j,k;
	long iseed=-9;
	particle *pos;
	TPtlStruct *grv;
	float *acc1,*acc2;



	if(Make_Total_Memory()==0){
		fprintf(stderr,"Error in making memory space \n");
		exit(99);
	}

	if(initflag == 0){
		i_force_spline(nx,rsphere);
		initflag = 1;
	}

	npos = N*27;
	ngrv = N;
	pos = (particle *)Malloc(sizeof(particle)*npos,PPTR(pos));
	grv = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*ngrv,PPTR(grv));
	acc1 = (float *) Malloc(sizeof(float)*npos,PPTR(acc1));
	acc2 = (float *) Malloc(sizeof(float)*npos,PPTR(acc2));
	for(i=0;i<ngrv;i++){
		grv[i].r[0] = pos[i].x = ran3(&iseed)*4;
		grv[i].r[1] = pos[i].y = ran3(&iseed)*4;
		grv[i].r[2] = pos[i].z = ran3(&iseed)*4;
	}
	for(i=ngrv;i<npos;i++){
		pos[i].x = 4+ ran3(&iseed)*4;
		pos[i].y = 4+ ran3(&iseed)*4;
		pos[i].z = 4+ ran3(&iseed)*4;
	}
	TIMER_START(45);
	DirectSummation(ngrv,grv,npos,pos);
	TIMER_STOP(45);
	printf("CPU TIME for Pair N-body %g\n",ELAPSED_TIME(45));
	for(i=0;i<npos;i++){
		acc1[i] = SQRT3(pos[i].ax,pos[i].ay,pos[i].az);
		pos[i].ax = pos[i].ay = pos[i].az = 0.;
//		printf("-P%d has acc = %g %g %g\n",i,pos[i].ax,pos[i].ay,pos[i].az);
	}
//	Direct_Nbody(pos,npos,grv,ngrv,GpusPerNode,myid,nx);/* This is for initialization of the texture */
	TIMER_START(45);
	Direct_Nbody(pos,npos,grv,ngrv,GpusPerNode,myid,nx);
	TIMER_STOP(45);
	printf("GPU TIME for Pair N-body %g\n",ELAPSED_TIME(45));
	for(i=0;i<npos;i++){
		float amp;
		acc2[i] = SQRT3(pos[i].ax,pos[i].ay,pos[i].az);
		amp = (acc2[i]-acc1[i])/acc1[i]*100;
		if(isnan(amp) || isinf(amp) || acc2[i]==0) {
			printf("+P%d has GPU accel %g but has CPU acc = %g \n",i,acc2[i],acc1[i]);
		}
		else {
			printf("+P%d has acc = %g %g %g\n",i,acc2[i],acc1[i],amp);
		}
	}


}
void DirectSummation(long ngrv,TPtlStruct *grv,long npos,particle *pos){
    float tmpx,tmpy,tmpz,fplmf;
    float px,py,pz,dist2,ptlmass,dist;
	float ax,ay,az;
    long i,j,ntmp;

    for(i=0;i<npos;i++){
        px = pos[i].x; py = pos[i].y; pz = pos[i].z;
		ax = ay = az = 0;
        for(j=0;j<ngrv;j++){
            tmpx = px - grv[j].r[0];
            tmpy = py - grv[j].r[1];
            tmpz = pz - grv[j].r[2];
            dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
            if(dist2 <= rspheresq){
                dist = SQRT(dist2);
                ntmp = dist*invran2nran;
                fplmf = forcecorrectdiff(ntmp,0) + forcecorrectslope(ntmp,0) * (dist-ntmp*ran2nran);
				ax += tmpx*fplmf;
				ay += tmpy*fplmf;
				az += tmpz*fplmf;
            }
        }
        pos[i].ax = ax;
        pos[i].ay = ay;
        pos[i].az = az;
    }
}
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(long *idum)
{
    static int inext,inextp;
    static long ma[56];
    static int iff=0;
    long mj,mk;
    int i,ii,k;

    if (*idum < 0 || iff == 0) {
        iff=1;
        mj=MSEED-(*idum < 0 ? -*idum : *idum);
        mj %= MBIG;
        ma[55]=mj;
        mk=1;
        for (i=1;i<=54;i++) {
            ii=(21*i) % 55;
            ma[ii]=mk;
            mk=mj-mk;
            if (mk < MZ) mk += MBIG;
            mj=ma[ii];
        }
        for (k=1;k<=4;k++)
            for (i=1;i<=55;i++) {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
        inext=0;
        inextp=31;
        *idum=1;
    }
    if (++inext == 56) inext=1;
    if (++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext]=mj;
    return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */

