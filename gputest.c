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


int N;
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
	long iseed=-29;
	particle *pos;
	TPtlStruct *grv;
	float *acc1,*acc2,*acc3,*acc4;
	float theta=0.5;



	if(Make_Total_Memory()==0){
		fprintf(stderr,"Error in making memory space \n");
		exit(99);
	}

	if(initflag == 0){
		i_force_spline(nx,rsphere);
		initflag = 1;
	}
	N = atoi(argv[1]);
	theta = atof(argv[2]);

	npos = N*27;
	ngrv = N;
	pos = (particle *)Malloc(sizeof(particle)*npos,PPTR(pos));
	grv = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*ngrv,PPTR(grv));
	acc1 = (float *) Malloc(sizeof(float)*npos,PPTR(acc1));
	acc2 = (float *) Malloc(sizeof(float)*npos,PPTR(acc2));
	acc3 = (float *) Malloc(sizeof(float)*npos,PPTR(acc3));
	acc4 = (float *) Malloc(sizeof(float)*npos,PPTR(acc4));
	for(i=0;i<ngrv;i++){
		grv[i].r[0] = pos[i].x = ran3(&iseed)*4;
		grv[i].r[1] = pos[i].y = ran3(&iseed)*4;
		grv[i].r[2] = pos[i].z = ran3(&iseed)*4;
		grv[i].type = TYPE_PTL;
		grv[i].mass = 1;

	}
	i = ngrv;
	while(i<npos){
		do {
			pos[i].x = (ran3(&iseed)*12-4);
		} while(pos[i].x > 0 && pos[i].x < 4);
		do {
			pos[i].y = (ran3(&iseed)*12-4);
		} while(pos[i].y > 0 && pos[i].y < 4);
		do {
			pos[i].z = (ran3(&iseed)*12-4);
		} while(pos[i].z > 0 && pos[i].z < 4);
		i++;
	}

	{ /* This is the CPU direct Nbody */
		TIMER_START(45);
		DirectSummation(ngrv,grv,npos,pos);
		TIMER_STOP(45);
		printf("CPU TIME for Pair N-body %g\n",ELAPSED_TIME(45));
		for(i=0;i<npos;i++){
			acc1[i] = SQRT3(pos[i].ax,pos[i].ay,pos[i].az);
		}
	}
	Direct_Nbody(pos,npos,grv,ngrv,GpusPerNode,myid,nx);/* This is for initialization of the texture */
	{ /* This the GPU direct N-body */
		TIMER_START(45);
		Direct_Nbody(pos,npos,grv,ngrv,GpusPerNode,myid,nx);
		TIMER_STOP(45);
		printf("GPU TIME for Pair N-body %g\n",ELAPSED_TIME(45));
		for(i=0;i<npos;i++) acc2[i] = SQRT3(pos[i].ax,pos[i].ay,pos[i].az);
		for(i=0;i<npos;i++) pos[i].ax = pos[i].ay = pos[i].az = 0;
	}
	{/* This is the GPU cell force */
		int ncell;
		Box box;
		TStruct *TREECELL;
		Box CheckBoundingBox(long, TPtlStruct *);
		void gputreeforce(particle *, int ,TStruct *, int, TPtlStruct *, int, float,int,int,int);
		void gputreeforce2(particle *, int ,TStruct *, int, TPtlStruct *, int, float,int,int,int);
		TREECELL = (TStruct *) Malloc(sizeof(TStruct)*ngrv,PPTR(TREECELL));

		TIMER_START(45);
		box.x = box.y = box.z = 0; box.width = 4.;
		ncell = Make_Tree(TREECELL,grv,ngrv,box,theta);
		gputreeforce(pos,npos,TREECELL,ncell,grv,ngrv,rspheresq,1,myid,nx);
		TIMER_STOP(45);
		printf("GPU TIME for Tree N-body 1: %g\n",ELAPSED_TIME(45));

		TIMER_START(45);
		box.x = box.y = box.z = 0; box.width = 4.;
		ncell = Make_Tree(TREECELL,grv,ngrv,box,theta);
		gputreeforce2(pos,npos,TREECELL,ncell,grv,ngrv,rspheresq,1,myid,nx);
		TIMER_STOP(45);
		printf("GPU TIME for Tree N-body 2: %g\n",ELAPSED_TIME(45));


		Free(TREECELL);
		for(i=0;i<npos;i++) acc3[i] = SQRT3(pos[i].ax,pos[i].ay,pos[i].az);
	}

	{
		Box box;
		TStruct *TREECELL;
		TREECELL = (TStruct *) Malloc(sizeof(TStruct)*ngrv,PPTR(TREECELL));
		float theta2 = theta*theta;
		TIMER_START(45);
		box.x = box.y = box.z = 0;
		box.width = 4.;
		Make_Tree(TREECELL,grv,ngrv,box,theta);
		for(i=0;i<npos;i++){
			treeforce(pos+i,theta2,TREECELL,grv,rspheresq);
		}
		TIMER_STOP(45);
		for(i=0;i<npos;i++) acc4[i] = SQRT3(pos[i].ax,pos[i].ay,pos[i].az);
		printf("CPU TIME for Tree N-body %g\n",ELAPSED_TIME(45));
		for(i=0;i<npos;i++) acc4[i] = SQRT3(pos[i].ax,pos[i].ay,pos[i].az);

	}
	for(i=0;i<npos;i++){
		float amp1,amp2,amp3;
		amp1 = (acc2[i]-acc1[i])/acc1[i]*100;
		amp2 = (acc3[i]-acc1[i])/acc1[i]*100;
		amp3 = (acc4[i]-acc1[i])/acc1[i]*100;
		if(isnan(amp2) || isinf(amp2) || acc2[i]==0) {
			printf("+P%d has GPU accel %g but has CPU acc = %g \n",i,acc2[i],acc1[i]);
		}
		else {
			printf("+P%d has acc = %g %g || %g %g : amp = %g %g | %g\n",i,acc3[i],acc2[i],acc4[i],acc1[i],amp2,amp1,amp3);
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

#define min(a,b) (a) > (b) ? (b): (a)
#define max(a,b) (a) < (b) ? (b): (a)
Box CheckBoundingBox(long ngrv,TPtlStruct *grv){
    float xmin,xmax,ymin,ymax,zmin,zmax;
    float width;
    long i;
    int si;
    Box box;
    xmin = ymin = zmin = 1.e25;
    xmax = ymax = zmax = -1.e25;
    for(i=0;i<ngrv;i++){
        xmin = min(xmin,grv[i].r[0]);
        ymin = min(ymin,grv[i].r[1]);
        zmin = min(zmin,grv[i].r[2]);
        xmax = max(xmax,grv[i].r[0]);
        ymax = max(ymax,grv[i].r[1]);
        zmax = max(zmax,grv[i].r[2]);
    }
    width = xmax-xmin;
    width = max(width,ymax-ymin);
    width = max(width,zmax-zmin);
    if(width !=0){
        box.x = xmin-width*1.e-6;
        box.y = ymin-width*1.e-6;
        box.z = zmin-width*1.e-6;
        box.width = width*(1.+2.e-6);
    }
    else {
        box.x = xmin-1.e-5;
        box.y = ymin-1.e-5;
        box.z = zmin-1.e-5;
        box.width = 2.e-5;
    }
    return box;
}
