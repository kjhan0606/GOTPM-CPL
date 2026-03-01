/*
 * 이 코드는 찾아낸 헤일로의 peak density 를 구하는 것이다.
 * 여기에 헤일로의 물리적 성질을 첨가할것이다.
 */
/*
 * 지금 버젼은 trim_pden.c  를 수행해서 under dense peak 을 갖는 헤일로들을
 * 미리 제거하고 난 후를 가정하고 계산 한다. 22/01/2003
 * */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "Memory.h"
typedef struct PtlPos{
	float x,y,z;
} PtlPos;
typedef struct PtlVel{
	float vx,vy,vz;
} PtlVel;
#define MIN(A,B) ((A)<(B) ? (A): (B))
#define MAX(A,B) ((A)>(B) ? (A): (B))
#define sqr(A) ((A)*(A))
#define pow2(A) ((A)*(A))
#define pow3(A) ((A)*(A)*(A))
#define sqr(A) ((A)*(A))
double com2real,real2com;
double pntmass,r1kineticfact,r2kineticfact;
double potentfact;
double onesolarmass=1.989E33L;
double  G= 6.672E-8L;
double pc= 3.08567802E18L;
double pi= 3.1415926535L;
double rotation;
int nspace;
float r1[3],r2[3],r3[3],radius;
float a;
float amax;
float redshift;
float size;
float rng;
float omep;
float omeplam;
float hubble;
float omepk,Hsub;
int nx,ny,nz;
float *den;
float q,s;
float d[3],v[3][3];
double para_spin;
float tenergy,tmass,pvx;
float cx,cy,cz,cvx,cvy,cvz,pvr;

int main(int argc, char **argv){
	FILE *fp,*wp,*paraid,*mergerwp;
	PtlPos *r;
	PtlVel *vr;
    float peak;
	float *x,*y,*z,*vx,*vy,*vz;
	int tnend,nend,i,j,k;
	int ncount;
	int *indx,id;
    int mpeak,np;
    float xmin,ymin,xmax,ymax,width;
	int anend[10000];
	int haloindex;
	int index,peakindex;
	char infile[100],parafile[100];
	char wfile[100],mergerfile[100];
	int START,inid;
	long nsize;
	int halomergerindex;
	int *nowindex,*beforeindex;

	float tt1,tt2,tt3,tt4;
	int nindex,nmpeak,npeakindex;
	int nhalomergerindex,mtmp;
	int nnowhaloindex,nbeforehaloindex;
	int beforehalomergerindex;

	printf("input rng ");
	scanf("%f",&rng);
	printf("input nspace ");
	scanf("%d",&nspace);
	printf("input omep ");
	scanf("%f",&omep);
	printf("input omeplam ");
	scanf("%f",&omeplam);
	printf("input amax ");
	scanf("%f",&amax);
	printf("input a ");
	scanf("%f",&a);
	printf("input h ");
	scanf("%f",&hubble);
	printf("input size ");
	scanf("%f",&size);
	nsize = rng*rng*rng/nspace/nspace/nspace;
	omepk = 1.-omep-omeplam;
   		Hsub = sqrt(omep*pow3(amax/a)+omeplam+omepk*pow2(amax/a));
   		r1kineticfact = size/hubble/rng/amax*a*100.E5*hubble*Hsub;
	    r2kineticfact = size/hubble/amax*a*a*100.E5*hubble*Hsub;
	   	com2real = size/hubble/rng/amax*a*1.E6*pc;
	    real2com = 1./com2real;
    	pntmass = 3./8./pi/G*sqr(100.E5*hubble)*1.E6*pc/onesolarmass;
	    pntmass = pntmass*pow3(size/hubble)/pow3(rng/nspace)*omep;
		printf("point mass is %g h^{-1}\n",pntmass*hubble);
		printf("halo(N=20) mass is %gh^{-1}\n",20.*pntmass*hubble);
		printf("mean density is %g h^2\n",pntmass*pow3(rng/size)*hubble);
}
