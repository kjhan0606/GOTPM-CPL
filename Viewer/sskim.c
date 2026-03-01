/*
 * 이 코드는 찾아낸 헤일로의 peak density 를 구하는 것이다.
 * 여기에 헤일로의 물리적 성질을 첨가할것이다.
 */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<string.h>
int NMIN=10000;
float rscale;
typedef struct PtlPos{
	float x,y,z;
} PtlPos;
typedef struct PtlVel{
	float vx,vy,vz;
} PtlVel;
typedef struct Particle{
	float x,y,z,vx,vy,vz;
	long indx;
} Particle;
Particle *particle;
typedef struct Halo{
	size_t np;
	float x,y,z,vx,vy,vz;
} Halo;
Halo halo;
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
float a;
float amax;
float redshift;
float size;
float rng;
float omep;
float omeplam,omepb;
float hubble;
float omepk,Hsub;
int nx,ny,nz;
float *den;
void findden(int , PtlPos *, int ,int ,int );
int shape(float , float , float, float *,float *,float *,
		float *,float *,float *);
void triaxialshape_(int *,float *,float *,float *,float *,
		float *,float *,float *);
float wc(float );
float find_r(PtlPos *, int, float ,float, float *, float *,float *);
void plot_();
void sm_init_();
void haloinit_();
void contourinit_();
void sm_draww_();
void sm_boxx_();
void contourplot_();

double ang_mom[3];
float q,s;
float d[3],v[3][3];
double para_spin;
float tenergy,tmass,pvx;
double cx,cy,cz,cvx,cvy,cvz,pvr;
double gcx,gcy,gcz;
size_t nump;
float rr;
float r1[3],r2[3],r3[3];
size_t Fread(void *buff,size_t size,size_t nmem,FILE *fp){
	char a,b,c,d;
	char *data;
	int i,j;
	size_t result;
	result = fread(buff,size,nmem,fp);
	data = (char *) buff;
	for(i=0;i<nmem*size;i+=4){
		a = data[i]; b = data[i+1]; c = data[i+2]; d = data[i+3];
		data[i] = d;
		data[i+1] = c;
		data[i+2] = b;
		data[i+3] = a;
	}
	return result;
}

int main(int argc, char **argv){
	FILE *fp,*wp,*fp2,*paraid,*smp;
	PtlPos *r,*or;
	PtlVel *vr;
    float peak;
	float *x,*y,*z,*vx,*vy,*vz;
	int tnend,nend,i,j,k;
	int ncount;
	long *indx,id;
    int mpeak,np;
    float xmin,ymin,xmax,ymax,width;
	int anend[10000];
	int contournend[10000];
	int index,peakindex;
	char infile[100],parafile[100];
	char infile2[100];
	int inid;

	if(argc != 3) {
		fprintf(stderr,"Error input : plot infile_number outfilename\n");
		exit(99);
	}
	inid = atoi(argv[1]);
	sprintf(infile,"FoF_halo_cat.%.5d",inid);
	sprintf(infile2,"FoF_member_particle.%.5d",inid);
	if((fp = fopen(infile,"rb")) == NULL||
			(fp2 = fopen(infile2,"rb")) == NULL){
		fprintf(stderr,"error opening file %s && %s \n",
				infile,infile2);
		exit(98);
	}
	wp = fopen(argv[2],"w");
    {
	    float npower,omepl,bias,astep,anow;
	    fread(&size,sizeof(float),1,fp);
		fread(&hubble,sizeof(float),1,fp);
		fread(&npower,sizeof(float),1,fp);
		fread(&omep,sizeof(float),1,fp);
		fread(&omepb,sizeof(float),1,fp);
		fread(&omepl,sizeof(float),1,fp);
		fread(&bias,sizeof(float),1,fp);
		fread(&nx,sizeof(int),1,fp);
		fread(&nspace,sizeof(int),1,fp);
		fread(&amax,sizeof(int),1,fp);
		fread(&astep,sizeof(int),1,fp);
		fread(&anow,sizeof(int),1,fp);
		printf("L=%g h=%g omep=%g omepl=%g\n",size,hubble,omep,omepl);
		rscale = size/nx;
    }



	omepk = 1.-omep-omeplam;
    Hsub = sqrt(omep*pow3(amax/a)+omeplam+omepk*pow2(amax/a));
    r1kineticfact = size/hubble/rng/amax*a*100.E5*hubble*Hsub;
    r2kineticfact = size/hubble/amax*a*a*100.E5*hubble*Hsub;
    com2real = size/hubble/rng/amax*a*1.E6*pc;
    real2com = 1./com2real;
    pntmass = 3./8./pi/G*sqr(100.E5*hubble)*1.E6*pc;
    pntmass = pntmass*pow3(size/hubble)/pow3(rng/nspace)*omep;
	potentfact = G*pntmass/com2real;
	index = 0;
   	while(fread(&halo,sizeof(Halo),1,fp) == 1){
		double mass;

		nend=halo.np;
		particle=(Particle*)malloc(sizeof(Particle)*nend);
		fread(particle,sizeof(Particle),nend,fp2);
		{
			double cx,cy,cz;
			cx = cy = cz =0;
			for(i=0;i<nend;i++){
				cx += particle[i].x;
				cy += particle[i].y;
				cz += particle[i].z;
			}
			cx = cx/nend;
			cy = cy/nend;
			cz = cz/nend;
			if(cx > 20 && cx < 60 && cy > 480 && cy < 520
					&& cz > 420 && cz<460){
				fwrite(particle,sizeof(Particle),nend,wp);
			}
		}
		free(particle);
	}
	fclose(fp);
	fclose(wp);
}
