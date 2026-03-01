#define NFILTER 9
#define NBIN 10000
#define NSED 5
#define NGAL 5
/*
#define ML_B_E 2000
#define ML_B_Sc 540
*/
#define ML_B_E (2000.*10./4.)
#define ML_B_Sc (540.*10./4.)
#define MLPOWER 0.75L  /* gamma_b = A * M^(n/3) if n = 1 Bahcall's paper*/
#define R_MAX 0.2
enum {CIRRUS = 0, SBG = 1, A220 = 2, AGN = 3, E = 4};
typedef struct Halo{
	float mass,x,y,z;
	int indx,mpeak,peakindx;
	float spin,pv,ang;
	float r,q,s;
	float a[3],b[3],c[3];
	float rot[3];
	float M2Lratio;
	float L;
	float density;
	float alpha;
	float fE;
	char type[20];
	int sending,receiving;
} Halo;

typedef struct galaxy_sed {
	int n;
	float *lambda,*flux;
} galaxy_sed;
typedef struct filterrespos {
	int n;
	float *lambda,*qe;
} filterrespos;
float red_flux(float, galaxy_sed, filterrespos,int );
galaxy_sed GAL_SED_INIT(char*);
filterrespos FILTER_INIT(char*);
float gasdev(long *);
float ran1(long *);
void lint(float *, float *,int,float, float *);
void dlint(double *, double *,int,double, double *);
void Mass2Light(char *, float *, float *);
#define MIN(A,B) ((A)<(B) ? (A): (B))
#define MAX(A,B) ((A)>(B) ? (A): (B))
#define VEGA_MAG 0.03
/*
struct {
	struct {
		float U=0.366, B=0.438, V=0.545, R=0.641, I=0.798;
		float J=1.22,H=1.63,K=2.19,Kp=2.12,L=3.45,Lstar=3.80;
	} efflambda;
	struct {
		float U=417.5,B=632,V=363.1,R=217.7,I=112.6;
		float J=31.47,H=11.38,K=3.961,Kp=4.479,L=0.708,Lstar=0.489;
	} flambda;
	struct {
		float U=0.770,B=-0.120,V=0.000,R=0.186,I=0.444;
		float J=0.899,H=1.379,K=1.886,Kp=1.826,L=2.765,Lstar=2.961;
	} zp;
} lum2mag;
*/
