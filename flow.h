typedef struct HALF{
	int first,second;
}HALF;
typedef struct EvolFact{
	float afact,bfact,pfact;
	float fact1,fact2;
	float fact1_push,fact2_push;
	float fact1_pull,fact2_pull;
} EvolFact;

#ifdef DEFINE_SIM_FLOW
#define EXTERN
#else
#define EXTERN extern
#endif
EXTERN HALF halfstep;
EXTERN int  maxsubT;
#undef EXTERN


#define Teta 0.5
#define GetNSub(dist) ((dist)/(EPSILON*Teta))

#define isnowstep(isub,msub) (( isub & (( 1<< (msub-isub))-1)) ==0 ? 1 : 0)
int flagpsmeasure(float , float , float ,int );
int flagBinnedDen(float , float , float );
int flagPreFoF(float , float , float , int );
int flagsyncpdata(float , float , float );

int flagwholeden(float , float , float );
int flagsuddenstop(int );

void Savexzslice(float *);
void SaveWholeDen(float *);
void SaveBinnedData(float *);
