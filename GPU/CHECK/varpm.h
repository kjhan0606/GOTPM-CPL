typedef struct{
	/* floor and ceil must be in front of nfloor, nceil */
	int nfloor,nceil;
} NZBound;
void VarPM_tsc(float *,int ,int ,int , float, float, pmparticletype *, 
		int ,int , int , float , int , int , float *, int , int);
void VarPM_tsc2(float *,int ,int ,int , float, float, pmparticletype *, 
		int ,int , int , float , float *, int , int);
void VarPM_veltsc2(float *,int ,int ,int , float, float, pmparticletype *, 
		int ,int , int , float , float *, int , int, int);
void VarPM_fda4_mesh( float *, int , int , float *, float , float , int ,int ,
		int );
void adVarPM_fda4_mesh( float *, int , int , float *, float , float , int ,int ,
		int );
void paddingpotentialmesh(float *,int , int , int ,int , int ,int, int ,int);
void paddingonepotentialmesh(float *,int , int , int ,int , int ,int, int, int);
#define UpperSending 1
#define BelowSending -1
