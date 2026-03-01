#define First 0
#define Second 1
#define Last 2
#define None -999999
#define CLOSE -999
#define FindZofParticle(n,p) p[p[n].u.sort-1].z
#define zsortparticles(n,p) indexx(n,p-1)
typedef struct Pair{
    int dest;
    int flag;
}Pair;
typedef struct CommPair{
    Pair up;
    Pair down;
}CommPair;
typedef struct NSlice{
    int up,down;
}NSlice;
/********************************************************************/
#define NELZBOX 4
typedef struct ZBox{
	float zstart,zfinal,zwidth;
} ZBox;

#ifdef EQNBox
MPI_Aint indices[NELZBOX];
int blocklens[] = {1,1,1,1};
MPI_Datatype MPI_ZBox_type;
MPI_Datatype old_types[]={MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_UB};
MPI_Datatype stype;
MPI_Op myop;
void Zboxsum(ZBox *in, ZBox *inout, int *len, MPI_Datatype *dptr){
	int i;
	ZBox c;
	for(i=0;i<*len;i++){
		c.zstart = in->zstart + inout->zstart;
		c.zfinal = in->zfinal + inout->zfinal;
		c.zwidth = in->zwidth + inout->zwidth;
		*inout = c;
		in++; inout++;
	}
}
#define Delete_ZBox_MPI_Op {\
	MPI_Op_free(&myop);\
}

#define Define_ZBox_MPI_Type {\
	indices[0] = 0;\
    indices[1] = (MPI_Aint)offsetof(ZBox,zfinal);\
    indices[2] = (MPI_Aint)offsetof(ZBox,zwidth);\
    indices[3] = sizeof(ZBox);\
	MPI_Type_contiguous(3,MPI_FLOAT,&stype);\
	MPI_Type_commit(&stype);\
	MPI_Op_create((MPI_User_function *)Zboxsum,0,&myop);\
}

#endif
/********************************************************************/

void DetermineProcid();
void FindLocalDenwidth(float,float,NSlice);
void FFTW2PotSlice(int,int,float,float,float*,int,int,int,float*);
void Den2FFTWSlice(int,int,int,int,float*,int,int,int,float*);
void n2fft_mesh(float *, int ,int , float *, int , int , int ,int ,int );
void EQNZBox(int, pmparticletype *, ZBox *, float, float,int);
void tsc_fftw_mesh(float *, int ,int , float *, int , int , int ,int ,int ,int);
