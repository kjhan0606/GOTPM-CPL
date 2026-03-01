#define sqr(A) ((A)*(A))
#define NODE 0
#define PARTICLE 1
#define NODE_HAVE_PARTICLE 4
#define EPSILON (0.1)
enum boolean {YES=01, NO=02};
/*
#define YES '0'
#define NO '1'
*/
typedef struct Box{
#ifdef XYZDBL
	double x,y,z;
#else
	float x,y,z;
#endif
	float width;
} Box;
/*      */
#define GENERAL_TYPE unsigned int type: 1
#define GENERALTPtlPOINTER GENERAL_TYPE;void *sibling
/*      */
enum {TYPE_TREE = 0,TYPE_PTL = 1};
typedef struct TYPE {
	GENERAL_TYPE;
} TYPE;
typedef struct GENERAL_TPtl_POINTER {
	GENERALTPtlPOINTER;
} GENERAL_TPtl_POINTER;
typedef struct TPtlStruct{
	GENERALTPtlPOINTER;
#ifdef XYZDBL
	double r[3];
#else
	float r[3];
#endif
	float mass;
	int indx;
} TPtlStruct;
typedef struct TStruct{
	GENERALTPtlPOINTER;
	void *daughter;
	float L;
	float dist2;
#ifdef XYZDBL
	double r0[3];
#else
	float r0[3];
#endif
	float mass;
	int Nparticle;
	float mono[3];
	float mrr;
	float quad[6];
} TStruct;
typedef struct BeginEndTree{
	TStruct *start;
	void *EndTree;
	void *EndPtl;
} BeginEndTree;
/* for Fof */
typedef struct FoFTPtlStruct{
	GENERALTPtlPOINTER;
#ifdef XYZDBL
	double r[3];
#else
	float r[3];
#endif
	int indx;
	TPtlStruct *pointer;
	enum boolean included;
/*
	char included;
*/
} FoFTPtlStruct;
typedef struct FoFTStruct{
	GENERALTPtlPOINTER;
	void *daughter;
	float L;
	float dist2;
	float dist;
#ifdef XYZDBL
	double r0[3];
#else
	float r0[3];
#endif
	float mono[3];
	int Nparticle;
} FoFTStruct;
typedef struct FoFBeginEndTree{
	FoFTStruct *start;
	void *EndTree;
	void *EndPtl;
} FoFBeginEndTree;

/*      */
typedef struct DMParticle{
#ifdef XYZDBL
	double x,y,z;
#else
	float x,y,z;
#endif
	float vx,vy,vz;
} DMParticle;
typedef struct particle {
#ifdef XYZDBL
	double x,y,z;
#else
	float x,y,z;
#endif
} particle;
typedef struct pforce {
#ifdef XYZDBL
	double x,y,z;
#else
	float x,y,z;
#endif
} pforce;

#define EraseFromTree(optr,ptr,nptr) {\
	switch(((TYPE*)optr)->type) {\
		case TYPE_TREE:\
			if(((FoFTStruct*)optr)->daughter == ptr) \
				((FoFTStruct*)optr)->daughter == nptr;\
			else ((FoFTStruct*)optr)->sibling = nptr;\
			break;\
		default :\
			((FoFTPtlStruct*)ptr)->sibling = nptr;\
	}\
}\
	
void treeforce(particle*,float, TStruct *,TPtlStruct *,pforce *);
void Make_Tree(TStruct *,TPtlStruct *,int ,Box );
float treeplumpotential(particle*,float, TStruct *,TPtlStruct *);
BeginEndTree divide_node(TStruct *,TStruct *, TPtlStruct *, Box ,TStruct *); 
FoFBeginEndTree FoF_divide_node(FoFTStruct *,FoFTStruct *, FoFTPtlStruct *, Box ,FoFTStruct *); 
void FoF_Make_Tree(FoFTStruct *,FoFTPtlStruct *,int ,Box );
void new_fof_link(particle*,float, FoFTStruct *, FoFTPtlStruct *,particle *);
int Find_Near(particle *,int ,TStruct *, TPtlStruct *,float *,float);
Box findbox(TPtlStruct *, int );
