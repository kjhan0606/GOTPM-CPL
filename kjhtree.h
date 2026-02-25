/* The calling program also includes pmheader.h */
#define sqr(A) ((A)*(A))
#define NODE 0
#define PARTICLE 1
#define NODE_HAVE_PARTICLE 16
#define MINWIDTH 1.e-6
#define EPSILON (0.1L)
#define MINFOFNP 20



#define DIRECTSUM 16
/* This value should be determined from the gputest.c */
#ifdef USE_GPU
#define GPUCPUDIRECTSUM (25)
/*
#define GPUDIRECTSUM (250)
*/
int GPUDIRECTSUM;
#endif




typedef float TREEPRECISION;


enum boolean {YES=01, NO=00};
typedef struct Box{
	TREEPRECISION x,y,z,width;
} Box;
#define SQRT(A) sqrt(A)


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

/******************************************************
 *
 *    For the Tree force correction 
 *
 ********************************************************/
typedef struct TPtlStruct{
	GENERALTPtlPOINTER;
	TREEPRECISION r[3];
	TREEPRECISION mass;
	TREEPRECISION vx,vy,vz;
	indextype indx;
} TPtlStruct;
typedef struct TStruct{
	GENERALTPtlPOINTER;
	void *daughter;
	TREEPRECISION L;
	TREEPRECISION dist2;
	TREEPRECISION dist_over_thetasq;
	TREEPRECISION mono[3];
	TREEPRECISION mass;
	int Nparticle;
	TREEPRECISION mrr;
	TREEPRECISION quad[6];
	TREEPRECISION trQ;
} TStruct;
typedef struct BeginEndTree{
	TStruct *start;
	void *EndTree;
	void *EndPtl;
} BeginEndTree;
/******************************************************
 *
 *    For the FoF 
 *
 ********************************************************/
typedef struct FoFTPtlStruct{
	GENERALTPtlPOINTER;
	TREEPRECISION r[3];
	treeparticletype *bp;
	enum boolean included;
	size_t haloindx;
} FoFTPtlStruct;
typedef struct FoFTStruct{
	GENERALTPtlPOINTER;
	void *daughter;
	TREEPRECISION L;
	TREEPRECISION dist2;
	TREEPRECISION dist;
	TREEPRECISION mono[3];
	int Nparticle;
} FoFTStruct;
typedef struct FoFBeginEndTree{
	FoFTStruct *start;
	void *EndTree;
	void *EndPtl;
} FoFBeginEndTree;


typedef struct FoFBox{
	TREEPRECISION xi,yi,zi;
	TREEPRECISION xf,yf,zf,maxwidth;

} FoFBox;
typedef struct HaloBound{
	size_t nmem;
	float xmin,xmax;
	float ymin,ymax;
	float zmin,zmax;
	int dumpflag;
	FoFTPtlStruct *sibling;
}HaloBound;



/*      */
typedef struct DMParticle{
	TREEPRECISION x,y,z;
	TREEPRECISION vx,vy,vz;
} DMParticle;
typedef struct particle {
	TREEPRECISION x,y,z;
	TREEPRECISION ax,ay,az;
	treeparticletype *bp;
} particle;
typedef struct pforce {
	TREEPRECISION x,y,z;
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
	
void treeforce(particle*,float, TStruct *,TPtlStruct *,float);
int Make_Tree(TStruct *,TPtlStruct *,long ,Box ,float);
float treeplumpotential(particle*,float, TStruct *,TPtlStruct *);
BeginEndTree divide_node(TStruct *,TStruct *, TPtlStruct *, Box ,TStruct *); 
FoFBeginEndTree FoF_divide_node(FoFTStruct *,FoFTStruct *, FoFTPtlStruct *,
		FoFBox ,FoFTStruct *);
void FoF_Make_Tree(FoFTStruct *,FoFTPtlStruct *,long ,FoFBox );
void new_fof_link(particle*,float, FoFTStruct *, FoFTPtlStruct *,particle *);
int Find_Near(particle *,int ,TStruct *, TPtlStruct *,float *,float);
Box findbox(TPtlStruct *, long );
int pnew_fof_link(particle *,float ,FoFTStruct *, FoFTPtlStruct *, particle *,int);


