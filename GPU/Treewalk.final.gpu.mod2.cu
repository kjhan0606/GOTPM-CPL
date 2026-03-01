/* This version seems ok. 21/02/2008 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include "cutil.h"
#include "pmheader.h"
#include "kjhtree.h"
#include "force_spline.h"
//#include "Memory.h"
#include "Memorygpu.h"
#include "nvgpu.h"
/* __constant__  for input, maximum:  16,384 elements */
/* __shared__ for output, maximum: 4096 elements*/


#define MIN(a,b) ( (a) < (b) ? (a): (b))
#define MAX(a,b) ( (a) > (b) ? (a): (b))

#define YES 1
#define NO 0
#define END 0

#define NPCMAX 5000

static cudaDeviceProp mydev;
static int mydevid;
#define GetBlockSize1D(npos,threadSize) \
	(npos/threadSize < mydev.maxThreadsDim[0]?  (npos/threadSize+((npos)%(threadSize)>0? 1:0)) : 0)

double  aa[8] = {
	0.000015258789062500L, /* == 1/256/256 */
	0.78224E+00L,
	0.37971E-06L,
	-0.60338E+00L,
	-0.34419E-07L,
	-0.39741E+01L,
	-0.10607E+01L,
	-0.38145E+00L
};

dim3 GetBlockSize3D(int np, int threadsize){
	dim3 pGrvBlockSize;
	int nx,ny,nz;
	int nblock;
	nx = mydev.maxGridSize[0];
	ny = mydev.maxGridSize[1];
	nz = mydev.maxGridSize[2];
	nblock = GetBlockSize1D(np,threadsize);
	pGrvBlockSize.z = nblock/(nx*ny);
	if(pGrvBlockSize.z > nz){
		fprintf(stderr,"Error occurred in the blocksize %d max. allowable= %g\n",
				pGrvBlockSize.z,nz);
		exit(124);
	}
	else if(pGrvBlockSize.z > 1) {
		pGrvBlockSize.y = ny;
		pGrvBlockSize.x = nx;
	}
	else if(nblock%nx ==0){
		pGrvBlockSize.y = (nblock-pGrvBlockSize.z*nx*ny)/nx;
		pGrvBlockSize.x = nx;
	}
	else{
		pGrvBlockSize.y = (nblock-pGrvBlockSize.z*nx*ny)/nx + 1;
		if(pGrvBlockSize.y ==0) {
			pGrvBlockSize.x = nblock;
		}
		else{
			pGrvBlockSize.x = nx;
		}
	}
	return pGrvBlockSize;
}

typedef struct __align__(16) ContactCell{
	int ipos;
	int icell;
	float ax,ay,az;
} ContactCell;

typedef struct __align__(16) ContactGrav{
	int ipos;
	int igrv;
	float ax,ay,az;
} ContactGrav;

typedef int oContactCell;
typedef int oContactGrav;

double  pmcfd(double r,int nx){
	double fr,tmp1,LL,rr;
	if(r<=0.L) return 0.L;
	LL = nx*nx;
	tmp1 = cosh(aa[1]*r);
	rr = r*r;
	fr = r/pow(rr+DEPSILON*DEPSILON,1.5L)
		-(
				(1.L/(rr)*tanh(aa[1]*r)-aa[1]/(tmp1*tmp1*r)
				 -2.L*aa[2]/aa[0]*r*exp(aa[3]*rr)
				 -2.L*aa[2]*aa[3]/aa[0]*rr*r*exp(aa[3]*rr)
				 +aa[4]/aa[0]*(1.L+aa[5]*rr+aa[6]*rr*rr)*exp(aa[7]*rr))
		 );
	fr = fr/LL;
	return fr;
}
#define EPS 1.E-1L;
float4 pmcf4(float sr, int nx){
	double r,dr,rm1,rm2,rp1,rp2;
	double fp,fpp,fr;
	float4 res;
	r = (double)sr;
	dr = r*EPS;
	rm1 = r - dr; rm2 = rm1 - dr;
	rp1 = r + dr; rp2 = rp1 + dr;
	fp = (pmcfd(rp1,nx)-pmcfd(rm1,nx))/(dr+dr);
	fpp = (pmcfd(rp2,nx)-2.L*pmcfd(r,nx)+pmcfd(rm2,nx))/(4.L*dr*dr);
	fr = pmcfd(r,nx);
	res.x = -(fr/r);
	res.y = -(0.5L*(fp*r-fr)/r/r/r);
	res.z = -(0.5L*(fpp*r*r-3.L*fp*r+3.L*fr)/r/r/r/r/r);
	return res;
}
#undef EPS

texture<float,1,cudaReadModeNormalizedFloat> tex1,tex2,tex3;
static cudaArray  *cuArray1, *cuArray2, *cuArray3, *cuArray4;
float *force1,*force2,*force3;


__global__ void kernel_PtclGravityMeasure(particle *dp,int npos,
		ContactGrav *outptcl, int tnptcl,TPtlStruct *dpgrv,
		float rspheresq,int nx){
	int idx = (gridDim.x*blockIdx.y+blockIdx.x)*threadSize+ threadIdx.x; 
	if(idx< npos){
		float dist2;
		float tmpx,tmpy,tmpz,fplmf;
		int you,me;
		particle mypos; 
		TPtlStruct urpos;
		ContactGrav contactgrav;

		contactgrav = outptcl[idx];
		you = contactgrav.igrv;
		me = contactgrav.ipos; 
		urpos = dpgrv[you]; /* This is getting data from the global memory space */
		mypos = dp[me]; /* This is getting data from the global memory space */

		tmpx = mypos.x - urpos.r[0]; 
		tmpy = mypos.y - urpos.r[1];
		tmpz = mypos.z - urpos.r[2];

		dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
		float rr = __fdividef(0.5f*__log10f(dist2) - LOG10RMIN,LOG10RMAXmLOG10RMIN);
		if(rr<1. && rr > 0.){
			fplmf = tex1D(tex1,rr);
			fplmf = fplmf*urpos.mass;
			outptcl[idx].ax = tmpx*fplmf;
			outptcl[idx].ay = tmpy*fplmf;
			outptcl[idx].az = tmpz*fplmf;
		}

	}
}
__global__ void kernel_CellGravityMeasure(particle *dp,int npos,ContactCell *outcell,
		int tnpcell, TStruct *dpcell, float rspheresq, int nx){
	int idx = (gridDim.x*blockIdx.y+blockIdx.x)*threadSize+ threadIdx.x; 
	if(idx< npos){
		float tmpx,tmpy,tmpz,dist2;
		int you,me;
		TStruct urpos;

		ContactCell contactcell = outcell[idx];
		you = contactcell.icell;
		me = contactcell.ipos;
		urpos = dpcell[you]; /* This is getting data from the global memory space */
		particle mypos = dp[me]; /* This is getting data from the global memory space */
		tmpx = mypos.x - urpos.mono[0];
		tmpy = mypos.y - urpos.mono[1];
		tmpz = mypos.z - urpos.mono[2];
		dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		float rr = __fdividef(0.5f*__log10f(dist2) - LOG10RMIN, LOG10RMAXmLOG10RMIN);
		if(rr <1. && rr > 0.){
			float Qxx = urpos.quad[0];
			float Qyy = urpos.quad[1];
			float Qzz = urpos.quad[2];
			float Qxy = urpos.quad[3];
			float Qxz = urpos.quad[4];
			float Qyz = urpos.quad[5];
			float trQ = (Qxx+Qyy+Qzz);
			float qxx1 = Qxx*tmpx+Qxy*tmpy+Qxz*tmpz;
			float qxx2 = Qxy*tmpx+Qyy*tmpy+Qyz*tmpz;
			float qxx3 = Qxz*tmpx+Qyz*tmpy+Qzz*tmpz;
			float qxy = Qxx*tmpx*tmpx+Qyy*tmpy*tmpy+Qzz*tmpz*tmpz+
				2.*(Qxy*tmpx*tmpy+Qxz*tmpx*tmpz+Qyz*tmpy*tmpz);
			float4 fplmf;
			fplmf.x = tex1D(tex1,rr);
			fplmf.y = tex1D(tex2,rr);
			fplmf.z = tex1D(tex3,rr);
			float tmptmp = urpos.mass*fplmf.x + trQ*fplmf.y + qxy*fplmf.z;
			float twofplmf2 = 2.*fplmf.y;
			outcell[idx].ax = tmpx*tmptmp+qxx1*twofplmf2;
			outcell[idx].ay = tmpy*tmptmp+qxx2*twofplmf2;
			outcell[idx].az = tmpz*tmptmp+qxx3*twofplmf2;
		}

	}
}
#define cellopen2(dist2,cell) (dist2<cell->dist_over_thetasq ? YES:NO)
#define cellopen(dist2,ptr) (dist2<((TStruct*)ptr)->dist_over_thetasq ? YES:NO)
__device__ void gpu_treeforce(particle *p, int mytidx, TStruct *tree, 
		TPtlStruct *grv, oContactGrav *ptcl, oContactCell *cell, int *nc,int *np ){
	void *ptr;
	float tmpx,tmpy,tmpz; 
	float dist2; 
	int mp,mc;
//	p->ax = p->ay = p->az = 0.;
	ptr = (void*)tree;
	mp = mc  = 0;
	while(ptr != NULL){
		if(((TYPE*)ptr)->type == TYPE_TREE){
			tmpx = p->x-((TStruct*)ptr)->mono[0];
			tmpy = p->y-((TStruct*)ptr)->mono[1];
			tmpz = p->z-((TStruct*)ptr)->mono[2];
			dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
			if(cellopen(dist2,ptr)==YES) 
				ptr = (void*)(((TStruct*)ptr)->daughter);
			else {
				cell[mc++] = (int)((TStruct*)ptr-tree);
				ptr = (void*)(((TStruct*)ptr)->sibling);
			}
		}
		else{
			ptcl[mp++] = (int)((TPtlStruct*)ptr-grv);
			ptr = (void*)(((TPtlStruct*)ptr)->sibling);
		}
	}
	*nc = mc; *np = mp;
}



__global__ void intsum1Dblock_kernel(int *a, int *out){
	extern __shared__ int share_a[];
	const int idx= blockIdx.x*blockDim.x + threadIdx.x;
	share_a[threadIdx.x] = a[idx];
	__syncthreads();
	for(int distance=blockDim.x/2;distance>0;distance/=2){
		if(threadIdx.x<distance) share_a[threadIdx.x] += share_a[threadIdx.x+distance];
		__syncthreads();
	}
	if(threadIdx.x==0) out[blockIdx.x] = share_a[0];
}
void intsum_kernel1D(int *in, int *out, int *blocks,int NumofBlocks){
	intsum1Dblock_kernel<<<NumofBlocks,threadSize>>>(in,blocks);
	intsum1Dblock_kernel<<<1,NumofBlocks>>>(in,blocks);
	*out = blocks[0];
}
//
__global__ void forcesum1Dblock_kernel(ContactCell *a, int npos,float4 *out){
	extern __shared__ float4 share_b[];
	const int idx= blockIdx.x*blockDim.x + threadIdx.x;
	if(idx < npos) {
		share_b[threadIdx.x].x = a[idx].ax;
		share_b[threadIdx.x].y = a[idx].ay;
		share_b[threadIdx.x].z = a[idx].az;
	}
	else {
		share_b[threadIdx.x] = make_float4(0.f,0.f,0.f,0.f);
	}
	__syncthreads();
	for(int distance=blockDim.x/2;distance>0;distance/=2){
		if(threadIdx.x<distance) {
			share_b[threadIdx.x].x = share_b[threadIdx.x].x+share_b[threadIdx.x+distance].x;
			share_b[threadIdx.x].y = share_b[threadIdx.x].y+share_b[threadIdx.x+distance].y;
			share_b[threadIdx.x].z = share_b[threadIdx.x].z+share_b[threadIdx.x+distance].z;
		}
		__syncthreads();
	}
	if(threadIdx.x==0) out[blockIdx.x] = share_b[0];
}
inline void forcesum_kernel1D(ContactCell *in, int npos,float4 *out, float4 *blocks,int NumofBlocks){
	forcesum1Dblock_kernel<<<NumofBlocks,threadSize>>>(in,npos,blocks);
	forcesum1Dblock_kernel<<<1,NumofBlocks>>>(in,NumofBlocks,blocks);
	*out = blocks[0];
}
#define PARTICLEFORCE(p,ptr,rspheresq){\
	point = (TPtlStruct*)ptr;\
	tmpx = p.x - point->r[0];\
	tmpy = p.y - point->r[1];\
	tmpz = p.z - point->r[2];\
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;\
	float rr = __fdividef(0.5f*__log10f(dist2)-LOG10RMIN, LOG10RMAXmLOG10RMIN);\
	if(rr>0. && rr < 1.) \
	{\
		float fplmf = point->mass * tex1D(tex1,rr);\
		accel.x += tmpx*fplmf;\
		accel.y += tmpy*fplmf;\
		accel.z += tmpz*fplmf;\
	}\
}
#define CELLFORCE2(p,cell,rspheresq){\
	float rr = __fdividef(0.5f*__log10f(dist2)-LOG10RMIN, LOG10RMAXmLOG10RMIN);\
	if(rr>0. && rr < 1.) \
	{\
		float4 fplmf;\
		fplmf.x = tex1D(tex1,rr);\
		fplmf.y = tex1D(tex2,rr);\
		fplmf.z = tex1D(tex3,rr);\
		float Qxx = cell->quad[0];\
		float Qyy = cell->quad[1];\
		float Qzz = cell->quad[2];\
		float Qxy = cell->quad[3];\
		float Qxz = cell->quad[4];\
		float Qyz = cell->quad[5];\
		float trQ = (Qxx+Qyy+Qzz);\
		float qxx1 = Qxx*tmpx+Qxy*tmpy+Qxz*tmpz;\
		float qxx2 = Qxy*tmpx+Qyy*tmpy+Qyz*tmpz;\
		float qxx3 = Qxz*tmpx+Qyz*tmpy+Qzz*tmpz;\
		float qxy = Qxx*tmpx*tmpx+Qyy*tmpy*tmpy+Qzz*tmpz*tmpz+ \
			2.*(Qxy*tmpx*tmpy+Qxz*tmpx*tmpz+Qyz*tmpy*tmpz);\
		float tmptmp = cell->mass*fplmf.x + trQ*fplmf.y + qxy*fplmf.z;\
		float twofplmf2 = 2.*fplmf.y;\
		accel.x += tmpx*tmptmp+qxx1*twofplmf2;\
		accel.y += tmpy*tmptmp+qxx2*twofplmf2;\
		accel.z += tmpz*tmptmp+qxx3*twofplmf2;\
	}\
}
__global__ void All_gpu_treeforce(particle *pp, int npos, TStruct *tree,TPtlStruct *grv,
		float rspheresq){
	void *ptr,*optr;
	float tmpx,tmpy,tmpz; 
	float dist2; 
	TPtlStruct *point;
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if(idx<npos){
		float4 accel;
		accel = make_float4( 0.f,0.f,0.f,0.f);
		particle p = pp[idx];
		ptr = (void*)tree;
		while(ptr != NULL){
			if(((TYPE*)ptr)->type == TYPE_TREE){
				TStruct *cell;
				cell = (TStruct*)ptr;
				tmpx = p.x-cell->mono[0];
				tmpy = p.y-cell->mono[1];
				tmpz = p.z-cell->mono[2];
				dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
				if(cellopen2(dist2,cell)==YES) 
					ptr = (void*)(cell->daughter);
				else {
					ptr = (void*)(cell->sibling);
					CELLFORCE2(p,cell,rspheresq);
				}
			}
			else{
				PARTICLEFORCE(p,ptr,rspheresq);
				ptr = (void*)(((TPtlStruct*)ptr)->sibling);
			}
		}
		pp[idx].ax = accel.x; pp[idx].ay = accel.y; pp[idx].az = accel.z;
	}
}


__global__ void kernel_GetContactList(particle *dp, int npos, int sidx,TStruct *dtree,
		TPtlStruct *dpgrv, oContactGrav *listptcl,oContactCell *listcell, int *nc, int *np){
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int mytidx;
	if(idx < npos) {
		mytidx = sidx+idx;
		gpu_treeforce(dp+idx,mytidx,dtree,dpgrv,listptcl+idx*NPCMAX,
				listcell+idx*NPCMAX,nc+idx,np+idx);
	}
}
__global__ void	kernel_StackLeftWardCell(oContactCell *list,int *np,ContactCell *out){
	const int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int ileft=0;
	for(int i = 0;i<idx;i++) ileft += np[i];
	for(int i = 0;i<np[idx];i++){
		out[ileft].ipos = idx;
		out[ileft++].icell = list[i+idx*NPCMAX];
	}
}
__global__ void	kernel_StackLeftWardGrav(oContactGrav *list,int *np,ContactGrav *out){
	const int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int ileft=0;
	for(int i = 0;i<idx;i++) ileft += np[i];
	for(int i = 0;i<np[idx];i++){
		out[ileft].ipos = idx;
		out[ileft++].igrv = list[i+idx*NPCMAX];
	}
}
/*
#define GTREE_RESET 0x00
#define GTREE_SIB_PTCL 0x01
#define GTREE_SIB_CELL 0x02
#define GTREE_DAU_PTCL 0x04
#define GTREE_DAU_CELL 0x08
#define RESET_WHOLE_FLAGS(x,i) (x[i].gpuflag = GTREE_RESET)

#define SET_SIB_PTCL(x,i) (x[i].gpuflag |= GTREE_SIB_PTCL)
#define SET_SIB_CELL(x,i) (x[i].gpuflag |= GTREE_SIB_CELL)
#define SET_DAU_PTCL(x,i) (x[i].gpuflag |= GTREE_DAU_PTCL)
#define SET_DAU_CELL(x,i) (x[i].gpuflag |= GTREE_DAU_CELL)

#define UNSET_SIB_PTCL(x,i) (x[i].gpuflag &= (~GTREE_SIB_PTCL))
#define UNSET_SIB_CELL(x,i) (x[i].gpuflag &= (~GTREE_SIB_CELL))
#define UNSET_DAU_PTCL(x,i) (x[i].gpuflag &= (~GTREE_DAU_PTCL))
#define UNSET_DAU_CELL(x,i) (x[i].gpuflag &= (~GTREE_DAU_CELL))

#define IS_SIB_PTCL(x,i) (x[i].gpuflag & GTREE_SIB_PTCL)
#define IS_SIB_CELL(x,i) (x[i].gpuflag & GTREE_SIB_CELL)
#define IS_DAU_PTCL(x,i) (x[i].gpuflag & GTREE_DAU_PTCL)
#define IS_DAU_CELL(x,i) (x[i].gpuflag & GTREE_DAU_CELL)

__global__ void kernel_initCELLPointer(TStruct *dtree,int ncell,TStruct *TREECELL,
		TPtlStruct *dpgrv, int ngrv, TPtlStruct *grv){
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if(idx<ncell){
		if(IS_SIB_PTCL(dtree,idx))
			dtree[idx].sibling = (void*)(((TPtlStruct *)(dtree[idx].sibling)-grv)+dpgrv);
		else if(IS_SIB_CELL(dtree,idx))
			dtree[idx].sibling = (void*)(((TStruct *)(dtree[idx].sibling)-TREECELL)+dtree);
		else
			dtree[idx].sibling = NULL;
		if(IS_DAU_PTCL(dtree,idx))
			dtree[idx].daughter = (void*)(((TPtlStruct *)(dtree[idx].daughter)-grv)+dpgrv);
		else if(IS_SIB_CELL(dtree,idx))
			dtree[idx].daughter = (void*)(((TStruct *)(dtree[idx].daughter)-TREECELL)+dtree);
		else
			dtree[idx].daughter = NULL;
	}
}
__global__ void kernel_initPTCLPointer(TStruct *dtree,int ncell,TStruct *TREECELL,
		TPtlStruct *dpgrv, int ngrv, TPtlStruct *grv){
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if(idx<ngrv){
		if(IS_SIB_PTCL(dpgrv,idx))
			dpgrv[idx].sibling = (void*)(((TPtlStruct *)(dpgrv[idx].sibling)-grv)+dpgrv);
		else if(IS_SIB_CELL(dpgrv,idx))
			dpgrv[idx].sibling = (void*)(((TStruct *)(dpgrv[idx].sibling)-TREECELL)+dtree);
		else
			dpgrv[idx].sibling = NULL;
	}
}
*/
__global__ void	Crr_CellPointer(TStruct *dtree,TStruct *TREECELL,
		int ncell,TPtlStruct *dpgrv,TPtlStruct *grv,int ngrv){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx<ncell){
		char *a;
		int offset,cellwidth;
		cellwidth = ncell*sizeof(TStruct);

		a = (char *)(dtree[idx].sibling);
		offset = (int)(a - (char*) TREECELL);
		if(a == NULL){}
		else if(offset < cellwidth && offset >= 0){
			dtree[idx].sibling = dtree + offset/(sizeof(TStruct)/sizeof(char));
		}
		else {
			offset = (int)(a-(char*) grv);
			dtree[idx].sibling = dpgrv + offset/(sizeof(TPtlStruct)/sizeof(char));
		}

		a = (char *)(dtree[idx].daughter);
		offset = (int)(a - (char*) TREECELL);
		if(a == NULL){}
		else if(offset < cellwidth && offset >= 0){
			dtree[idx].daughter = dtree + offset/(sizeof(TStruct)/sizeof(char));
		}
		else {
			offset = (int)(a-(char*) grv);
			dtree[idx].daughter = dpgrv + offset/(sizeof(TPtlStruct)/sizeof(char));
		}
	}
}
__global__ void	Crr_GrvPointer(TStruct *dtree,TStruct *TREECELL,
		int ncell,TPtlStruct *dpgrv,TPtlStruct *grv,int ngrv){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if(idx<ngrv){
		char *a;
		int offset,cellwidth;
		cellwidth = ncell*sizeof(TStruct);

		a = (char *)(dpgrv[idx].sibling);
		offset = (int)(a - (char*) TREECELL);
		if(a == NULL){}
		else if(offset < cellwidth && offset >= 0){
			dpgrv[idx].sibling = dtree + offset/(sizeof(TStruct)/sizeof(char));
		}
		else {
			offset = (int)(a-(char*) grv);
			dpgrv[idx].sibling = dpgrv + offset/(sizeof(TPtlStruct)/sizeof(char));
		}
	}
}


particle *dpos;
TStruct *dtree;
TPtlStruct *dpgrv;
int *ncll,*nptcl;
#define FORCETEXTURETYPE float4

static int initialflag=1;

void CUDA_TERMINATE(int myid){
	fprintf(stdout,"P%d Terminating Cuda Device\n",myid);
	cudaFreeArray(cuArray1);
	cudaFreeArray(cuArray2);
	cudaFreeArray(cuArray3);
	cudaFreeArray(cuArray4);
	cudaUnbindTexture(tex1);
	cudaUnbindTexture(tex2);
	cudaUnbindTexture(tex3);
}

extern "C" void gputreeforce(particle *pos,int npos,
		TStruct *TREECELL, int ncell,TPtlStruct *grv,int ngrv,
		float rspheresq, int GpusPerNode, int myid, int nx){
	unsigned int blockSize;
	int ndev;



	if(initialflag){
		float *hforce1,*hforce2,*hforce3;
		mydevid = myid%GpusPerNode;
		cudaSetDevice(mydevid);
		cudaGetDeviceProperties(&mydev,mydevid);

		cudaGetDeviceCount(&ndev);
		if(mydevid >= ndev){
			fflush(stderr);
			fprintf(stderr,"P%d has error allocating GPU device to this thread",myid);
			fprintf(stderr," : mydevid=%d ndev=%d\n",mydevid,ndev);
			exit(999);
		}
		else{
			fflush(stderr);
			if(myid==0) fprintf(stderr,"###########################################################\n");
			fprintf(stderr,"P%d has now been initializing the CPU setting and textures\n",myid);
			fprintf(stderr,"P%d has gpu id=%d with shared global memsize= %d\n",myid,mydevid,mydev.totalGlobalMem);
			fprintf(stderr,"P%d has gpu Rev. %d.%d Clock %d khz\n",myid,mydev.major,mydev.minor,mydev.clockRate);
		}
		if(GPUMake_Total_Memory((INT8)((mydev.totalGlobalMem)*0.7L/(double)GpusPerNode)) == 0 ){
			fprintf(stderr,"Can't initialize memory - aborting job\n");
			exit(3);
		}
		hforce1 = (float*)malloc(sizeof(float)*NTextureSize);
		hforce2 = (float*)malloc(sizeof(float)*NTextureSize);
		hforce3 = (float*)malloc(sizeof(float)*NTextureSize);
		CUDA_SAFE_CALL( cudaMalloc(PPTR(force1),sizeof(float)*NTextureSize));
		CUDA_SAFE_CALL( cudaMalloc(PPTR(force2),sizeof(float)*NTextureSize));
		CUDA_SAFE_CALL( cudaMalloc(PPTR(force3),sizeof(float)*NTextureSize));
		for(int i=0;i<NTextureSize;i++){
			float r = LOG10RMIN+(LOG10RMAX-LOG10RMIN)*(i+0.5)/(float) NTextureSize;
			float4 hforce4;
			r = powf(10.,r);
			hforce4 = pmcf4(r,nx);
			hforce1[i] = hforce4.x;
			hforce2[i] = hforce4.y;
			hforce3[i] = hforce4.z;
		}
		CUDA_SAFE_CALL( cudaMemcpy(force1,hforce1,sizeof(float)*NTextureSize,
					cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL( cudaMemcpy(force2,hforce2,sizeof(float)*NTextureSize,
					cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL( cudaMemcpy(force3,hforce3,sizeof(float)*NTextureSize,
					cudaMemcpyHostToDevice));

		free(hforce3); free(hforce2); free(hforce1);

		tex1.normalized= tex2.normalized= tex3.normalized= YES;
		tex1.filterMode = tex2.filterMode = tex3.filterMode 
			= cudaFilterModeLinear;
		tex1.addressMode[0] = tex2.addressMode[0] 
			= tex3.addressMode[0] = cudaAddressModeClamp;

		CUDA_SAFE_CALL( cudaMallocArray(&cuArray1,&tex1.channelDesc,NTextureSize,1));
		CUDA_SAFE_CALL( cudaMallocArray(&cuArray2,&tex2.channelDesc,NTextureSize,1));
		CUDA_SAFE_CALL( cudaMallocArray(&cuArray3,&tex3.channelDesc,NTextureSize,1));
		CUDA_SAFE_CALL( cudaBindTextureToArray(tex1,cuArray1));
		CUDA_SAFE_CALL( cudaBindTextureToArray(tex2,cuArray2));
		CUDA_SAFE_CALL( cudaBindTextureToArray(tex3,cuArray3));


		CUDA_SAFE_CALL( cudaMemcpyToArray(cuArray1,0,0,force1,
					sizeof(float)*NTextureSize,cudaMemcpyDeviceToDevice));
		CUDA_SAFE_CALL( cudaMemcpyToArray(cuArray2,0,0,force2,
					sizeof(float)*NTextureSize,cudaMemcpyDeviceToDevice));
		CUDA_SAFE_CALL( cudaMemcpyToArray(cuArray3,0,0,force3,
					sizeof(float)*NTextureSize,cudaMemcpyDeviceToDevice));

		CUDA_SAFE_CALL( cudaFree(force3));
		CUDA_SAFE_CALL( cudaFree(force2));
		CUDA_SAFE_CALL( cudaFree(force1));

		INT8 size = GPUCheckAvailableMemory();
		printf("P%d has %d bytes gpu memory pool\n",myid,size);

		initialflag = 0;
	}



	{
		/*
		CUDA_SAFE_CALL( cudaMalloc(PPTR(dtree),sizeof(TStruct)*ncell));
		CUDA_SAFE_CALL( cudaMalloc(PPTR(dpgrv),sizeof(TPtlStruct)*ngrv));
		*/

		/*
		unsigned int timer = 0;
		CUDA_SAFE_CALL( cudaThreadSynchronize() );
		CUT_SAFE_CALL( cutCreateTimer( &timer));
		CUT_SAFE_CALL( cutStartTimer( timer));
		*/

		dtree = (TStruct *)GPUMalloc(sizeof(TStruct)*ncell,PPTR(dtree));
		dpgrv = (TPtlStruct *)GPUMalloc(sizeof(TPtlStruct)*ngrv,PPTR(dpgrv));


		if(1){
			for(int i=0;i<ncell;i++){
				if(TREECELL[i].sibling == NULL){}
				else if((((TYPE*)(TREECELL[i].sibling))->type)==TYPE_TREE)
   	             TREECELL[i].sibling = (void*)(dtree +((TStruct*)(TREECELL[i].sibling)-TREECELL));
				else 
					TREECELL[i].sibling = (void*)(dpgrv +((TPtlStruct*)(TREECELL[i].sibling)-grv));
	
				if(TREECELL[i].daughter == NULL){}
				else if((((TYPE*)(TREECELL[i].daughter))->type)==TYPE_TREE)
					TREECELL[i].daughter = (void*)(dtree +((TStruct*)(TREECELL[i].daughter)-TREECELL));
				else 
					TREECELL[i].daughter = (void*)(dpgrv +((TPtlStruct*)(TREECELL[i].daughter)-grv));
			}
			for(int i=0;i<ngrv;i++){
				if(grv[i].sibling == NULL){}
				else if((((TYPE*)(grv[i].sibling))->type)==TYPE_TREE)
					grv[i].sibling = (void*)(dtree+((TStruct*)(grv[i].sibling)-TREECELL));
				else 
					grv[i].sibling = (void*)(dpgrv+((TPtlStruct*)(grv[i].sibling)-grv));
			}
			CUDA_SAFE_CALL( cudaMemcpy(dtree,TREECELL,sizeof(TStruct)*ncell,cudaMemcpyHostToDevice));
			CUDA_SAFE_CALL( cudaMemcpy(dpgrv,grv,sizeof(TPtlStruct)*ngrv,cudaMemcpyHostToDevice));
		}
		else{
			CUDA_SAFE_CALL( cudaMemcpy(dtree,TREECELL,sizeof(TStruct)*ncell,cudaMemcpyHostToDevice));
			CUDA_SAFE_CALL( cudaMemcpy(dpgrv,grv,sizeof(TPtlStruct)*ngrv,cudaMemcpyHostToDevice));

			blockSize = GetBlockSize1D(ncell,threadSize);
			Crr_CellPointer<<<blockSize,threadSize>>>(dtree,TREECELL,ncell,
					dpgrv,grv,ngrv);

			blockSize = GetBlockSize1D(ngrv,threadSize);
			Crr_GrvPointer<<<blockSize,threadSize>>>(dtree,TREECELL,ncell,
					dpgrv,grv,ngrv);
		}
		/*

		*/



		/*
		CUDA_SAFE_CALL( cudaMalloc(PPTR(dpos),sizeof(particle)*npos));
		*/
		dpos = (particle *)GPUMalloc(sizeof(particle)*npos,PPTR(dpos));
		CUDA_SAFE_CALL( cudaMemcpy(dpos,pos,sizeof(particle)*npos,cudaMemcpyHostToDevice));



		blockSize = GetBlockSize1D(npos,threadSize);
		All_gpu_treeforce<<<blockSize,threadSize>>>(dpos,npos,dtree,dpgrv,rspheresq);




		CUDA_SAFE_CALL( cudaMemcpy(pos,dpos,sizeof(particle)*npos,cudaMemcpyDeviceToHost));


		/*
		CUDA_SAFE_CALL( cudaFree(dtree));
		CUDA_SAFE_CALL( cudaFree(dpgrv));
		CUDA_SAFE_CALL( cudaFree(dpos));
		*/
		GPUFree(dpos);
		GPUFree(dpgrv);
		GPUFree(dtree);

		/*
		CUDA_SAFE_CALL( cudaThreadSynchronize() );
		CUT_SAFE_CALL( cutStopTimer( timer));
		float TimeInMs = cutGetTimerValue( timer);
		printf( "Processing time per batch: %f (ms)\n", TimeInMs);
		CUT_SAFE_CALL( cutDeleteTimer( timer));
		CUDA_TERMINATE(myid);exit(999);
		*/


	}
	/*
	if(myid==0){
		FILE *wp;
		wp = fopen("CHECK.DAT","w");
		for(int i=0;i<npos;i++){
			fprintf(wp,"%d %g %g %g : \n",i,pos[i].ax,pos[i].ay,pos[i].az);
		}
		fclose(wp);
	}
	*/

	return;




	int ipos,localnpos,stepnpos;

	stepnpos = mydev.totalGlobalMem/(2.*sizeof(oContactGrav)*NPCMAX)/GpusPerNode*0.7;

	// This is for the force array
	float4 *force1,*force2;
	float4 *hforce1,*hforce2,*blockforce;
	blockSize = GetBlockSize1D(NPCMAX,threadSize);
	cudaMallocHost(PPTR(hforce1),stepnpos*sizeof(float4));
	cudaMallocHost(PPTR(hforce2),stepnpos*sizeof(float4));
	cudaMalloc(PPTR(force1),stepnpos*sizeof(float4));
	cudaMalloc(PPTR(force2),stepnpos*sizeof(float4));
	cudaMalloc(PPTR(blockforce),blockSize*sizeof(float4));

	for(ipos=0;ipos<npos;ipos+=stepnpos){
		localnpos = MIN(stepnpos,npos-ipos);

		oContactCell *listcell;
		oContactGrav *listptcl;
		CUDA_SAFE_CALL( cudaMalloc(PPTR(listcell),localnpos*NPCMAX*sizeof(oContactCell)));
		CUDA_SAFE_CALL( cudaMalloc(PPTR(listptcl),localnpos*NPCMAX*sizeof(oContactGrav)));

		CUDA_SAFE_CALL( cudaMalloc(PPTR(ncll),localnpos*sizeof(int)));
		CUDA_SAFE_CALL( cudaMalloc(PPTR(nptcl),localnpos*sizeof(int)));
		/* This is to extract the contact list of positions dp of number npos 
		 * and nc, np is the contact number of cells and point mass */
		blockSize = GetBlockSize1D(localnpos,threadSize);
		kernel_GetContactList<<<blockSize,threadSize>>>(dpos+ipos,localnpos,ipos,dtree,dpgrv,
				listptcl,listcell,ncll,nptcl);

		int tncell,tnptcl; /* Get the total number of contact cells and particles */
		int *blocksum;
		cudaMalloc(PPTR(blocksum),blockSize*sizeof(int));
		intsum_kernel1D(ncll,&tncell,blocksum,blockSize);
		intsum_kernel1D(nptcl,&tnptcl,blocksum,blockSize);
		CUDA_SAFE_CALL( cudaFree(blocksum));

		ContactCell *outcell;
		ContactGrav *outptcl;
		CUDA_SAFE_CALL( cudaMalloc(PPTR(outcell),tncell*sizeof(ContactCell)));
		kernel_StackLeftWardCell<<<blockSize,threadSize>>>(listcell,ncll,outcell);

		CUDA_SAFE_CALL( cudaMalloc(PPTR(outptcl),tnptcl*sizeof(ContactGrav)));
		kernel_StackLeftWardGrav<<<blockSize,threadSize>>>(listptcl,nptcl,outptcl);
		CUDA_SAFE_CALL( cudaFree(listcell));
		CUDA_SAFE_CALL( cudaFree(listptcl));

		dim3 blockSize3D;
		blockSize3D = GetBlockSize3D(tnptcl,threadSize);
		kernel_PtclGravityMeasure<<<blockSize3D,threadSize>>>
			(dpos+ipos,localnpos,outptcl,tnptcl,dpgrv,rspheresq,nx);

		blockSize3D = GetBlockSize3D(tncell,threadSize);
		kernel_CellGravityMeasure<<<blockSize3D,threadSize>>>
			(dpos+ipos,localnpos,outcell,tncell,dtree,rspheresq,nx);

		{
			for(int i=0;i<localnpos;i++){
				int localblockSize = GetBlockSize1D(ncll[i],threadSize);
				forcesum_kernel1D(outcell,ncll[i],force1+i,blockforce,localblockSize);

				localblockSize = GetBlockSize1D(nptcl[i],threadSize);
				forcesum_kernel1D((ContactCell*)outptcl,nptcl[i],force2+i,blockforce,localblockSize);
			}
			CUDA_SAFE_CALL( cudaMemcpy(hforce1,force1,sizeof(float4)*localnpos,cudaMemcpyDeviceToHost));
			CUDA_SAFE_CALL( cudaMemcpy(hforce2,force2,sizeof(float4)*localnpos,cudaMemcpyDeviceToHost));
			for(int i=0;i<localnpos;i++){
				pos[ipos+i].ax = hforce1[i].x + hforce2[i].x;
				pos[ipos+i].ay = hforce1[i].y + hforce2[i].y;
				pos[ipos+i].az = hforce1[i].z + hforce2[i].z;
			}
		}

		CUDA_SAFE_CALL( cudaFree(ncll));
		CUDA_SAFE_CALL( cudaFree(nptcl));
		CUDA_SAFE_CALL( cudaFree(outcell));
		CUDA_SAFE_CALL( cudaFree(outptcl));

	}
	// This is to free the temporary memories 
	CUDA_SAFE_CALL( cudaFree(blockforce));
	CUDA_SAFE_CALL( cudaFree(force2));
	CUDA_SAFE_CALL( cudaFree(force1));
	CUDA_SAFE_CALL( cudaFree(hforce2));
	CUDA_SAFE_CALL( cudaFree(hforce1));

	CUDA_SAFE_CALL( cudaFree(dpos));
	CUDA_SAFE_CALL( cudaFree(dpgrv)); 
	CUDA_SAFE_CALL( cudaFree(dtree)); 
}
