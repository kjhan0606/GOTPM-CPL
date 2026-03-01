/* This version seems ok. 21/02/2008 */
/* This code adopts a constant memory space for TStruct & TPtlStruct */
// This is the fastest version in 27/02/2008.
// GpusPerNode should be modified because it is a temporary thing.
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include "cutil.h"
#include "cuda.h"
#include "pmheader.h"
#include "kjhtree.h"
#include "force_spline.h"
#include "Memorycpp.h"
#include "Memorygpu.h"
#include "nvgpu.h"
/* __constant__  for input, maximum:  16,384 elements */
/* __shared__ for output, maximum: 4096 elements*/




#define NEEDLESS 0 


#define MIN(a,b) ( (a) < (b) ? (a): (b))
#define MAX(a,b) ( (a) > (b) ? (a): (b))



#define YES 1
#define NO 0
#define END 0


cudaDeviceProp mydev;
static int mydevid;
/*
#define GetNumberofBlock(npos,BlockSize) \
	(npos/BlockSize < mydev.maxGridSize[0]?  (npos/BlockSize+((npos)%(BlockSize)>0? 1:0)) : 0)
	*/
#define GetNumberofBlock(npos,BlockSize) ((npos+BlockSize-1)/BlockSize)

#define SetTexture(a,b,c,d,e) {\
	a.channelDesc = cudaCreateChannelDesc<b>();\
	a.normalized = c;\
	a.filterMode = d;\
	a.addressMode[0] = e;\
}

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


double  pmcfd(double r,int nx){
	double fr,tmp1,LL,rr;
	if(r<=0.L) return 0.L;
	LL = nx*nx;
	tmp1 = cosh(aa[1]*r);
	rr = r*r;
	fr = r/pow(rr+EPSILON*EPSILON,1.5L)
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

/*
texture<float,1,cudaReadModeNormalizedFloat> tex1;
texture<float4,1,cudaReadModeNormalizedFloat> tex4;
*/
texture<float,1,cudaReadModeElementType> tex1;
texture<float4,1,cudaReadModeElementType> tex4;
static cudaArray  *cuArray1, *cuArray4;
float *force1;
float4 *hforce4;


/* This is texture in the global memory */
texture<int2  ,1,cudaReadModeElementType> tex_cellsibdau;
texture<int   ,1,cudaReadModeElementType> tex_ptclsib;
texture<float4,1,cudaReadModeElementType> tex_cellbasic;
texture<float4,1,cudaReadModeElementType> tex_cellquad1;
texture<float4,1,cudaReadModeElementType> tex_cellquad2;
texture<float4,1,cudaReadModeElementType> tex_ptclbasic;

/* This is texture in the texture memory */
texture<int2  ,1,cudaReadModeElementType> texcu_cellsibdau;
texture<int   ,1,cudaReadModeElementType> texcu_ptclsib;
texture<float4,1,cudaReadModeElementType> texcu_cellbasic;
texture<float4,1,cudaReadModeElementType> texcu_cellquad1;
texture<float4,1,cudaReadModeElementType> texcu_cellquad2;
texture<float4,1,cudaReadModeElementType> texcu_ptclbasic;

#define cellopen(dist2,dist_over_thetasq) (dist2<dist_over_thetasq ? YES:NO)

/* The actual arguments are p, accel, and now */
#define PARTICLEFORCE(p,now){\
	float4 point;\
	point = cuptcl ? tex1D(texcu_ptclbasic,now): tex1Dfetch(tex_ptclbasic,now);\
	tmpx = p.x - point.x;\
	tmpy = p.y - point.y;\
	tmpz = p.z - point.z;\
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;\
	float rr = (0.5f*__log10f(dist2)-LOG10RMIN)* invLOG10RMAXmLOG10RMIN;\
	if(rr < 1.) \
	{ /* point.w === particle mass */\
		float fplmf = point.w * tex1D(tex1,rr);\
		accel.x += tmpx*fplmf;\
		accel.y += tmpy*fplmf;\
		accel.z += tmpz*fplmf;\
	}\
}
/* The actual arguments are now, dist2, tmpx,tmpy,tmpz, accel */
#define CELLFORCE(p,now){\
	float rr = (0.5f*__log10f(dist2)-LOG10RMIN)* invLOG10RMAXmLOG10RMIN;\
	if( rr < 1.) \
	{\
		float4 fplmf;\
		fplmf = tex1D(tex4,rr);\
		float4 quad1,quad2;\
		quad1 = cucell ? tex1D(texcu_cellquad1,now): tex1Dfetch(tex_cellquad1,now);\
		quad2 = cucell ? tex1D(texcu_cellquad2,now): tex1Dfetch(tex_cellquad2,now);\
		float qxx1 = quad1.x*tmpx+quad1.w*tmpy+quad2.x*tmpz;\
		float qxx2 = quad1.w*tmpx+quad1.y*tmpy+quad2.y*tmpz;\
		float qxx3 = quad2.x*tmpx+quad2.y*tmpy+quad1.z*tmpz;\
		float qxy = quad1.x*tmpx*tmpx+quad1.y*tmpy*tmpy+quad1.z*tmpz*tmpz+ \
			2.*(quad1.w*tmpx*tmpy+quad2.x*tmpx*tmpz+quad2.y*tmpy*tmpz);\
		float tmptmp = quad2.z *fplmf.x + quad2.w*fplmf.y + qxy*fplmf.z;\
		float twofplmf2 = fplmf.y + fplmf.y;\
		accel.x += tmpx*tmptmp+qxx1*twofplmf2;\
		accel.y += tmpy*tmptmp+qxx2*twofplmf2;\
		accel.z += tmpz*tmptmp+qxx3*twofplmf2;\
	}\
}

#define ADDCELLFORCE(now,nlistcell) {work[idx*MaxNList+nlistcell] = now; nlistcell++;}
#define ADDPARTICLEFORCE(now,nlistcell) {work[idx*MaxNList + MaxNContactCell+ nlistptl] = now; nlistptl++;}

template<bool cucell,bool cuptcl>
__global__ void Adv_All_gpu_treeforce(particle *pp, int npos, float rspheresq,int *work){
	int now;
	int nlistcell,nlistptl;
	float tmpx,tmpy,tmpz; 
	float dist2; 
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	nlistcell = nlistptl= 0 ;
	if(idx<npos){
		particle p = pp[idx];
		now = 1;
		while(now != 0){
			if(now >0 ){
				float4 mono;
				int2 Nowxy;
				mono = cucell ? tex1D(texcu_cellbasic,now): tex1Dfetch(tex_cellbasic,now);
				tmpx = p.x-mono.x;
				tmpy = p.y-mono.y;
				tmpz = p.z-mono.z;
				dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
				Nowxy = cucell ? tex1D(texcu_cellsibdau,now): tex1Dfetch(tex_cellsibdau,now);
				if(cellopen(dist2,mono.w)==YES) 
					now = Nowxy.y;
				else {
					ADDCELLFORCE(now,nlistcell);
					now = Nowxy.x;
				}
			}
			else{
				now = -now;
				ADDPARTICLEFORCE(now,nlistptl);
				now = cuptcl ? tex1D(texcu_ptclsib,now): tex1Dfetch(tex_ptclsib,now);
			}
		}
		float4 accel = { 0.,0.,0.,0.};
		for(int i=0;i<nlistcell;i++){
			float tmpx,tmpy,tmpz;
			float4 mono;
			now = work[idx*MaxNList+i];
			mono = cucell ? tex1D(texcu_cellbasic,now): tex1Dfetch(tex_cellbasic,now);
			tmpx = p.x-mono.x;
			tmpy = p.y-mono.y;
			tmpz = p.z-mono.z;
			dist2 = tmpx*tmpx+ tmpy*tmpy + tmpz*tmpz;
			CELLFORCE(p,now);
		}
		for(int i=0;i<nlistptl;i++){
			now = work[idx*MaxNList+MaxNContactCell+i];
			PARTICLEFORCE(p,now);
		}
		pp[idx].ax = accel.x; pp[idx].ay = accel.y; pp[idx].az = accel.z;
	}
#ifndef __DEVICE_EMULATION__
	__syncthreads();
#endif
}
template<bool cucell,bool cuptcl>
__global__ void All_gpu_treeforce(particle *pp, int npos, float rspheresq){
	int now;
	float tmpx,tmpy,tmpz; 
	float dist2; 
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if(idx<npos){
		float3 accel = { 0.,0.,0.};
		particle p = pp[idx];
		now = 1;
		while(now != 0){
			if(now >0 ){
				float4 mono;
				int2 Nowxy;
				mono = cucell ? tex1D(texcu_cellbasic,now): tex1Dfetch(tex_cellbasic,now);
				tmpx = p.x-mono.x;
				tmpy = p.y-mono.y;
				tmpz = p.z-mono.z;
				dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
				Nowxy = cucell ? tex1D(texcu_cellsibdau,now): tex1Dfetch(tex_cellsibdau,now);
				if(cellopen(dist2,mono.w)==YES) 
					now = Nowxy.y;
				else {
					CELLFORCE(p,now);
					now = Nowxy.x;
				}
			}
			else{
				now = -now;
				PARTICLEFORCE(p,now);
				now = cuptcl ? tex1D(texcu_ptclsib,now): tex1Dfetch(tex_ptclsib,now);
			}
		}
		pp[idx].ax = accel.x; pp[idx].ay = accel.y; pp[idx].az = accel.z;
	}
#ifndef __DEVICE_EMULATION__
	__syncthreads();
#endif
}

#define TEXADDRESSMODECLAMP {\
	tex_cellsibdau.addressmode[0] = \
	tex_ptclsib.addressmode[0] = \
	tex_cellbasic.addressmode[0] = \
	tex_cellquad1.addressmode[0] = \
	tex_cellquad2.addressmode[0] = \
	tex_ptclbasic.addressmode[0] = cudaAddressModeClamp;\
}

#define TEXFILTERMODEPOINT {\ tex_cellsibdau.addressmode[0] = \
	tex_ptclsib.addressmode[0] = \
	tex_cellbasic.addressmode[0] = \
	tex_cellquad1.addressmode[0] = \
	tex_cellquad2.addressmode[0] = \
	tex_ptclbasic.addressmode[0] = cudaFilterModePoint;}


static int initial_flag=1;
static int cuda_initial_flag=1;

void CUDA_TERMINATE(int myid){
	fprintf(stdout,"P%d Terminating Cuda Device\n",myid);
	cudaFreeArray(cuArray1);
	cudaFreeArray(cuArray4);
	cudaUnbindTexture(tex1);
	cudaUnbindTexture(tex4);
}
static cudaArray  *cu_cellsibdau,*cu_ptclsib,*cu_cellbasic,*cu_cellquad1,*cu_cellquad2,*cu_ptclbasic;
static int2 *de_cellsibdau; int *de_ptclsib;
static float4 *de_cellbasic,*de_cellquad1,*de_cellquad2,*de_ptclbasic;

static int2 *h_cellsibdau;int *h_ptclsib;
static float4 *h_cellbasic,*h_cellquad1,*h_cellquad2,*h_ptclbasic;
static particle *dpos;

static int nde_cell=NTextureBuffer, nde_ptcl=NTextureBuffer;
static int nde_pos=NPosBuffer;

/*
#define CUDA_SAFE_CALL(a) a
*/

void cuda_setting(int myid, int GpuPerNode){
	int ndev;
	if(cuda_initial_flag ==1){
		mydevid = myid%GpuPerNode;
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
			if(myid==0) {
				fprintf(stderr,"###########################################################\n");
				fprintf(stderr,"P%d has now been initializing the CPU setting and textures\n",myid);
				fprintf(stderr,"P%d has gpu id=%d with global memsize= %lu Mbytes\n",myid,mydevid,mydev.totalGlobalMem/1048576);
				fprintf(stderr,"P%d has gpu  with shared memsize per Block= %lu bytes\n",myid,mydev.sharedMemPerBlock);
				fprintf(stderr,"P%d has gpu  with register 32-bit per Block= %lu E/A\n",myid,mydev.regsPerBlock);
				fprintf(stderr,"P%d has gpu  with maximum Grid size = %lu x %lu x %lu\n",
					myid,mydev.maxGridSize[0],mydev.maxGridSize[1],mydev.maxGridSize[2]);
				fprintf(stderr,"P%d has gpu Rev. %d.%d Clock %d khz\n",myid,mydev.major,mydev.minor,mydev.clockRate);
				fflush(stderr);
			}
		}
		/*
		if(GPUMake_Total_Memory((INT8)((mydev.totalGlobalMem)*0.7L/(double)GpuPerCpu)) == 0 ){
			fprintf(stderr,"Can't initialize memory - aborting job\n");
			exit(3);
		}
		*/
		cuda_initial_flag = 0;
	}
}



void INITIALIZE(int myid, int GpuPerNode, int nx){
	float *hforce1;
	{
		cuda_setting(myid, GpuPerNode);
		hforce1 = (float*)Malloc(sizeof(float)*NForceTexture,PPTR(hforce1));
		hforce4 = (float4*)Malloc(sizeof(float4)*NForceTexture,PPTR(hforce4));

		for(int i=0;i<NForceTexture;i++){
			float r = LOG10RMIN+(LOG10RMAX-LOG10RMIN)*(i+0.5)/(float) NForceTexture;
			float4 aforce4;
			r = powf(10.,r);
			aforce4 = pmcf4(r,nx);
			hforce1[i] = aforce4.x;
#ifdef DEBUG
			printf("+P%d has forcediff = %g\n",i,hforce1[i]);
#endif
		}
		for(int i=0;i<NForceTexture;i++){
			float r = LOG10RMIN+(LOG10RMAX-LOG10RMIN)*(i+0.5)/(float) NForceTexture;
			r = powf(10.,r);
			hforce4[i] = pmcf4(r,nx);
		}

		{
			/* Texture bindings for the Direct N-body */
			SetTexture(tex1,float,true,cudaFilterModeLinear,cudaAddressModeClamp);
			SetTexture(tex4,float4,true,cudaFilterModeLinear,cudaAddressModeClamp);
			CUDA_SAFE_CALL( cudaMallocArray(&cuArray1,&tex1.channelDesc,NForceTexture,1));
			CUDA_SAFE_CALL( cudaMallocArray(&cuArray4,&tex4.channelDesc,NForceTexture,1));
	
			CUDA_SAFE_CALL( cudaMemcpyToArray(cuArray1,0,0,hforce1,
						sizeof(float)*NForceTexture,cudaMemcpyHostToDevice));
			CUDA_SAFE_CALL( cudaMemcpyToArray(cuArray4,0,0,hforce4,
						sizeof(float4)*NForceTexture,cudaMemcpyHostToDevice));
	
	
			CUDA_SAFE_CALL( cudaBindTextureToArray(tex1,cuArray1));
			CUDA_SAFE_CALL( cudaBindTextureToArray(tex4,cuArray4));
		}
	


		Free(hforce4); Free(hforce1);

		cudaMalloc(PPTR(dpos),sizeof(particle)*nde_pos);

		if(1){ /* Texture bindings for the global memory space */
			SetTexture(tex_ptclsib,   int,   false,cudaFilterModePoint,cudaAddressModeClamp);
			SetTexture(tex_ptclbasic, float4,false,cudaFilterModePoint,cudaAddressModeClamp);
			SetTexture(tex_cellsibdau,int2,  false,cudaFilterModePoint,cudaAddressModeClamp);
			SetTexture(tex_cellbasic, float4,false,cudaFilterModePoint,cudaAddressModeClamp);
			SetTexture(tex_cellquad1, float4,false,cudaFilterModePoint,cudaAddressModeClamp);
			SetTexture(tex_cellquad2, float4,false,cudaFilterModePoint,cudaAddressModeClamp);

			cudaMalloc(PPTR(de_cellsibdau),sizeof(int2)*nde_cell);
			cudaMalloc(PPTR(de_cellbasic),sizeof(float4)*nde_cell);
			cudaMalloc(PPTR(de_cellquad1),sizeof(float4)*nde_cell);
			cudaMalloc(PPTR(de_cellquad2),sizeof(float4)*nde_cell);
			cudaMalloc(PPTR(de_ptclsib),sizeof(int)*nde_ptcl);
			cudaMalloc(PPTR(de_ptclbasic),sizeof(float4)*nde_ptcl);


			CUDA_SAFE_CALL(cudaBindTexture(0,tex_ptclsib,de_ptclsib,sizeof(int)*nde_ptcl));
			CUDA_SAFE_CALL(cudaBindTexture(0,tex_ptclbasic,de_ptclbasic,sizeof(float4)*nde_ptcl));

			CUDA_SAFE_CALL(cudaBindTexture(0,tex_cellsibdau,de_cellsibdau,sizeof(int2)*nde_cell));
			CUDA_SAFE_CALL(cudaBindTexture(0,tex_cellbasic,de_cellbasic,sizeof(float4)*nde_cell));
			CUDA_SAFE_CALL(cudaBindTexture(0,tex_cellquad1,de_cellquad1,sizeof(float4)*nde_cell));
			CUDA_SAFE_CALL(cudaBindTexture(0,tex_cellquad2,de_cellquad2,sizeof(float4)*nde_cell));
		}
		if(1){ /* Texture bindings for the texture memory space */
			SetTexture(texcu_ptclsib,   int,   false,cudaFilterModePoint,cudaAddressModeClamp);
			SetTexture(texcu_ptclbasic, float4,false,cudaFilterModePoint,cudaAddressModeClamp);

			SetTexture(texcu_cellsibdau,int2,  false,cudaFilterModePoint,cudaAddressModeClamp);
			SetTexture(texcu_cellbasic, float4,false,cudaFilterModePoint,cudaAddressModeClamp);
			SetTexture(texcu_cellquad1, float4,false,cudaFilterModePoint,cudaAddressModeClamp);
			SetTexture(texcu_cellquad2, float4,false,cudaFilterModePoint,cudaAddressModeClamp);

			cudaMallocArray(&cu_cellsibdau,&texcu_cellsibdau.channelDesc,NTextureSize4,1);
			cudaMallocArray(&cu_cellbasic, &texcu_cellbasic.channelDesc,NTextureSize4,1);
			cudaMallocArray(&cu_cellquad1, &texcu_cellquad1.channelDesc,NTextureSize4,1);
			cudaMallocArray(&cu_cellquad2, &texcu_cellquad2.channelDesc,NTextureSize4,1);

			cudaMallocArray(&cu_ptclsib,&texcu_ptclsib.channelDesc,NTextureSize,1);
			cudaMallocArray(&cu_ptclbasic,&texcu_ptclbasic.channelDesc,NTextureSize,1);

			cudaBindTextureToArray(texcu_cellsibdau,cu_cellsibdau);
			cudaBindTextureToArray(texcu_cellbasic, cu_cellbasic);
			cudaBindTextureToArray(texcu_cellquad1, cu_cellquad1);
			cudaBindTextureToArray(texcu_cellquad2, cu_cellquad2);

			cudaBindTextureToArray(texcu_ptclsib,cu_ptclsib);
			cudaBindTextureToArray(texcu_ptclbasic,cu_ptclbasic);
		}
		initial_flag = 0;
	}
}

extern "C" void gputreeforce(particle *pos,int npos,
		TStruct *TREECELL, int ncell,TPtlStruct *grv,int ngrv,
		float rspheresq, int GpusPerNode, int myid, int nx){
	bool bcell,bptcl;

	if(initial_flag)INITIALIZE(myid,GpusPerNode,nx);

	{

		int ncellp,ngrvp;
		ncellp = ncell+1;
		ngrvp = ngrv+1;

		h_cellsibdau = (int2*) Malloc(sizeof(int2)*ncellp,PPTR(h_cellsibdau));
		h_cellbasic = (float4*)Malloc(sizeof(float4)*ncellp,PPTR(h_cellbasic));
		h_cellquad1 = (float4*)Malloc(sizeof(float4)*ncellp,PPTR(h_cellquad1));
		h_cellquad2 = (float4*)Malloc(sizeof(float4)*ncellp,PPTR(h_cellquad2));
		h_ptclsib = (int*)Malloc(sizeof(int)*ngrvp,PPTR(h_ptclsib));
		h_ptclbasic = (float4*)Malloc(sizeof(float4)*ngrvp,PPTR(h_ptclbasic));


		int ip;
		/* Building tree structure in the gpu memory */
		for(int i=0;i<ncell;i++){
			ip = i + 1;
			if(TREECELL[i].sibling == NULL)
				h_cellsibdau[ip].x = 0;
			else if((((TYPE*)(TREECELL[i].sibling))->type)==TYPE_TREE)
               	h_cellsibdau[ip].x = ((TStruct*)(TREECELL[i].sibling)-TREECELL)+1;
			else 
				h_cellsibdau[ip].x = -1*(int)(((TPtlStruct*)(TREECELL[i].sibling)-grv)+1);
			
			if(TREECELL[i].daughter == NULL)
				h_cellsibdau[ip].y = 0;
			else if((((TYPE*)(TREECELL[i].daughter))->type)==TYPE_TREE)
               	h_cellsibdau[ip].y = (int)((TStruct*)(TREECELL[i].daughter)-TREECELL)+1;
			else 
				h_cellsibdau[ip].y = -1*(int)(((TPtlStruct*)(TREECELL[i].daughter)-grv)+1);
			h_cellbasic[ip].x = TREECELL[i].mono[0];
			h_cellbasic[ip].y = TREECELL[i].mono[1];
			h_cellbasic[ip].z = TREECELL[i].mono[2];
			h_cellbasic[ip].w = TREECELL[i].dist_over_thetasq;
			h_cellquad1[ip].x = TREECELL[i].quad[0];
			h_cellquad1[ip].y = TREECELL[i].quad[1];
			h_cellquad1[ip].z = TREECELL[i].quad[2];
			h_cellquad1[ip].w = TREECELL[i].quad[3];
			h_cellquad2[ip].x = TREECELL[i].quad[4];
			h_cellquad2[ip].y = TREECELL[i].quad[5];
			h_cellquad2[ip].z = TREECELL[i].mass;
			h_cellquad2[ip].w = TREECELL[i].trQ;
		}


		for(int i=0;i<ngrv;i++){
			ip = i + 1;
			if(grv[i].sibling == NULL)
				h_ptclsib[ip] = 0;
			else if((((TYPE*)(grv[i].sibling))->type)==TYPE_TREE)
               	h_ptclsib[ip] = (int)(((TStruct*)(grv[i].sibling)-TREECELL)+1);
			else 
				h_ptclsib[ip] = -1*(int)(((TPtlStruct*)(grv[i].sibling)-grv)+1);
			h_ptclbasic[ip].x = grv[i].r[0];
			h_ptclbasic[ip].y = grv[i].r[1];
			h_ptclbasic[ip].z = grv[i].r[2];
			h_ptclbasic[ip].w = grv[i].mass;
		}


		if(ncellp <= NTextureSize4){
			bcell = true;
			cudaMemcpyToArray(cu_cellsibdau,0,0,h_cellsibdau,sizeof(int2)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpyToArray(cu_cellbasic,0,0,h_cellbasic,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpyToArray(cu_cellquad1,0,0,h_cellquad1,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpyToArray(cu_cellquad2,0,0,h_cellquad2,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);

		}
		else {
			bcell = false;
			if(ncellp > nde_cell){
				nde_cell = ncellp;
				cudaFree(de_cellsibdau);
				cudaFree(de_cellbasic);
				cudaFree(de_cellquad1);
				cudaFree(de_cellquad2);
				cudaMalloc(PPTR(de_cellsibdau),sizeof(int2)*nde_cell);
				cudaMalloc(PPTR(de_cellbasic),sizeof(float4)*nde_cell);
				cudaMalloc(PPTR(de_cellquad1),sizeof(float4)*nde_cell);
				cudaMalloc(PPTR(de_cellquad2),sizeof(float4)*nde_cell);

				CUDA_SAFE_CALL(cudaBindTexture(0,tex_cellsibdau,de_cellsibdau,sizeof(int2)*nde_cell));
				CUDA_SAFE_CALL(cudaBindTexture(0,tex_cellbasic,de_cellbasic,sizeof(float4)*nde_cell));
				CUDA_SAFE_CALL(cudaBindTexture(0,tex_cellquad1,de_cellquad1,sizeof(float4)*nde_cell));
				CUDA_SAFE_CALL(cudaBindTexture(0,tex_cellquad2,de_cellquad2,sizeof(float4)*nde_cell));
			}

			cudaMemcpy(de_cellsibdau,h_cellsibdau,sizeof(int2)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpy(de_cellbasic,h_cellbasic,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpy(de_cellquad1,h_cellquad1,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpy(de_cellquad2,h_cellquad2,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
		}

		if(ngrvp <= NTextureSize){
			bptcl = true;
			cudaMemcpyToArray(cu_ptclsib,0,0,h_ptclsib,sizeof(int)*ngrvp,cudaMemcpyHostToDevice);
			cudaMemcpyToArray(cu_ptclbasic,0,0,h_ptclbasic,sizeof(float4)*ngrvp,cudaMemcpyHostToDevice);

		}
		else {
			bptcl = false;
			if(ngrvp > nde_ptcl){
				nde_ptcl = ngrvp;
				cudaFree(de_ptclsib);
				cudaFree(de_ptclbasic);
				cudaMalloc(PPTR(de_ptclsib),sizeof(int)*nde_ptcl);
				cudaMalloc(PPTR(de_ptclbasic),sizeof(float4)*nde_ptcl);
				CUDA_SAFE_CALL(cudaBindTexture(0,tex_ptclsib,(void*)de_ptclsib,sizeof(int)*nde_ptcl));
				CUDA_SAFE_CALL(cudaBindTexture(0,tex_ptclbasic,de_ptclbasic,sizeof(float4)*nde_ptcl));
			}
			cudaMemcpy(de_ptclsib,h_ptclsib,sizeof(int)*ngrvp,cudaMemcpyHostToDevice);
			cudaMemcpy(de_ptclbasic,h_ptclbasic,sizeof(float4)*ngrvp,cudaMemcpyHostToDevice);
		}
		Free(h_ptclbasic);
		Free(h_ptclsib);
		Free(h_cellquad2);
		Free(h_cellquad1);
		Free(h_cellbasic);
		Free(h_cellsibdau);


		if(npos> nde_pos){
			nde_pos = npos;
			cudaFree(dpos);
			cudaMalloc(PPTR(dpos),sizeof(particle)*nde_pos);
		}

		CUDA_SAFE_CALL( cudaMemcpy(dpos,pos,sizeof(particle)*npos,cudaMemcpyHostToDevice));


		/*
		int numthreads = (mydev.maxGridSize[0]/2)*BlockSize;
		int nsplit = (npos+numthreads-1)/numthreads;

		for(int i=0;i<nsplit;i++){
			int mp = (numthreads<(npos-i*numthreads) ? numthreads: (npos-i*numthreads));
			int nblock  = (mp+BlockSize-1)/BlockSize;
			if(bcell == true && bptcl == true) {
				All_gpu_treeforce<true,true><<<nblock,BlockSize>>>(&(dpos[i*numthreads]),
						mp,rspheresq);
			}
			else if(bcell == true && bptcl == false){
				All_gpu_treeforce<true,false><<<nblock,BlockSize>>>(&(dpos[i*numthreads]),
						mp,rspheresq);
			}
			else if(bcell == false && bptcl == true){
				All_gpu_treeforce<false,true><<<nblock,BlockSize>>>(&(dpos[i*numthreads]),
						mp,rspheresq);
			}
			else {
				All_gpu_treeforce<false,false><<<nblock,BlockSize>>>(&(dpos[i*numthreads]),
						mp,rspheresq);
			}

		}
		*/
		int nblock = GetNumberofBlock(npos,BlockSize);
		if(0){ 
			int *work;
			cudaMalloc(PPTR(work),sizeof(int)*npos*MaxNList);
			if(bcell == true && bptcl == true) {
				Adv_All_gpu_treeforce<true,true><<<nblock,BlockSize>>>(dpos,npos,rspheresq,work);
			}
			else if(bcell == true && bptcl == false){
				Adv_All_gpu_treeforce<true,false><<<nblock,BlockSize>>>(dpos,npos,rspheresq,work);
			}
			else if(bcell == false && bptcl == true){
				Adv_All_gpu_treeforce<false,true><<<nblock,BlockSize>>>(dpos,npos,rspheresq,work);
			}
			else {
				Adv_All_gpu_treeforce<false,false><<<nblock,BlockSize>>>(dpos,npos,rspheresq,work);
			}

			cudaFree(work);
		}
		else {
			if(bcell == true && bptcl == true) {
				All_gpu_treeforce<true,true><<<nblock,BlockSize>>>(dpos,npos,rspheresq);
			}
			else if(bcell == true && bptcl == false){
				All_gpu_treeforce<true,false><<<nblock,BlockSize>>>(dpos,npos,rspheresq);
			}
			else if(bcell == false && bptcl == true){
				All_gpu_treeforce<false,true><<<nblock,BlockSize>>>(dpos,npos,rspheresq);
			}
			else {
				All_gpu_treeforce<false,false><<<nblock,BlockSize>>>(dpos,npos,rspheresq);
			}
		}
		CUDA_SAFE_CALL(cudaMemcpy(pos,dpos,sizeof(particle)*npos,cudaMemcpyDeviceToHost));

		if(NEEDLESS){
			cudaUnbindTexture(tex_cellsibdau);
			cudaUnbindTexture(tex_cellbasic);
			cudaUnbindTexture(tex_cellquad1);
			cudaUnbindTexture(tex_cellquad2);
			cudaUnbindTexture(tex_ptclsib);
			cudaUnbindTexture(tex_ptclbasic);
		}


	}

	return;
}

/* Put tiles over the 2-dimensional workspace on the npos x nrgv array */
__global__ void GGet_Gravity(int ngrv, float4 *grv, int npos, float4 *Pos, float3 *acc){
	__shared__ float4 shGrv[BlockSize];
	int gid = blockIdx.x*blockDim.x + threadIdx.x;
	int ntile,itile,gjd;
	float3 myAcc;
	float4 myPos = Pos[gid];

	if(gid<npos) myPos = Pos[gid];
	else myPos = make_float4(0,0,0,0);
	acc[gid] = make_float3(0,0,0);
	ntile = (ngrv+blockDim.x-1)/blockDim.x;
	for(itile=0;itile<ntile;itile++){
		gjd = itile*blockDim.x + threadIdx.x;
		if(gjd<ngrv) shGrv[threadIdx.x] = grv[gjd];
		else shGrv[threadIdx.x].w = 0;
#ifndef __DEVICE_EMULATION__
		__syncthreads();
#endif
		myAcc = make_float3(0,0,0);
 		for(int i=0;i<blockDim.x;i++){
			int j = (i+threadIdx.x)%blockDim.x; /* for memory coalescing */
			{
				float3 r;
				r.x = myPos.x-shGrv[j].x;
				r.y = myPos.y-shGrv[j].y;
				r.z = myPos.z-shGrv[j].z;
				float dist2 = r.x*r.x+r.y*r.y+r.z*r.z;
				float rr = (0.5f*log10(dist2)-LOG10RMIN)* invLOG10RMAXmLOG10RMIN;
				if(rr<1.){
					float fplmf = shGrv[j].w * tex1D(tex1,rr);
					myAcc.x += r.x*fplmf;
					myAcc.y += r.y*fplmf;
					myAcc.z += r.z*fplmf;
				}
			}

		}
		if(gid<npos) {
			acc[gid].x = acc[gid].x + myAcc.x;
			acc[gid].y = acc[gid].y + myAcc.y;
			acc[gid].z = acc[gid].z + myAcc.z;
		}
#ifndef __DEVICE_EMULATION__
		/*
		__syncthreads();
		*/
#endif
	}
}
__global__ void Get_Gravity(int ngrv, float4 *grv, int npos,float4 *Pos,float3 *acc){
    __shared__ float4 sharedPos[BlockSize];
    int gid = blockIdx.x*blockDim.x + threadIdx.x;
    int ibstart,ibfinal;
    int nblock = (ngrv+blockDim.x-1)/blockDim.x;
    float3 myAcc = {0,0,0};
    float4 myPos;
	if(gid < npos) myPos = Pos[gid];
	else myPos = make_float4(0,0,0,0);
   	ibstart = blockIdx.x;
   	ibfinal = ibstart + nblock;
   	for(int i=ibstart;i<ibfinal;i++){
   	    int gjd = (i%nblock)*blockDim.x + threadIdx.x;
   	    if(gjd < ngrv) sharedPos[threadIdx.x] = grv[gjd];
   	    else sharedPos[threadIdx.x] = make_float4(0,0,0,0);
#ifdef __DEVICE_EMULATION__
   	    __syncthreads();
#endif
   	    for(int j =0;j<blockDim.x;j++){
   	       float3 r;
   	       r.x = myPos.x - sharedPos[j].x;
   	       r.y = myPos.y - sharedPos[j].y;
   	       r.z = myPos.z - sharedPos[j].z;
   	       float dist2 = r.x*r.x+r.y*r.y+r.z*r.z;
   	       float rr = (0.5f*__log10f(dist2)-LOG10RMIN)* invLOG10RMAXmLOG10RMIN;
   	       if(rr < 1.){
   	       		float fplmf = sharedPos[j].w * tex1D(tex1,rr);
   	             myAcc.x += r.x*fplmf;
	             myAcc.y += r.y*fplmf;
   	           	 myAcc.z += r.z*fplmf;
   	   		}
   		}
 	}
	if(gid<npos) acc[gid] = myAcc;
}

static float4 *position, *nbody;
static float3 *accel;
static float4 *devposition, *devnbody;
static float3 *devaccel;
static int nposbuffer=NPosBuffer,nbodybuffer=NBodyBuffer;
static int dn_initial_flag = 1;
void DIRECTINTIALIZE(){
	cudaMalloc(PPTR(devposition),sizeof(float4)*nposbuffer);
	cudaMalloc(PPTR(devaccel),sizeof(float3)*nposbuffer);
	cudaMalloc(PPTR(devnbody),sizeof(float4)*nbodybuffer);
	dn_initial_flag = 0;
}

extern "C" void Direct_Nbody( particle *pos,int npos, TPtlStruct *grv,int ngrv,
		int GpusPerNode, int myid, int nx){
	if(initial_flag) INITIALIZE(myid,GpusPerNode,nx);
	if(dn_initial_flag) DIRECTINTIALIZE();
	int nblocks = (npos+BlockSize-1)/BlockSize;


	if(npos> nposbuffer){
		// There is no realloc in the CUDA
		nposbuffer = npos;
		cudaFree(devposition);
		cudaFree(devaccel);
		cudaMalloc(PPTR(devposition),sizeof(float4)*npos);/* Make space in the device memory */
		cudaMalloc(PPTR(devaccel),sizeof(float3)*npos);
	}
	if(ngrv > nbodybuffer){
		nbodybuffer = ngrv;
		cudaFree(devnbody);
		cudaMalloc(PPTR(devnbody),sizeof(float4)*ngrv);
	}
	position = (float4*)Malloc(sizeof(float4)*npos,PPTR(position));
	nbody = (float4*)Malloc(sizeof(float4)*ngrv,PPTR(nbody));
	accel = (float3*)Malloc(sizeof(float3)*npos,PPTR(accel)); /* There is no performance gain in reading
	the page locked memory */

	for(int i=0;i<npos;i++){
		position[i].x = pos[i].x;
		position[i].y = pos[i].y;
		position[i].z = pos[i].z;
	}
	for(int i=0;i<ngrv;i++){
		nbody[i].x = grv[i].r[0];
		nbody[i].y = grv[i].r[1];
		nbody[i].z = grv[i].r[2];
		nbody[i].w = 1;
	}
	CUDA_SAFE_CALL(cudaMemcpy(devposition,position,sizeof(float4)*npos,cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(devnbody,nbody,sizeof(float4)*ngrv,cudaMemcpyHostToDevice));

  	Get_Gravity<<<nblocks,BlockSize>>>(ngrv,devnbody,npos,devposition,devaccel);
	CUDA_SAFE_CALL(cudaMemcpy(accel,devaccel,sizeof(float3)*npos,cudaMemcpyDeviceToHost));

	for(int i=0;i<npos;i++){
		pos[i].ax = accel[i].x;
		pos[i].ay = accel[i].y;
		pos[i].az = accel[i].z;
	}
	Free(accel); 
	Free(nbody);Free(position);
}

/* The actual arguments are p, accel, and now */
#define MA_PARTICLEFORCE(p,now){\
	float4 point;\
	point = ptclbasic[now];\
	tmpx = p.x - point.x;\
	tmpy = p.y - point.y;\
	tmpz = p.z - point.z;\
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;\
	float rr = (0.5f*__log10f(dist2)-LOG10RMIN)* invLOG10RMAXmLOG10RMIN;\
	if(rr < 1.) \
	{\
		float fplmf = point.w * tex1D(tex1,rr);\
		accel.x += tmpx*fplmf;\
		accel.y += tmpy*fplmf;\
		accel.z += tmpz*fplmf;\
	}\
}
/* The actual arguments are now, dist2, tmpx,tmpy,tmpz, accel */
#define MA_CELLFORCE(p,now){\
	float rr = (0.5f*__log10f(dist2)-LOG10RMIN)* invLOG10RMAXmLOG10RMIN;\
	if( rr < 1.) \
	{\
		float4 fplmf;\
		fplmf = tex1D(tex4,rr);\
		float4 quad1,quad2;\
		quad1 = cellquad1[now];\
		quad2 = cellquad2[now];\
		/*\
		float trQ = (quad1.x+quad1.y+quad1.z);\
		*/\
		float trQ = quad2.w;\
		float qxx1 = quad1.x*tmpx+quad1.w*tmpy+quad2.x*tmpz;\
		float qxx2 = quad1.w*tmpx+quad1.y*tmpy+quad2.y*tmpz;\
		float qxx3 = quad2.x*tmpx+quad2.y*tmpy+quad1.z*tmpz;\
		float qxy = quad1.x*tmpx*tmpx+quad1.y*tmpy*tmpy+quad1.z*tmpz*tmpz+ \
			2.*(quad1.w*tmpx*tmpy+quad2.x*tmpx*tmpz+quad2.y*tmpy*tmpz);\
		float tmptmp = quad2.z *fplmf.x + trQ*fplmf.y + qxy*fplmf.z;\
		float twofplmf2 = fplmf.y + fplmf.y;\
		accel.x += tmpx*tmptmp+qxx1*twofplmf2;\
		accel.y += tmpy*tmptmp+qxx2*twofplmf2;\
		accel.z += tmpz*tmptmp+qxx3*twofplmf2;\
	}\
}

template<bool cucell,bool cuptcl>
__global__ void More_Adv_All_gpu_treeforce(particle *pp, int npos, float rspheresq,
		float4 *cellbasic, int2 *cellsibdau, int *ptclsib, float4 *ptclbasic,
		float4 *cellquad1, float4 *cellquad2){
	int now;
	float tmpx,tmpy,tmpz; 
	float dist2; 
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if(idx<npos){
		float4 accel = { 0.,0.,0.,0.};
		particle p = pp[idx];
		now = 1;
		while(now != 0){
			if(now >0 ){
				float4 mono;
				int2 Nowxy;
				mono = cellbasic[now];
				tmpx = p.x-mono.x;
				tmpy = p.y-mono.y;
				tmpz = p.z-mono.z;
				dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
				Nowxy = cellsibdau[now];
				if(cellopen(dist2,mono.w)==YES) 
					now = Nowxy.y;
				else {
					MA_CELLFORCE(p,now);
					now = Nowxy.x;
				}
			}
			else{
				now = -now;
				MA_PARTICLEFORCE(p,now);
				now = ptclsib[now];
			}
		}
		pp[idx].ax = accel.x; pp[idx].ay = accel.y; pp[idx].az = accel.z;
	}
}

extern "C" void gputreeforce2(particle *pos,int npos,
		TStruct *TREECELL, int ncell,TPtlStruct *grv,int ngrv,
		float rspheresq, int GpusPerNode, int myid, int nx){
	unsigned int nblock;

	if(initial_flag)INITIALIZE(myid,GpusPerNode,nx);

	{

		int ncellp,ngrvp;
		ncellp = ncell+1;
		ngrvp = ngrv+1;

		h_cellsibdau = (int2*) Malloc(sizeof(int2)*ncellp,PPTR(h_cellsibdau));
		h_cellbasic = (float4*)Malloc(sizeof(float4)*ncellp,PPTR(h_cellbasic));
		h_cellquad1 = (float4*)Malloc(sizeof(float4)*ncellp,PPTR(h_cellquad1));
		h_cellquad2 = (float4*)Malloc(sizeof(float4)*ncellp,PPTR(h_cellquad2));
		h_ptclsib = (int*)Malloc(sizeof(int)*ngrvp,PPTR(h_ptclsib));
		h_ptclbasic = (float4*)Malloc(sizeof(float4)*ngrvp,PPTR(h_ptclbasic));


		int ip;
		/* Building tree structure in the gpu memory */
		for(int i=0;i<ncell;i++){
			ip = i + 1;
			if(TREECELL[i].sibling == NULL)
				h_cellsibdau[ip].x = 0;
			else if((((TYPE*)(TREECELL[i].sibling))->type)==TYPE_TREE)
               	h_cellsibdau[ip].x = ((TStruct*)(TREECELL[i].sibling)-TREECELL)+1;
			else 
				h_cellsibdau[ip].x = -1*(int)(((TPtlStruct*)(TREECELL[i].sibling)-grv)+1);
			
			if(TREECELL[i].daughter == NULL)
				h_cellsibdau[ip].y = 0;
			else if((((TYPE*)(TREECELL[i].daughter))->type)==TYPE_TREE)
               	h_cellsibdau[ip].y = (int)((TStruct*)(TREECELL[i].daughter)-TREECELL)+1;
			else 
				h_cellsibdau[ip].y = -1*(int)(((TPtlStruct*)(TREECELL[i].daughter)-grv)+1);
			h_cellbasic[ip].x = TREECELL[i].mono[0];
			h_cellbasic[ip].y = TREECELL[i].mono[1];
			h_cellbasic[ip].z = TREECELL[i].mono[2];
			h_cellbasic[ip].w = TREECELL[i].dist_over_thetasq;
			h_cellquad1[ip].x = TREECELL[i].quad[0];
			h_cellquad1[ip].y = TREECELL[i].quad[1];
			h_cellquad1[ip].z = TREECELL[i].quad[2];
			h_cellquad1[ip].w = TREECELL[i].quad[3];
			h_cellquad2[ip].x = TREECELL[i].quad[4];
			h_cellquad2[ip].y = TREECELL[i].quad[5];
			h_cellquad2[ip].z = TREECELL[i].mass;
			h_cellquad2[ip].w = TREECELL[i].trQ;
			/*
			h_cellquad2[ip].w = 0.f;
			*/
		}


		for(int i=0;i<ngrv;i++){
			ip = i + 1;
			if(grv[i].sibling == NULL)
				h_ptclsib[ip] = 0;
			else if((((TYPE*)(grv[i].sibling))->type)==TYPE_TREE)
               	h_ptclsib[ip] = (int)(((TStruct*)(grv[i].sibling)-TREECELL)+1);
			else 
				h_ptclsib[ip] = -1*(int)(((TPtlStruct*)(grv[i].sibling)-grv)+1);
			h_ptclbasic[ip].x = grv[i].r[0];
			h_ptclbasic[ip].y = grv[i].r[1];
			h_ptclbasic[ip].z = grv[i].r[2];
			h_ptclbasic[ip].w = grv[i].mass;
		}


		{
			if(ncellp > nde_cell){
				nde_cell = ncellp;
				cudaFree(de_cellsibdau);
				cudaFree(de_cellbasic);
				cudaFree(de_cellquad1);
				cudaFree(de_cellquad2);
				cudaMalloc(PPTR(de_cellsibdau),sizeof(int2)*nde_cell);
				cudaMalloc(PPTR(de_cellbasic),sizeof(float4)*nde_cell);
				cudaMalloc(PPTR(de_cellquad1),sizeof(float4)*nde_cell);
				cudaMalloc(PPTR(de_cellquad2),sizeof(float4)*nde_cell);
				CUDA_SAFE_CALL(cudaBindTexture(0,tex_cellsibdau,de_cellsibdau,sizeof(int2)*nde_cell));
				CUDA_SAFE_CALL(cudaBindTexture(0,tex_cellbasic,de_cellbasic,sizeof(float4)*nde_cell));
				CUDA_SAFE_CALL(cudaBindTexture(0,tex_cellquad1,de_cellquad1,sizeof(float4)*nde_cell));
				CUDA_SAFE_CALL(cudaBindTexture(0,tex_cellquad2,de_cellquad2,sizeof(float4)*nde_cell));
			}

			cudaMemcpy(de_cellsibdau,h_cellsibdau,sizeof(int2)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpy(de_cellbasic,h_cellbasic,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpy(de_cellquad1,h_cellquad1,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpy(de_cellquad2,h_cellquad2,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
		}

		{
			if(ngrvp > nde_ptcl){
				nde_ptcl = ngrvp;
				cudaFree(de_ptclsib);
				cudaFree(de_ptclbasic);
				cudaMalloc(PPTR(de_ptclsib),sizeof(int)*nde_ptcl);
				cudaMalloc(PPTR(de_ptclbasic),sizeof(float4)*nde_ptcl);
				CUDA_SAFE_CALL(cudaBindTexture(0,tex_ptclsib,(void*)de_ptclsib,sizeof(int)*nde_ptcl));
				CUDA_SAFE_CALL(cudaBindTexture(0,tex_ptclbasic,de_ptclbasic,sizeof(float4)*nde_ptcl));
			}
			cudaMemcpy(de_ptclsib,h_ptclsib,sizeof(int)*ngrvp,cudaMemcpyHostToDevice);
			cudaMemcpy(de_ptclbasic,h_ptclbasic,sizeof(float4)*ngrvp,cudaMemcpyHostToDevice);
		}
		Free(h_ptclbasic);
		Free(h_ptclsib);
		Free(h_cellquad2);
		Free(h_cellquad1);
		Free(h_cellbasic);
		Free(h_cellsibdau);


		if(npos> nde_pos){
			nde_pos = npos;
			cudaFree(dpos);
			cudaMalloc(PPTR(dpos),sizeof(particle)*nde_pos);
		}

		CUDA_SAFE_CALL( cudaMemcpy(dpos,pos,sizeof(particle)*npos,cudaMemcpyHostToDevice));

		nblock = GetNumberofBlock(npos,BlockSize);
		{ 
			More_Adv_All_gpu_treeforce<true,true><<<nblock,BlockSize>>>(dpos,npos,rspheresq,
				de_cellbasic, de_cellsibdau, de_ptclsib, de_ptclbasic, de_cellquad1, de_cellquad2);

		}


		CUDA_SAFE_CALL(cudaMemcpy(pos,dpos,sizeof(particle)*npos,cudaMemcpyDeviceToHost));



		if(NEEDLESS){
			cudaUnbindTexture(tex_cellsibdau);
			cudaUnbindTexture(tex_cellbasic);
			cudaUnbindTexture(tex_cellquad1);
			cudaUnbindTexture(tex_cellquad2);
			cudaUnbindTexture(tex_ptclsib);
			cudaUnbindTexture(tex_ptclbasic);
		}


	}

	return;
}
/*
__global__ void Get_Gravity(int ngrv, float4 *grv, int npos,float4 *Pos,float4 *acc){
    __shared__ float4 sharedPos[BlockSize];
    int gid = blockIdx.x*blockDim.x + threadIdx.x;
    int ibstart,ibfinal;
    int nblock = (ngrv+blockDim.x-1)/blockDim.x;
    float4 myAcc = {0,0,0,0};
    float4 myPos;
    myPos = Pos[gid];
    ibstart = blockIdx.x;
    ibfinal = ibstart + nblock;
    for(int i=ibstart;i<ibfinal;i++){
        int gjd = (i%nblock)*blockDim.x + threadIdx.x;
        if(gjd < ngrv) sharedPos[threadIdx.x] = grv[gjd];
        else sharedPos[threadIdx.x] = make_float4(0,0,0,0);
#ifdef __DEVICE_EMULATION__
        __syncthreads();
#endif
        for(int j =0;j<blockDim.x;j++){
            float3 r;
            r.x = myPos.x - sharedPos[j].x;
            r.y = myPos.y - sharedPos[j].y;
            r.z = myPos.z - sharedPos[j].z;
            float dist2 = r.x*r.x+r.y*r.y+r.z*r.z;
            float rr = (0.5f*__log10f(dist2)-LOG10RMIN)* invLOG10RMAXmLOG10RMIN;
            if(rr < 1.){
                float fplmf = sharedPos[j].w * tex1D(tex1,rr);
                myAcc.x += r.x*fplmf;
                myAcc.y += r.y*fplmf;
                myAcc.z += r.z*fplmf;
            }
        }
    }
    acc[gid] = myAcc;
}
*/
