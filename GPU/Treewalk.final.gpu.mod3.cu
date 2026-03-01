/* This version seems ok. 21/02/2008 */
/* This code adopts a constant memory space for TStruct & TPtlStruct */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include "cutil.h"
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

texture<float,1,cudaReadModeNormalizedFloat> tex1;
texture<float4,1,cudaReadModeNormalizedFloat> tex4;
static cudaArray  *cuArray1, *cuArray4;
float *force1;
float4 *hforce4;


texture<int2  ,1,cudaReadModeElementType> tex_cellsibdau;
texture<int   ,1,cudaReadModeElementType> tex_ptclsib;
texture<float4,1,cudaReadModeElementType> tex_cellbasic;
texture<float4,1,cudaReadModeElementType> tex_cellquad1;
texture<float4,1,cudaReadModeElementType> tex_cellquad2;
texture<float4,1,cudaReadModeElementType> tex_ptclbasic;

#define cellopen(dist2,dist_over_thetasq) (dist2<dist_over_thetasq ? YES:NO)

#define PARTICLEFORCE(p,now,rspheresq){\
	float4 point;\
	point = tex1Dfetch(tex_ptclbasic,now);\
	tmpx = p.x - point.x;\
	tmpy = p.y - point.y;\
	tmpz = p.z - point.z;\
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;\
	float rr = __fdividef(0.5f*__log10f(dist2)-LOG10RMIN, LOG10RMAXmLOG10RMIN);\
	if(rr>0. && rr < 1.) \
	{\
		float fplmf = point.w * tex1D(tex1,rr);\
		accel.x += tmpx*fplmf;\
		accel.y += tmpy*fplmf;\
		accel.z += tmpz*fplmf;\
	}\
}
#define CELLFORCE(p,now,rspheresq){\
	float rr = __fdividef(0.5f*__log10f(dist2)-LOG10RMIN, LOG10RMAXmLOG10RMIN);\
	if(rr>0. && rr < 1.) \
	{\
		float4 fplmf;\
		fplmf = tex1D(tex4,rr);\
		float4 quad1,quad2;\
		quad1 = tex1Dfetch(tex_cellquad1,now);\
		quad2 = tex1Dfetch(tex_cellquad2,now);\
		float Qxx = quad1.x;\
		float Qyy = quad1.y;\
		float Qzz = quad1.z;\
		float Qxy = quad1.w;\
		float Qxz = quad2.x;\
		float Qyz = quad2.y;\
		float trQ = (Qxx+Qyy+Qzz);\
		float qxx1 = Qxx*tmpx+Qxy*tmpy+Qxz*tmpz;\
		float qxx2 = Qxy*tmpx+Qyy*tmpy+Qyz*tmpz;\
		float qxx3 = Qxz*tmpx+Qyz*tmpy+Qzz*tmpz;\
		float qxy = Qxx*tmpx*tmpx+Qyy*tmpy*tmpy+Qzz*tmpz*tmpz+ \
			2.*(Qxy*tmpx*tmpy+Qxz*tmpx*tmpz+Qyz*tmpy*tmpz);\
		float tmptmp = quad2.z *fplmf.x + trQ*fplmf.y + qxy*fplmf.z;\
		float twofplmf2 = 2.*fplmf.y;\
		accel.x += tmpx*tmptmp+qxx1*twofplmf2;\
		accel.y += tmpy*tmptmp+qxx2*twofplmf2;\
		accel.z += tmpz*tmptmp+qxx3*twofplmf2;\
	}\
}
//template<bool cellcuar,bool ptclcuar>
__global__ void All_gpu_treeforce11(particle *pp, int npos, float rspheresq){
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
				mono = tex1Dfetch(tex_cellbasic,now);
				tmpx = p.x-mono.x;
				tmpy = p.y-mono.y;
				tmpz = p.z-mono.z;
				dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
				Nowxy = tex1Dfetch(tex_cellsibdau,now);
				if(cellopen(dist2,mono.w)==YES) 
					now = Nowxy.y;
				else {
					CELLFORCE(p,now,rspheresq);
					now = Nowxy.x;
				}
			}
			else{
				now = -now;
				PARTICLEFORCE(p,now,rspheresq);
				now = tex1Dfetch(tex_ptclsib,now);
			}
		}
		pp[idx].ax = accel.x; pp[idx].ay = accel.y; pp[idx].az = accel.z;
	}
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

particle *dpos;

static int initialflag=1;

void CUDA_TERMINATE(int myid){
	fprintf(stdout,"P%d Terminating Cuda Device\n",myid);
	cudaFreeArray(cuArray1);
	cudaFreeArray(cuArray4);
	cudaUnbindTexture(tex1);
	cudaUnbindTexture(tex4);
}

extern "C" void gputreeforce(particle *pos,int npos,
		TStruct *TREECELL, int ncell,TPtlStruct *grv,int ngrv,
		float rspheresq, int GpusPerNode, int myid, int nx){
	unsigned int blockSize;
	int ndev;



	if(initialflag){
		float *hforce1;
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
		hforce1 = (float*)Malloc(sizeof(float)*NTextureSize,PPTR(hforce1));
		hforce4 = (float4*)Malloc(sizeof(float4)*NTextureSize4,PPTR(hforce4));

		for(int i=0;i<NTextureSize;i++){
			float r = LOG10RMIN+(LOG10RMAX-LOG10RMIN)*(i+0.5)/(float) NTextureSize;
			float4 aforce4;
			r = powf(10.,r);
			aforce4 = pmcf4(r,nx);
			hforce1[i] = aforce4.x;
		}
		for(int i=0;i<NTextureSize4;i++){
			float r = LOG10RMIN+(LOG10RMAX-LOG10RMIN)*(i+0.5)/(float) NTextureSize4;
			r = powf(10.,r);
			hforce4[i] = pmcf4(r,nx);
		}

		tex4.normalized=tex1.normalized= YES;
		tex4.filterMode = tex1.filterMode = cudaFilterModeLinear;
		tex4.addressMode[0] = tex1.addressMode[0] = cudaAddressModeClamp;

		CUDA_SAFE_CALL( cudaMallocArray(&cuArray1,&tex1.channelDesc,NTextureSize,1));
		CUDA_SAFE_CALL( cudaMallocArray(&cuArray4,&tex4.channelDesc,NTextureSize4,1));

		CUDA_SAFE_CALL( cudaBindTextureToArray(tex1,cuArray1));
		CUDA_SAFE_CALL( cudaBindTextureToArray(tex4,cuArray4));


		CUDA_SAFE_CALL( cudaMemcpyToArray(cuArray1,0,0,hforce1,
					sizeof(float)*NTextureSize,cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL( cudaMemcpyToArray(cuArray4,0,0,hforce4,
					sizeof(float4)*NTextureSize4,cudaMemcpyHostToDevice));

		Free(hforce4); Free(hforce1);

		INT8 size = GPUCheckAvailableMemory();
		printf("P%d has %d bytes gpu memory pool\n",myid,size);

		initialflag = 0;
	}



	{
		cudaArray  *cu_cellsibdau,*cu_ptclsib,*cu_cellbasic,*cu_cellquad1,*cu_cellquad2,*cu_ptclbasic;
		int2 *de_cellsibdau; int *de_ptclsib;
		float4 *de_cellbasic,*de_cellquad1,*de_cellquad2,*de_ptclbasic;

		int2 *h_cellsibdau;int *h_ptclsib;
		float4 *h_cellbasic,*h_cellquad1,*h_cellquad2,*h_ptclbasic;
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
			h_cellquad2[ip].w = 0.f;
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


		/*
		if(ncell <= NTextureSize4){
		*/
		if(0){
			cudaMallocArray(&cu_cellsibdau,&tex_cellsibdau.channelDesc,ncellp,1);
			cudaMallocArray(&cu_cellbasic, &tex_cellbasic.channelDesc,ncellp,1);
			cudaMallocArray(&cu_cellquad1, &tex_cellquad1.channelDesc,ncellp,1);
			cudaMallocArray(&cu_cellquad2, &tex_cellquad2.channelDesc,ncellp,1);
			cudaBindTextureToArray(tex_cellsibdau,cu_cellsibdau);
			cudaBindTextureToArray(tex_cellbasic, cu_cellbasic);
			cudaBindTextureToArray(tex_cellquad1, cu_cellquad1);
			cudaBindTextureToArray(tex_cellquad2, cu_cellquad2);
			cudaMemcpyToArray(cu_cellsibdau,0,0,h_cellsibdau,sizeof(int2)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpyToArray(cu_cellbasic,0,0,h_cellbasic,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpyToArray(cu_cellquad1,0,0,h_cellquad1,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpyToArray(cu_cellquad2,0,0,h_cellquad2,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
		}
		else {
			/*
			cudaMalloc(PPTR(de_cellsibdau),sizeof(int2)*ncellp);
			cudaMalloc(PPTR(de_cellbasic),sizeof(float4)*ncellp);
			cudaMalloc(PPTR(de_cellquad1),sizeof(float4)*ncellp);
			cudaMalloc(PPTR(de_cellquad2),sizeof(float4)*ncellp);
			*/
			de_cellsibdau = (int2*)GPUMalloc(sizeof(int2)*ncellp,PPTR(de_cellsibdau));
			de_cellbasic = (float4*)GPUMalloc(sizeof(float4)*ncellp,PPTR(de_cellbasic));
			de_cellquad1 = (float4*)GPUMalloc(sizeof(float4)*ncellp,PPTR(de_cellquad1));
			de_cellquad2 = (float4*)GPUMalloc(sizeof(float4)*ncellp,PPTR(de_cellquad2));

			cudaMemcpy(de_cellsibdau,h_cellsibdau,sizeof(int2)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpy(de_cellbasic,h_cellbasic,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpy(de_cellquad1,h_cellquad1,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
			cudaMemcpy(de_cellquad2,h_cellquad2,sizeof(float4)*ncellp,cudaMemcpyHostToDevice);
			cudaBindTexture(0,tex_cellsibdau,de_cellsibdau,sizeof(int2)*ncellp);
			cudaBindTexture(0,tex_cellbasic,de_cellbasic,sizeof(float4)*ncellp);
			cudaBindTexture(0,tex_cellquad1,de_cellquad1,sizeof(float4)*ncellp);
			cudaBindTexture(0,tex_cellquad2,de_cellquad2,sizeof(float4)*ncellp);
		}

		/*
		if(ngrv <= NTextureSize4){
		*/
		if(0){
			cudaMallocArray(&cu_ptclsib,&tex_ptclsib.channelDesc,ngrvp,1);
			cudaMallocArray(&cu_ptclbasic,&tex_ptclbasic.channelDesc,ngrvp,1);
			cudaBindTextureToArray(tex_ptclsib,cu_ptclsib);
			cudaBindTextureToArray(tex_ptclbasic,cu_ptclbasic);
			cudaMemcpyToArray(cu_ptclsib,0,0,h_ptclsib,sizeof(int)*ngrvp,cudaMemcpyHostToDevice);
			cudaMemcpyToArray(cu_ptclbasic,0,0,h_ptclbasic,sizeof(float4)*ngrvp,cudaMemcpyHostToDevice);
		}
		else {
			/*
			cudaMalloc(PPTR(de_ptclsib),sizeof(int)*ngrvp);
			cudaMalloc(PPTR(de_ptclbasic),sizeof(float4)*ngrvp);
			*/
			de_ptclsib = (int*)GPUMalloc(sizeof(int)*ngrvp,PPTR(de_ptclsib));
			de_ptclbasic = (float4*)GPUMalloc(sizeof(float4)*ngrvp,PPTR(de_ptclbasic));
			cudaMemcpy(de_ptclsib,h_ptclsib,sizeof(int)*ngrvp,cudaMemcpyHostToDevice);
			cudaMemcpy(de_ptclbasic,h_ptclbasic,sizeof(float4)*ngrvp,cudaMemcpyHostToDevice);
			cudaBindTexture(0,tex_ptclsib,de_ptclsib,sizeof(int)*ngrvp);
			cudaBindTexture(0,tex_ptclbasic,de_ptclbasic,sizeof(float4)*ngrvp);
		}
		Free(h_ptclbasic);
		Free(h_ptclsib);
		Free(h_cellquad2);
		Free(h_cellquad1);
		Free(h_cellbasic);
		Free(h_cellsibdau);


		/*
		CUDA_SAFE_CALL( cudaMalloc(PPTR(dpos),sizeof(particle)*npos));
		*/
		dpos = (particle *)GPUMalloc(sizeof(particle)*npos,PPTR(dpos));
		CUDA_SAFE_CALL( cudaMemcpy(dpos,pos,sizeof(particle)*npos,cudaMemcpyHostToDevice));

		blockSize = GetBlockSize1D(npos,threadSize);

		All_gpu_treeforce11<<<blockSize,threadSize>>>(dpos,npos,rspheresq);


		CUDA_SAFE_CALL( cudaMemcpy(pos,dpos,sizeof(particle)*npos,cudaMemcpyDeviceToHost));


		GPUFree(dpos);
		/*
		cudaFree(dpos);
		*/


		if(NEEDLESS){
			cudaUnbindTexture(tex_cellsibdau);
			cudaUnbindTexture(tex_cellbasic);
			cudaUnbindTexture(tex_cellquad1);
			cudaUnbindTexture(tex_cellquad2);
			cudaUnbindTexture(tex_ptclsib);
			cudaUnbindTexture(tex_ptclbasic);
		}

		/*
		cudaFree(de_cellbasic);
		cudaFree(de_cellsibdau);
		cudaFree(de_cellquad1);
		cudaFree(de_cellquad2);
		cudaFree(de_ptclsib);
		cudaFree(de_ptclbasic);
		*/
		GPUFree(de_ptclbasic);
		GPUFree(de_ptclsib);
		GPUFree(de_cellquad2);
		GPUFree(de_cellquad1);
		GPUFree(de_cellbasic);
		GPUFree(de_cellsibdau);



		/*
		cudaFreeArray(cu_cellbasic);
		cudaFreeArray(cu_cellsibdau);
		cudaFreeArray(cu_cellquad1);
		cudaFreeArray(cu_cellquad2);
		cudaFreeArray(cu_ptclsib);
		cudaFreeArray(cu_ptclbasic);
		*/

	}

	return;
}
