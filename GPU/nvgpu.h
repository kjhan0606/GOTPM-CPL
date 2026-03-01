/*
#include "cutil.h"
*/
#define threadSize 256
#define MaxthreadSize 512
#define GPUSPERNODE 2
#define CPUSPERGPU 2
#define  GPU_Sync CUDA_SAFE_CALL(cudaThreadSynchronize())
#define NCONST 16384
#define cuNSPLINE (1<<13)
#define NTextureSize (1<<13)
#define NTextureSize4 (1<<(13-1))
#define NTextureBuffer (1<<20)
#define NPosBuffer (1<<20)
#define NBodyBuffer (1<<20)
