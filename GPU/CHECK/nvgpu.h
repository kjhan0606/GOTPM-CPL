/*
#include "cutil.h"
*/
#define threadSize 32
/* Maximum number of cells that are gravitationally contacting the position */
#define MaxNContactCell 1000
/* Maximum number of particles that are gravitationally contacting the position */
#define MaxNContactptl  2000
#define MaxNList (MaxNContactCell + MaxNContactptl)
#define GPUSPERNODE 4
#define CPUSPERGPU 1
#define  GPU_Sync CUDA_SAFE_CALL(cudaThreadSynchronize())
#define NCONST 16384
#define cuNSPLINE (1<<13)


#define NForceTexture (1<<12)

#define NTextureSize (1<<12)
#define NTextureSize4 (1<<12)

#define NTextureBuffer (1<<20)
#define NPosBuffer (1<<20)
#define NBodyBuffer (1<<13)

#define BlockSize threadSize



