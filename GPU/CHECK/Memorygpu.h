#define MAX_MALLOC 100000
#define PPTR(A) ((void **)(&A))
#define MemAlign 256 /* bytes */
/* for particle num > 2000000000 */
#define MPI_INT8	MPI_LONG_LONG
typedef  long long INT8;
INT8 GPUMake_Total_Memory(INT8);
INT8 GPUCheckAvailableMemory();
void *GPUCalloc(INT8,INT8,void **);
void *GPUMalloc(INT8, void **);
void *GPURealloc(void *, INT8,void **);
INT8 GPUfreespace();
void GPUfreelast(void *);
INT8 GPUptrsize(void *);
void *GPUresizelast(void *,INT8);
void GPUNumMemStack();
void GPUFreeRightNumMemStack();
void GPULastSwitchPointer(void **);
void GPUMemSwitchPointer(void **,void **);
INT8 GPUCurMemStack();
void GPUInitialOldMemStack(INT8);
void GPUStackPosition(void *a);
void GPUFree(void *a);
#define MEGABYTE 1048576L
#ifndef NMEG
#define NMEG 100L
#endif
#ifndef INFINITY 
#define INFINITY -1L
#endif
/*
#define INFINITY -1L
*/

#ifndef MEMMAIN
#define Realloc(A,B) Realloc(A,B,PPTR(A))
#define MemSwitchPointer(A,B) MemSwitchPointer(PPTR(A),PPTR(B))
#endif
