#define MAX_MALLOC 100000
#define PPTR(A) ((void **)(&A))
/* for particle num > 2000000000 */
#define MPI_INT8	MPI_LONG_LONG
typedef  long long INT8;
extern "C" INT8 Make_Total_Memory();
extern "C" INT8 CheckAvailableMemory();
extern "C" void *Calloc(INT8,INT8,void **);
extern "C" void *Malloc(INT8, void **);
extern "C" void *Realloc(void *, INT8,void **);
extern "C" INT8 freespace();
extern "C" void freelast(void *);
extern "C" INT8 ptrsize(void *);
extern "C" void *resizelast(void *,INT8);
extern "C" void NumMemStack();
extern "C" void FreeRightNumMemStack();
extern "C" void LastSwitchPointer(void **);
extern "C" void MemSwitchPointer(void **,void **);
extern "C" INT8 CurMemStack();
extern "C" void InitialOldMemStack(INT8);
extern "C" void StackPosition(void *a);
extern "C" void Free(void *a);
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
