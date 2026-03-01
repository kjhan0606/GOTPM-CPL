/* 64 bits supported by Juhan Kim 30/08/02 */
#define MAX_MALLOC 100000
#define PPTR(A) ((void **)(&A))
long Make_Total_Memory(long);
void *Calloc(long,long,void **);
void *Malloc(long, void **);
void *Realloc(void *, long,void **);
long freespace();
void freelast(void *);
long ptrsize(void *);
void *resizelast(void *,long);
void NumMemStack();
void FreeRightNumMemStack();
void LastSwitchPointer(void **);
void MemSwitchPointer(void **,void **);
long CurMemStack();
void InitialOldMemStack(long);
#define MEGABYTE 1048576L
#ifndef NMEG
#define NMEG 1000L
#endif
#define INFINITY -1

#ifndef MEMMAIN
#define Realloc(A,B) Realloc(A,B,PPTR(A))
#define MemSwitchPointer(A,B) MemSwitchPointer(PPTR(A),PPTR(B))
#endif
