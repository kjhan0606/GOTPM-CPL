#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#include<math.h>
#include "cutil.h"
#define MEMMAIN
#include "Memorygpu.h"
typedef struct memorystrcut {
	INT8 Size;
	void *Starting;
	void **PtrToVariable;
} GPUmemorystruct;
static GPUmemorystruct Memory[MAX_MALLOC];
static void *GPUFREE=NULL;
static INT8 GPUtsize=0;
static INT8 GPUCurrent_Stack=0,GPUCurrent_Stack_org=0;
static void *GPUSTART_TOTAL_MEMORY=NULL;
static void *GPULastEnd=NULL;
INT8 GPUCheckAvailableMemory(){
	INT8 i;
	INT8 size;
	size = 0;
	for(i=0;i<GPUCurrent_Stack;i++) size += Memory[i].Size;
	size = GPUtsize - size;
	return size;
}
INT8 GPUMake_Total_Memory(INT8 size){
	INT8 i;
	GPUCurrent_Stack = 0;
	GPUFREE = NULL;
	for(i=0;i<MAX_MALLOC;i++){
		Memory[i].Starting = NULL;
		Memory[i].Size = 0;
		Memory[i].PtrToVariable = NULL;
	}
	cudaMalloc((void**)&GPUFREE,size);
	GPULastEnd = (void*)((char *)GPUFREE + size);
	GPUtsize = size;
	GPUSTART_TOTAL_MEMORY = GPUFREE;
	if(GPUFREE == NULL){
		size /= MEGABYTE;
		fprintf(stderr,"Error when initializing %d Mbytes memory space\n",size);
		return 0;
	}
	else {
		size /= MEGABYTE;
		fprintf(stderr,"Initializing %d Mbytes memory space\n",size);
		return 1;
	}

}
INT8 GPUptrsize(void *p){
	INT8 i;
	for(i=0;i<GPUCurrent_Stack;i++){
		if(Memory[i].Starting == p) return Memory[i].Size;
	}
	return NULL;
}
INT8 GPUGPUfreespace(){
	INT8 GPUfree;
	if(GPUCurrent_Stack == 0){
		return GPUtsize;
	}
	else{
		GPUfree= GPUtsize - ((INT8)((char *)Memory[GPUCurrent_Stack-1].Starting
				-(char *)GPUSTART_TOTAL_MEMORY) + Memory[GPUCurrent_Stack-1].Size);
		return GPUfree;
	}
}
void *GPUMalloc(INT8 size, void **src){
	void *value;
	// This is to enforce the GPU memory alignment.
	// Please refer to pg. 51
//	size = (size/MemAlign + 1) *MemAlign;
	size = (int)ceil((double)size/(double)MemAlign)*MemAlign;

	if(size == 0) {
		printf("Warning : size %ld is being allocated\n",size);
		return NULL;
	}
	if(size == INFINITY){
		if(GPUCurrent_Stack == 0){
			size = GPUtsize;
		}
		else {
			size = GPUtsize - ( (INT8)((char *)Memory[GPUCurrent_Stack-1].Starting-
					(char *)GPUSTART_TOTAL_MEMORY) + Memory[GPUCurrent_Stack-1].Size);
		}
		printf("change allocated memory size from infinity to %ld\n",size);
	}
	Memory[GPUCurrent_Stack].Size = size;
	Memory[GPUCurrent_Stack].Starting = GPUFREE;
	Memory[GPUCurrent_Stack].PtrToVariable = src;
	value = GPUFREE;
	GPUFREE = (void *)((char *)GPUFREE + size);
	if(GPUFREE > (void *)((char *)GPUSTART_TOTAL_MEMORY + GPUtsize)) {
		INT8 tmpGPUtsize;
		GPUFREE = (void *)((char *)GPUFREE-size); /* restore the original mem */
		fprintf(stderr,"Error allocating memory %d bytes in GPUMalloc()\n",size);
		fprintf(stderr,"need more %d bytes \n",(int)(size - 
					(GPUtsize-((char *)GPUFREE-(char *)GPUSTART_TOTAL_MEMORY))));
        tmpGPUtsize = GPUtsize/NMEG;
		fprintf(stderr,"total memory is %d Mbytes\n",tmpGPUtsize);
		fprintf(stderr," GPUFREE=%p GPUSTART_TOTAL_MEMORY=%p GPUtsize=%ld\n",GPUFREE,
				GPUSTART_TOTAL_MEMORY,GPUtsize);
		fprintf(stderr," GPUCurrent_Stack=%ld\n",GPUCurrent_Stack);



		value = NULL;
		cudaMalloc(src,size);
		value = src;
		if(src == NULL) {
			fprintf(stderr,"Cannot allocate memory any more\n");
			exit(0);
		}
		else {
			return value;
		}
	}
	GPUCurrent_Stack++;
	return value;
}
__global__ void initialize(char *a,INT8 size){
	int idx = threadIdx.x;
	if(idx==0)
		for(INT8 i=0;i<size;i++) a[i] =0;
}
__global__ void GPUmemmove(void *t, void *s, INT8 size){
	int idx = threadIdx.x;
	char *target, *src;
	target = (char *) t;
	src = (char *) s;
	if(idx==0)
		for(INT8 i=0;i<size;i++)
			*(target++) = *(src++);
}
void *GPUCalloc(INT8 num, INT8 size,void **src){
	char *value;
//	size = (size/MemAlign + 1) *MemAlign;
	size = (int)ceil((double)size/(double)MemAlign)*MemAlign;

	size = size *num;
	value = (char *)GPUMalloc(size,src);
	/*
	for(i=0;i<size;i++){
		value[i] = 0;
	}
	*/
	initialize<<<1,32>>>(value,size);

	return (void *) value;
}
void GPUEraseMemory(INT8 srcn){
	INT8 size=0,esize;
	INT8 i;
	void *from;
	GPUFREE = (void *)((char *)GPUFREE - Memory[srcn].Size);
	*(Memory[srcn].PtrToVariable) = NULL;
	esize = Memory[srcn].Size;
	if(srcn == GPUCurrent_Stack-1){
		GPUCurrent_Stack--;
		return;
	}
	else {
		from=Memory[srcn+1].Starting;
		for(i=srcn+1;i<GPUCurrent_Stack;i++){
			Memory[i-1].Starting = (void *)((char *)Memory[i].Starting-esize);
			Memory[i-1].Size = Memory[i].Size;
			Memory[i-1].PtrToVariable = Memory[i].PtrToVariable;
			*(Memory[i-1].PtrToVariable) = Memory[i-1].Starting;
			size += Memory[i].Size;
		}
		GPUCurrent_Stack--;
		GPUmemmove<<<1,32>>>(Memory[srcn].Starting,from,size);
	}
}
void GPUFree(void *a){
	INT8 i;
	if(a == NULL){
	}
	else {
		for(i=0;i<GPUCurrent_Stack;i++){
			if((char *)a == Memory[i].Starting){
				GPUEraseMemory(i);
				return;
			}
		}
		if(i==GPUCurrent_Stack) {
			cudaFree(a);
			a = NULL;
		}
	}
}
__global__ void cudamemcpy(void *dest, void *src, INT8 size){
	char *d,*s;
	int idx = threadIdx.x;
	if(idx==0){
		d = (char*)dest;
		s = (char*)src;
		for(INT8 i=0;i<size;i++) d[i] = s[i];
	}
}
__global__ void cudaleftpush(void *dest, void *src, INT8 size){
	char *d,*s;
	int idx=threadIdx.x;
	if(idx==0){
		d = (char*)dest;
		s = (char*)src;
		for(INT8 i=0;i<size;i++){
			*(d--) = *(s--);
		}
	}
}
void *GPURealloc(void *a, INT8 size,void **src){
	INT8 i,j,diffsize,nmove;
	void *from,*to;
	void *ptr2;
	void *value;
	char *tmpptr;
	/* freeing the memory if size == 0 */
	if(size <0) {
		fprintf(stderr,"Error,,, minus size reallocation \n");
		fflush(stderr);
		exit(9999);
	}
	if(size == 0) {
		GPUFree(a);
		return NULL;
	}
	if(a == NULL){
		printf("Strange in Realloc ");
		printf("%d  size of %d\n",src,size);fflush(stdout);
		a = (void *)GPUMalloc(size,src);
		return a;
	}
	for(i=0;i<GPUCurrent_Stack;i++){
		if((char *)a == Memory[i].Starting){
			break;
		}
	}
	/* If this memory was allocated using malloc */
	if(i == GPUCurrent_Stack) {
		fprintf(stderr,"Error reallocating device memory for CUDA\n");
		exit(9999);
		return (void *) realloc(a,size);
	}

//	size = (size/MemAlign + 1) *MemAlign;
	size = (int)ceil((double)size/(double)MemAlign)*MemAlign;

	diffsize = size-Memory[i].Size;
	if(diffsize == 0) return a;
	/* if there is no more free space */
	/*
	if(GPUtsize <  (INT8)((char *)SMalloc[i]-
		(char *)GPUSTART_TOTAL_MEMORY + size)) {
	*/
	if(GPUtsize < (INT8)((char *)GPUFREE-(char *)GPUSTART_TOTAL_MEMORY + diffsize)) {
		fprintf(stderr,"Error allocating memory %d bytes in Realloc()\n",size);
		fprintf(stderr,"need more %d bytes \n",(int)(size - 
					(GPUtsize-((char *)GPUFREE-(char *)GPUSTART_TOTAL_MEMORY))));
		fflush(stderr);
		value = NULL;
		cudaMalloc(src,size);
		value = *src;
		if(value == NULL) {
			fprintf(stderr,"Cannot reallocate any more memory\n");
			fprintf(stderr,"Finish this program !!!!!!!!!!!!!!!! \n");
			exit(0);
		}
		ptr2 = Memory[i].Starting;

		/*
		memcpy(value,ptr2,Memory[i].Size);
		*/
		cudamemcpy<<<1,32>>>(value,ptr2,Memory[i].Size);

		/*
		GPUCurrent_Stack--;
		*/
		GPUFree(a);
		return value;
	}
	if(i != GPUCurrent_Stack-1){
		if(size > Memory[i].Size) {
			from = (void *)((char *)GPUFREE-1);
			to = Memory[i+1].Starting;
			j = 0;
			tmpptr = (char *)Memory[GPUCurrent_Stack-1].Starting+
				Memory[GPUCurrent_Stack-1].Size +diffsize;

			/*
			for(ptr=from;ptr>=to;ptr=(void *)((char *)ptr-1)){
				*(tmpptr+j-1) = *((char *)ptr);
				j--;
			}
			*/
			{
				INT8 np = (INT8)((char*)from-(char*)to);
				cudaleftpush<<<1,32>>>(tmpptr,from,np);
			}

		}
		else if( size < Memory[i].Size) {
			nmove = 0;
			for(j=i+1;j<GPUCurrent_Stack;j++){
				nmove += Memory[j].Size;
			}
			GPUmemmove<<<1,32>>>((void *)((char *)Memory[i].Starting+size),
					Memory[i+1].Starting,nmove);
		}
		for(j=i+1;j<GPUCurrent_Stack;j++){
			*(Memory[j].PtrToVariable) = 
				(void **) ((char *)*(Memory[j].PtrToVariable) + diffsize);
			/*
			SMalloc[j] += diffsize;
			*/
			Memory[j].Starting = (void *)((char *) Memory[j].Starting+diffsize);
		}
	}
	GPUFREE = (void *)((char *)GPUFREE + diffsize);
	Memory[i].Size = size;
	return a;
}
void *GPUresizelast(void *p, INT8 size){
	p = (void *)GPURealloc(p,size,Memory[GPUCurrent_Stack-1].PtrToVariable);
	return p;
}
void *GPUReallocLast(void *p,INT8 size){
	INT8 diffsize;
//	size = (size/MemAlign + 1) *MemAlign;
	size = (int)ceil((double)size/(double)MemAlign)*MemAlign;

	diffsize = size - Memory[GPUCurrent_Stack-1].Size;
	GPUFREE = (void *)((char *)GPUFREE + diffsize);
	if((char *)GPULastEnd <=(char *) GPUFREE) 
		return NULL;
	else {
		Memory[GPUCurrent_Stack-1].Size = size; 
		return Memory[GPUCurrent_Stack-1].PtrToVariable;
	}
}
void GPUfreelast(void *p){
	if(GPUCurrent_Stack == 0){
		fprintf(stderr,"No current allocated space\n");
	}
	else {
		/*
		GPGPUUFree(*P2Array[GPUCurrent_Stack-1]);
		*/
		GPUFree(p);
	}
	return;
}
void GPUdumpptr(){
	return;
}
void GPUNumMemStack(){
	GPUCurrent_Stack_org = GPUCurrent_Stack;
	return ;
}
void GPUFreeRightNumMemStack(){
	INT8 i;
	for(i=GPUCurrent_Stack_org;i<GPUCurrent_Stack;i++){
		*(Memory[i].PtrToVariable) = NULL;
	}
	GPUFREE = (void *)((char *) Memory[GPUCurrent_Stack_org-1].Starting
			+Memory[GPUCurrent_Stack_org-1].Size);
	GPUCurrent_Stack = GPUCurrent_Stack_org;
	return;
}
INT8 GPUCurMemStack(){
	return GPUCurrent_Stack;
}
void GPUInitialOldMemStack(INT8 oldstacknum){
	INT8 i;
	for(i=oldstacknum;i<GPUCurrent_Stack;i++){
		*(Memory[i].PtrToVariable) = NULL;
	}
	GPUFREE = (void *)((char *) Memory[oldstacknum-1].Starting+
			Memory[oldstacknum-1].Size);
	GPUCurrent_Stack = oldstacknum;
	return;
}
void GPULastSwitchPointer(void **den){
	Memory[GPUCurrent_Stack-1].PtrToVariable = den;
}
void GPUMemSwitchPointer(void **den, void **den2){
	INT8 i;
	for(i=0;i<GPUCurrent_Stack;i++){
		if(Memory[i].PtrToVariable == den){
			Memory[i].PtrToVariable = den2;
			return;
		}
	}
}
void GPUStackPosition(void *a){
	INT8 i;
	for(i=0;i<GPUCurrent_Stack;i++){
		if((char *)a == Memory[i].Starting){
			printf("Now at stack # %ld in %ld\n",i,GPUCurrent_Stack);fflush(stdout);
		}
	}
}
