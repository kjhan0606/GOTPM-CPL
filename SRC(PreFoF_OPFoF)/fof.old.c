#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<omp.h>
#include "pmheader.h"
#include "mpi.h"
#include "fof.h"
#define IMOD(A,B) ((A) - ((A)/(B))*(B))
#define MIN(A,B) ((A)<(B) ? (A):(B))
#define MAX(A,B) ((A)>(B) ? (A):(B))

void Omp_FoF_Make_Tree(FoFTStruct *TREE_START,FoFTPtlStruct *ptl,size_t np,Box box,
		int nx, long long ntreemax){
	int nthreads;
#pragma omp parallel 
	{
#pragma omp master
		{
			nthreads = omp_get_num_threads();
		}
	}
	FoFTPtlStruct **Grids1D = (FoFTPtlStruct **)malloc(sizeof(FoFTPtlStruct*)*nthreads);
	double nsize = (double ) nx/(double ) nthreads;
	{
		size_t i;
		for(i=0;i<nthreads;i++){
			Grids1D[i] = NULL;
		}
	}
#pragma omp parallel 
	{
		size_t i;
		int idthread = omp_get_thread_num();
		for(i=0;i<np;i++){
			int ix = ptl[i].r[0]/nsize;
			if(ix >= nthreads) ix = nthreads-1;
			else if(ix < 0) ix = 0;
			if(ix == idthread){
				FoFTPtlStruct *tp = Grids1D[ix];
				Grids1D[ix] = ptl + i;
				ptl[i].sibling = tp;
				ptl[i].included = NO;
			}
		}
		long long treeSize = ntreemax/nthreads;

		FoFTStruct *thread_TREE_START = TREE_START + idthread*treeSize;
		((GENERAL_TPtl_POINTER*)thread_TREE_START)->sibling = (thread_TREE_START + treeSize);
		if(idthread == nthreads-1) ((GENERAL_TPtl_POINTER*)thread_TREE_START)->sibling = NULL;

		FoFBeginEndTree FoF_divide_node(FoFTStruct *,FoFTStruct *, FoFTPtlStruct *, Box ,FoFTStruct *); 
		FoFBeginEndTree beginend = FoF_divide_node(thread_TREE_START,thread_TREE_START+1,
				Grids1D[idthread], box,thread_TREE_START);

		if(beginend.start - thread_TREE_START >= treeSize){
			printf("Erorr overflowing the Tree structure memory: %ld : %ld\n", 
					(beginend.start-thread_TREE_START), treeSize);
			exit(99);
		}
	}

	free(Grids1D);
}
void FoF_Make_Tree(FoFTStruct *TREE_START,FoFTPtlStruct *ptl,size_t np,Box box){
	FoFBeginEndTree beginend;
	FoFTStruct *NewTree;
	FoFTPtlStruct *ptr;
	ptr = ptl;
	while(ptr != NULL){
		ptr->included = NO;
		ptr = ptr->sibling;
	}
	TREE_START->sibling = NULL;
	FoFBeginEndTree FoF_divide_node(FoFTStruct *,FoFTStruct *, FoFTPtlStruct *, Box ,FoFTStruct *); 
	beginend = FoF_divide_node(TREE_START,TREE_START+1,ptl,box,TREE_START);
}
FoFBeginEndTree FoF_divide_node(FoFTStruct *TREE_START,FoFTStruct *NewTree, 
		FoFTPtlStruct *ptl, Box box,FoFTStruct *ThisTree){ 
	FoFBeginEndTree beginend;
	FoFTStruct *p2tree,tmpnode[8];
	FoFTStruct *NowCal;
	void *from_sibling;
	FoFTPtlStruct *p2ptl,*tmpptr,*tmpptr2;
	Box tmpbox[8];
	int i,j,k,mnode,mx,my,mz;
	POSTYPE x0,y0,z0,inv_halfw,halfw;
	POSTYPE tmpx,tmpy,tmpz,tmpdist2,distmax;
	POSTYPE ptlmass;
	int count;
	ThisTree->type = TYPE_TREE;
	/*
	ThisTree->sibling = NULL;
	*/
	ThisTree->daughter = NULL;
//ThisTree->L = box.width;
	/*
	ThisTree->r0[0] = box.x;
	ThisTree->r0[1] = box.y;
	ThisTree->r0[2] = box.z;
	*/
	ThisTree->Nparticle = 0;
	ThisTree->mono[0]= ThisTree->mono[1]= ThisTree->mono[2]= 0.;

	p2ptl = ptl;
	while(p2ptl != NULL){
		ThisTree->Nparticle ++;
		ThisTree->mono[0] += p2ptl->r[0];
		ThisTree->mono[1] += p2ptl->r[1];
		ThisTree->mono[2] += p2ptl->r[2];
		p2ptl = p2ptl->sibling;
	}
	ThisTree->mono[0] /= ThisTree->Nparticle;
	ThisTree->mono[1] /= ThisTree->Nparticle;
	ThisTree->mono[2] /= ThisTree->Nparticle;
	distmax = -1.E20;
	p2ptl = ptl;
	while(p2ptl != NULL){
		tmpx = p2ptl->r[0] - ThisTree->mono[0];
		tmpy = p2ptl->r[1] - ThisTree->mono[1];
		tmpz = p2ptl->r[2] - ThisTree->mono[2];
		tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		distmax = MAX(distmax,tmpdist2);
		p2ptl = p2ptl->sibling;
	}
	ThisTree->dist2 = distmax;
	ThisTree->dist = sqrt(distmax);
	x0 = box.x;
	y0 = box.y;
	z0 = box.z;
	halfw = box.width*0.5L;
	inv_halfw = 1.L/halfw;
	/* initialize temporary tree array */
	for(i=0;i<8;i++) {
		tmpnode[i].sibling = tmpnode[i].daughter = NULL;
		tmpnode[i].Nparticle = 0;
		tmpbox[i].width = halfw;
		/*
		tmpbox[i].x = x0+IMOD(i,2)*halfw;
		tmpbox[i].y = y0+(IMOD(i,4)/2)*halfw;
		*/
		tmpbox[i].x = x0+(i%2)*halfw;
		tmpbox[i].y = y0+((i%4)/2)*halfw;
		tmpbox[i].z = z0+(i/4)*halfw;
	}
	p2ptl = ptl;
	while(p2ptl != NULL){
		mx = (int)((p2ptl->r[0] - x0)*inv_halfw);
		my = (int)((p2ptl->r[1] - y0)*inv_halfw);
		mz = (int)((p2ptl->r[2] - z0)*inv_halfw);
		/*
		mx = my = mz = 0;
		if(p2ptl->r[0] >= ThisTree->mono[0]) mx= 1;
		if(p2ptl->r[1] >= ThisTree->mono[1]) my= 1;
		if(p2ptl->r[2] >= ThisTree->mono[2]) mz= 1;
		*/
		mx = MIN(mx,1);
		mx = MAX(mx,0);
		my = MIN(my,1);
		my = MAX(my,0);
		mz = MIN(mz,1);
		mz = MAX(mz,0);

		mnode = mx + 2*my + 4*mz;
		tmpnode[mnode].Nparticle ++; 
		tmpptr = tmpnode[mnode].daughter;
		tmpptr2 = p2ptl->sibling;
		tmpnode[mnode].daughter = p2ptl;
		p2ptl->sibling = tmpptr;
		p2ptl = tmpptr2;
	}
	/* Making link from Mother Node */
    for(i=0;i<8;i++){
        if(tmpnode[i].Nparticle > 0) break;
    }
    if(tmpnode[i].Nparticle >= MIN_CELL_PARTICLE_NUM && ThisTree->dist>0.5*MINCELLWIDTH){
        ThisTree->daughter = (void *)(NewTree);
    }
    else {
        ThisTree->daughter = (void *) tmpnode[i].daughter ;
    }
//  count = 0;
    NowCal = NewTree;
    from_sibling = NULL;
    for(i=0;i<8;i++){
        if(tmpnode[i].Nparticle >= MIN_CELL_PARTICLE_NUM && ThisTree->dist>0.5*MINCELLWIDTH){
            NewTree->daughter = tmpnode[i].daughter;
            if(from_sibling != NULL)
                ((GENERAL_TPtl_POINTER*)from_sibling)->sibling = NewTree;
            from_sibling = NewTree;
            NewTree++;
//          count ++;
        }
        else if(tmpnode[i].Nparticle > 0 ){
            tmpptr = tmpnode[i].daughter;
            if(from_sibling != NULL)
                ((GENERAL_TPtl_POINTER*)from_sibling)->sibling = tmpptr;
            while(tmpptr != NULL){
                from_sibling = tmpptr;
                tmpptr = tmpptr->sibling;
            }
        }
    }
    ((GENERAL_TPtl_POINTER*)from_sibling)->sibling = ThisTree->sibling;
    for(i=0;i<8;i++){
        if(tmpnode[i].Nparticle >= MIN_CELL_PARTICLE_NUM && ThisTree->dist >0.5*MINCELLWIDTH){
            beginend = FoF_divide_node(TREE_START,NewTree,
                    tmpnode[i].daughter, tmpbox[i],NowCal);
            NewTree = beginend.start;
            NowCal ++;
        }
    }
    beginend.start = NewTree;
    return beginend;
}
