#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<float.h>
#include<unistd.h>
#ifdef _OPENMP
#include<omp.h>
#endif

#include "pmheader.h"
#include "fof.h"


int nullfct0() { return 0;}
int nullfct1() { return 1;}
#ifndef _OPENMP
#define omp_get_thread_num() nullfct0()
#define omp_get_num_threads() nullfct1()
#endif



#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

/*
#define YES 1
#define NO 0
*/


#define SubCellDicision(a,b) ((a)>(b)? 1:0)
#define DivideNode(ThisNode,nparticles) ((ThisNode->dist > 0.5*MINCELLWIDTH ? YES:NO) & (nparticles>=MIN_CELL_PARTICLE_NUM ? YES:NO))

/*
#define EraseFromTree(optr,ptr,nptr) {\
	switch(((TYPE*)optr)->type) {\
		case TYPE_TREE:\
			if(((FoFTStruct*)optr)->daughter == ptr) ((FoFTStruct*)optr)->daughter == nptr;\
			else ((FoFTStruct*)optr)->sibling = nptr;\
			break;\
		default :\
			((FoFTPtlStruct*)ptr)->sibling = nptr;\
	}\
}\
*/




FoFTStruct *Omp2_FoF_divide_node(FoFTStruct *ThisNode, FoFTStruct *SpareNode, Box box,
		int recursiveflag)
{
	FoFTStruct *p2tree, tmpnode[8];
	Box tmpbox[8];
	void *LeftSibling;
	FoFTPtlStruct *p2ptl, *tmpptr, *tmpptr2, *nodeparticles;
	size_t i, j, k, mnode, mx,my,mz;
	POSTYPE tmpx,tmpy,tmpz,tmpdist2, distmax;
	float ptlmass;

	nodeparticles = ThisNode->daughter;
	ThisNode->type = TYPE_TREE;
	ThisNode->Nparticle = 0;
	ThisNode->mono[0] = ThisNode->mono[1] = ThisNode->mono[2] = 0;

	for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
		ThisNode->Nparticle += 1;
		ThisNode->mono[0] += p2ptl->r[0];
		ThisNode->mono[1] += p2ptl->r[1];
		ThisNode->mono[2] += p2ptl->r[2];
	}
	ThisNode->mono[0] /= ThisNode->Nparticle;
	ThisNode->mono[1] /= ThisNode->Nparticle;
	ThisNode->mono[2] /= ThisNode->Nparticle;
	distmax = -1.E20;
	/*
	float q0,q1,q2,q3,q4,q5;
	q0=q1=q2=q3=q4=q5=0;
	*/
	for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
		tmpx = p2ptl->r[0] - ThisNode->mono[0];
		tmpy = p2ptl->r[1] - ThisNode->mono[1];
		tmpz = p2ptl->r[2] - ThisNode->mono[2];
		tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		distmax = MAX(distmax, tmpdist2);
	}
	ThisNode->dist = sqrt(distmax);
	ThisNode->dist2 = (distmax);

	for(i=0;i<8;i++) {
		tmpnode[i].sibling = tmpnode[i].daughter = NULL;
		tmpnode[i].Nparticle = 0;
	}
	POSTYPE hw = box.width*0.5L;
	POSTYPE invhw = 1.L/hw;
	POSTYPE x0,y0,z0;
	x0 = box.x;
	y0 = box.y;
	z0 = box.z;
	for(i=0;i<8;i++) {
		tmpbox[i].x = box.x +(i%2)*hw;
		tmpbox[i].y = box.y +((i%4)/2)*hw;
		tmpbox[i].z = box.z +(i/4)*hw;
	}

	for(p2ptl=nodeparticles;p2ptl;) {

		/*
		mx = SubCellDicision(p2ptl->r[0],ThisNode->mono[0]);
		my = SubCellDicision(p2ptl->r[1],ThisNode->mono[1]);
		mz = SubCellDicision(p2ptl->r[2],ThisNode->mono[2]);
		*/
		mx = (int)(p2ptl->r[0]-x0)*invhw;
		my = (int)(p2ptl->r[1]-y0)*invhw;
		mz = (int)(p2ptl->r[2]-z0)*invhw;
		mx = MIN(mx,1);
		my = MIN(my,1);
		mz = MIN(mz,1);
		mx = MAX(mx,0);
		my = MAX(my,0);
		mz = MAX(mz,0);

		mnode = mx + 2*my + 4*mz;
		tmpnode[mnode].Nparticle ++;
		tmpptr = tmpnode[mnode].daughter;
		tmpptr2 = p2ptl->sibling;
		tmpnode[mnode].daughter = p2ptl;
		p2ptl->sibling = tmpptr;
		p2ptl = tmpptr2;
	}
	for(i=0;i<8;i++) if(tmpnode[i].Nparticle >0) break;

	FoFTStruct *FirstDaughter,*NowDaughter;
	FirstDaughter = NowDaughter =  SpareNode;
	/* Link to the first daughter or first particle */
	if( DivideNode(ThisNode, tmpnode[i].Nparticle) )
		ThisNode->daughter = (void*)FirstDaughter;
	else 
		ThisNode->daughter = (void*) tmpnode[i].daughter;

	LeftSibling = NULL;

	for(i=0;i<8;i++) {
		if( DivideNode(ThisNode, tmpnode[i].Nparticle) ){
			NowDaughter->daughter = tmpnode[i].daughter;
			NowDaughter->Nparticle = tmpnode[i].Nparticle;
			if(LeftSibling) ((GENERAL_TPtl_POINTER*) LeftSibling)->sibling = NowDaughter;
			LeftSibling = NowDaughter;
			NowDaughter++;
		}
		else if(tmpnode[i].Nparticle >0) {
			tmpptr = tmpnode[i].daughter;
			if(LeftSibling) ((GENERAL_TPtl_POINTER*)LeftSibling)->sibling = tmpptr;
			for(;tmpptr;tmpptr=tmpptr->sibling) LeftSibling = tmpptr;
		}
	}
	((GENERAL_TPtl_POINTER *)LeftSibling)->sibling = ThisNode->sibling;
	SpareNode = NowDaughter;


	if(recursiveflag == RECURSIVE){
		FoFTStruct *NextJobNode = FirstDaughter;
		for(i=0;i<8;i++){
			if(DivideNode(ThisNode, tmpnode[i].Nparticle) ){
				SpareNode = Omp2_FoF_divide_node(NextJobNode,SpareNode, tmpbox[i], recursiveflag);
				NextJobNode ++;
			}
		}
	}
	return SpareNode;
}

void Omp2_FoF_Make_Tree(FoFTStruct *TREE_START, size_t nnode, FoFTPtlStruct *ptl, size_t np, 
		Box box,
		int recursiveflag){
	size_t i;
	FoFTStruct *SpareNode = TREE_START+1;
	TREE_START->sibling = NULL;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<np;i++) {
		ptl[i].sibling = ptl + (i +1);
		ptl[i].type = TYPE_PTL;
		ptl[i].included = NO;
	}
	ptl[np-1].sibling = NULL;
	TREE_START->daughter = &(ptl[0]);
	TREE_START->Nparticle = np;




	if(recursiveflag == RECURSIVE) SpareNode = Omp2_FoF_divide_node(TREE_START, SpareNode, box,
			recursiveflag);
	else if(recursiveflag == SERIALIZED){
		FoFTStruct *work;
		for(work=TREE_START;SpareNode-work >0; work++){
			SpareNode = Omp2_FoF_divide_node(work, SpareNode, box, SERIALIZED);
		}
	}
	else if(recursiveflag == PTHREAD) 
	{

		int npid=1;

#ifdef _OPENMP
#pragma omp parallel
		{
#pragma omp master
			{
				npid= omp_get_num_threads();
			}
		}
#endif
		FoFTStruct *work = TREE_START;
		do{
			SpareNode = Omp2_FoF_divide_node(work, SpareNode, box, PTHREAD);
			work ++;
		}
		while( (SpareNode-work) <= (npid*3));
		size_t twork = (SpareNode-work);
		printf("total twork for omp tree building: %ld\n", twork);

#ifdef _OPENMP
#pragma omp parallel 
#endif
		{
			int pid = omp_get_thread_num();
			size_t mys, myf;
			size_t worksize = (twork + npid - 1)/npid;
			size_t j;
			mys = worksize *pid;
			myf = MIN(worksize *(pid+1), twork);
			FoFTStruct *threadSpareNode = SpareNode + pid * ( (nnode - (SpareNode-TREE_START))/ npid);
			FoFTStruct *threadNodeStart = threadSpareNode;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
			for(j=0;j<twork;j++) {
//			for(j=mys;j<myf;j++) {
			//	printf("T%d takes the job of j= %ld/%ld\n", pid, j, twork);
				threadSpareNode = Omp2_FoF_divide_node(work+j,threadSpareNode,box, RECURSIVE);
			//	printf("T%d finishes the job of j= %ld/%ld\n", pid, j, twork);
			}
			printf("p%d has made tree nodes : %ld / %ld\n", 
					pid, (threadSpareNode-threadNodeStart), (nnode - (SpareNode-TREE_START))/ npid);
		}
	}
	printf("Completed the tree buidling\n");
}


