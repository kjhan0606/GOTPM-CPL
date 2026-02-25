/* This makes Tree structure and walks along tree structures.
 * For performance issue it is needed to check which one is better,
 * 2*a or a+a.
 * */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include "pmheader.h"
#include "kjhtree.h"
#include "force_spline.h"
#include "common_simparam.h"
double pmcf(double);
#define MIN(A,B) ((A)<(B) ? (A):(B))
#define MAX(A,B) ((A)>(B) ? (A):(B))
enum boolean fof_open(particle p,FoFTStruct *tree, float fof_link){
	TREEPRECISION tmpx,tmpy,tmpz; 
	TREEPRECISION dist2,dist,r,diff; 
	TREEPRECISION ratio; 
	tmpx = p.x-tree->mono[0]; 
	tmpy = p.y-tree->mono[1]; 
	tmpz = p.z-tree->mono[2]; 
	r = tree->dist;
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz; 
	dist = SQRT(dist2);
	diff = dist - r;
	if(diff <= fof_link) return YES;
	else return NO;
}
enum boolean cellopen(particle *p,TStruct *tree){
	TREEPRECISION tmpx,tmpy,tmpz; 
	TREEPRECISION dist2; 
	TREEPRECISION ratio; 
	tmpx = p->x-tree->mono[0]; 
	tmpy = p->y-tree->mono[1]; 
	tmpz = p->z-tree->mono[2]; 
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz; 
	if(dist2 < tree->dist_over_thetasq) return YES;
	else return NO;
}
#define CELLFORCE(p,ptr,rspheresq) \
{\
	pointer = (TStruct *) ptr;   \
	tmpx = (p->x)-(pointer->mono[0]);\
	tmpy = (p->y)-(pointer->mono[1]);\
	tmpz = (p->z)-(pointer->mono[2]);\
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;   \
	if(dist2 <= rspheresq){\
		dist=SQRT(dist2);\
		ntmp = dist*invran2nran;   \
		Qxx = pointer->quad[0];\
		Qyy = pointer->quad[1];\
		Qzz = pointer->quad[2];\
		Qxy = pointer->quad[3];\
		Qxz = pointer->quad[4];\
		Qyz = pointer->quad[5];\
		/*\
		trQ = (Qxx + Qyy + Qzz);   \
		*/\
		trQ = pointer->trQ;\
		qxx1 = Qxx*tmpx+Qxy*tmpy+Qxz*tmpz;   \
		qxx2 = Qxy*tmpx+Qyy*tmpy+Qyz*tmpz;   \
		qxx3 = Qxz*tmpx+Qyz*tmpy+Qzz*tmpz;   \
		qxy = Qxx*tmpx*tmpx+Qyy*tmpy*tmpy+Qzz*tmpz*tmpz+   \
			2*(Qxy*tmpx*tmpy+Qxz*tmpx*tmpz+Qyz*tmpy*tmpz); \
		tmpdist = (dist-ntmp*ran2nran);   \
		ntmp = MIN(ntmp,NSPLINE-1);\
		fplmf1 = forcecorrectdiff(ntmp,0) + forcecorrectslope(ntmp,0)*tmpdist;   \
		fplmf2 = forcecorrectdiff(ntmp,1) + forcecorrectslope(ntmp,1)*tmpdist;   \
		fplmf3 = forcecorrectdiff(ntmp,2) + forcecorrectslope(ntmp,2)*tmpdist;   \
		tmptmp = pointer->mass*fplmf1 + trQ*fplmf2 + qxy*fplmf3;\
		twofplmf2 = fplmf2+fplmf2;\
		float pcorr = getpcorr(dist);\
		p->ax += tmpx*tmptmp+qxx1*twofplmf2;\
		p->ay += tmpy*tmptmp+qxx2*twofplmf2;\
		p->az += tmpz*tmptmp+qxx3*twofplmf2;\
	}\
}

// 2024_06_21 minseong position pcorr deleted 
#define PARTICLEFORCE(p,ptr,rspheresq) \
{ \
	ppointer = (TPtlStruct *) ptr; \
	ptlmass = ppointer->mass; \
	tmpx = p->x-ppointer->r[0];\
	tmpy = p->y-ppointer->r[1];\
	tmpz = p->z-ppointer->r[2];\
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;  \
	if(dist2 <= rspheresq){\
		dist = SQRT(dist2);\
		/* \
		ntmp = dist/ran2nran;  \
		*/\
		ntmp = dist*invran2nran;  \
		fplmf = forcecorrectdiff(ntmp,0) + forcecorrectslope(ntmp,0) * (dist-ntmp*ran2nran);  \
		fplmf *= ptlmass; \
		float pcorr = getpcorr(dist);\
		p->ax +=  tmpx * fplmf;\
		p->ay +=  tmpy * fplmf;\
		p->az +=  tmpz * fplmf;\
	}\
}


int iorder[8] = { 0, 7, 1, 6, 2, 5, 3, 4};


void new_fof_link(particle *p,float fof_link,FoFTStruct *tree,
		FoFTPtlStruct *ptl,particle *linked){
	int ncount, now;
	void *ptr,*optr,*nptr;
	TREEPRECISION tmpx,tmpy,tmpz;
	TREEPRECISION fof_link2,dist2;
	particle point;
	fof_link2 = fof_link*fof_link;
	ncount = now = 0;
	point.x = p->x;
	point.y = p->y;
	point.z = p->z;
	do {
		optr = (void *) tree;
		ptr = (void*) tree;
		while(ptr != NULL){
			switch(((TYPE*)ptr)->type){
				case TYPE_TREE:
					if(((FoFTStruct *)ptr)->sibling ==
							((FoFTStruct *)ptr)->daughter){
						EraseFromTree(optr,ptr,((FoFTStruct *)ptr)->sibling);
						ptr = ((FoFTStruct *)ptr)->sibling;
					}
					else
					switch(fof_open(point,ptr,fof_link)){
						case YES:
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->daughter);
							break;
						default :
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->sibling);
					}
					break;
				default :
					if(((FoFTPtlStruct*)ptr)->included == YES){
						nptr = ((FoFTPtlStruct *)ptr)->sibling;
						EraseFromTree(optr,ptr,nptr);
						ptr = nptr;
					}
					else {
						tmpx = point.x - ((FoFTPtlStruct*)ptr)->r[0];
						tmpy = point.y - ((FoFTPtlStruct*)ptr)->r[1];
						tmpz = point.z - ((FoFTPtlStruct*)ptr)->r[2];
						dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
						if(dist2 <= fof_link2){
							linked[ncount].x = ((FoFTPtlStruct*)ptr)->r[0];
							linked[ncount].y = ((FoFTPtlStruct*)ptr)->r[1];
							linked[ncount].z = ((FoFTPtlStruct*)ptr)->r[2];
							((FoFTPtlStruct*)ptr)->included = YES;
							ncount ++;
							nptr = ((FoFTPtlStruct *)ptr)->sibling;
							EraseFromTree(optr,ptr,nptr);
						}
						else optr = ptr;
						ptr = (void*)(((FoFTPtlStruct*)ptr)->sibling);
					}
			}
		}
		point = linked[now];
		now ++;
	} while( now < ncount);
	return;
}
void treeforce(particle *p,float theta2,TStruct *tree, TPtlStruct *ptl,
		float rspheresq){
	void *ptr,*optr;
	int ntmp;   
	TREEPRECISION tmpx,tmpy,tmpz,dist2,dist3,dist5,dist;   
	TREEPRECISION Qxx,Qyy,Qzz,Qxy,Qxz,Qyz,trQ; 
	TREEPRECISION qxx1,qxx2,qxx3,qxy,tmpdist;   
	TREEPRECISION fplmf1,fplmf2,twofplmf2,fplmf3,fplmf4,tmptmp;   
	TREEPRECISION fplmf;  
	TStruct *pointer;
	TPtlStruct *ppointer;   
	TREEPRECISION ptlmass;
	ptr = (void *)tree;
	p->ax=p->ay=p->az = 0.;
	while(ptr != NULL){
		switch(((TYPE*)ptr)->type ) {
			case TYPE_TREE:
				switch(cellopen(p,ptr)){
					case YES: 
						ptr = (void *)(((TStruct*)ptr)->daughter);
						break;
					default :
						CELLFORCE(p,ptr,rspheresq);
						ptr = (void *)(((TStruct*)ptr)->sibling);
				}
				break;
			default: 
				PARTICLEFORCE(p,ptr,rspheresq);
				ptr = (void*)(((TPtlStruct*)ptr)->sibling);
		}
	}
	return;
}
void FoF_Make_Tree(FoFTStruct *TREE_START,FoFTPtlStruct *ptl,long np,
		FoFBox box){
	FoFBeginEndTree beginend;
	FoFTStruct *NewTree;
	FoFTPtlStruct *ptr;
	TREE_START->sibling = NULL;
	int i;
	for(i=0;i<np;i++) {
		ptl[i].sibling = &ptl[i+1];
		ptl[i].included = NO;
	}
	ptl[np-1].sibling = NULL;
	beginend = FoF_divide_node(TREE_START,TREE_START+1,ptl,box,TREE_START);
}
FoFBeginEndTree FoF_divide_node(FoFTStruct *TREE_START,FoFTStruct *NewTree, 
		FoFTPtlStruct *ptl, FoFBox box,FoFTStruct *ThisTree){ 
	FoFBeginEndTree beginend;
	FoFTStruct *p2tree,tmpnode[8];
	FoFTStruct *NowCal;
	void *from_sibling;
	FoFTPtlStruct *p2ptl,*tmpptr,*tmpptr2;
	FoFBox tmpbox[8];
	int i,j,k,mnode,mx,my,mz;
	TREEPRECISION x0,y0,z0,inv_halfw,halfw;
	TREEPRECISION tmpx,tmpy,tmpz,tmpdist2,distmax;
	TREEPRECISION ptlmass;
	int count;
	ThisTree->type = TYPE_TREE;
	ThisTree->daughter = NULL;
	ThisTree->L = box.maxwidth;
	ThisTree->Nparticle = 0;
	{
		double monx,mony,monz;
		monx = mony = monz = 0.L;
		p2ptl = ptl;
		while(p2ptl != NULL){
			ThisTree->Nparticle ++;
			monx += p2ptl->r[0];
			mony += p2ptl->r[1];
			monz += p2ptl->r[2];

			p2ptl = p2ptl->sibling;
		}
		ThisTree->mono[0] = monx/ ThisTree->Nparticle;
		ThisTree->mono[1] = mony/ ThisTree->Nparticle;
		ThisTree->mono[2] = monz/ ThisTree->Nparticle;
	}
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
	ThisTree->dist = SQRT(distmax);
	x0 = box.xi;
	y0 = box.yi;
	z0 = box.zi;
	halfw = box.maxwidth*0.5;
	inv_halfw = 1./halfw;
	for(i=0;i<8;i++) {
		tmpnode[i].sibling = tmpnode[i].daughter = NULL;
		tmpnode[i].Nparticle = 0;
		tmpbox[i].maxwidth = halfw;
		tmpbox[i].xi = x0+(i%2)*halfw;
		tmpbox[i].yi = y0+((i%4)/2)*halfw;
		tmpbox[i].zi = z0+(i/4)*halfw;
	}
	p2ptl = ptl;
	while(p2ptl != NULL){
		mx = (int)((p2ptl->r[0] - x0)*inv_halfw);
		my = (int)((p2ptl->r[1] - y0)*inv_halfw);
		mz = (int)((p2ptl->r[2] - z0)*inv_halfw);
		mx = MIN(mx,1); my = MIN(my,1); mz = MIN(mz,1);
		mnode = mx + 2*my + 4*mz;
		tmpnode[mnode].Nparticle ++; 
		tmpptr = tmpnode[mnode].daughter;
		tmpptr2 = p2ptl->sibling;
		tmpnode[mnode].daughter = p2ptl;
		p2ptl->sibling = tmpptr;
		p2ptl = tmpptr2;
	}
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle > 0) break;
	}
	if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE 
			&& tmpbox[i].maxwidth >= MINWIDTH){
		ThisTree->daughter = (void *)(NewTree);
	}
	else {
		ThisTree->daughter = (void *) tmpnode[i].daughter ;
	}
	count = 0;
	NowCal = NewTree;
	from_sibling = NULL;
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE 
				&& tmpbox[i].maxwidth >= MINWIDTH){
			NewTree->daughter = tmpnode[i].daughter;
			if(from_sibling != NULL) 
				((GENERAL_TPtl_POINTER*)from_sibling)->sibling = NewTree;
			from_sibling = NewTree;
			NewTree++;
			count ++;
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
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE
				&& tmpbox[i].maxwidth >= MINWIDTH){
			beginend = FoF_divide_node(TREE_START,NewTree,
					tmpnode[i].daughter, tmpbox[i],NowCal);
			NewTree = beginend.start;
			NowCal ++;
		}
	}
	beginend.start = NewTree;
	return beginend;
}

/*
extern float theta;
*/
TREEPRECISION thetasq;
int Make_Tree(TStruct *TREE_START,TPtlStruct *ptl,long np,Box box, float theta){
	BeginEndTree beginend;
	int ncell;
	long i;
	TStruct *NewTree;
	thetasq = theta*theta;
	TREE_START->sibling = NULL;
	for(i=0;i<np;i++) ptl[i].sibling = &ptl[i+1];
	ptl[np-1].sibling = NULL;
	beginend = divide_node(TREE_START,TREE_START+1,ptl,box,TREE_START);
 	ncell = (int)(beginend.start - TREE_START);
	return ncell;
}
BeginEndTree divide_node(TStruct *TREE_START,TStruct *NewTree, 
		TPtlStruct *ptl, Box box,TStruct *ThisTree){ 
	BeginEndTree beginend;
	TStruct *p2tree, tmpnode[8];
	TStruct *NowCal;
	void *from_sibling;
	TPtlStruct *p2ptl,*tmpptr,*tmpptr2;
	Box tmpbox[8];
	int i,j,k,mnode,mx,my,mz;
	TREEPRECISION x0,y0,z0,inv_halfw,halfw;
	TREEPRECISION tmpx,tmpy,tmpz,tmpdist2,distmax;
	TREEPRECISION ptlmass;
	long ncount;

	ThisTree->type = TYPE_TREE;
	ThisTree->daughter = NULL;
	ThisTree->L = box.width;
	/*
	ThisTree->r0[0] = box.x;
	ThisTree->r0[1] = box.y;
	ThisTree->r0[2] = box.z;
	*/
	ThisTree->mass = 0;
	ThisTree->Nparticle = 0;
	/*
	ThisTree->mono[0]= ThisTree->mono[1]= ThisTree->mono[2]= 0.;
	ThisTree->quad[0] = ThisTree->quad[1] = ThisTree->quad[2] = 0.;
	ThisTree->quad[3] = ThisTree->quad[4] = ThisTree->quad[5] = 0.;
	ThisTree->mrr = 0.;
	*/

	{
		double monx,mony,monz;
		monx = mony = monz = 0L;
		p2ptl = ptl;
		while(p2ptl != NULL){
			ptlmass = p2ptl->mass;
			ThisTree->mass += ptlmass;
			ThisTree->Nparticle ++;
			/*
			ThisTree->mono[0] += ptlmass * p2ptl->r[0];
			ThisTree->mono[1] += ptlmass * p2ptl->r[1];
			ThisTree->mono[2] += ptlmass * p2ptl->r[2];
			*/
			monx += ptlmass* p2ptl->r[0];
			mony += ptlmass* p2ptl->r[1];
			monz += ptlmass* p2ptl->r[2];
			p2ptl = p2ptl->sibling;
		}
		/*
		ThisTree->mono[0] /= ThisTree->mass;
		ThisTree->mono[1] /= ThisTree->mass;
		ThisTree->mono[2] /= ThisTree->mass;
		*/
		ThisTree->mono[0] = monx/ThisTree->mass;
		ThisTree->mono[1] = mony/ThisTree->mass;
		ThisTree->mono[2] = monz/ThisTree->mass;
	}
	{
		/*
		float mrr=0;
		*/
		double mrr,quad0,quad1,quad2,quad3,quad4,quad5;
		mrr = quad0=quad1=quad2=quad3=quad4=quad5 = 0;
		distmax = -1.E20;
		p2ptl = ptl;
		while(p2ptl != NULL){
			ptlmass = p2ptl->mass;
			tmpx = p2ptl->r[0] - ThisTree->mono[0];
			tmpy = p2ptl->r[1] - ThisTree->mono[1];
			tmpz = p2ptl->r[2] - ThisTree->mono[2];
			tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
			distmax = MAX(distmax,tmpdist2);
			/*
			ThisTree->mrr += ptlmass * tmpdist2;
			ThisTree->quad[0] += ptlmass * tmpx * tmpx;
			ThisTree->quad[1] += ptlmass * tmpy * tmpy;
			ThisTree->quad[2] += ptlmass * tmpz * tmpz;
			ThisTree->quad[3] += ptlmass * tmpx * tmpy;
			ThisTree->quad[4] += ptlmass * tmpx * tmpz;
			ThisTree->quad[5] += ptlmass * tmpy * tmpz;
			*/
			mrr += ptlmass * tmpdist2;
			quad0 += ptlmass * tmpx * tmpx;
			quad1 += ptlmass * tmpy * tmpy;
			quad2 += ptlmass * tmpz * tmpz;
			quad3 += ptlmass * tmpx * tmpy;
			quad4 += ptlmass * tmpx * tmpz;
			quad5 += ptlmass * tmpy * tmpz;

			p2ptl = p2ptl->sibling;
		}
		ThisTree->dist2 = distmax;
		ThisTree->dist_over_thetasq = (distmax)/thetasq;
		ThisTree->mrr = mrr;
		ThisTree->quad[0] = quad0;
		ThisTree->quad[1] = quad1;
		ThisTree->quad[2] = quad2;
		ThisTree->quad[3] = quad3;
		ThisTree->quad[4] = quad4;
		ThisTree->quad[5] = quad5;
		ThisTree->trQ = quad0+quad1+quad2;
		/*
		ThisTree->trQ = ThisTree->quad[0]+ThisTree->quad[1]+ThisTree->quad[2];
		*/

	}
	x0 = box.x;
	y0 = box.y;
	z0 = box.z;
	halfw = box.width*0.5;
	inv_halfw = 1./halfw;
	for(i=0;i<8;i++) {
		tmpnode[i].sibling = tmpnode[i].daughter = NULL;
		tmpnode[i].Nparticle = 0;
		tmpbox[i].width = halfw;
		tmpbox[i].x = x0+(i%2)*halfw;
		tmpbox[i].y = y0+((i%4)/2)*halfw;
		tmpbox[i].z = z0+(i/4)*halfw;
	}

	p2ptl = ptl;
	while(p2ptl != NULL){
		mx = (int)((p2ptl->r[0] - x0)*inv_halfw);
		my = (int)((p2ptl->r[1] - y0)*inv_halfw);
		mz = (int)((p2ptl->r[2] - z0)*inv_halfw);
		mx = MIN(mx,1); my = MIN(my,1); mz = MIN(mz,1);
		mnode = mx + 2*my + 4*mz;
		tmpnode[mnode].Nparticle ++; 
		tmpptr = tmpnode[mnode].daughter;
		tmpptr2 = p2ptl->sibling;
		tmpnode[mnode].daughter = p2ptl;
		p2ptl->sibling = tmpptr;
		p2ptl = tmpptr2;
	}
	for(j=0;j<8;j++){ i = iorder[j];
		if(tmpnode[i].Nparticle > 0) break;
	}
	if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE
			&& tmpbox[i].width >= MINWIDTH){
		ThisTree->daughter = (void *)(NewTree);
	}
	else {
		ThisTree->daughter = (void *) tmpnode[i].daughter ;
	}
	ncount = 0;
	NowCal = NewTree;
	from_sibling = NULL;
	for(j=0;j<8;j++){ i = iorder[j];
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE
				&& tmpbox[i].width >= MINWIDTH){
			NewTree->daughter = tmpnode[i].daughter;
			if(from_sibling != NULL) 
				((GENERAL_TPtl_POINTER*)from_sibling)->sibling = NewTree;
			from_sibling = NewTree;
			NewTree++;
			ncount ++;
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
	for(j=0;j<8;j++){ i = iorder[j];
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE
				&& tmpbox[i].width >= MINWIDTH){
			beginend = divide_node(TREE_START,NewTree,
					tmpnode[i].daughter, tmpbox[i],NowCal);
			NewTree = beginend.start;
			NowCal ++;
		}
	}
	beginend.start = NewTree;
	return beginend;
}
/*
#define OFFSET 0.0001
Box findbox(TPtlStruct *ptl, size_t nend){
	float xmin,ymin,zmin,xmax,ymax,zmax;
	float width;
	Box box;
	size_t i;
	box.x = box.y = box.z = 0.;
	xmax=ymax=zmax = -1.E25;
	xmin=ymin=zmin = 1.E25;
	for(i=0;i<nend;i++){
		xmin = MIN(xmin,ptl[i].r[0]);
		ymin = MIN(ymin,ptl[i].r[1]);
		zmin = MIN(zmin,ptl[i].r[2]);
		xmax = MAX(xmax,ptl[i].r[0]);
		ymax = MAX(ymax,ptl[i].r[1]);
		zmax = MAX(zmax,ptl[i].r[2]);
	}
	width = MAX(width,xmax-xmin);
	width = MAX(width,ymax-ymin);
	width = MAX(width,zmax-zmin);
	box.x = xmin;
	box.y = ymin;
	box.z = zmin;
	box.width = width;
	box.x = box.x -box.width*OFFSET;
	box.y = box.y -box.width*OFFSET;
	box.z = box.z -box.width*OFFSET;
	box.width = box.width + box.width*OFFSET*2.;
	return box;
}
#undef OFFSET
*/
enum boolean pfof_open(particle p,FoFTStruct *tree, float fof_link){
    TREEPRECISION tmpx,tmpy,tmpz;
    TREEPRECISION dist2,dist,r,diff;
    TREEPRECISION ratio;
    tmpx = fabs(p.x-tree->mono[0]);
    tmpy = fabs(p.y-tree->mono[1]);
    tmpz = fabs(p.z-tree->mono[2]);
    r = tree->dist;
    dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
    dist = sqrt(dist2);
    diff = dist - r;
    if(diff <= fof_link) return YES;
    else return NO;
}
int pnew_fof_link(particle *p,float fof_link,FoFTStruct *tree,
		FoFTPtlStruct *ptl,particle *linked,int nhalo){
	int ncount, now;
	void *ptr,*optr,*nptr;
	TREEPRECISION tmpx,tmpy,tmpz;
	TREEPRECISION fof_link2,dist2;
	TREEPRECISION Lpx,Lpy,Lpz;
	particle point;

	fof_link2 = fof_link*fof_link;
	ncount = now = 0;
	point.x = p->x;
	point.y = p->y;
	point.z = p->z;
	do {
		optr = ptr = (void*) tree;
		while(ptr != NULL){
			switch(((TYPE*)ptr)->type){
				case TYPE_TREE:
					if(((FoFTStruct *)ptr)->sibling ==
							((FoFTStruct *)ptr)->daughter){
						EraseFromTree(optr,ptr,((FoFTStruct *)ptr)->sibling);
						ptr = ((FoFTStruct *)ptr)->sibling;
					}
					else
					switch(pfof_open(point,ptr,fof_link)){
						case YES:
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->daughter);
							break;
						default :
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->sibling);
					}
					break;
				default :
					if(((FoFTPtlStruct*)ptr)->included == YES){
						nptr = ((FoFTPtlStruct *)ptr)->sibling;
						EraseFromTree(optr,ptr,nptr);
						ptr = nptr;
					}
					else {
						tmpx = fabs(point.x - ((FoFTPtlStruct*)ptr)->r[0]);
						tmpy = fabs(point.y - ((FoFTPtlStruct*)ptr)->r[1]);
						tmpz = fabs(point.z - ((FoFTPtlStruct*)ptr)->r[2]);

						dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
						dist2 = sqrt(dist2);
						if(dist2 <= fof_link){
							linked[ncount].x = ((FoFTPtlStruct*)ptr)->r[0];
							linked[ncount].y = ((FoFTPtlStruct*)ptr)->r[1];
							linked[ncount].z = ((FoFTPtlStruct*)ptr)->r[2];
							((FoFTPtlStruct*)ptr)->haloindx = nhalo;
							((FoFTPtlStruct*)ptr)->included = YES;
							ncount ++;
							nptr = ((FoFTPtlStruct *)ptr)->sibling;
							EraseFromTree(optr,ptr,nptr);
						}
						else optr = ptr;
						ptr = (void*)(((FoFTPtlStruct*)ptr)->sibling);
					}
			}
		}
		point = linked[now];
		now ++;
	} while( now <= ncount);
	return (ncount);
}
