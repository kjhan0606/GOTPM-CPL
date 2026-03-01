/* This makes Tree structure and walks along tree structures.
 *
 * */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "tree.h"
#include "force_spline.h"
/*
#define IMOD(A,B) ((A) - ((A)/(B))*(B))
*/
#define EPS 5.E-6
#define IMOD(A,B) (A%B)
#define MIN(A,B) ((A)<(B) ? (A):(B))
#define MAX(A,B) ((A)>(B) ? (A):(B))
/*
static float maxarg1,maxarg2;
static float minarg1,minarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define MIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ? (minarg1) : (minarg2))
*/
extern float diff[NSPLINE][3],slope[NSPLINE][3],ran2nran;
enum boolean fof_open(particle p,FoFTStruct *tree, float fof_link){
	float tmpx,tmpy,tmpz; 
	float dist2,dist,r,diff; 
	float ratio; 
	tmpx = p.x-tree->mono[0]; 
	tmpy = p.y-tree->mono[1]; 
	tmpz = p.z-tree->mono[2]; 
	r = tree->dist;
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz; 
	dist = sqrt(dist2);
	diff = dist - r;
	if(diff <= fof_link) return YES;
	else return NO;
}
enum boolean nopen(particle *p,TStruct *tree, float theta2){
	float tmpx,tmpy,tmpz; 
	float dist2; 
	float ratio; 
	tmpx = p->x-tree->mono[0]; 
	tmpy = p->y-tree->mono[1]; 
	tmpz = p->z-tree->mono[2]; 
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz; 
	ratio = tree->dist2/theta2; 
	if(dist2 < ratio) return YES;
	else return NO;
}
static int ntmp;   
static float tmpx,tmpy,tmpz,dist2;   
static float xx,yy,zz,xy,xz,yz,tmpxx; 
static float qxx1,qxx2,qxx3,qxy,tmpdist;   
static float fplmf1,fplmf2,fplmf3,fplmf4,tmptmp;   
static float fplmf;  
static TStruct *pointer;   
static TPtlStruct *ppointer;   
static float ptlmass;
#define CELLFORCE(p,ptr,force,ran2nran) \
{\
	pointer = (TStruct *) ptr;   \
	tmpx = pointer->mono[0] - p->x;   \
	tmpy = pointer->mono[1] - p->y; \
	tmpz = pointer->mono[2] - p->z;   \
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;   \
	xx = pointer->quad[0];   \
	yy = pointer->quad[1];   \
	zz = pointer->quad[2];   \
	xy = pointer->quad[3];   \
	xz = pointer->quad[4];   \
	yz = pointer->quad[5];   \
	tmpxx = xx + yy + zz;   \
	qxx1 = xx*tmpx+xy*tmpy+xz*tmpz;   \
	qxx2 = xy*tmpx+yy*tmpy+yz*tmpz;   \
	qxx3 = xz*tmpx+yz*tmpy+zz*tmpz;   \
	qxy = xx*tmpx*tmpx+yy*tmpy*tmpy+zz*tmpz*tmpz+   \
		2.*(xy*tmpx*tmpy+xz*tmpx*tmpz+yz*tmpy*tmpz); \
	ntmp = dist2/ran2nran;   \
	tmpdist = (dist2-ntmp*ran2nran);   \
	ntmp = MIN(ntmp,NSPLINE-1);\
	fplmf1 = diff[ntmp][0] + slope[ntmp][0]*tmpdist;   \
	fplmf2 = diff[ntmp][1] + slope[ntmp][1]*tmpdist;   \
	fplmf3 = diff[ntmp][2] + slope[ntmp][2]*tmpdist;   \
	tmptmp = pointer->mass * fplmf1 + tmpxx*fplmf2 +   \
		qxy*fplmf3;   \
	force->x += tmpx*tmptmp+2.*qxx1*fplmf2; \
	force->y += tmpy*tmptmp+2.*qxx2*fplmf2; \
	force->z += tmpz*tmptmp+2.*qxx3*fplmf2; \
}
#define PARTICLEFORCE(p,ptr,force,ran2nran) \
{ \
	ppointer = (TPtlStruct *) ptr; \
	ptlmass = ppointer->mass; \
	tmpx = ppointer->r[0] - p->x;  \
	tmpy = ppointer->r[1] - p->y;  \
	tmpz = ppointer->r[2] - p->z;  \
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;  \
	ntmp = dist2/ran2nran;  \
	ntmp = MIN(ntmp,NSPLINE-1);\
	fplmf = diff[ntmp][0] + slope[ntmp][0] * (dist2-ntmp*ran2nran);  \
	fplmf *= ptlmass; \
	force->x +=  tmpx * fplmf;  \
	force->y +=  tmpy * fplmf;  \
	force->z +=  tmpz * fplmf;  \
}

static float idist2,sqrtdist2,isqrtdist2; 
#define CELLPLUMPOTENT(p,ptr,potent,epsilon2) \
{\
	pointer = (TStruct *) ptr;   \
	tmpx = p->x - pointer->mono[0];   \
	tmpy = p->y - pointer->mono[1];   \
	tmpz = p->z - pointer->mono[2];   \
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;   \
	xx = pointer->quad[0];   \
	yy = pointer->quad[1];   \
	zz = pointer->quad[2];   \
	xy = pointer->quad[3];   \
	xz = pointer->quad[4];   \
	yz = pointer->quad[5];   \
	tmpxx = xx + yy + zz;   \
	qxy = xx*tmpx*tmpx+yy*tmpy*tmpy+zz*tmpz*tmpz+   \
		2.*(xy*tmpx*tmpy+xz*tmpx*tmpz+yz*tmpy*tmpz); \
	ntmp = dist2/ran2nran;   \
	tmpdist = (dist2-ntmp*ran2nran);   \
	ntmp = MIN(ntmp,NSPLINE-1);\
	fplmf1 = diff[ntmp][0] + slope[ntmp][0]*tmpdist;   \
	fplmf2 = diff[ntmp][1] + slope[ntmp][1]*tmpdist;   \
	fplmf3 = diff[ntmp][2] + slope[ntmp][2]*tmpdist;   \
	potent += pointer->mass * fplmf1 + tmpxx*fplmf2 + \
		qxy*fplmf3; \
}

#define PARTICLEPLUMPOTENT(p,ptr,potent,epsilon2) \
{\
	ppointer = (TPtlStruct *) ptr; \
	ptlmass = ppointer->mass; \
	tmpx = p->x -  ppointer->r[0]; \
	tmpy = p->y -  ppointer->r[1]; \
	tmpz = p->z -  ppointer->r[2]; \
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;  \
	ntmp = dist2/ran2nran;  \
	ntmp = MIN(ntmp,NSPLINE-1);\
	fplmf = diff[ntmp][0] + slope[ntmp][0] * (dist2-ntmp*ran2nran);  \
	fplmf *= ptlmass; \
	potent += fplmf; \
}


/* 
 * *p is the position at which you want to calculate force using tree
 * theta2 is the opening angle for tree walk
 * *tree is the tree structure.
 */

void new_fof_link(particle *p,float fof_link,FoFTStruct *tree,
		FoFTPtlStruct *ptl,particle *linked){
	int ncount, now;
	void *ptr,*optr,*nptr;
	float tmpx,tmpy,tmpz;
	float fof_link2,dist2;
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
		pforce *force){
	void *ptr,*optr;
	ptr = (void *)tree;
	force->x = force->y = force->z = 0.;
	while(ptr != NULL){
		switch(((TYPE*)ptr)->type ) {
			case TYPE_TREE:
				switch(nopen(p,ptr,theta2)){
					case YES: 
						ptr = (void *)(((TStruct*)ptr)->daughter);
						break;
					default :
						CELLFORCE(p,ptr,force,ran2nran);
						ptr = (void *)(((TStruct*)ptr)->sibling);
				}
				break;
			default: 
				PARTICLEFORCE(p,ptr,force,ran2nran);
				ptr = (void*)(((TPtlStruct*)ptr)->sibling);
		}
	}
	return;
}
float treeplumpotential(particle *p,float theta2,TStruct *tree,
		TPtlStruct *ptl){
	void *ptr;
	float potent=0.;
	float epsilon2;
	epsilon2 = EPSILON * EPSILON;
	ptr = (void *) tree;
	while(ptr != NULL){
		switch(((TYPE*)ptr)->type ) {
			case TYPE_TREE:
				switch(nopen(p,ptr,theta2)){
					case YES: 
						ptr = (void *)(((TStruct*)ptr)->daughter);
						break;
					default :
						CELLPLUMPOTENT(p,ptr,potent,epsilon2);
						ptr = (void *)(((TStruct*)ptr)->sibling);
				}
				break;
			default: 
				PARTICLEPLUMPOTENT(p,ptr,potent,epsilon2);
				ptr = (void*)(((TPtlStruct*)ptr)->sibling);
		}
	}
	return potent;
}
void FoF_Make_Tree(FoFTStruct *TREE_START,FoFTPtlStruct *ptl,int np,Box box){
	FoFBeginEndTree beginend;
	FoFTStruct *NewTree;
	FoFTPtlStruct *ptr;
	ptr = ptl;
	while(ptr != NULL){
		ptr->included = NO;
		ptr = ptr->sibling;
	}
	TREE_START->sibling = NULL;
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
	float x0,y0,z0,inv_halfw,halfw;
	float tmpx,tmpy,tmpz,tmpdist2,distmax;
	float ptlmass;
	int count;
	ThisTree->type = TYPE_TREE;
	/*
	ThisTree->sibling = NULL;
	*/
	ThisTree->daughter = NULL;
	ThisTree->L = box.width;
	ThisTree->r0[0] = box.x;
	ThisTree->r0[1] = box.y;
	ThisTree->r0[2] = box.z;
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
	halfw = box.width*0.5;
	inv_halfw = 1./halfw;
	/* initialize temporary tree array */
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
	if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
		ThisTree->daughter = (void *)(NewTree);
	}
	else {
		ThisTree->daughter = (void *) tmpnode[i].daughter ;
	}
	count = 0;
	NowCal = NewTree;
	from_sibling = NULL;
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
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
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
			beginend = FoF_divide_node(TREE_START,NewTree,
					tmpnode[i].daughter, tmpbox[i],NowCal);
			NewTree = beginend.start;
			NowCal ++;
		}
	}
	beginend.start = NewTree;
	return beginend;
}

void Make_Tree(TStruct *TREE_START,TPtlStruct *ptl,int np,Box box){
	BeginEndTree beginend;
	int i;
	TStruct *NewTree;
	TREE_START->sibling = NULL;
	for(i=0;i<np;i++) ptl[i].sibling = &ptl[i+1];
	ptl[np-1].sibling = NULL;
	beginend = divide_node(TREE_START,TREE_START+1,ptl,box,TREE_START);
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
	float x0,y0,z0,inv_halfw,halfw;
	float tmpx,tmpy,tmpz,tmpdist2,distmax;
	float ptlmass;
	int count;

	/*   */
	ThisTree->type = TYPE_TREE;
	/*
	ThisTree->sibling = NULL;
	*/
	ThisTree->daughter = NULL;
	ThisTree->L = box.width;
	ThisTree->r0[0] = box.x;
	ThisTree->r0[1] = box.y;
	ThisTree->r0[2] = box.z;
	ThisTree->mass = 0;
	ThisTree->Nparticle = 0;
	ThisTree->mono[0]= ThisTree->mono[1]= ThisTree->mono[2]= 0.;
	ThisTree->quad[0] = ThisTree->quad[1] = ThisTree->quad[2] = 0.;
	ThisTree->quad[3] = ThisTree->quad[4] = ThisTree->quad[5] = 0.;
	ThisTree->mrr = 0.;

	p2ptl = ptl;
	while(p2ptl != NULL){
		ptlmass = p2ptl->mass;
		ThisTree->mass += ptlmass;
		ThisTree->Nparticle ++;
		ThisTree->mono[0] += ptlmass * p2ptl->r[0];
		ThisTree->mono[1] += ptlmass * p2ptl->r[1];
		ThisTree->mono[2] += ptlmass * p2ptl->r[2];
		p2ptl = p2ptl->sibling;
	}
	ThisTree->mono[0] /= ThisTree->mass;
	ThisTree->mono[1] /= ThisTree->mass;
	ThisTree->mono[2] /= ThisTree->mass;
	distmax = -1.E20;
	p2ptl = ptl;
	while(p2ptl != NULL){
		ptlmass = p2ptl->mass;
		tmpx = p2ptl->r[0] - ThisTree->mono[0];
		tmpy = p2ptl->r[1] - ThisTree->mono[1];
		tmpz = p2ptl->r[2] - ThisTree->mono[2];
		tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		distmax = MAX(distmax,tmpdist2);
		ThisTree->mrr += ptlmass * tmpdist2;
		ThisTree->quad[0] += ptlmass * tmpx * tmpx;
		ThisTree->quad[1] += ptlmass * tmpy * tmpy;
		ThisTree->quad[2] += ptlmass * tmpz * tmpz;
		ThisTree->quad[3] += ptlmass * tmpx * tmpy;
		ThisTree->quad[4] += ptlmass * tmpx * tmpz;
		ThisTree->quad[5] += ptlmass * tmpy * tmpz;
		p2ptl = p2ptl->sibling;
	}
	ThisTree->dist2 = distmax;
	p2ptl = ptl;
	x0 = box.x;
	y0 = box.y;
	z0 = box.z;
	halfw = box.width*0.5;
	inv_halfw = 1./halfw;
	/* initialize temporary tree array */
	for(i=0;i<8;i++) {
		tmpnode[i].sibling = tmpnode[i].daughter = NULL;
		tmpnode[i].Nparticle = 0;
		tmpbox[i].width = halfw;
		tmpbox[i].x = x0+(i%2)*halfw;
		tmpbox[i].y = y0+((i%4)/2)*halfw;
		tmpbox[i].z = z0+(i/4)*halfw;
	}

	while(p2ptl != NULL){
		mx = (int)((p2ptl->r[0] - x0)*inv_halfw);
		my = (int)((p2ptl->r[1] - y0)*inv_halfw);
		mz = (int)((p2ptl->r[2] - z0)*inv_halfw);
		mx = MIN(mx,1);
		my = MIN(my,1);
		mz = MIN(mz,1);
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
	if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
		ThisTree->daughter = (void *)(NewTree);
	}
	else {
		ThisTree->daughter = (void *) tmpnode[i].daughter ;
	}
	count = 0;
	NowCal = NewTree;
	from_sibling = NULL;
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
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
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
			beginend = divide_node(TREE_START,NewTree,
					tmpnode[i].daughter, tmpbox[i],NowCal);
			NewTree = beginend.start;
			NowCal ++;
		}
	}
	beginend.start = NewTree;
	return beginend;
}
#define OFFSET 0.001
Box findbox(TPtlStruct *ptl, int nend){
	static float xmin,ymin,zmin,xmax,ymax,zmax;
	static float width;
	static Box box;
	static int i;
	box.x = box.y = box.z = 0.;
	xmax=ymax=zmax = -1.E25;
	xmin=ymin=zmin = 1.E25;
	for(i=0;i<nend;i++){
		/*
		xmin = MIN(xmin,ptl[i].r[0]);
		ymin = MIN(ymin,ptl[i].r[1]);
		zmin = MIN(zmin,ptl[i].r[2]);
		xmax = MAX(xmax,ptl[i].r[0]);
		ymax = MAX(ymax,ptl[i].r[1]);
		zmax = MAX(zmax,ptl[i].r[2]);
		*/
		if(xmin > ptl[i].r[0]) xmin = ptl[i].r[0];
		if(ymin > ptl[i].r[1]) ymin = ptl[i].r[1];
		if(zmin > ptl[i].r[2]) zmin = ptl[i].r[2];
		if(xmax < ptl[i].r[0]) xmax = ptl[i].r[0];
		if(ymax < ptl[i].r[1]) ymax = ptl[i].r[1];
		if(zmax < ptl[i].r[2]) zmax = ptl[i].r[2];
	}

	width = MAX(width,(xmax-xmin));
	width = MAX(width,(ymax-ymin));
	width = MAX(width,(zmax-zmin));
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
