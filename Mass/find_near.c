#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include "near.h"
#define MAX(A,B) ((A)>(B) ? (A):(B))
#define MIN(A,B) ((A)<(B) ? (A):(B))
#define MAX_NUM_NEAR 1000
typedef struct nearestneighbor{
	float dist2;
	int index;
} Neighbor;
Neighbor neighbor[MAX_NUM_NEAR];
int NUM_NEIGHBOR;
enum boolean NEAR_OPEN(particle p,TStruct *tree,int now,float maxdist){
	float tmpx,tmpy,tmpz;
	float dist2;
	float r2;
	float dist,r,sortdist;
	if(now < NUM_NEIGHBOR) return YES;
	tmpx = p.x - tree->mono[0];
	tmpy = p.y - tree->mono[1];
	tmpz = p.z - tree->mono[2];
	r2 = tree->dist2;
	r = sqrt(r2);
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
	dist = sqrt(dist2);
	sortdist = (maxdist);
	if(dist-r > sortdist) return NO;
	else return YES;
}
#define INSERT(dist2,maxdist,maxdist2) \
{\
	for(i=0;i<now;i++) if(neighbor[i].dist2>dist2) break;\
	for(j=MIN(now,NUM_NEIGHBOR-1)-1;j>=i;j--) neighbor[j+1] = neighbor[j];\
	neighbor[i].dist2 = dist2;\
	neighbor[i].index = (TPtlStruct*)ptr - ptl;\
	now ++;\
	now=MIN(NUM_NEIGHBOR,now);\
	maxdist = sqrt(maxdist2=neighbor[now-1].dist2);\
}
/* Num_neighbor --> 찾고자하는 주변 입자의 수
 * result       --> Num_neighbor번째 입자까지의 거리
 */
int Find_Near(particle *p,int Num_neighbor,TStruct *tree,
		TPtlStruct *ptl, float *result, int *nearindex,float *DIST2){
	int i,j,k;
	float dist2,tmpx,tmpy,tmpz;
	float maxdist,maxdist2;
	void *ptr;
	int now;
	particle point;

    NUM_NEIGHBOR = Num_neighbor;

	point.x = p->x;
	point.y = p->y;
	point.z = p->z;

	now = 0;
	maxdist = maxdist2 = 1.2E24;
	ptr = (void *) tree;
	while(ptr != NULL){
		switch(((TYPE*)ptr)->type){
			case TYPE_TREE:
				switch(NEAR_OPEN(point,ptr,now,maxdist)){
					case YES:
						ptr = (void *)(((TStruct*)ptr)->daughter);
						break;
					default :
						ptr = (void *)(((TStruct*)ptr)->sibling);
				}
				break;
			default :
				tmpx = point.x - ((TPtlStruct*)ptr)->r[0];
				tmpy = point.y - ((TPtlStruct*)ptr)->r[1];
				tmpz = point.z - ((TPtlStruct*)ptr)->r[2];
				dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
				if(now < NUM_NEIGHBOR || dist2 < maxdist2) INSERT(dist2,maxdist,maxdist2);
				ptr = (void*)(((TPtlStruct*)ptr)->sibling);
		}

	}
	/*
	*result = 0.5*(sqrt(neighbor[now-1].dist2)+sqrt(neighbor[now-2].dist2));
	*/
	*result = (sqrt(neighbor[now-1].dist2));
	for(i=0;i<now;i++){
		nearindex[i] = neighbor[i].index;
		DIST2[i] = neighbor[i].dist2;
	}
	/*
	{
		for(i=0;i<now;i++){
			printf("%g ",sortdist2[i]);
		}
		printf("\n");
	}
	*/
	return now;
}
