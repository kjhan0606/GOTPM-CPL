/* This code outputs the major merging tree of the most massive halos at the input time step */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<string.h>

#include "merger.h"

#define NSTEPS 76
int nsteps[NSTEPS] ={107, 136, 165, 183, 205, 234, 250, 270, 292, 318, 
	321, 353, 385, 415, 444, 473, 501,  555, 632, 657, 681, 
	706, 730, 754, 801, 824, 916, 938, 961, 983, 1006, 1028, 1051, 
	1073, 1095, 1118, 1140, 1162, 1184, 1207, 1229, 1244, 1252, 1274, 1297, 1319, 
	1328, 1342, 1365, 1388, 1410, 1434, 1457, 1480, 1503, 1527, 1535, 1550, 
	1574, 1598, 1622, 1646, 1670, 1694, 1738, 1744, 1769, 1794, 1819, 1844, 
	1870, 1896, 1922, 1948, 1974, 2001};


int main(int argc, char **argv){
	IDTYPE i,j,k;
	int fstep,istep,istart;
	IDTYPE np,ip;
	TrHalo *nmbp,imbp;
	MbpRV imbprv;
	FILE *fp1;


	printf("###########################################################################\n");
	printf("This program first find the most massive halo at the given time step\n");
	printf("Then, find the merging tree upward.\n");
	printf("###########################################################################\n");

	printf("please input the last step to get the major merging & aid\n");
	scanf("%d",&fstep);

	for(i=0;i<NSTEPS;i++){
		if(nsteps[i] == fstep){
			istart = i;
			break;
		}
	}

	char infile[190],infile1[190];
	sprintf(infile,"HaloLinkedList.%.5d",nsteps[istart]);
	FILE *fp = fopen(infile,"rb");
	float nredshift,iredshift;
	fread(&nredshift,sizeof(float),1,fp);
	fread(&np,sizeof(IDTYPE),1,fp);
	nmbp = (TrHalo*)malloc(sizeof(TrHalo)*np);
	fread(nmbp,sizeof(TrHalo),np,fp); 
	fclose(fp);

	IDTYPE aid;
	float mass = -1E20;
	for(i=0;i<np;i++){
		if(nmbp[i].mass > mass) {
			aid = i;
			mass = nmbp[i].mass;
		}

	}
	scanf("%ld",&aid);
	if(aid >= np) return 0;

	printf("The most massive halo or mbp is selected: mbp = %ld m= %g  upaid= %ld in np= %ld nexthid = %ld\n",
			nmbp[aid].mbp,nmbp[aid].mass, nmbp[aid].upaid, np, nmbp[aid].nexthid);
	imbp = nmbp[aid];
	printf("%ld nstep=%d : z= %g , hid=%ld aid= %ld downaid = %ld upaid = %ld mass/hmass= %g %gx/y/z= %g %g %g vx/y/z= %g %g %g\n",
			imbp.mbp,nsteps[istart],
			nredshift,nmbp[aid].nowhid,aid,imbp.downaid,imbp.upaid,imbp.mass, imbp.fofhmass,
			imbprv.x,imbprv.y,imbprv.z,imbprv.vx,imbprv.vy,imbprv.vz);

	/* nmbp is already sorted with aid. */
	IDTYPE nextaid = nmbp[aid].upaid;
	if(nextaid == TERMINATE) return 0;


	for(i=istart +1; i<NSTEPS  ;i++){
		IDTYPE inp;
		sprintf(infile,"HaloLinkedList.%.5d",nsteps[i]);
		sprintf(infile1,"HaloRVList.%.5d",nsteps[i]);
		fp = fopen(infile,"rb");
		if(fp == NULL) break;
		fread(&iredshift,sizeof(float),1,fp);
		fread(&inp,sizeof(IDTYPE),1,fp);
		fseek(fp,sizeof(TrHalo)*nextaid,SEEK_CUR);
		fread(&imbp,sizeof(TrHalo),1,fp); 
		fclose(fp);
		fp1 = fopen(infile1,"rb");
		fseek(fp1,sizeof(MbpRV)*nextaid,SEEK_SET);
		fread(&imbprv,sizeof(MbpRV),1,fp1);
		fclose(fp1);
		printf("%ld nstep=%d : z= %g , hid=%ld aid= %ld downaid = %ld upaid = %ld mass/hmass= %g %g x/y/z= %g %g %g vx/y/z= %g %g %g\n",
				imbp.mbp,nsteps[i],
				iredshift,imbp.nowhid,nextaid,imbp.downaid,imbp.upaid,imbp.mass, imbp.fofhmass,
				imbprv.x,imbprv.y,imbprv.z,imbprv.vx,imbprv.vy,imbprv.vz);
		nextaid = imbp.upaid;
		if(nextaid  == TERMINATE) break;
	}
	return 0;

}

