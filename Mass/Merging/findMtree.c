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
	printf("and find all the most bound particles (mbp's) in the halo\n");
	printf("And, then, trace them down to higher redshift printing the mass attached to the mbp\n");
	printf("###########################################################################\n");

	printf("please input the last step to get the major merging\n");
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

	/* aid is the id number of most bound particle at the redshift */
	IDTYPE aid;
	IDTYPE hid;
	float mass = -1E20;
	/* hid is the halo id number */
	for(i=0;i<np;i++){
		if(nmbp[i].mass > mass) {
			aid = i;
			mass = nmbp[i].mass;
			hid = nmbp[i].nowhid;
		}

	}
	printf("The most massive mbp is selected: m= %g\n",mass);
	/* Now find the mbp's in the most massive halos */
	IDTYPE *downaid = (IDTYPE*)malloc(sizeof(IDTYPE)*10000000);
	IDTYPE ndownaid = 0;
	for(i=0;i<np;i++){
		if(nmbp[i].nowhid == hid){
			downaid[ndownaid] = i;
			ndownaid ++;
		}
	}

	/* nmbp is already sorted with aid. */
	IDTYPE nextaid = nmbp[aid].downaid;


	for(i=istart; i>=0 ;i--){
		IDTYPE inp;
		sprintf(infile,"HaloLinkedList.%.5d",nsteps[i]);
		sprintf(infile1,"HaloRVList.%.5d",nsteps[i]);
		fp = fopen(infile,"rb");
		if(fp==NULL) break;
		fp1 = fopen(infile1,"rb");
		fread(&iredshift,sizeof(float),1,fp);
		fread(&inp,sizeof(IDTYPE),1,fp);
		fclose(fp);
		printf("nstep = %d redshift = %g\n",nsteps[i],iredshift);



		fp = fopen(infile,"r");
		for(j=0;j<ndownaid;j++){
			char type;
			IDTYPE id = downaid[j];
			if(id == TERMINATE) continue;
			IDTYPE offset = sizeof(float)+sizeof(IDTYPE)+sizeof(TrHalo)*id;
			fseek(fp,offset,SEEK_SET);
			fread(&imbp,sizeof(TrHalo),1,fp);
			if(IS_ALONE_SUB((&imbp),0,statusflag)) type = 'A'; /* Alone halo? */
			else if(IS_NEW_BORN((&imbp),0,statusflag)) type = 'N'; /* New born halo? */
			else if(IS_MAJOR_SATELLITE((&imbp),0,statusflag)) type = 'J'; /* Majro satellite? */
			else if(IS_MINOR_SATELLITE((&imbp),0,statusflag)) type = 'I'; /* Minor satellite? */
			else type = 'U';

			offset = sizeof(MbpRV)*id;
			fseek(fp1,offset,SEEK_SET);
			fread(&imbprv,sizeof(MbpRV),1,fp1);
			printf("  %c  hid= %d mbp= %ld aid= %d mass/hmass= %g/%g x/y/z= %g %g %g vx/y/z= %g %g %g\n",type,
					imbp.nowhid,
					imbp.mbp, id,imbp.mass,imbp.fofhmass,
					imbprv.x,imbprv.y,imbprv.z,imbprv.vx,imbprv.vy,imbprv.vz);
			downaid[j] = imbp.downaid;
		}
		if(fp!= NULL) fclose(fp);
		if(fp1!= NULL) fclose(fp1);
	}
	return 0;

}

