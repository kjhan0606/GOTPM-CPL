#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include "merger.h"



#define NO 0x00
#define YES 0x01
#define WasMB(a,n) ((a>>n) & YES)


int main(int argc, char **argv)
{
	FILE *fp;
	long nsize = 1200L*1200L*1200L/8L;

	fp = fopen(argv[1],"r");
	char *mask = (char *)calloc(nsize,sizeof(char));
	float red;
	IDTYPE nnowmbplist;


	fread(&red, sizeof(float),1,fp);
	fread(&nnowmbplist, sizeof(IDTYPE),1,fp);
	TrHalo buff[1000000];

	IDTYPE maxhid=-1;
	IDTYPE totnum = 0;
	IDTYPE ndup = 0;

	IDTYPE np,mp=0;
	while((np = fread(buff,sizeof(TrHalo),1000000,fp))>0)
	{
		int i;
		for(i=0;i<np;i++)
		{
			IDTYPE hid = buff[i].nowhid;
			IDTYPE ii = (hid >> 3);
			IDTYPE ioff = (hid & 7L);
			if(hid > maxhid) maxhid = hid;
			if(buff[i].mass == buff[i].fofhmass) {
				if(! WasMB(mask[ii],ioff)) {
					mask[ii] |= (1<<ioff);
					totnum ++;
				}
				else {
					printf("p%ld has duplication hid= %ld : mass %g %g\n",i,hid,buff[i].mass,buff[i].fofhmass);
					ndup ++;
				}
			}
		}
		mp += np;
	}
	fclose(fp);
	printf("%s Total Trhalo array = %ld :: total halo num %ld and maximum hid = %ld\n",argv[1],mp,totnum, maxhid);
	printf("duplication number is %ld\n",ndup);
	
}
