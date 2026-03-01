/* Thi code is to make the linked merging tree */
/* It is made by Juhan Kim and tested on 06/03/2014. */

#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<limits.h>



#define DEFINE_SIM_PARA
#include "pmheader.h"
#undef DEFINE_SIM_PARA



#include "merger.h"


#define HaloLinkFile "HaloLinkedList"
#define HaloRVFile   "HaloRVList"



#define BitsOfChar 8L
#define MIN(a,b) ((a)<(b)? (a):(b))
#define MAX(a,b) ((a)>(b)? (a):(b))


int *p2halo;

TrHalo *nowmbplist,*nextmbplist;
MbpRV *mbprv;
#define NO 0x00
#define YES 0x01
#define WasMB(a,n) ((a>>n) & YES)

IDTYPE nnowmbplist,nnextmbplist,nnextsubhalo;

int sortNowHidTrHalo(const void *a, const void *b){
	TrHalo *aa, *bb;
	aa = (TrHalo*)a;
	bb = (TrHalo*)b;
	if(aa->nowhid < bb->nowhid) return -1;
	else if(aa->nowhid > bb->nowhid) return 1;
	else return 0;
}
int sortAidTrHalo(const void *a, const void *b){
	TrHalo *aa, *bb;
	aa = (TrHalo*)a;
	bb = (TrHalo*)b;
	if(aa->aid < bb->aid) return -1;
	else if(aa->aid > bb->aid) return 1;
	else return 0;
}

int sortNextHidLink(const void *a, const void *b){
	UpLink *aa, *bb;
	aa = (UpLink*)a;
	bb = (UpLink*)b;
	if(aa->nexthid < bb->nexthid) return -1;
	else if(aa->nexthid > bb->nexthid) return 1;
	else return 0;
}

int sortMbTrHalo(const void *a, const void *b){
	TrHalo *aa, *bb;
	aa = (TrHalo*)a;
	bb = (TrHalo*)b;
	if(aa->mbp < bb->mbp) return -1;
	else if(aa->mbp > bb->mbp) return 1;
	else return 0;
}

int sortMbUpLink(const void *a, const void *b){
	UpLink *aa, *bb;
	aa = (UpLink*)a;
	bb = (UpLink*)b;
	if(aa->mbp < bb->mbp) return -1;
	else if(aa->mbp > bb->mbp) return 1;
	else return 0;
}


float boxsize,hubble,npower,omep,omepb,omeplam;
float bias;
int nx,ny,nz,nspace;
float amax,astep,anow;
float redshift1,redshift2,redshift;


void GetHead(char *infile, FILE **fp){
	fread(&boxsize,sizeof(float),1,*fp);
	fread(&hubble,sizeof(float),1,*fp);
	fread(&npower,sizeof(float),1,*fp);
	fread(&omep,sizeof(float),1,*fp);
	fread(&omepb,sizeof(float),1,*fp);
	fread(&omeplam,sizeof(float),1,*fp);
	fread(&bias,sizeof(float),1,*fp);
	fread(&nx,sizeof(int),1,*fp);
	fread(&nspace,sizeof(int),1,*fp);
	fread(&amax,sizeof(float),1,*fp);
	fread(&astep,sizeof(float),1,*fp);
	fread(&anow,sizeof(float),1,*fp);
}


#define NPSTEP 10000000
void *GetHalo(char *infile, IDTYPE *mhalo){
	FILE *fp = fopen(infile,"r");
	HaloQ halo[1000000];
	TrHalo *trhalo;
	IDTYPE nhalo = 0;
	IDTYPE maxnhalo = NPSTEP;
	int np;
	trhalo = (TrHalo*)malloc(sizeof(TrHalo)*maxnhalo);
	fread(&redshift,sizeof(float),1,fp);
	while((np = fread(halo,sizeof(HaloQ),1000000,fp))>0){
		IDTYPE i;
		for(i=0;i<np;i++){
			trhalo[nhalo+i].mbp = halo[i].mbp;
			trhalo[nhalo+i].nowhid = halo[i].hid;
		}
		nhalo += np;
		if(nhalo + NPSTEP > maxnhalo){
			maxnhalo += NPSTEP;
			trhalo = (TrHalo *)realloc(trhalo,sizeof(TrHalo)*maxnhalo);
		}
	}
	fclose(fp);
	*mhalo = nhalo;
	return (void*)trhalo;
}


void *GetHalo2(int nstep, IDTYPE *mhalo){
	char infile1[190],infile2[190];
	FILE *fp1,*fp2;
	IDTYPE np;
	TrHalo *a = NULL;

	*mhalo = 0;

	sprintf(infile1,HaloLinkFile".%.5d",nstep);
	printf("%s is being opened \n",infile1);
	if((fp1 = fopen(infile1,"r")) != NULL){
		fread(&redshift,sizeof(float),1,fp1);
		fread(&np,sizeof(IDTYPE),1,fp1);
		a = (TrHalo*)malloc(sizeof(TrHalo)*np);
		fread(a,sizeof(TrHalo),np,fp1);
		fclose(fp1);
		*mhalo += np;
		printf("%ld is read\n",np);
	}
	return a;
}
void GetMask(int nstep, char *mask, IDTYPE *mhalo){
	char infile1[190],infile2[190];
	FILE *fp1,*fp2;
	IDTYPE np;
	TrHalo *a;
	int NNN = 3000000;


	a = (TrHalo*)malloc(sizeof(TrHalo)*NNN);

	*mhalo = 0;

	sprintf(infile1,HaloLinkFile".%.5d",nstep);
	printf("%s is being opened \n",infile1);
	if((fp1 = fopen(infile1,"r")) != NULL){
		fread(&redshift,sizeof(float),1,fp1);
		fread(&np,sizeof(IDTYPE),1,fp1);
		while((np = fread(a,sizeof(TrHalo),NNN,fp1))>0){
			int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for(i=0;i<np;i++){
				INDXTYPE ii = (a[i].mbp>>3);
				int ioff =  (int)(a[i].mbp & 7L);
				mask[ii] |= (1<<ioff);
			}
			*mhalo += np;
		}
		fclose(fp1);
		printf("mask is set for np = %ld\n",*mhalo);
	}
	free(a);
}



#define npstep 50000000


int sorthid(const void *a, const void *b){
	IDTYPE *aa,*bb;
	aa = (IDTYPE*) a;
	bb = (IDTYPE*) b;
	if(*aa < *bb) return -1;
	else if(*aa>*bb) return 1;
	else return 0;
}

void FindMajorMergermbp(UpLink *uplink,IDTYPE nuplink,TrHalo *nextmbplist,IDTYPE nnextmbplist){
	IDTYPE i,j,k,ilink,kk;
	if(nuplink ==0) return;
	j = 0;
	for(i=0;i<nnextmbplist;){
		if(j>=nuplink) break;
		if(uplink[j].nexthid < nextmbplist[i].nowhid ) {
			while(uplink[j].nexthid < nextmbplist[i].nowhid && j < nuplink) j++;
			continue;
		}
		else if(uplink[j].nexthid > nextmbplist[i].nowhid){
			while(uplink[j].nexthid > nextmbplist[i].nowhid && i < nnextmbplist) i++;
			continue;
		}
		if(j>=nuplink || i >=nnextmbplist) break;

		/*
		if(nextmbplist[i].nowhid == 320786){
			printf("here\n");
		}
		*/

		IDTYPE  nexthid = nextmbplist[i].nowhid;
		int jlink = 0;
		while(nexthid == uplink[j+jlink].nexthid){
			jlink ++;
			if(j+jlink == nuplink) break;
		}
		if(jlink > 1){
			float maxmass = uplink[j].mass;
			IDTYPE jmajor = j;
			for(k=j;k<j+jlink;k++){
				if(uplink[k].mass > maxmass){
					maxmass = uplink[k].mass;
					jmajor = k;
				}
			}
			for(k=j;k<j+jlink;k++){
				if(k==jmajor) SET_MAJOR_LINE(uplink,k,linkflag);
				else SET_MINOR_LINE(uplink,k,linkflag);
			}
			for(k=i;k<nnextmbplist;k++){
				if(nextmbplist[k].nowhid != nexthid){
					ilink = k;
					break;
				}
			}
			for(k=i;k<ilink;k++){
				if(nextmbplist[k].mbp == uplink[jmajor].mbp){
					SET_MAJOR_SATELLITE(nextmbplist,k,statusflag);
					nextmbplist[k].downaid = uplink[jmajor].aid;
					uplink[jmajor].upaid = nextmbplist[k].aid;
				}
				else {
					SET_MINOR_SATELLITE(nextmbplist,k,statusflag);
					/* */
					for(kk = j;kk<j+jlink;kk++){
						if(uplink[kk].mbp == nextmbplist[k].mbp){
							nextmbplist[k].mass = uplink[kk].mass;
							nextmbplist[k].downaid = uplink[kk].aid;
							uplink[kk].upaid = nextmbplist[k].aid;
						}
					}
				}
			}
			i = ilink;
		}
		else if(jlink ==1){
			SET_MAJOR_LINE(uplink,j,linkflag);
			nextmbplist[i].downaid = uplink[j].aid;
			uplink[j].upaid = nextmbplist[i].aid;
			i ++;
		}
		else {
			i++;
		}
		j += jlink;
	}
}

double rscale,vscale;


int main(int argc, char **argv){
	char infile1[190],infile2[100];
	char infile3[190],infile4[100];
	char infile5[190],infile6[100];
	char inmerge[190];
	char outfile[190];
	FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6;
	int nstep1,nstep2;
	float pmass;
	IDTYPE i,j;





	nstep1 = atoi(argv[1]);
	nstep2 = atoi(argv[2]);
	if(nstep1 >= nstep2){
		fprintf(stderr,"Error large nstep1 than nstep2\n");
		exit(99);
	}

	sprintf(infile1,"FoFHaloQuantities%.5d.bin",nstep1);
	sprintf(infile2,"FoFHaloQuantities%.5d.bin",nstep2);

	sprintf(infile3,"FoF_halo_cat.%.5d",nstep2);
	sprintf(infile4,"FoF_member_particle.%.5d",nstep2);
	sprintf(infile5,"FoF_halo_cat.%.5d",nstep1);
	sprintf(infile6,"FoF_member_particle.%.5d",nstep1);


	/*
	nowmbplist = GetHalo2(nstep1,&nnowmbplist);
	*/

	fp5 = fopen(infile5,"r");
	GetHead(infile5,&fp5);
	fclose(fp5);
	redshift1 = amax/anow-1;


	ny = nz = nx;

	pmass = 2.7755e11*omep*pow((double)boxsize/(double)nx,3.L);



	INDXTYPE nsize = ((long)nx*(long)ny*(long)nz-1L)/BitsOfChar + 1L;
	char *mask = (char *)calloc(nsize,sizeof(char));
	
	/* Mask for the most bound member particle */
	/*
	for(i=0;i<nnowmbplist;i++){
		INDXTYPE ii = (nowmbplist[i].mbp>>3);
		int ioff =  (nowmbplist[i].mbp & 7L);
		mask[ii] |= (1<<ioff);
	}
	*/
	GetMask(nstep1,mask,&nnowmbplist);
	if(0){
		int ioff,icount = 0,jcount = 0;
		for(i=0;i<nsize;i++){	
			if(mask[i] !=0) {
				jcount ++;
			}
			for(ioff=0;ioff<8*sizeof(char);ioff++) {
				if(WasMB(mask[i],ioff)) icount ++;
			}
		}
		printf("it has mask out pixels %d %d\n",icount,jcount);fflush(stdout);
	}

	IDTYPE iuplink=0;
	IDTYPE nexthid=0;
	UpLink *uplink = (UpLink*) malloc(sizeof(UpLink)*nnowmbplist);

	for(i=0;i<nnowmbplist;i++) RESET_WHOLE_FLAGS_MP(uplink,i,linkflag);



	fp2 = fopen(infile2,"r");

	IDTYPE nmaxnextmbplist = 10000000 + nnowmbplist;
	nextmbplist = (TrHalo*)malloc(sizeof(TrHalo)*nmaxnextmbplist);
	mbprv = (MbpRV*)malloc(sizeof(MbpRV)*nmaxnextmbplist);
	for(i=0;i<nmaxnextmbplist;i++) {
		RESET_WHOLE_FLAGS_MP(nextmbplist,i,statusflag);
		RESET_WHOLE_FLAGS_MP(nextmbplist,i,linkflag);
	}
	nnextmbplist = 0;


	fp3 = fopen(infile3,"r");
	GetHead(infile3,&fp3);
	redshift2 = amax/anow -1;
	fp4 = fopen(infile4,"r");

	double red2p1 = redshift2+1;

	rscale = boxsize/nx;
	double Hsub = sqrt(omep*red2p1*red2p1*red2p1 + omeplam +(1-omep-omeplam)*red2p1*red2p1);
	vscale = 100.L*boxsize/hubble/amax*anow*anow*hubble*Hsub; /* In unit of km/second */


	simpar.mx = simpar.nx = simpar.ny = simpar.nz = nx;
	simpar.nspace = 1;
	simpar.mxmy  = simpar.nx*simpar.ny;


	printf("Now making the linking \n");fflush(stdout);


	IDTYPE nmaxbp = 10000000;
	basicparticletype *bp = (basicparticletype*)malloc(nmaxbp*sizeof(basicparticletype));
	IDTYPE aid = 0;

	Halo halo;
	HaloQ haloq;

	while(fread(&halo,sizeof(Halo),1,fp3)==1){
		fread(&haloq,sizeof(HaloQ),1,fp2);
		if(halo.np > nmaxbp ) {
			nmaxbp += 1000000;
			bp = (basicparticletype*)realloc(bp,sizeof(basicparticletype)*nmaxbp);
		}
		fread(bp,sizeof(basicparticletype),halo.np,fp4);

		if(halo.np ==0 || halo.np < NMINMEMBERS) continue;

		int nsub = 0;

		float halomass = pmass*halo.np;
		int imajor = 0;

		for(i=0;i<halo.np;i++){
			INDXTYPE ii = (bp[i].indx >>3);
			INDXTYPE ioff = bp[i].indx & 7L;
			if(WasMB(mask[ii],ioff)){
				SET_ONE_OF_MULTI_SUB(nextmbplist,nnextmbplist,statusflag);
				uplink[iuplink].nexthid = nexthid;
				uplink[iuplink].mbp = bp[i].indx;
				uplink[iuplink].upaid = TERMINATE;

				SET_MULTI_LINE(uplink,iuplink,linkflag);
				nextmbplist[nnextmbplist].aid = aid++;
				nextmbplist[nnextmbplist].nowhid = nexthid;
				mbprv[nnextmbplist].x = XofP(bp+i)*rscale;
				mbprv[nnextmbplist].y = YofP(bp+i)*rscale;
				mbprv[nnextmbplist].z = ZofP(bp+i)*rscale;
				mbprv[nnextmbplist].vx = bp[i].vx*vscale;
				mbprv[nnextmbplist].vy = bp[i].vy*vscale;
				mbprv[nnextmbplist].vz = bp[i].vz*vscale;
				nextmbplist[nnextmbplist].mass = haloq.mass;
				nextmbplist[nnextmbplist].fofhmass = haloq.mass;
				nextmbplist[nnextmbplist].mbp = bp[i].indx;
				nextmbplist[nnextmbplist].downaid = TERMINATE;
				nextmbplist[nnextmbplist].upaid = TERMINATE;
				nnextmbplist ++;
				nsub ++;
				iuplink ++;
			}
		}
		/*
                if(nexthid == 6307269L) {
                     printf("Now look: %d nsub\n",nsub);
                }
				*/

		if(nsub < 2){
			if(nsub == 0) {
				SET_NEW_BORN(nextmbplist,nnextmbplist,statusflag);
				IDTYPE ii = -1;
				for(i=0;i<halo.np;i++){
					if(bp[i].indx == haloq.mbp) {
						ii = i;
						break;
					}
				}
				nextmbplist[nnextmbplist].aid = aid++;
				nextmbplist[nnextmbplist].nowhid = nexthid;
				nextmbplist[nnextmbplist].mbp = haloq.mbp;
				mbprv[nnextmbplist].x = XofP(bp+ii)*rscale;
				mbprv[nnextmbplist].y = YofP(bp+ii)*rscale;
				mbprv[nnextmbplist].z = ZofP(bp+ii)*rscale;
				mbprv[nnextmbplist].vx = bp[ii].vx*vscale;
				mbprv[nnextmbplist].vy = bp[ii].vy*vscale;
				mbprv[nnextmbplist].vz = bp[ii].vz*vscale;
				nextmbplist[nnextmbplist].mass = haloq.mass;
				nextmbplist[nnextmbplist].fofhmass = haloq.mass;
				nextmbplist[nnextmbplist].downaid = TERMINATE;
				nextmbplist[nnextmbplist].upaid = TERMINATE;
				nnextmbplist++;
			}
			else {

				SET_ALONE_SUB(nextmbplist,nnextmbplist-1,statusflag);
				UNSET_ONE_OF_MULTI_SUB(nextmbplist,nnextmbplist-1,statusflag);

				SET_SINGLE_LINE(uplink,iuplink-1,linkflag);
				UNSET_MULTI_LINE(uplink,iuplink-1,linkflag);
				
			}
		}
		if(nmaxnextmbplist - nnextmbplist < 100000){
			nmaxnextmbplist += 1000000;
			nextmbplist = (TrHalo*)realloc(nextmbplist,sizeof(TrHalo)*nmaxnextmbplist);
			mbprv = (MbpRV*)realloc(mbprv,sizeof(MbpRV)*nmaxnextmbplist);
			for(i=nnextmbplist;i<nmaxnextmbplist;i++) {
				RESET_WHOLE_FLAGS_MP(nextmbplist,i,statusflag);
				RESET_WHOLE_FLAGS_MP(nextmbplist,i,linkflag);
			}
		}
		nexthid ++;
	}
	fclose(fp2);fclose(fp4);fclose(fp3);
	free(bp);free(mask);
	printf("Well made the link nnextmbplist/iuplink = %ld %ld\n",nnextmbplist,iuplink);fflush(stdout);





	char out2[190];
	FILE *wp;
	sprintf(out2,HaloRVFile".%.5d",nstep2);
	if(nnextmbplist>0){
		FILE *wp2 = fopen(out2,"w");
		fwrite(mbprv,sizeof(MbpRV),nnextmbplist,wp2);
		fclose(wp2);
	}
	free(mbprv);

	/* reading nowmbplist */
	nowmbplist = GetHalo2(nstep1,&nnowmbplist);

	IDTYPE nuplink = iuplink;
	if(nuplink != nnowmbplist){
		fprintf(stderr,"missing mbp between nowmbplist and uplink %ld\n",nnowmbplist-nuplink);
	}

	/*
	for(i=0;i<nnowmbplist;i++) if(nowmbplist[i].aid==24) printf("nowmbplist[24] %ld %ld %ld\n",nowmbplist[i].aid, 
			nowmbplist[i].upaid,nowmbplist[i].downaid);
			*/


	printf("in the massive sorting\n");
	/* dump to nowmbplist link flag */
	qsort(nowmbplist,nnowmbplist,sizeof(TrHalo),sortMbTrHalo);
	qsort(uplink,nuplink,sizeof(UpLink),sortMbUpLink);
	j = 0;
	for(i=0;i<nnowmbplist;i++){
		int icount = 0;
		for(;j<nuplink;j++) {
			if(nowmbplist[i].mbp == uplink[j].mbp) {
				uplink[j].mass = nowmbplist[i].mass;
				uplink[j].nowhid = nowmbplist[i].nowhid;
				uplink[j].aid = nowmbplist[i].aid;
				nowmbplist[i].linkflag = uplink[j].linkflag;
				icount ++;
			}
			else if(uplink[j].mbp > nowmbplist[i].mbp) break;
		}
		if(icount ==0) SET_MISSING_LINE(nowmbplist,i,linkflag);
	}

	/* Find the major merger line */
	/*
	qsort(nextmbplist,nnextmbplist,sizeof(TrHalo),sortAidTrHalo);
	*/
	qsort(nextmbplist,nnextmbplist,sizeof(TrHalo),sortNowHidTrHalo);
	qsort(uplink,nuplink,sizeof(UpLink),sortNextHidLink);
	printf("in identifying the major satellite \n");
	FindMajorMergermbp(uplink,nuplink,nextmbplist,nnextmbplist);

	qsort(uplink,nuplink,sizeof(UpLink),sortMbUpLink);
	j = 0;
	for(i=0;i<nnowmbplist;i++){
		int icount = 0;
		for(;j<nuplink;j++) {
			if(nowmbplist[i].mbp == uplink[j].mbp) {
				nowmbplist[i].linkflag = uplink[j].linkflag;
				nowmbplist[i].upaid = uplink[j].upaid;
				icount ++;
			}
			else if(uplink[j].mbp > nowmbplist[i].mbp) break;
		}
	}
	free(uplink);



	qsort(nowmbplist,nnowmbplist,sizeof(TrHalo),sortAidTrHalo);

	printf("in saving data 1 \n");
	if(nnowmbplist >0){
		sprintf(infile1,HaloLinkFile".%.5d",nstep1);
		FILE *wp2 = fopen(infile1,"w");
		fwrite(&redshift1,sizeof(float),1,wp2);
		fwrite(&nnowmbplist,sizeof(IDTYPE),1,wp2);
		fwrite(nowmbplist,sizeof(TrHalo),nnowmbplist,wp2);
		fclose(wp2);
	}


	printf("in saving data 2 \n");
	qsort(nextmbplist,nnextmbplist,sizeof(TrHalo),sortAidTrHalo);
	if(nnextmbplist >0){
		sprintf(out2,HaloLinkFile".%.5d",nstep2);
		FILE *wp2 = fopen(out2,"w");
		fwrite(&redshift2,sizeof(float),1,wp2);
		fwrite(&nnextmbplist,sizeof(IDTYPE),1,wp2);
		fwrite(nextmbplist,sizeof(TrHalo),nnextmbplist,wp2);
		fclose(wp2);
	}


	printf("saved to %s with %ld mbp's\n",out2,nnextmbplist);
	return 0;
}
