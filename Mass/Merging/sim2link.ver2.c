/* Thi code is to make the linked merging tree */
/* It is made by Juhan Kim and tested on 06/03/2014. */
/* It is changed to host an extensive pthreads */
/* Un confirmed .Please do not use this one. It is not as good as seemed. */

#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<limits.h>
#include<omp.h>


int nullfct0(){ return 0;}
int nullfct1(){ return 1;}


#ifndef _OPENMP
#define omp_get_thread_num() nullfct0()
#define omp_get_num_threads() nullfct1()
#endif



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
#pragma omp parallel for
			for(i=0;i<np;i++){
				INDXTYPE ii = (a[i].mbp>>3);
				int ioff =  (a[i].mbp & 7L);
				mask[ii] |= (1<<ioff);
			}
			*mhalo += np;
		}
		fclose(fp1);
		printf("mask is set\n");
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
	IDTYPE i,j,k,jj,kk;
	if(nuplink ==0) return;
	j = 0;
	for(i=0;i<nnextmbplist;){
		HIDTYPE  nexthid = nextmbplist[i].nowhid;
		int ilink = 0;
		while(nexthid == uplink[j+ilink].nexthid){
			ilink ++;
		}
		if(ilink > 1){
			float maxmass = uplink[j].mass;
			IDTYPE ii = j;
			for(k=j;k<j+ilink;k++){
				if(uplink[k].mass > maxmass){
					maxmass = uplink[k].mass;
					ii = k;
				}
			}
			for(k=j;k<j+ilink;k++){
				if(k==ii) SET_MAJOR_LINE(uplink,k,linkflag);
				else SET_MINOR_LINE(uplink,k,linkflag);
			}
			for(k=i;k<nnextmbplist;k++){
				if(nextmbplist[k].nowhid != nexthid){
					jj = k;
					break;
				}
			}
			for(k=i;k<jj;k++){
				if(nextmbplist[k].mbp == uplink[ii].mbp){
					SET_MAJOR_SATELLITE(nextmbplist,k,statusflag);
					nextmbplist[k].downaid = uplink[ii].aid;
					uplink[ii].upaid = nextmbplist[k].aid;
				}
				else {
					SET_MINOR_SATELLITE(nextmbplist,k,statusflag);
					/* */
					for(kk = j;kk<j+ilink;kk++){
						if(uplink[kk].mbp == nextmbplist[k].mbp){
							nextmbplist[k].mass = uplink[kk].mass;
							nextmbplist[k].downaid = uplink[kk].aid;
							uplink[kk].upaid = nextmbplist[k].aid;
						}
					}
				}
			}
		}
		else if(ilink ==1){
			SET_MAJOR_LINE(uplink,j,linkflag);
			nextmbplist[i].downaid = uplink[j].aid;
			uplink[j].upaid = nextmbplist[i].aid;
			i ++;
		}
		else {
			i++;
		}
		j += ilink;
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
	IDTYPE i,j,k;





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



	INDXTYPE nsize = (long)nx*(long)ny*(long)nz/BitsOfChar + 1L;
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

	IDTYPE iuplink=0;
	HIDTYPE nexthid=0;
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

	Halo halo[1000];
	HaloQ haloq[1000];
	int runningnhalo[1000];
	int nhalo;

	while((nhalo = fread(&halo,sizeof(Halo),1000,fp3))>0){
		fread(haloq,sizeof(HaloQ),1000,fp2);
		int tnhalo = 0;
		for(j=0;j<nhalo;j++) {
			runningnhalo[j] = tnhalo;
			tnhalo += halo[j].np;
		}
		if(tnhalo > nmaxbp){
			nmaxbp += 1000000;
			bp = (basicparticletype*)realloc(bp,sizeof(basicparticletype)*nmaxbp);
		}
		fread(bp,sizeof(basicparticletype),tnhalo,fp4);
		i = 0;
		int kk = 0;
		for(j=0;j<nhalo;j++){
			if(halo[j].np > 30) {
				halo[i] = halo[j];
				haloq[i] = haloq[j];
				for(k=0;k<halo[j].np;k++){
					bp[kk+k] = bp[runningnhalo[j]+k];
				}
				runningnhalo[i] = kk;
				kk += halo[j].np;
				i++;
			}
		}
		nhalo = i;

		int *runningiuplink;
		int *runningnnextmbplist;
#ifdef _OPENMP
#pragma omp parallel 
#endif
		{
			if(omp_get_thread_num() == 0){
				runningiuplink = (int *) malloc(sizeof(int)*omp_get_num_threads());
				runningnnextmbplist = (int *)malloc(sizeof(int)*omp_get_num_threads());
			}
			UpLink *p_uplink;
			TrHalo *p_nextmbplist;
			MbpRV  *p_mbprv;
			IDTYPE p_iuplink,p_nnextmbplist;
			p_uplink = (UpLink*)malloc(sizeof(UpLink)*tnhalo);
			p_nextmbplist = (TrHalo *)malloc(sizeof(TrHalo )*tnhalo);
			p_mbprv = (MbpRV *)malloc(sizeof(MbpRV )*tnhalo);

			int trID,nthr;
			trID = omp_get_thread_num();
			nthr = omp_get_num_threads();
			int chunk = (nhalo+nthr-1)/nthr;
			int jstart = chunk*trID;
			int jend = MIN(jstart+chunk,nhalo);
			int jhalo;


			p_iuplink = p_nnextmbplist = 0;


			for(jhalo=jstart;j<jend;j++){
				int nsub = 0;
		
				float halomass = pmass*halo[jhalo].np;
				int imajor = 0;
		
				for(i=0;i<halo[jhalo].np;i++){
					INDXTYPE ii = (bp[runningnhalo[jhalo]+i].indx >>3);
					INDXTYPE ioff = bp[runningnhalo[jhalo]+i].indx & 7L;
					if(WasMB(mask[ii],ioff)){
						SET_ONE_OF_MULTI_SUB(p_nextmbplist,p_nnextmbplist,statusflag);
						p_uplink[p_iuplink].nexthid = jhalo+nexthid;
						p_uplink[p_iuplink].mbp = bp[runningnhalo[jhalo]+i].indx;
						p_uplink[p_iuplink].upaid = TERMINATE;
		
						SET_MULTI_LINE(p_uplink,p_iuplink,linkflag);
						p_nextmbplist[p_nnextmbplist].aid = aid++;
						p_nextmbplist[p_nnextmbplist].nowhid = jhalo+nexthid;
						p_mbprv[nnextmbplist].x = XofP(bp+runningnhalo[jhalo]+i)*rscale;
						p_mbprv[nnextmbplist].y = YofP(bp+runningnhalo[jhalo]+i)*rscale;
						p_mbprv[nnextmbplist].z = ZofP(bp+runningnhalo[jhalo]+i)*rscale;
						p_mbprv[nnextmbplist].vx = bp[runningnhalo[jhalo]+i].vx*vscale;
						p_mbprv[nnextmbplist].vy = bp[runningnhalo[jhalo]+i].vy*vscale;
						p_mbprv[nnextmbplist].vz = bp[runningnhalo[jhalo]+i].vz*vscale;
						p_nextmbplist[nnextmbplist].mass = haloq[jhalo].mass;
						p_nextmbplist[nnextmbplist].mbp = bp[runningnhalo[jhalo]+i].indx;
						p_nextmbplist[nnextmbplist].downaid = TERMINATE;
						p_nnextmbplist ++;
						nsub ++;
						p_iuplink ++;
					}
				}
				if(nsub < 2){
					if(nsub == 0) {
						SET_NEW_BORN(p_nextmbplist,p_nnextmbplist,statusflag);
						IDTYPE ii = -1;
						for(i=0;i<halo[jhalo].np;i++){
							if(bp[runningnhalo[jhalo]+i].indx == haloq[jhalo].mbp) {
								ii = i;
								break;
							}
						}
						p_nextmbplist[p_nnextmbplist].aid = aid++;
						p_nextmbplist[p_nnextmbplist].nowhid = jhalo+nexthid;
						p_nextmbplist[p_nnextmbplist].mbp = haloq[jhalo].mbp;
						p_mbprv[nnextmbplist].x = XofP(bp+runningnhalo[jhalo]+ii)*rscale;
						p_mbprv[nnextmbplist].y = YofP(bp+runningnhalo[jhalo]+ii)*rscale;
						p_mbprv[nnextmbplist].z = ZofP(bp+runningnhalo[jhalo]+ii)*rscale;
						p_mbprv[nnextmbplist].vx = bp[runningnhalo[jhalo]+ii].vx*vscale;
						p_mbprv[nnextmbplist].vy = bp[runningnhalo[jhalo]+ii].vy*vscale;
						p_mbprv[nnextmbplist].vz = bp[runningnhalo[jhalo]+ii].vz*vscale;
						p_nextmbplist[nnextmbplist].mass = haloq[jhalo].mass;
						p_nextmbplist[nnextmbplist].downaid = TERMINATE;
						p_nnextmbplist++;
					}
					else {
		
						SET_ALONE_SUB(p_nextmbplist,p_nnextmbplist-1,statusflag);
						UNSET_ONE_OF_MULTI_SUB(p_nextmbplist,p_nnextmbplist-1,statusflag);
		
						SET_SINGLE_LINE(p_uplink,p_iuplink-1,linkflag);
						UNSET_MULTI_LINE(p_uplink,p_iuplink-1,linkflag);
						
					}
				}
			}
			runningiuplink[trID] = p_iuplink;
			runningnnextmbplist[trID] = p_nnextmbplist;
#pragma omp barrier
			if(trID ==0){
				for(i=1;i<trID;i++){
					runningiuplink[trID] += runningiuplink[trID-1];
					runningnnextmbplist[trID] += runningnnextmbplist[trID-1];
				}
				if(nmaxnextmbplist - nnextmbplist < runningnnextmbplist[nthr-1]){
					nmaxnextmbplist += (nnextmbplist +runningnnextmbplist[nthr-1]+10000);
					nextmbplist = (TrHalo*)realloc(nextmbplist,sizeof(TrHalo)*nmaxnextmbplist);
					mbprv = (MbpRV*)realloc(mbprv,sizeof(MbpRV)*nmaxnextmbplist);
					for(i=nnextmbplist;i<nmaxnextmbplist;i++) {
						RESET_WHOLE_FLAGS_MP(nextmbplist,i,statusflag);
						RESET_WHOLE_FLAGS_MP(nextmbplist,i,linkflag);
					}
				}
			}
#pragma omp barrier
			for(i=0;i<p_iuplink;i++){
				int istart ;
				if(trID==0) istart = 0;
				else istart = runningiuplink[trID-1];
				uplink[iuplink + istart +i] = p_uplink[i];
			}
			for(i=0;i<p_nnextmbplist;i++){
				int istart ;
				if(trID==0) istart = 0;
				else istart = runningnnextmbplist[trID-1];
				nextmbplist[nnextmbplist + istart +i] = p_nextmbplist[i];
				mbprv[nnextmbplist + istart +i] = p_mbprv[i];
			}
			if(trID==0){
				for(i=nnextmbplist;i<nnextmbplist+runningnnextmbplist[nthr-1];i++){
					nextmbplist[i].aid = i;
				}
				nexthid += nhalo;
				nnextmbplist += runningnnextmbplist[nthr-1];
				iuplink += runningiuplink[nthr-1];
			}
			free(p_uplink);
			free(p_nextmbplist);
			free(p_mbprv);
			if(omp_get_thread_num() == 0){
				free(runningiuplink);
				free(runningnnextmbplist);
			}

		}
	}
	fclose(fp2);fclose(fp4);fclose(fp3);
	free(bp);free(mask);
	printf("Well made the link \n");fflush(stdout);



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
	qsort(nextmbplist,nnextmbplist,sizeof(TrHalo),sortAidTrHalo);
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
