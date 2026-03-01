/* Not Yet Confirmed !!!! */
/* Thi code is to make the linked merging tree */
/* It is made by Juhan Kim and tested on 06/03/2014. */

#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
/*
#undef __GNU_C__
#define __USE_ISOC99
#include<limits.h>
#undef __USE_ISOC99
*/
#include<limits.h>
#include<features.h>
#include<omp.h>
#include<mpi.h>
#define DEFINE_SIM_PARA
#include "pmheader.h"
#undef DEFINE_SIM_PARA

#include "merger.h"


#include "pqsort.h"



/*
GenerateParallelComp(TrHalo,INDXTYPE,mbp);
GenerateParallelComp(TrHalo,INDXTYPE,nowhid);
GenerateParallelComp(TrHalo,INDXTYPE,nexthid);
GenerateParallelComp(UpLink,INDXTYPE,mbp);
*/


int nid, myid;


IDTYPE maxhid,tmaxhid,minhid,tminhid;


#define READY 0
#define WRITING 1
#define NP_TAG 11
#define R_TAG 456
int ready = READY,writing=WRITING;
IDTYPE finish=-99999;
int nstatus;
int motherrank=0;





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

IDTYPE nnowmbplist,nnextmbplist,nnextsubhalo,nuplink;
IDTYPE tnnowmbplist, tnnextmbplist, tnuplink;

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
int sortNowhidTrHalo(const void *a, const void *b){
	TrHalo *aa, *bb;
	aa = (TrHalo*)a;
	bb = (TrHalo*)b;
	if(aa->nowhid < bb->nowhid) return -1;
	else if(aa->nowhid > bb->nowhid) return 1;
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
	if(myid==0) {
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
	MPI_Bcast(&boxsize,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&hubble,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&npower,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&omep,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&omepb,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&omeplam,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&bias,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nspace,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&amax,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&astep,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&anow,1,MPI_FLOAT,0,MPI_COMM_WORLD);
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

#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))

void *GetHalo2(int nstep, IDTYPE *mhalo){
	char infile1[190],infile2[190];
	FILE *fp1,*fp2;
	IDTYPE np,mp;
	TrHalo *a = NULL;
	int i,j,k, WGROUPSIZE=10;
	IDTYPE isend,iget=0;
	isend = iget;
	int src,tgt,itag=9;
	MPI_Status status;
	src = myid-1;
	tgt = myid+1;


	*mhalo = 0;

	sprintf(infile1,HaloLinkFile".%.5d",nstep);
	if(myid==0) printf("%s is being opened \n",infile1);
	if(myid !=0) {
		MPI_Recv(&iget,sizeof(IDTYPE),MPI_BYTE,src,itag,MPI_COMM_WORLD,&status);
	}
	if((fp1 = fopen(infile1,"r")) != NULL){

		fread(&redshift,sizeof(float),1,fp1);
		fread(&np,sizeof(IDTYPE),1,fp1);
		mp = (np-1)/nid+1;
		iget = mp* myid;

		if(mp+iget > np) mp = np - iget;

		a = (TrHalo*)malloc(sizeof(TrHalo)*mp);
		fseek(fp1,sizeof(TrHalo)*iget,SEEK_CUR);
		fread(a,sizeof(TrHalo),mp,fp1);

		fclose(fp1);
		*mhalo = mp;
		printf("P%d : %ld is read\n",myid,mp);
	}
	if(myid != nid-1) {
		isend = iget + mp;
		MPI_Send(&isend,sizeof(IDTYPE),MPI_BYTE,tgt,itag,MPI_COMM_WORLD);
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
		printf("mask is set for np = %ld\n", *mhalo);
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
		while(uplink[j].nexthid < nextmbplist[i].nowhid && j < nuplink) j++;
		if(nextmbplist[i].mbp == 53693757289L){
			printf("here i am\n");
		}
		IDTYPE  nexthid = nextmbplist[i].nowhid;
		int jlink = 0;
		while(nexthid == uplink[j+jlink].nexthid){
			jlink ++;
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
			i= ilink;
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


void debug(TrHalo *next, IDTYPE np){
	return;
	IDTYPE i;
	for(i=0;i<np;i++){
		if(next[i].mbp == 30478370506L){
			printf("P%d has mbp = %ld at %d\n",myid,next[i].mbp,i);
		}
	}
	
}


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




	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);

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


	if(myid==0) fp5 = fopen(infile5,"r");
	GetHead(infile5,&fp5);
	if(myid==0) fclose(fp5);
	redshift1 = amax/anow-1;

	ny = nz = nx;

	pmass = 2.7755e11*omep*pow((double)boxsize/(double)nx,3.L);



	INDXTYPE nsize = ((long)nx*(long)ny*(long)nz/BitsOfChar-1L) + 1L;
	char *mask = (char *)calloc(nsize,sizeof(char));
	
	if(myid==0) fp2 = fopen(infile2,"r");
	if(myid==0) GetMask(nstep1,mask,&nnowmbplist);
	MPI_Bcast(&nnowmbplist,sizeof(INDXTYPE),MPI_BYTE,0,MPI_COMM_WORLD);
	{
		if(myid==0) printf("now before broadcasting\n");fflush(stdout);
		IDTYPE nchunk = INT_MAX;
		for(i=0;i<nsize;i+=nchunk) {
			nchunk = MIN(nchunk,nsize-i);
			MPI_Bcast(mask+i,nchunk,MPI_CHAR,0,MPI_COMM_WORLD);
			if(myid==0) printf(" ..");fflush(stdout);
		}
		if(myid==0) printf("now completing broadcasting\n");fflush(stdout);
	}
	if(0){
		int ioff;
		int icount = 0;
/*
		for(i=0;i<nsize;i++){
			for(ioff=0;ioff<8*sizeof(char);ioff++) if(WasMB(mask[i],ioff)) icount++;
		}
*/
		for(i=0;i<nsize;i++){
			if(mask[i] != 0) {
				icount ++;
			}
		}
		printf("P%d has %d mask out pixel\n",myid,icount);fflush(stdout);
/*
		MPI_Finalize();exit(0);
*/
	}




	if(myid==0) fp3 = fopen(infile3,"r");
	GetHead(infile3,&fp3);
	redshift2 = amax/anow -1;
	if(myid==0) fp4 = fopen(infile4,"r");

	double red2p1 = redshift2+1;

	rscale = boxsize/nx;
	double Hsub = sqrt(omep*red2p1*red2p1*red2p1 + omeplam +(1-omep-omeplam)*red2p1*red2p1);
	vscale = 100.L*boxsize/hubble/amax*anow*anow*hubble*Hsub; /* In unit of km/second */


	simpar.mx = simpar.nx = simpar.ny = simpar.nz = nx;
	simpar.nspace = 1;
	simpar.mxmy  = simpar.nx*simpar.ny;


	if(myid==0) printf("Now making the linking \n");fflush(stdout);


	IDTYPE nmaxbp = 10000000;
	basicparticletype *bp ;
	if(myid==0) bp= (basicparticletype*)malloc(nmaxbp*sizeof(basicparticletype));

	IDTYPE iuplink=0;
	IDTYPE nexthid;
	UpLink *uplink;

	if(myid ==motherrank) { 
		MPI_Status mstatus,cstatus;
		nexthid=0;
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

			do {
				MPI_Probe(MPI_ANY_SOURCE, READY,MPI_COMM_WORLD,&mstatus);
				int src = mstatus.MPI_SOURCE;
				int dest = mstatus. MPI_SOURCE;
				MPI_Recv(&ready,1,MPI_INT,src,READY,MPI_COMM_WORLD,&cstatus);
				if(ready == READY){
					MPI_Send(&nexthid,sizeof(IDTYPE),MPI_BYTE,dest,NP_TAG,MPI_COMM_WORLD);
					MPI_Send(&halo,sizeof(Halo),MPI_BYTE,dest,R_TAG,MPI_COMM_WORLD);
					MPI_Send(&haloq,sizeof(HaloQ),MPI_BYTE,dest,R_TAG,MPI_COMM_WORLD);
					MPI_Send(bp,halo.np*sizeof(basicparticletype),MPI_BYTE,dest,R_TAG,MPI_COMM_WORLD);
				}
			}while(0);  /* while(ready != READY);*/
			nexthid ++;
		}
		for(i=1;i<nid;i++){
			MPI_Probe(MPI_ANY_SOURCE,READY,MPI_COMM_WORLD,&mstatus);
			MPI_Recv(&ready,1,MPI_INT,mstatus.MPI_SOURCE,READY,MPI_COMM_WORLD,&cstatus);
			MPI_Send(&finish,sizeof(IDTYPE),MPI_BYTE,mstatus.MPI_SOURCE,NP_TAG,MPI_COMM_WORLD);
		}
		printf("P%d is now finishing the distribution of halo data \n",myid);
		fclose(fp2);fclose(fp4);fclose(fp3);
		free(bp);
		/* Now gathering RV file for save */
		char out2[100];
		sprintf(out2,HaloRVFile".%.5d",nstep2);
		FILE *wp2;

		wp2 = fopen(out2,"w");
		for(i=0;i<nid;i++){
			if(i==motherrank) continue;
			IDTYPE nlocal;
			MPI_Recv(&nlocal,sizeof(IDTYPE),MPI_BYTE,i,R_TAG,MPI_COMM_WORLD,&cstatus);
			mbprv = (MbpRV*)malloc(sizeof(MbpRV)*nlocal);
			int nchunk = INT_MAX;
			nsize = nlocal*sizeof(MbpRV);
			for(j=0;j<nsize;j+=nchunk){
				if(nsize-j < nchunk) nchunk = nsize-j;
				MPI_Recv(((char*)mbprv)+j,nchunk,MPI_BYTE,i,R_TAG,MPI_COMM_WORLD,&cstatus);
			}
			fwrite(mbprv,sizeof(MbpRV),nlocal,wp2);
			printf("getting rvfile from P%d\n",i);
			free(mbprv);
		}
		fclose(wp2);

		nuplink = nnextmbplist = 0;
		nextmbplist = (TrHalo*)malloc(sizeof(TrHalo));

		uplink = (UpLink*)malloc(sizeof(UpLink));

		/* Now gathering nextmbplist file for anlaysis */
		/*
		nextmbplist = (TrHalo*)malloc(sizeof(TrHalo));
		nnextmbplist = 0;
		for(i=0;i<nid;i++){
			if(i==motherrank) continue;
			IDTYPE nlocal;
			MPI_Recv(&nlocal,sizeof(IDTYPE),MPI_BYTE,i,R_TAG,MPI_COMM_WORLD,&cstatus);
			nextmbplist = (TrHalo*)realloc(nextmbplist,sizeof(TrHalo)*(nnextmbplist+nlocal));
			MPI_Recv(nextmbplist+nnextmbplist,nlocal*sizeof(TrHalo),MPI_BYTE,i,R_TAG,MPI_COMM_WORLD,&cstatus);
			nnextmbplist += nlocal;
		}
		for(i=0;i<nnextmbplist;i++) nextmbplist[i].aid = i;

		uplink = (UpLink*)malloc(sizeof(UpLink));
		iuplink = 0;
		for(i=0;i<nid;i++){
			if(i==motherrank) continue;
			IDTYPE nlocal;
			MPI_Recv(&nlocal,sizeof(IDTYPE),MPI_BYTE,i,R_TAG,MPI_COMM_WORLD,&cstatus);
			uplink = (UpLink*)reall(uplink,sizeof(UpLink)*(iuplink+nlocal));
			MPI_Recv(uplink+iuplink,nlocal*sizeof(UpLink),MPI_BYTE,i,R_TAG,MPI_COMM_WORLD,&cstatus);
			iuplink += nlocal;
		}
		*/

	}
	else {
		Halo halo; 
		HaloQ haloq; 
		MPI_Status cstatus,mstatus;
		IDTYPE nmaxnextmbplist = 1000000 + nnowmbplist/(nid-0);
		IDTYPE nmaxnuplink = 1000000 + nnowmbplist/(nid-1);
		uplink = (UpLink*) malloc(sizeof(UpLink)*nmaxnuplink);
		for(i=0;i<nmaxnextmbplist;i++) RESET_WHOLE_FLAGS_MP(uplink,i,linkflag);
		nextmbplist = (TrHalo*)malloc(sizeof(TrHalo)*nmaxnextmbplist);
		mbprv = (MbpRV*)malloc(sizeof(MbpRV)*nmaxnextmbplist);
		for(i=0;i<nmaxnextmbplist;i++) {
			RESET_WHOLE_FLAGS_MP(nextmbplist,i,statusflag);
			RESET_WHOLE_FLAGS_MP(nextmbplist,i,linkflag);
		}
		iuplink = 0;
		nnextmbplist = 0;
		ready = READY;
		MPI_Send(&ready,1,MPI_INT,motherrank,READY,MPI_COMM_WORLD);
		MPI_Recv(&nexthid,sizeof(IDTYPE),MPI_BYTE,motherrank,NP_TAG,MPI_COMM_WORLD,&cstatus);
		while(nexthid != finish){
			MPI_Recv(&halo,sizeof(Halo),MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
			MPI_Recv(&haloq,sizeof(HaloQ),MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
			if(nmaxbp < halo.np ) {
				nmaxbp += halo.np;
				bp = (basicparticletype*)realloc(bp,sizeof(basicparticletype)*nmaxbp);
			}
			MPI_Recv(bp,halo.np*sizeof(basicparticletype),MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
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
/*
					nextmbplist[nnextmbplist].aid = aid++;
*/
					nextmbplist[nnextmbplist].nowhid = nexthid;
					mbprv[nnextmbplist].x = XofP(bp+i)*rscale;
					mbprv[nnextmbplist].y = YofP(bp+i)*rscale;
					mbprv[nnextmbplist].z = ZofP(bp+i)*rscale;
					mbprv[nnextmbplist].vx = bp[i].vx*vscale;
					mbprv[nnextmbplist].vy = bp[i].vy*vscale;
					mbprv[nnextmbplist].vz = bp[i].vz*vscale;
					nextmbplist[nnextmbplist].mass = haloq.mass;
					nextmbplist[nnextmbplist].mbp = bp[i].indx;
					nextmbplist[nnextmbplist].downaid = TERMINATE;
					nnextmbplist ++;
					nsub ++;
					iuplink ++;
				}
			}
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
/*
					nextmbplist[nnextmbplist].aid = aid++;
*/
					nextmbplist[nnextmbplist].nowhid = nexthid;
					mbprv[nnextmbplist].x = XofP(bp+ii)*rscale;
					mbprv[nnextmbplist].y = YofP(bp+ii)*rscale;
					mbprv[nnextmbplist].z = ZofP(bp+ii)*rscale;
					mbprv[nnextmbplist].vx = bp[ii].vx*vscale;
					mbprv[nnextmbplist].vy = bp[ii].vy*vscale;
					mbprv[nnextmbplist].vz = bp[ii].vz*vscale;
					nextmbplist[nnextmbplist].mass = haloq.mass;
					nextmbplist[nnextmbplist].mbp = haloq.mbp;
					nextmbplist[nnextmbplist].downaid = TERMINATE;
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
			if(nmaxnuplink - iuplink < 10000000){
				nmaxnuplink += 10000000;
				uplink = (UpLink*)realloc(uplink,sizeof(TrHalo)*nmaxnuplink);
				for(i=iuplink;i<nmaxnuplink;i++) RESET_WHOLE_FLAGS_MP(uplink,i,linkflag);
			}
			if(nnextmbplist%100000 ==0) printf("P%d has now passed %ld'th \n",myid,nnextmbplist);
			free(bp);
			ready = READY;
			MPI_Send(&ready,1,MPI_INT,motherrank,READY,MPI_COMM_WORLD);
			MPI_Recv(&nexthid,sizeof(IDTYPE),MPI_BYTE,motherrank,NP_TAG,MPI_COMM_WORLD,&cstatus);
		}
		MPI_Send(&nnextmbplist,sizeof(IDTYPE),MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD);
		int nchunk = INT_MAX;
		nsize = sizeof(MbpRV)*nnextmbplist;
		for(i=0;i<nsize;i+= nchunk){
			if(nsize -i < nchunk) nchunk = nsize-i;
			MPI_Send(((char*)mbprv)+i,nchunk,MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD);
		}
		free(mbprv);
		/*
		MPI_Send(&nnextmbplist,sizeof(IDTYPE),MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD);
		MPI_Send(nextmbplist,nnextmbplist*sizeof(TrHalo),MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD);
		MPI_Send(&iuplink,sizeof(IDTYPE),MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD);
		MPI_Send(uplink,iuplink*sizeof(UpLink),MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD);
		
		free(nextmbplist);free(uplink);
		*/
		IDTYPE ncumul,mcumul,icumul,jcumul;
		mcumul = jcumul = 0;
		int dest = (myid+1+nid)%nid;
		int src = (myid-1+nid)%nid;
		if(myid!= 1){
			MPI_Recv(&mcumul,sizeof(IDTYPE),MPI_BYTE,src,1,MPI_COMM_WORLD,&cstatus);
			MPI_Recv(&jcumul,sizeof(IDTYPE),MPI_BYTE,src,0,MPI_COMM_WORLD,&cstatus);
		}
		ncumul = mcumul + nnextmbplist;
		icumul = jcumul + iuplink;
		if(myid != nid-1){
			MPI_Send(&ncumul,sizeof(IDTYPE),MPI_BYTE,dest,1,MPI_COMM_WORLD);
			MPI_Send(&icumul,sizeof(IDTYPE),MPI_BYTE,dest,0,MPI_COMM_WORLD);
		}
		for(i=0;i<nnextmbplist;i++){
			nextmbplist[i].aid = i + mcumul;
		}
		for(i=0;i<iuplink;i++){
			uplink[i].aid = i + jcumul;
		}
		printf("P%d has nnextmbplist/nuplink %ld %ld\n",myid,nnextmbplist, iuplink);
		nuplink = iuplink;
		/*
		tnnextmbplist = mcumul + nnextmbplist;
		tnuplink = jcumul + iuplink;
		MPI_Bcast(&tnnextmbplist, sizeof(IDTYPE),MPI_BYTE, nid-1,MPI_COMM_WORLD);
		MPI_Bcast(&tnuplink, sizeof(IDTYPE),MPI_BYTE, nid-1,MPI_COMM_WORLD);
		*/


	}
	free(mask);
	printf("Well made the link \n");fflush(stdout);


	debug(nextmbplist,nnextmbplist);



	char out2[190];
	FILE *wp;
	/* reading nowmbplist */

	nowmbplist = GetHalo2(nstep1,&nnowmbplist);

	if(myid==0){
		if(nuplink != nnowmbplist){
			fprintf(stderr,"missing mbp between nowmbplist and uplink %ld\n",nnowmbplist-nuplink);
		}
	}



	{
		if(myid==0) printf("in the massive sorting\n");
		/* dump to nowmbplist link flag */
		/*
		qsort(nowmbplist,nnowmbplist,sizeof(TrHalo),sortMbTrHalo);
		qsort(uplink,nuplink,sizeof(UpLink),sortMbUpLink);
		*/
		if(1){ 
			INDXTYPE min,max;
			min  = 0;
			max = (INDXTYPE)simpar.nx*(INDXTYPE)simpar.ny*(INDXTYPE)simpar.nz;
			MyConstRangePqsort(nowmbplist,nnowmbplist,sizeof(TrHalo),TrHalo,INDXTYPE,mbp,min,max,MPI_COMM_WORLD);
			MyConstRangePqsort(uplink,nuplink,sizeof(UpLink),UpLink, INDXTYPE, mbp,min,max,MPI_COMM_WORLD);
		}
		else {
			MyPqsort(nowmbplist,nnowmbplist,sizeof(TrHalo),TrHalo,INDXTYPE,mbp,MPI_COMM_WORLD);
			MyPqsort(uplink,nuplink,sizeof(UpLink),UpLink, INDXTYPE, mbp,MPI_COMM_WORLD);
		}

		debug(nextmbplist,nnextmbplist);
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
	}
	debug(nextmbplist,nnextmbplist);
	
	/* Find the major merger line */
	if(1){
		IDTYPE min,max;
		min = INT_MAX;
		max = INT_MIN;
		maxhid = INT_MIN;
		minhid = INT_MAX;
		for(i=0;i<nnextmbplist;i++){
			maxhid = MAX(maxhid, nextmbplist[i].nowhid);
			minhid = MIN(minhid, nextmbplist[i].nowhid);
		}
		for(i=0;i<nuplink;i++){
			maxhid = MAX(maxhid, uplink[i].nexthid);
			minhid = MIN(minhid, uplink[i].nexthid);
		}
		if(sizeof(IDTYPE) != sizeof(long long)) {
			fprintf(stderr,"Error in matching the size of IDTYPE and long long\n");
			MPI_Finalize();exit(999);
		}
		MPI_Reduce(&maxhid,&tmaxhid, 1, MPI_LONG_LONG, MPI_MAX,0, MPI_COMM_WORLD);
		maxhid = tmaxhid;
		MPI_Bcast(&maxhid,1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
		max = maxhid+1;

		MPI_Reduce(&minhid,&tminhid, 1, MPI_LONG_LONG, MPI_MIN,0, MPI_COMM_WORLD);
		minhid = tminhid;
		MPI_Bcast(&minhid,1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
		min = minhid;

		MyConstRangePqsort(uplink,nuplink, sizeof(UpLink), UpLink, IDTYPE,nexthid, min,max,MPI_COMM_WORLD);
		MyConstRangePqsort(nextmbplist,nnextmbplist, sizeof(TrHalo), TrHalo, IDTYPE,nowhid, min,max,MPI_COMM_WORLD);

	}
	else{
		MyPqsort(uplink,nuplink, sizeof(UpLink), UpLink, IDTYPE,nexthid, MPI_COMM_WORLD);
		MyPqsort(nextmbplist,nnextmbplist, sizeof(TrHalo), TrHalo, IDTYPE,nowhid, MPI_COMM_WORLD);
	}
	debug(nextmbplist,nnextmbplist);
	printf("in identifying the major satellite \n");
	FindMajorMergermbp(uplink,nuplink,nextmbplist,nnextmbplist);
	debug(nextmbplist,nnextmbplist);

	if(1){ 
		INDXTYPE min,max;
		min = 0;
		max = (INDXTYPE)nx*(INDXTYPE)ny*(INDXTYPE)nz;
		MyConstRangePqsort(uplink,nuplink,sizeof(UpLink),UpLink, INDXTYPE, mbp,min,max,MPI_COMM_WORLD);
	}
	else {
		MyPqsort(uplink,nuplink,sizeof(UpLink),UpLink, INDXTYPE, mbp,MPI_COMM_WORLD);
	}
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

	IDTYPE tnnowmbplist,tnnextmbplist;
	MPI_Reduce(&nnowmbplist,&tnnowmbplist,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&nnextmbplist,&tnnextmbplist,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&tnnowmbplist,1, MPI_LONG_LONG, 0,MPI_COMM_WORLD);
	MPI_Bcast(&tnnextmbplist,1, MPI_LONG_LONG, 0,MPI_COMM_WORLD);

	{ 
		MyPqsort(nextmbplist,nnextmbplist,sizeof(TrHalo),TrHalo, IDTYPE, aid,MPI_COMM_WORLD);
		MyPqsort(nowmbplist,nnowmbplist,sizeof(TrHalo),TrHalo, IDTYPE, aid,MPI_COMM_WORLD);
	}

#define AHaloLinkFile "aaaaaaaa"

	if(myid==0) printf("in saving data 1 \n");
	for(i=0;i<nid;i++){
		if(nnowmbplist >0 && myid == i){
			sprintf(infile1,HaloLinkFile".%.5d",nstep1);
			FILE *wp2;
			if(myid==0) wp2 = fopen(infile1,"w");
			else wp2 = fopen(infile1,"a");
			if(myid==0){
				fwrite(&redshift1,sizeof(float),1,wp2);
				fwrite(&tnnowmbplist,sizeof(IDTYPE),1,wp2);
			}
			fwrite(nowmbplist,sizeof(TrHalo),nnowmbplist,wp2);
			fclose(wp2);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}


	debug(nextmbplist,nnextmbplist);
	if(myid==0) printf("in saving data 2 \n");


	for(i=0;i<nid;i++){
		if(nnextmbplist >0 && myid == i){
			sprintf(out2,HaloLinkFile".%.5d",nstep2);
			FILE *wp2;
			if(myid==0) wp2 = fopen(out2,"w");
			else wp2 = fopen(out2,"a");
			if(myid==0) {
				fwrite(&redshift2,sizeof(float),1,wp2);
				fwrite(&tnnextmbplist,sizeof(IDTYPE),1,wp2);
			}
			fwrite(nextmbplist,sizeof(TrHalo),nnextmbplist,wp2);
			fclose(wp2);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}


	if(myid==0) printf("saved to %s with %ld mbp's\n",out2,nnextmbplist);
	debug(nextmbplist,nnextmbplist);
	MPI_Finalize();
	return 0;
}
