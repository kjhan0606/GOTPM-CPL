#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<sys/times.h>
#include<mpi.h>
#include "pmheader.h"
#include "Memory.h"
#define EQNBox
#include "slicemigrate.h"
#undef EQNBox



int nprocsflag=1;
CommPair *commpair;
int *nstridenum;
int *tmpnstridenum;
/*
CommPair commpair[NPROCS];
int nstridenum[NPROCS];
int tmpnstridenum[NPROCS];
*/
NSlice *nfftwslice,*ndenslice,*npotslice;
static int llocal_z_start,llocal_nz;
static int nx,ny,nz;
static int myid,nid;
long numpixelslice;

void DetermineProcid(){
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
}
void DetermineStrideComm(int istride){
	int i,j,k;
	int dest,src,jj,iflag;
	/* Initialize the communication table */

	for(j=0;j<nid;j++) commpair[j].up.flag = commpair[j].down.flag 
		= None;
	for(j=0;j<nid;j++){
		if(commpair[j].up.flag != None) continue;
		commpair[j].down.flag = Second;
		jj = j;
		dest = (jj + istride + nid) % nid;
		iflag = 0;
		while(commpair[dest].down.flag == None){
			commpair[jj].up.flag = commpair[dest].down.flag = 
				iflag % 2;
			commpair[jj].up.dest = dest;
			commpair[dest].down.dest = jj;
			jj = dest;
			dest = (jj + istride + nid) % nid;
			iflag ++;
		}
		if(commpair[dest].down.flag != ((iflag) %2)){
			commpair[dest].down.flag = commpair[jj].up.flag = Last;
		}
		else commpair[jj].up.flag = commpair[dest].down.flag;
		commpair[jj].up.dest = dest;
		commpair[dest].down.dest = jj;
	}
}
void determine_comm(int nid, CommPair *commpair){
	int istride,j,k;
	for(istride=1;istride<nid;istride++){
		/* Initialize the communication table */
		printf("################################\n");
		for(j=0;j<nid;j++) commpair[j].up.flag = commpair[j].down.flag 
			= None;
		if(myid == 0 ){
			int dest,src,jj,iflag;
			for(j=0;j<nid;j++){
				if(commpair[j].up.flag != None) continue;
				commpair[j].down.flag = Second;
				printf("%d ",j);
				jj = j;
				dest = (jj + istride + nid) % nid;
				iflag = 0;
				while(commpair[dest].down.flag == None){
					commpair[jj].up.flag = commpair[dest].down.flag = 
						iflag % 2;
					commpair[jj].up.dest = dest;
					commpair[dest].down.dest = jj;
					printf(" ---> %d",dest);
					jj = dest;
					dest = (jj + istride + nid) % nid;
					iflag ++;
				}
				if(commpair[dest].down.flag != ((iflag) %2)){
					commpair[dest].down.flag = commpair[jj].up.flag = Last;
					printf("Last");
				}
				else commpair[jj].up.flag = commpair[dest].down.flag;
				commpair[jj].up.dest = dest;
				commpair[dest].down.dest = jj;
				printf("\n");
			}
		}
	}
}
/* */
/*
int nowslices[NPROCS][NNZ];
int afterslices[NPROCS][NNZ];
*/
int *nowslices,*afterslices;
#define nowslices(i,j) (nowslices[(i)*nz+(j)])
#define afterslices(i,j) (afterslices[(i)*nz+(j)])
/* this routine checks the range of slices in current and in after mesh 
 * Multiple processors can have slice ranges that overlap each other.
 * */
void InitNslices(int nowzstart,int nowzfinal,int afterzstart,int afterzfinal){
	int i,j,k;
	MPI_Status mpistatus;
	nowslices = (int*)Malloc(sizeof(int)*nid*nz,PPTR(nowslices));
	afterslices = (int*)Malloc(sizeof(int)*nid*nz,PPTR(afterslices));
	for(i=0;i<nz;i++){
		nowslices(myid,i) = None;
		afterslices(myid,i) = None;
	}
	k = 0;
	for(i=nowzstart;i<=nowzfinal;i++){
		j = (i+nz)% nz;
		nowslices(myid,j) = k;
		k++;
	}
	k = 0;
	for(i=afterzstart;i<=afterzfinal;i++){
		j = (i+nz)% nz;
		afterslices(myid,j) = k;
		k++;
	}
	if(myid==0) printf("P%d made mapping\n",myid);

#ifdef MPI_SLOW
	if(myid ==0){
		for(i=1;i<nid;i++){
			MPI_Recv(&nowslices(i,0),nz,MPI_INT,i,0,MPI_COMM_WORLD,&mpistatus);
		}
	}
	else MPI_Send(&nowslices(myid,0),nz,MPI_INT,0,0,MPI_COMM_WORLD);
#else
	if(myid ==0){
		MPI_Status rstatus;
		for(i=1;i<nid;i++){
			MPI_Probe(MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&rstatus);
			MPI_Recv(&nowslices(rstatus.MPI_SOURCE,0),nz,MPI_INT,rstatus.MPI_SOURCE,0,MPI_COMM_WORLD,&mpistatus);
		}
	}
	else MPI_Send(&nowslices(myid,0),nz,MPI_INT,0,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(myid==0) printf("P%d before broadcasting ",myid);
	MPI_Bcast(nowslices,nz*nid,MPI_INT,0,MPI_COMM_WORLD);
	if(myid==0) printf("after broadcasting \n");
	/*
	for(i=0;i<nid;i++){
		MPI_Bcast(&nowslices(i,0),nz,MPI_INT,0,MPI_COMM_WORLD);
	}
	*/
#ifdef MPI_SLOW
	if(myid ==0){
		for(i=1;i<nid;i++){
			MPI_Recv(&afterslices(i,0),nz,MPI_INT,i,0,MPI_COMM_WORLD,
					&mpistatus);
		}
	}
	else MPI_Send(&afterslices(myid,0),nz,MPI_INT,0,0,MPI_COMM_WORLD);
#else
	if(myid ==0){
		MPI_Status rstatus;
		for(i=1;i<nid;i++){
			MPI_Probe(MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&rstatus);
			MPI_Recv(&afterslices(rstatus.MPI_SOURCE,0),nz,MPI_INT,rstatus.MPI_SOURCE,0,MPI_COMM_WORLD,&mpistatus);
		}
	}
	else MPI_Send(&afterslices(myid,0),nz,MPI_INT,0,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	if(myid==0) printf("P%d before broadcasting ",myid);
	MPI_Bcast(afterslices,nid*nz,MPI_INT,0,MPI_COMM_WORLD);
	if(myid==0) printf("after broadcasting \n");
	/*
	for(i=0;i<nid;i++){
		MPI_Bcast(&afterslices(i,0),nz,MPI_INT,0,MPI_COMM_WORLD);
	}
	*/
}
void FinalizeSlices(){
	Free(afterslices);
	Free(nowslices);
}
/* Copy noewslices to afterslices if they are in the same processors 
 * and return meaningful nstridenum value
 * if interprocessor communication is needed.
 * */
void CheckAndCopyStride(float *nowslice,float *afterslice,
		int nowzstart,int nowzfinal,int afterzstart,int afterzfinal){
	int i,j,k;
	long ii;
	int nowstart,afterstart;
	int otherid;
	float *now,*after;
	DetermineProcid();
	for(i=0;i<nid;i++) nstridenum[i] = 0;
	for(i=0;i<nz;i++){
		if(nowslices(myid,i)!= None){
			for(otherid=0;otherid<nid;otherid++){
				if(afterslices(otherid,i) !=  None){
					k = (otherid-myid+nid)%nid;
					nstridenum[k]++;
					if(myid==otherid) {
						now = nowslice+(nowslices(myid,i))*numpixelslice;
						after = afterslice+(afterslices(myid,i))*numpixelslice;
						for(ii=0;ii<numpixelslice;ii++){
							*after = *now;now++;after++;
						}
					}
				}
			}
		}
	}
	(void) MPI_Reduce(nstridenum,tmpnstridenum,nid,MPI_INT,MPI_SUM,
					  0,MPI_COMM_WORLD);
	if(myid==0) {
		for(i=0;i<nid;i++) {
			nstridenum[i] = tmpnstridenum[i];
		}
	}
	MPI_Bcast(nstridenum,nid,MPI_INT,0,MPI_COMM_WORLD);
}
void SliceSend(int dest, float *den){
	int i,j,k,iflag;
	long nowid,afterid;
	MPI_Status status;
	MPI_Request request;
	for(i=0;i<nz;i++){
		if((nowid=nowslices(myid,i))!=None &&
				(afterid=afterslices(dest,i))!=None){
			MPI_Send(den+nowid*numpixelslice,numpixelslice,MPI_FLOAT,dest,0,
					MPI_COMM_WORLD);
			MPI_Recv(&iflag,1,MPI_INT,dest,1,MPI_COMM_WORLD,&status);
			/* This is a gabage process to confirm the synchronized communications */
			if(iflag < 1){
				if(den[0]<-1.E20) fprintf(stderr,"ERRORR in networking %d\n",iflag);
				break;
			}
		}
	}
}
void SliceRecvAdd(int src, float *fftwmesh, float *tmpslice){
	long ii;
	int i,iflag=1;
	long nowid,afterid;
	MPI_Status status;
	MPI_Request request;
	for(i=0;i<nz;i++){
		if((nowid=nowslices(src,i))!=None &&
				(afterid=afterslices(myid,i))!=None){
			MPI_Recv(tmpslice,numpixelslice,MPI_FLOAT,src,0,
					MPI_COMM_WORLD,&status);
			for(ii=0;ii<numpixelslice;ii++){
				fftwmesh[afterid*numpixelslice+ii] += tmpslice[ii];
			}
			MPI_Send(&iflag,1,MPI_INT,src,1,MPI_COMM_WORLD);
			iflag ++;
		}
	}
}
void SliceRecv(int src, float *den){
	int i,j,k,iflag=1;
	long nowid,afterid;
	MPI_Status status;
	MPI_Request request;
	for(i=0;i<nz;i++){
		if((nowid=nowslices(src,i))!=None &&
				(afterid=afterslices(myid,i))!=None){
			MPI_Irecv(den+afterid*numpixelslice,numpixelslice,MPI_FLOAT,src,0,
					MPI_COMM_WORLD,&request);
			MPI_Wait(&request,&status);
			MPI_Send(&iflag,1,MPI_INT,src,1,MPI_COMM_WORLD);
			iflag ++;
		}
	}
}
void RearrangePot(int Nx,int Ny, int Nz, int zwidth,float *den){
	long i,j,k;
	long numpixelslice,numtmp;
	float *tmpden,*a,*b;

	numpixelslice = Ny*(2*(Nx/2+1));
	tmpden = (float *) Malloc(sizeof(float)*3*numpixelslice,PPTR(tmpden));
	a = tmpden;
	b = den;

	for(i=0;i<3*numpixelslice;i++){
		*a = *b; a++;b++;
	}
	numtmp = (zwidth-3)*numpixelslice;
	a = den;
	b = den+3*numpixelslice;
	for(i=0;i<numtmp;i++){
		*a = *b ; a++;b++;
	}
	b = tmpden;
	for(i=0;i<3*numpixelslice;i++){
		*a = *b; a++;b++;
	}
	Free(tmpden);
}
void tsc_fftw_mesh(float *den, int Zstart,int Zheight,
	float *fftwmesh, int LLocal_z_start, int LLocal_nz, int Nx,int Ny,int Nz,
	int mode){
	int i,j,k,nzstart,nzheight;
	int dest,src;
	int nowzstart,nowzfinal,afterzstart,afterzfinal;
	float *tmpslice;
	int zstart,zheight;

	DetermineProcid();
	/* From argument to Local Global Viriables */
	llocal_z_start=LLocal_z_start; llocal_nz=LLocal_nz; zstart=Zstart; 
	zheight=Zheight; 
	nowzstart = zstart; 
	nowzfinal = zstart + zheight-1; 
								
	afterzstart = llocal_z_start;
	afterzfinal = llocal_z_start + llocal_nz-1;
	nx = Nx; ny = Ny; nz = Nz;
	numpixelslice = 2*(nx/2+1)*ny;



	if(nprocsflag){
		nstridenum = (int *) Malloc(sizeof(int)*nid,PPTR(nstridenum));
		tmpnstridenum = (int *) Malloc(sizeof(int)*nid,PPTR(tmpnstridenum));
		commpair = (CommPair *)Malloc(sizeof(CommPair)*nid,PPTR(commpair));
		nprocsflag = 0;
	}

	if(myid==0) printf("P%d has zslices from %d to %d <---> %d %d (tsc-->fftw)\n",
			myid,nowzstart,nowzfinal,afterzstart,afterzfinal);
	InitNslices(nowzstart,nowzfinal,afterzstart,afterzfinal);
	CheckAndCopyStride(den,fftwmesh,nowzstart,nowzfinal,afterzstart,
			afterzfinal);
	tmpslice = (float *) Malloc(sizeof(float)*numpixelslice,PPTR(tmpslice));
	fflush(stdout);fflush(stderr);
	MPI_Barrier(MPI_COMM_WORLD);

	for(i=1;i<nid;i++){
		if(nstridenum[i] > 0){
			DetermineStrideComm(i);
			if(myid==nid-1) printf("Stride #=%d\n",i);fflush(stdout);
			dest = commpair[myid].up.dest;
			src = commpair[myid].down.dest;
			if(commpair[myid].up.flag == First){
				SliceSend(dest,den);
			}
			else if(commpair[myid].down.flag == First){
				SliceRecvAdd(src,fftwmesh,tmpslice);
			}
			if(commpair[myid].up.flag == Second){
				SliceSend(dest,den);
			}
			else if(commpair[myid].down.flag == Second){
				SliceRecvAdd(src,fftwmesh,tmpslice);
			}
			if(commpair[myid].up.flag == Last){
				SliceSend(dest,den);
			}
			else if(commpair[myid].down.flag == Last){
				SliceRecvAdd(src,fftwmesh,tmpslice);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	/*
	if(myid==0){
		for(i=0;i<nid;i++){
			if(nstridenum[i] != 0) printf("nstrid[%d] = %d\n",i,nstridenum[i]);
		}
	}

	*/
	Free(tmpslice);
/*
	Free(commpair);
	Free(tmpnstridenum);
	Free(nstridenum);
	{
		int ii;
		double sum;
		for(ii=0;ii<nid;ii++){
			if(ii == myid){
				for(k = 0;k<local_nz;k++){
					sum = 0;
				   for(i=0;i<2*(nx/2+1)*ny;i++) {
				     	sum += fftwmesh[i+k*2*(nx/2+1)*ny];
				   }
			 	   printf("P%d has density sum %f in  %d z\n",myid,sum,k);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

	}
	*/

	/*
	{
		int ii;
		double sum,tsum;
		sum = 0.;
	   for(i=0;i<2*(nx/2+1)*ny*local_nz;i++) {
	     	sum += fftwmesh[i];
	   }
 	   printf("P%d has local total density sum %f : %d \n",myid,sum);
	   MPI_Reduce(&sum,&tsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	   if(myid==0) printf("Total density sum is %f\n",tsum);
	   exit(99);
	}
	*/
	FinalizeSlices();
	if(mode ==0) for(i=0;i<2*(nx/2+1)*ny*llocal_nz;i++) fftwmesh[i] = fftwmesh[i]-1.;
	if(myid==0) printf("Finish the tsc_fftw_mesh \n");
}
void fftw_fda4_mesh(float *den, int Zstart,int Zheight,
	float *fftwmesh, int LLocal_z_start, int LLocal_nz, int Nx,int Ny,int Nz,
	int rearrangeflag){
	int i,j,k,nzstart,nzheight;
	int dest,src;
	int nowzstart,nowzfinal,afterzstart,afterzfinal;
	float *tmpslice;
	int zstart,zheight;


	DetermineProcid();
	if(nprocsflag){
		nstridenum = (int *) Malloc(sizeof(int)*nid,PPTR(nstridenum));
		tmpnstridenum = (int *) Malloc(sizeof(int)*nid,PPTR(tmpnstridenum));
		commpair = (CommPair *)Malloc(sizeof(CommPair)*nid,PPTR(commpair));
		nprocsflag = 0;
	}

	/* From argument to Local Global Viriables */
	llocal_z_start=LLocal_z_start; llocal_nz=LLocal_nz; zstart=Zstart; 
	zheight=Zheight; 
	nowzstart = zstart;
	nowzfinal = zstart + zheight-1;
	afterzstart = llocal_z_start;
	afterzfinal = llocal_z_start + llocal_nz-1;
	nx = Nx; ny = Ny; nz = Nz;
	numpixelslice = ny*(2*(nx/2+1));



	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0) printf("P%d has now been passing initnslices\n",myid);
	InitNslices(nowzstart,nowzfinal,afterzstart,afterzfinal);
	if(myid==0) printf("P%d has now passed initnslices\n",myid);
	CheckAndCopyStride(den,fftwmesh,nowzstart,nowzfinal,afterzstart,
			afterzfinal);
	if(myid==0) printf("P%d has zslices from %d to %d <---> %d %d (tsc<-->fftw)\n",
			myid,nowzstart,nowzfinal,afterzstart,afterzfinal);
	fflush(stdout);fflush(stderr);
	MPI_Barrier(MPI_COMM_WORLD);
	for(i=1;i<nid;i++){
		if(nstridenum[i] > 0){
			DetermineStrideComm(i);
			if(myid==nid-1) printf("Stride #=%d\n",i);fflush(stdout);
			dest = commpair[myid].up.dest;
			src = commpair[myid].down.dest;
			if(commpair[myid].up.flag == First){
				SliceSend(dest,den);
			}
			else if(commpair[myid].down.flag == First){
				SliceRecv(src,fftwmesh);
			}
			if(commpair[myid].up.flag == Second){
				SliceSend(dest,den);
			}
			else if(commpair[myid].down.flag == Second){
				SliceRecv(src,fftwmesh);
			}
			if(commpair[myid].up.flag == Last){
				SliceSend(dest,den);
			}
			else if(commpair[myid].down.flag == Last){
				SliceRecv(src,fftwmesh);
			}
		}
	}
	/*
	MPI_Barrier(MPI_COMM_WORLD);
	*/
	if(myid==0) printf("P%d has passed\n",myid);
	if(rearrangeflag)RearrangePot(nx,ny,nz,llocal_nz,fftwmesh);
	FinalizeSlices();
}


#define NRANSI
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50
#define NSTACKp1 51
#define NR_END 1

/*
void indexx(unsigned long n, pmparticletype *p) {
	unsigned long i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0;
	int istack[NSTACKp1];
	float a;

	for (j=1;j<=n;j++) p[j].u.sort=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=p[j].u.sort;
				a=p[indxt].z;
				for (i=j-1;i>=1;i--) {
					if (p[p[i].u.sort].z <= a) break;
					p[i+1].u.sort=p[i].u.sort;
				}
				p[i+1].u.sort=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(p[k].u.sort,p[l+1].u.sort);
			if (p[p[l+1].u.sort].z > p[p[ir].u.sort].z) {
				SWAP(p[l+1].u.sort,p[ir].u.sort)
			}
			if (p[p[l].u.sort].z > p[p[ir].u.sort].z) {
				SWAP(p[l].u.sort,p[ir].u.sort)
			}
			if (p[p[l+1].u.sort].z > p[p[l].u.sort].z) {
				SWAP(p[l+1].u.sort,p[l].u.sort)
			}
			i=l+1;
			j=ir;
			indxt=p[l].u.sort;
			a=p[indxt].z;
			for (;;) {
				do i++; while (p[p[i].u.sort].z < a);
				do j--; while (p[p[j].u.sort].z > a);
				if (j < i) break;
				SWAP(p[i].u.sort,p[j].u.sort)
			}
			p[l].u.sort=p[j].u.sort;
			p[j].u.sort=indxt;
			jstack += 2;
			if (jstack > NSTACK) fprintf(stderr,"NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}
#undef NR_END
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI
*/
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */

/*
void EQNZBox(int np,pmparticletype *p, ZBox *zbox, float zfloor,
		float zceil,int nz){
	int i,j,k;
	float cputime0[10],cputime1[10];
	float stime,ftime;
	indextype *numparticle;
	indextype *tmpnum;
	indextype totaln,localn;
	indextype mynstart,mynfinal;
	indextype ntmp,ntmp2;
	ZBox *tmpzbox,aaa;

	MPI_Comm_size(MPI_COMM_WORLD,&nid);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	numparticle = (indextype *) Malloc(sizeof(indextype)*nid,
			PPTR(numparticle));
	tmpnum = (indextype *) Malloc(sizeof(indextype)*nid,PPTR(tmpnum));
	for(i=0;i<nid;i++) tmpnum[i] = 0;
	tmpnum[myid] = np;
	MPI_Reduce(tmpnum,numparticle,nid,MPI_INDEX,MPI_SUM,0,
				MPI_COMM_WORLD);
	MPI_Bcast(numparticle,nid*sizeof(indextype),MPI_BYTE,0,MPI_COMM_WORLD);
	Free(tmpnum);

	totaln = 0;
	for(i=0;i<nid;i++) totaln += numparticle[i];
	localn = totaln / (indextype) nid;
	mynstart = 0;
	for(i=0;i<myid;i++) mynstart += numparticle[i];
	mynfinal = mynstart+numparticle[myid];

	Free(numparticle);

	for(i=0;i<nid;i++) 
		zbox[i].zstart = zbox[i].zfinal = zbox[i].zwidth = 0.;
	MPI_Barrier(MPI_COMM_WORLD);
	TIMER_START(9);
	zsortparticles(np,p);
	TIMER_STOP(9);
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0) fprintf(stdout,"The sorting time is %f\n",ELAPSED_TIME(9));
	for(i=0;i<nid;i++) {
		ntmp = i*localn;
		if(ntmp >=mynstart && ntmp < mynfinal){
			ntmp2 = ntmp-mynstart;
			zbox[i].zstart = FindZofParticle(ntmp2,p);
		}
	}
	{
		Define_ZBox_MPI_Type;
		MPI_Type_struct(NELZBOX,blocklens,indices,old_types,&MPI_ZBox_type);
		tmpzbox = (ZBox *)Malloc(sizeof(ZBox)*nid,PPTR(tmpzbox));
		for(i=0;i<nid;i++) tmpzbox[i] = zbox[i];
		MPI_Reduce(tmpzbox,zbox,nid,MPI_ZBox_type,myop,0,
				MPI_COMM_WORLD);
		Free(tmpzbox);
		if(myid==0) {
			zbox[0].zstart = zfloor;
			zbox[nid-1].zfinal = (float) nz;
			for(i=0;i<nid-1;i++) zbox[i].zfinal = zbox[i+1].zstart;
		}
		MPI_Bcast(zbox,nid,MPI_ZBox_type,0,MPI_COMM_WORLD);
		Delete_ZBox_MPI_Op;
	}
	return ;
}
*/
