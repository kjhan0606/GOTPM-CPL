#include "defs.h"
#include "common_simparam.h"
#include "pmheader.h"
#include "pm_common.h"
#include "params.h"
static int  myid,nid;
#define max(a,b) ((a) > (b) ? (a): (b))
#define min(a,b) ((a) < (b) ? (a): (b))

#include<sys/times.h>


/*
pmparticletype pblock[BLOCKSIZE];
*/
#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))

int jwrite(char *rvprefix, int stepcount) {
	FILE *fp;
	int i,j,k, p;
	int count;
	int ncount, nblock, nremain, ncheck;
	int dummy=0;
	float m;
	char *fileheader = FILEHEADER;
	char procname[3],number[3];
	char filename[80];
	MPI_Status status;
	size_t nwritecheck;
	int oldstepcount;
	int src,tgt,WGroupSize;
	int isend,iget,itag=1;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	sprintf(simpar.rvfilename,"%s.%.5d",rvprefix,stepcount);

	if( myid == 0 ) 
		fprintf(stderr," writing data at a=%g,z=%g \n",simpar.anow,simpar.amax/simpar.anow - 1.0);
#ifdef SLOW_RW
#ifndef OLD_RW
	for(i=0;i<WGROUPSIZE;i++){
		if(myid%WGROUPSIZE == i)
#else
	for(i=0;i<nid;i++){
		if(i == myid)
#endif
		{
			TIMER_START(5);
			sprintf(filename,"%s.%.5d%.5d",rvprefix,stepcount,myid);
			fprintf(stderr," writing file: %s %d ",filename,np);
			fp = fopen(filename,"w");
			write_head(fp,simpar);
			ncheck = np;
			count = 0;
			fwrite(pmparticles,sizeof(pmparticletype),np,fp);
			fclose(fp);
			TIMER_STOP(5);
			fprintf(stderr," written: %g \n",ELAPSED_TIME(5));

		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
#else
	isend = iget = 1;
	WGroupSize = WGROUPSIZE;
	src = myid-1;
	tgt = myid+1;
	if(RANKINGROUP(myid,WGroupSize) != 0 ) 
		MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);
	{
		TIMER_START(5);
		sprintf(filename,"%s.%.5d%.5d",rvprefix,stepcount,myid);
		fprintf(stderr," writing file: %s %ld ",filename,np);
		fp = fopen(filename,"w");
		write_head(fp,simpar);
		ncheck = np;
		count = 0;
		fwrite(pmparticles,sizeof(pmparticletype),np,fp);
		fclose(fp);
		TIMER_STOP(5);
		fprintf(stderr,"P%d written: %g \n",myid,ELAPSED_TIME(5));
	}
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) 
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);
#endif
/*########################################################### */
	if( myid == 0 ) {
		FILE *paramfile;
		char paramname[80], redpk[190];

		sprintf(paramname,"params.%.5d",stepcount);
		paramfile = fopen(paramname,"w");
		
		{
			void write_sim_parameter_file(FILE *, SimParameters );
			write_sim_parameter_file(paramfile,simpar);
		}

		fclose(paramfile);
	}
	return 0;
}

void jread(char *rvfilename)
{
	FILE *rvfile;
	int i, j, k, count;
	INT8 nmax;
	int nsend, nrecv;
	INT8 NTotal, NLocal, NLocal0, NTotalcheck;
        int Ntotal;
	INT8 nbytes;
	float m, mrv[6];
	char rvname[80],filename[80];
	MPI_Status status;
	float ksum, totalken;
	pmparticletype *p;
	int src,tgt,WGroupSize;
	int isend,iget,itag=1;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
/* passed the test already for the existence of this file */
	sprintf(filename,"%s%.5d",rvfilename,myid);
#ifdef SLOW_RW
#ifdef OLD_RW
	for(i=0;i<WGROUPSIZE;i++){
		if(myid%WGROUPSIZE == i){
#else
	for(i=0;i<nid;i++){
		if(i == myid){
#endif
			rvfile = fopen(filename,"r");
			simpar=read_head(rvfile);
			np = simpar.np;
			zstart = simpar.zmin;
			zheight = simpar.zmax-simpar.zmin;
			NLocal = np;
			nmax = ptrsize(pmparticles)/sizeof(pmparticletype);
			fprintf(stderr,"P%d nmax %ld\n",myid,nmax); fflush(stderr);
           	fprintf(stderr,"P%d starts reading data. %ld",myid,np);
			fflush(stderr);
			TIMER_START(5);
			fread(pmparticles,sizeof(pmparticletype),np,rvfile);
			fclose(rvfile)
			TIMER_STOP(5);
           	fprintf(stderr," --  P%d finishes reading data, %ld in %g sec. \n",myid, np,ELAPSED_TIME(5));
			fflush(stderr);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
#else
	isend = iget = 1;
	WGroupSize = WGROUPSIZE;
	src = myid-1;
	tgt = myid+1;
	if(RANKINGROUP(myid,WGroupSize) != 0 ) 
		MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);
	{
			rvfile = fopen(filename,"r");
			simpar=read_head(rvfile);
			np = simpar.np;
			zstart = simpar.zmin;
			zheight = simpar.zmax-simpar.zmin;
			NLocal = np;
			nmax = ptrsize(pmparticles)/sizeof(pmparticletype);
			fprintf(stderr,"P%d nmax %lld\n",myid,nmax); fflush(stderr);
           	fprintf(stderr,"P%d starts reading data. %ld",myid,np);
			fflush(stderr);
			TIMER_START(5);
			fread(pmparticles,sizeof(pmparticletype),np,rvfile);
			fclose(rvfile);
			TIMER_STOP(5);
           	fprintf(stderr," --  P%d finishes reading data, %ld in %g sec. \n",myid, np,ELAPSED_TIME(5));
			fflush(stderr);

	}
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) 
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);

#endif
/* Initialize kinetic energy term */
	ksum = 0;
	for(p=pmparticles; p<pmparticles+np; p++) {
#if defined(MULTIMASS)
		ksum += (p->mass)*((p->vx)*(p->vx) + (p->vy)*(p->vy) + (p->vz)*(p->vz));
#else
		ksum += ((p->vx)*(p->vx) + (p->vy)*(p->vy) + (p->vz)*(p->vz));
#endif
	}
	MPI_Allreduce(&ksum,&totalken,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	ken = totalken;
}
