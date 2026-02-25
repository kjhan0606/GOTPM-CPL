/*
 * PMParticleType Migration Function. 07/11/2005 by Juhan Kim
 * Adjusted for the low-usage of buffer space and the faster transfer rate.
 */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include "pmheader.h"
#include "mpi.h"
#include "Memory.h"
#define max(A,B) ( (A) > (B) ? (A) : (B) )
#define min(A,B) ( (A) < (B) ? (A) : (B) )
#define TICKS 100.
MPI_Datatype MPI_SIZE_T;
#define GETNUMCOMM(mp,mmp) {\
	MPI_Reduce(&mp,&mmp,1,MPI_SIZE_T,MPI_SUM,0,MPI_COMM_WORLD);\
	MPI_Bcast(&mmp,1,MPI_SIZE_T,0,MPI_COMM_WORLD);\
}
int nid,myid;
/*
#define MPI_SENDRECV(send,numa,tpa,dst,sf,recv,numb,tpb,src,rf,comm,sta){\
	for(i=0;i<nid;i++){\
		if(myid==i) MPI_Send(send,numa,tpa,dst,sf,comm);\
		else if(src==i) MPI_Recv(recv,numb,tpb,src,rf,comm,sta);\
		MPI_Barrier(comm);\
	}\
}
*/
#define MPI_SENDRECV(send,numa,tpa,dst,sf,recv,numb,tpb,src,rf,comm,sta){\
	MPI_Sendrecv(send,numa,tpa,dst,sf,recv,numb,tpb,src,rf,comm,sta);\
}

extern pmparticletype  *pmparticles;

size_t maxbuff,chunksize;
int nloop;


void migrate(size_t *mp,int nz,float zwidth, float zmin){
  	pmparticletype *p,tmp,*ptr,*bp,*sp;
	pmparticletype *sndp,*lp,*rp,*bp2send;
	MPI_Status status;
	int src,dest;
	size_t nrecv,nsend,commsize;
	size_t i,j,k,np,sndnp;
	float zmax;
#ifdef XYZDBL
	double zp,zup,zdown;
#else
	float zp,zup,zdown;
#endif
	size_t niter,tniter;


	if(sizeof(size_t)==sizeof(unsigned int)) {
		MPI_SIZE_T = MPI_UNSIGNED;
	}
	else if(sizeof(size_t) == sizeof(unsigned long)){
		MPI_SIZE_T = MPI_UNSIGNED_LONG;
	}
	else{
		fprintf(stderr,"Wrong size of size_t in adkjhmigrate.c\n");
		exit(99999);
	}

	np = *mp;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);


	/*
	zmax = zmin + zwidth;
	*/
	dest = (myid+nid-1)%nid;
	src = (myid+nid+1)%nid;
	MPI_SENDRECV(&zmin,1,MPI_FLOAT,dest,0,
			&zmax,1,MPI_FLOAT,src,0,MPI_COMM_WORLD,&status);
	if(myid == nid-1) zmax = (float) nz;


	commsize = 100;
	sp = pmparticles;


	dest = (myid+1+nid)%nid;
	src = (myid-1+nid)%nid;
	while(commsize != 0){
		rp = pmparticles + np;
		for(lp=sp;lp<rp;){
			zp = ZofP(lp);
			if( !(zp <zmax && zp >= zmin) ){
				zup = fabs(zmax - zp);
				zup = min(zup,nz-zup);
				zdown = fabs(zmin-zp);
				zdown = min(zdown,nz-zdown);
				if(zup <= zdown){
					rp --;
					tmp = *lp;
					*lp = *rp;
					*rp = tmp;
				}
				else lp ++;
			}
			else lp++;
		}
		/* rp is the right bound of the domain pmparticles 
		 * bp2send is the left bound of the pmparticles to send */
		bp2send = rp;
		sndnp = np-(bp2send-pmparticles);
		nsend = sndnp;


		MPI_SENDRECV(&nsend,1,MPI_SIZE_T,dest,0,
				&nrecv,1,MPI_SIZE_T,src,0,MPI_COMM_WORLD,&status);

		GETNUMCOMM(nsend,commsize);
		if(myid==0) printf("Total commsize = %ld",commsize);
		/* maxbuff: maximum available buffer size of local domain  
		 * = freespace - max(nrecv-nsend,0) after all migrations are completed.
		 * chunksize: the maximum allowable buffer size on whole domains 
		 * schunk: temporary size of buffer to send
		 * rchunk: temporary size of buffer to receive
		 * tniter: total number of iteration to send and receive chunks */
		maxbuff = freespace()/sizeof(pmparticletype)-max(nrecv-nsend,0)-10;
		MPI_Reduce(&maxbuff,&chunksize,1,MPI_SIZE_T,MPI_MIN,0,MPI_COMM_WORLD);
		MPI_Bcast(&chunksize,1,MPI_SIZE_T,0,MPI_COMM_WORLD);

		{
			size_t snchunk,rnchunk,schunk,rchunk;
			niter = (size_t)ceil((double)nsend/(double)chunksize);
			MPI_Reduce(&niter,&tniter,1,MPI_SIZE_T,MPI_MAX,0,MPI_COMM_WORLD);
			MPI_Bcast(&tniter,1,MPI_SIZE_T,0,MPI_COMM_WORLD);
			snchunk = rnchunk = 0;
			for(i=0;i<tniter;i++){
				/* Determine the chunk size to send */
				if(snchunk < nsend){
					schunk = min(chunksize,nsend-snchunk);
					sndp = (pmparticletype *)
						Malloc(sizeof(pmparticletype)*schunk,PPTR(sndp));
					for(p=bp2send+snchunk,ptr=sndp;p<bp2send+snchunk+schunk;
							p++,ptr++) *ptr = *p;
				}
				else schunk = 0;
				MPI_SENDRECV(&schunk,1,MPI_SIZE_T,dest,0,
						&rchunk,1,MPI_SIZE_T,src,0,MPI_COMM_WORLD,&status);
				/* If received chunk size is larger than the original 
				 * particle array, we have to expand it to host the
				 * newly incoming chunk*/
				if(rnchunk+rchunk > nsend){
					pmparticles = (pmparticletype *)Realloc(pmparticles,
						sizeof(pmparticletype)*(np-nsend+rnchunk+rchunk));
					bp2send = pmparticles + (np - nsend);
				}
				MPI_SENDRECV(sndp,schunk*sizeof(pmparticletype),
						MPI_BYTE,dest,0,bp2send+rnchunk,
						rchunk* sizeof(pmparticletype),MPI_BYTE,src,0,
						MPI_COMM_WORLD,&status);

				if(snchunk < nsend) Free(sndp);
				snchunk += schunk;
				rnchunk += rchunk;
			}
		}
		pmparticles = (pmparticletype *)Realloc(pmparticles,
				sizeof(pmparticletype)*(np+nrecv-nsend));
		sp = pmparticles+np-nsend;
		np = np+nrecv-nsend;
		if(myid==0) printf("+total com pmparticles =%ld in iter %ld\n",(long)commsize,(long)tniter);
		fflush(stdout);
	}

	commsize =100;
	sp = pmparticles;

	dest = (myid-1+nid)%nid;
	src = (myid+1+nid)%nid;
	/*
	j = 0;
	nloop = 0;
	*/
	while(commsize != 0){
		rp = pmparticles + np;
		for(lp=sp;lp<rp;){
			zp = ZofP(lp);
			if( !(zp >= zmin && zp <zmax) ){
				rp --;
				tmp = *lp;
				*lp = *rp;
				*rp = tmp;
			}
			else lp ++;
		}


		/* rp is the right bound of the domain pmparticles 
		 * bp2send is the left bound of the pmparticles to send */
		bp2send = rp;
		sndnp = np-(bp2send-pmparticles);
		nsend = sndnp;
		/*
		if(nloop > nid/2) {
			for(i=0;i<nsend;i++){
				printf("P%d %g %g ::%d : %g\n",myid,zmin,zmax,nloop,(bp2send+i)->z);
			}
			MPI_Finalize();
			exit(999);
		}
		*/
		MPI_SENDRECV(&nsend,1,MPI_SIZE_T,dest,0,
				&nrecv,1,MPI_SIZE_T,src,0,MPI_COMM_WORLD,&status);

		GETNUMCOMM(nsend,commsize);
		/* maxbuff: the maximum buffer size of local domain  
		 * = freespace - max(nrecv-nsend,0) after all migrations are completed.
		 * chunksize: the maximum allowable buffer size in whole domains 
		 * schunk: the temporary size of buffer to send
		 * rchunk: the temporary size of buffer to receive
		 * tniter: total number of iteration to send and receive chunks */
		maxbuff = freespace()/sizeof(pmparticletype)-max(nrecv-nsend,0)-10;
		MPI_Reduce(&maxbuff,&chunksize,1,MPI_SIZE_T,MPI_MIN,0,MPI_COMM_WORLD);
		MPI_Bcast(&chunksize,1,MPI_SIZE_T,0,MPI_COMM_WORLD);
		{
			size_t snchunk,rnchunk,schunk,rchunk;
			niter = (size_t)ceil((double)nsend/(double)chunksize);
			MPI_Reduce(&niter,&tniter,1,MPI_SIZE_T,MPI_MAX,0,MPI_COMM_WORLD);
			MPI_Bcast(&tniter,1,MPI_SIZE_T,0,MPI_COMM_WORLD);
			snchunk = rnchunk = 0;
			for(i=0;i<tniter;i++){
				/* Determine the chunk size to send */
				if(snchunk < nsend){
					schunk = min(chunksize,nsend-snchunk);
					sndp = (pmparticletype *)
						Malloc(sizeof(pmparticletype)*schunk,PPTR(sndp));
					for(p=bp2send+snchunk,ptr=sndp;p<bp2send+snchunk+schunk;
							p++,ptr++) *ptr = *p;
				}
				else schunk = 0;
				MPI_SENDRECV(&schunk,1,MPI_SIZE_T,dest,0,
						&rchunk,1,MPI_SIZE_T,src,0,MPI_COMM_WORLD,&status);
				/* If received chunk size is larger than the original 
				 * particle array, we have to expand it to host the
				 * newly incoming chunk*/
				if(rnchunk+rchunk > nsend){
					pmparticles = (pmparticletype *)Realloc(pmparticles,
						sizeof(pmparticletype)*(np-nsend+rnchunk+rchunk));
					bp2send = pmparticles + (np - nsend);
				}
				MPI_SENDRECV(sndp,schunk*sizeof(pmparticletype),
						MPI_BYTE,dest,0,bp2send+rnchunk,
						rchunk* sizeof(pmparticletype),MPI_BYTE,src,0,
						MPI_COMM_WORLD,&status);

				if(snchunk < nsend) Free(sndp);
				snchunk += schunk;
				rnchunk += rchunk;
			}
		}
		pmparticles = (pmparticletype *)Realloc(pmparticles,
				sizeof(pmparticletype)*(np+nrecv-nsend));
		sp = pmparticles+np-nsend;
		np = np+nrecv-nsend;
		if(myid==0) printf("-total com pmparticles =%ld in iter %ld\n",(long)commsize,(long)tniter);
		/*
		nloop ++;
		*/
	}
	*mp = np;
	simpar.np = np;
}
