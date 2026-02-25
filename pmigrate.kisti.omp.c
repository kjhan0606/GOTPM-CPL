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
#define GETNUMCOMM(mp,mmp) {\
	MPI_Reduce(&mp,&mmp,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);\
	MPI_Bcast(&mmp,1,MPI_LONG,0,MPI_COMM_WORLD);\
}
#define MPI_MAX_BYTE 200000000L
static int nid,myid;
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

long maxbuff,chunksize;
int nloop;

#define BitsOfChar 8L
#define NO_MASK 0x00
#define YES_MASK  0x01

#define IsMaskedBit(mask,ip) ( (mask[ip/BitsOfChar] >> ip%BitsOfChar) & YES_MASK)
#define MarkBit(mask,ip) do {\
	mask[ip/BitsOfChar] |= (1 << (ip%BitsOfChar));\
} while(0)
#define ClearBit(mask,ip) do {\
	mask[ip/BitsOfChar] &= ~(1 << (ip%BitsOfChar));\
} while(0)



void migrate(size_t *mp,int nz,float zwidth, float zmin){
  	pmparticletype *p,tmp,*ptr,*bp,*sp;
	pmparticletype *sndp,*lp,*rp,*bp2send;
	MPI_Status status;
	int src,dest;
	long nrecv,nsend,commsize;
	long i,j,k,np,sndnp;
	float zmax;
	/*
#ifdef XYZDBL
	double zp,zup,zdown;
#else
	float zp,zup,zdown;
#endif
*/
	long niter,tniter;


	/*
	if(sizeof(long)==sizeof(unsigned int)) {
		MPI_LONG = MPI_UNSIGNED;
	}
	else if(sizeof(long) == sizeof(unsigned long)){
		MPI_LONG = MPI_UNSIGNED_LONG;
	}
	else{
		fprintf(stderr,"Wrong size of long in adkjhmigrate.c\n");
		exit(99999);
	}
	*/

	TIMER_START(59);

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

		if(rp > sp) {
			size_t ntarget = rp-sp;
			long long nbytes = (ntarget + BitsOfChar-1)/BitsOfChar;
			char *mask =(char*)Malloc(sizeof(char)*nbytes, PPTR(mask));
			for(i=0;i<nbytes;i++) mask[i] = NO_MASK;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
			for(i=0;i<ntarget;i++){
#ifdef XYZDBL
				double zp, zup, zdown;
#else
				float zp, zup,zdown;
#endif
				zp = ZofP(sp+i);
				if( zp >= zmax || zp < zmin){
					zup = fabs(zmax - zp); 
					zup = min(zup,nz-zup); 
					zdown = fabs(zmin-zp); 
					zdown = min(zdown,nz-zdown);
					if(zup<=zdown ){
						MarkBit(mask,i);
					}
				}
			}
			for(lp=sp;lp<rp;){
				i = lp-sp;
				if(IsMaskedBit(mask,i)) {
					rp --;
					tmp = *lp;
					*lp = *rp;
					*rp = tmp;
					j = (rp-sp);
					if(IsMaskedBit(mask,j)) {
						MarkBit(mask,i); /* mask set */
					}
					else {
						ClearBit(mask,i); /* mask clear */
					}
					MarkBit(mask,j); /* mask set */
				}
				else lp++;
			}
			Free(mask);
		}

		/* rp is the right bound of the domain pmparticles 
		 * bp2send is the left bound of the pmparticles to send */
		bp2send = rp;
		sndnp = np-(bp2send-pmparticles);
		nsend = sndnp;


		MPI_SENDRECV(&nsend,1,MPI_LONG,dest,0,
				&nrecv,1,MPI_LONG,src,0,MPI_COMM_WORLD,&status);

		GETNUMCOMM(nsend,commsize);
		if(myid==0) printf("Total commsize = %ld",commsize);
		/* maxbuff: maximum available buffer size of local domain  
		 * = freespace - max(nrecv-nsend,0) after all migrations are completed.
		 * chunksize: the maximum allowable buffer size on whole domains 
		 * schunk: temporary size of buffer to send
		 * rchunk: temporary size of buffer to receive
		 * tniter: total number of iteration to send and receive chunks */
		maxbuff = freespace()/sizeof(pmparticletype)-max(nrecv-nsend,0)-10;
		maxbuff = min(maxbuff,MPI_MAX_BYTE/sizeof(pmparticletype));
		MPI_Reduce(&maxbuff,&chunksize,1,MPI_LONG,MPI_MIN,0,MPI_COMM_WORLD);
		MPI_Bcast(&chunksize,1,MPI_LONG,0,MPI_COMM_WORLD);

		{
			long snchunk,rnchunk,schunk,rchunk;
			niter = (long)ceil((double)nsend/(double)chunksize);
			MPI_Reduce(&niter,&tniter,1,MPI_LONG,MPI_MAX,0,MPI_COMM_WORLD);
			MPI_Bcast(&tniter,1,MPI_LONG,0,MPI_COMM_WORLD);
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
				MPI_SENDRECV(&schunk,1,MPI_LONG,dest,0,
						&rchunk,1,MPI_LONG,src,0,MPI_COMM_WORLD,&status);
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

		if(rp >sp){
			size_t ntarget = rp-sp;
			long long nbytes = (ntarget + BitsOfChar-1)/BitsOfChar;
			char *mask =(char*)Malloc(sizeof(char)*nbytes, PPTR(mask));
			for(i=0;i<nbytes;i++) mask[i] = NO_MASK;
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for(i=0;i<ntarget;i++){
#ifdef XYZDBL
				double zp;
#else
				float zp;
#endif
				zp = ZofP(sp+i);
				if( zp < zmin || zp >=zmax){
					MarkBit(mask,i);
				}
			}
			for(lp=sp;lp<rp;){
				i = lp-sp;
				if(IsMaskedBit(mask,i)) {
					rp --;
					tmp = *lp;
					*lp = *rp;
					*rp = tmp;
					j = (rp-sp);
					if(IsMaskedBit(mask,j)) {
						MarkBit(mask,i); /* mask set */
					}
					else {
						ClearBit(mask,i); /* mask clear */
					}
					MarkBit(mask,j); /* mask set */
				}
				else lp++;
			}
			Free(mask);
		}


		/* rp is the right bound of the domain pmparticles 
		 * bp2send is the left bound of the pmparticles to send */
		bp2send = rp;
		sndnp = np-(bp2send-pmparticles);
		nsend = sndnp;
		MPI_SENDRECV(&nsend,1,MPI_LONG,dest,0,
				&nrecv,1,MPI_LONG,src,0,MPI_COMM_WORLD,&status);

		GETNUMCOMM(nsend,commsize);
		/* maxbuff: the maximum buffer size of local domain  
		 * = freespace - max(nrecv-nsend,0) after all migrations are completed.
		 * chunksize: the maximum allowable buffer size in whole domains 
		 * schunk: the temporary size of buffer to send
		 * rchunk: the temporary size of buffer to receive
		 * tniter: total number of iteration to send and receive chunks */
		maxbuff = freespace()/sizeof(pmparticletype)-max(nrecv-nsend,0)-10;
		maxbuff = min(maxbuff,MPI_MAX_BYTE/sizeof(pmparticletype));
		MPI_Reduce(&maxbuff,&chunksize,1,MPI_LONG,MPI_MIN,0,MPI_COMM_WORLD);
		MPI_Bcast(&chunksize,1,MPI_LONG,0,MPI_COMM_WORLD);
		{
			long snchunk,rnchunk,schunk,rchunk;
			niter = (long)ceil((double)nsend/(double)chunksize);
			MPI_Reduce(&niter,&tniter,1,MPI_LONG,MPI_MAX,0,MPI_COMM_WORLD);
			MPI_Bcast(&tniter,1,MPI_LONG,0,MPI_COMM_WORLD);
			if(myid==0) printf("-Total com pmparticles =%ld in iter %ld\n",(long)commsize,(long)tniter);
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
				MPI_SENDRECV(&schunk,1,MPI_LONG,dest,0,
						&rchunk,1,MPI_LONG,src,0,MPI_COMM_WORLD,&status);
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
	TIMER_STOP(59);
	if(myid==0) fprintf(stdout,"CPU(migrate)   = %f\n",ELAPSED_TIME(59));
}
