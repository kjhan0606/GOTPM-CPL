/* Send and receive macros */
#ifdef TIMESEND
#define SEND(dest,tag,buf,size)  \
{ \
TIMER_START(33); \
MPI_Send((buf),(size),MPI_BYTE,(dest),(tag),MPI_COMM_WORLD); \
TIMER_STOP(33); \
if( myid == 0 ) { \
	if( size > 1024 ) { \
		totaltime += ELAPSED_TIME(33); \
fprintf(stdout,"Sent %d bytes in %g seconds - speed=%g MB/sec -- elapsed %g\n", \
	size, ELAPSED_TIME(33),size/ELAPSED_TIME(33)/1048576, totaltime); \
	} \
	} \
}
#else
#define SEND(dest,tag,buf,size)  \
MPI_Send((buf),(size),MPI_BYTE,(dest),(tag),MPI_COMM_WORLD)
#endif
#define RECEIVE(dest,tag,buf,size)               \
{                                                \
	MPI_Recv((buf),(size),MPI_BYTE,(dest),(tag), \
              MPI_COMM_WORLD,&receive_status);   \
}
#define SENDRECV(src, srctag, srcbuf, srcsize, dest, desttag, destbuf, destsize) \
{                                                            \
	MPI_Status sendrecv_status;                              \
	MPI_Sendrecv(destbuf, destsize, MPI_BYTE, dest, desttag, \
				 srcbuf,  srcsize,  MPI_BYTE,  src,  srctag, \
                 MPI_COMM_WORLD, &sendrecv_status);          \
}
#define SYNCHRONIZE() MPI_Barrier(MPI_COMM_WORLD)
#define ANY_TAG MPI_ANY_TAG
#define PROBE()

/* timer macros */
/*
#define TIMER_START(A)  ;
#define TIMER_STOP(A)   ;
#define ELAPSED_TIME(A) 0
#define BUSY_TIME(A) 0
#define IDLE_TIME(A) 0
*/
#define TIMER_START(A) cputime0[A] = gettime();
#define TIMER_STOP(A)  cputime1[A] = gettime();
#define ELAPSED_TIME(A) (cputime1[A] - cputime0[A])
#define BUSY_TIME(A)  0
#define IDLE_TIME(A)  0
#define EXIT(A) { MPI_Abort(MPI_COMM_WORLD,A); exit(A); }
#define CLEAN_EXIT(A) { MPI_Finalize(); exit(A); }
