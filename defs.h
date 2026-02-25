#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
/*
#include "struct.h"
*/
#include "macros.h"

/* Stride for observational mode */
#define ObservationStride 5 
#define ModulusOfObservationStride 1

float gettime();
/*
char *mymalloc(), *resizelast();
*/
#include "Memory.h"

/* Domain Decomposition define's in domdecomp.c */

#define WORKTOL 1.0e-3  /* precision for work division in decomposition */
#define MAXITER 30      /* maximum allowed number of iterations */
#define PSKIP   10      /* skip parameter for the work calculation */

/*
#define INFINITY -1
*/

/* ------- error codes --------- */
#define OPENERROR 1
#define MEMORY_ERROR 2
#define READ_ERROR 3
#define NPROCMAX 512                  /* maximum number of processors */

/* Axes definitions */
#define XAXIS 0
#define YAXIS 1
#define ZAXIS 2

/* pn which is the master */
#define MASTER_PN   0

/* User msg tags */
#define FREE_WORKER 0
#define WORK        1
#define ANSWER      2
#define STOP        3
#define HEADER_TAG 4
#define PARTICLE_TAG 5
#define DONE 6
#define WORK_TAG 7
#define SPLIT_TAG 8
#define NP_TAG 9
#define VEC_TAG 10
#define TREE_TAG 11
#define NNODES_TAG 12
#define RMIN_TAG 13
#define RMAX_TAG 14
#define CELL_TAG 15
#define RUN_STATS_TAG 16
#define FREE_MASTER 17
#define DOMAIN_STATS_TAG 18
#define READY_TAG 19
#define COFM_TAG 20
#define ACKNOWLEDGE 21
#define VISITOR 22
#define NUMRECV_TAG 23
#define PSEND_TAG 24
#define PRECV_TAG 25
#define PASSFLAG 33
#define COUNTFLAG 34
#define PARTICLEFLAG 35
#define VCOUNT 36
#define VPARTICLES 37
#define GOFLAG 38

#define BRANCH 0
#define LEAF   1

/* For SyncOrDie */
#define DIE  0
#define SYNC 1

#define P_BLOCK_SIZE 20000


#define XY 0
#define XZ 1
#define YZ 2


#define PARTICULE_TYPE 1
#define CELL_TYPE 2
#define NCHILDREN 8

#ifndef NEW_BARNES_HUT
#define NEW_BARNES_HUT
#endif

#define NBINS 4000000

/* Note NTREEMIN < GROUPNUM or there will be trouble (check SingleGroup in
group.c */

#ifndef NTREEMIN
#define NTREEMIN 32
#endif
#ifndef GROUPNUM
#define GROUPNUM 40
#endif
/* Factor of the softening length in the unit of mean particle separation */
#define FACTEPSSQ 0.1
/*
#ifndef EPSSQ
#define EPSSQ 4.0e-8
#define EPSSQ 5.96046e-10
#endif
*/

#define BLOCKSIZE 32768

#define FLOAT_TYPE float
