#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include <mpi.h>
#include "fof.h"
#define DEFINE_SIM_PARA
#include "pmheader.h"
#undef DEFINE_SIM_PARA
#include "params.h"
//#define MP 1000000
//READTYPE rp[MP];


#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))



#define MIN(a,b) (a)<(b)? (a):(b)
#define MAX(a,b) (a)>(b)? (a):(b)
POSTYPE zstart,zwidth;
size_t np;
size_t  Fread(void *a,size_t b,size_t c, FILE *fp){
        char *A;
        char t1,t2,t3,t4;
        size_t i,nmem;
        nmem = fread(a,b,c,fp);
        A = (char *)a;
        for(i=0;i<b*nmem;i+=4){
                t1 = A[i];
                t2 = A[i+1];
                t3 = A[i+2];
                t4 = A[i+3];
                A[i] =t4;
                A[i+1] =t3;
                A[i+2] =t2;
                A[i+3] =t1;
        }
        return nmem;
}

#define Fread(a,b,c,d) fread(a,b,c,d)


size_t readparticle(FoFTPtlStruct **Bp,size_t np,int nstep,int nfile,int nz,
		char *infile, int mp, MPI_Comm comm){
	FoFTPtlStruct *bp,*p;
	size_t j,k;
	int i;
	FILE *fp;
	char outfile[190];
	FILE *wp;
	int nid,myid;
	size_t npadd;
	MPI_Status status;
	int itag=1;
	int MP = 1000000;

	MPI_Comm_rank(comm,&myid);
	MPI_Comm_size(comm,&nid);

	/*
	int isend,iget;
	isend = iget = 1;
	int WGroupSize = WGROUPSIZE;
	int src = myid -1;
	int tgt = myid + 1;
	if(RANKINGROUP(myid,WGroupSize) !=0)
		MPI_Recv(&i,1,MPI_INT,src,itag,comm,&status);
		*/

	int ok = 2;
	int tgt = 0;
	int stag = 2;
	MPI_Send(&ok, 1, MPI_INT, tgt, stag, MPI_COMM_WORLD);
	printf("G%d is sending request to P0 to get permission to read\n", myid);
	MPI_Recv(&ok, 1, MPI_INT, tgt, stag, MPI_COMM_WORLD, &status);
	printf("G%d is getting permission to read\n", myid);


	npadd = 0;
	fp=fopen(infile,"r");
	if(fp){
			double zmin,zmax;
			simpar = read_head(fp);
			zstart = simpar.zmin;
			zwidth = simpar.zmax-simpar.zmin;
			mp = simpar.np;
				
			bp = *Bp;
			bp = (FoFTPtlStruct *)realloc(bp,sizeof(FoFTPtlStruct)*(np+mp));
			*Bp = bp;
			p = bp + np;

			printf("G%d Now opening %s with %d particles from np= %d",myid,infile,mp,(int)np);
			fflush(stdout);
			zmin = 1e20;
			zmax = -1e20;
			READTYPE *rp = (READTYPE*)malloc(sizeof(READTYPE)*MP);;
			while((mp=Fread(rp,sizeof(READTYPE),MP,fp))>0){
				for(j=0;j<mp;j++){
					p->r[0]=XofP(rp+j);
					p->r[1]=YofP(rp+j);
					p->r[2]=ZofP(rp+j);
					p->rv[0]=rp[j].vx;
					p->rv[1]=rp[j].vy;
					p->rv[2]=rp[j].vz;
					p->indx=rp[j].indx;
					zmin = MIN(zmin,p->r[2]);
					zmax = MAX(zmax,p->r[2]);
					npadd++;
					p++;
				}
			}
			free(rp);

			printf("Well read %ld particles zmin: zmax   %lg : %lg\n",npadd,zmin,zmax);
			fflush(stdout);
			fclose(fp);
	}
	else {
		printf("P%d is now allocated with a file = %s . However, maybe, it is out of the file boundary\n",myid, infile);
	}
	stag = 3;
	ok = 2;
	MPI_Send(&ok, 1, MPI_INT, tgt, stag, MPI_COMM_WORLD);
	printf("G%d is sending the  report of finishing the read\n", myid);
	/*
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid)
		MPI_Send(&i,1,MPI_INT,tgt,itag,comm);
		*/

	return npadd;
}
