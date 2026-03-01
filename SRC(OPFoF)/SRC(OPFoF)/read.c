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
#define MP 1000000

READTYPE rp[MP];
float zstart,zwidth;
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


size_t readparticle(FoFTPtlStruct **Bp,size_t np, char *infile, char pflag){
	FoFTPtlStruct *bp,*p;
	size_t j,k;
	int i,mp,mmp;
	int next,before;
	FILE *fp;
	char outfile[190];
	FILE *wp;
	int nid,myid;
	size_t npadd;
	MPI_Status status;
	int itag;

	printf("%s is being read with pflag =%c\n",infile,pflag);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);

	next = (myid+1+nid)%nid;
	before = (myid-1+nid)%nid;
	itag = 0;
	if(myid != 0) MPI_Recv(&i,1,MPI_INT,before, itag,MPI_COMM_WORLD,&status);
	{
			fp=fopen(infile,"r");
			if(pflag == 'N') {
				simpar = read_head(fp);
				zstart = simpar.zmin;
				zwidth = simpar.zmax-simpar.zmin;
				mp = simpar.np;
			}
			else {
				Fread(&zstart,sizeof(float),1,fp);
				Fread(&zwidth,sizeof(float),1,fp);
				Fread(&mp,sizeof(int),1,fp);
			}
				
			bp = *Bp;
			bp = (FoFTPtlStruct *)realloc(bp,sizeof(FoFTPtlStruct)*(np+mp));
			*Bp = bp;
			p = bp + np;

			npadd = 0;
			printf("Now opening %s with %d particles from np= %d\n",infile,mp,np);
			fflush(stdout);
			while((mmp=Fread(rp,sizeof(READTYPE),MP,fp))){
				for(j=0;j<mmp;j++){
					p->r[0]=XofP(rp+j);
					p->r[1]=YofP(rp+j);
					p->r[2]=ZofP(rp+j);
					p->rv[0]=rp[j].vx;
					p->rv[1]=rp[j].vy;
					p->rv[2]=rp[j].vz;
					p->indx=rp[j].indx;
					npadd++;
					if(npadd > mp){
						fprintf(stderr,"Error: # of simpar.np is less than the real #\n");
						exit(9999);
					}
					p++;
				}
			}
			printf("Well read %ld particles \n",npadd);
			fflush(stdout);
			fclose(fp);
	}
	if(myid != nid-1) MPI_Send(&i,1,MPI_INT,next, itag,MPI_COMM_WORLD);
	return npadd;
}
