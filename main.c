#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "pmheader.h"
#include "Memory.h"


int Zeldovich(int, char**);
int SecondOrderLPT(int, char**);
SimParameters read_sim_parameter_file(FILE*),tmppar;

#ifdef INTEL
int main(int argc, char *argv[]){
#elif PGCC
int MAIN_(int argc, char *argv[]){
#else
#error You must type the compiler for the fortran and C compatibility setting for main.
#endif
	int myid,nid;
	FILE *simfile;
#ifdef _OPENMP
	void initfftw4MpiOpenMP(int , char **);
    initfftw4MpiOpenMP(argc, argv);
#else
    MPI_Init(&argc,&argv);
    fftwf_mpi_init();
#endif

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);


	if( argc == 2 ) {
		int sflag=1;
		if(myid==0) simfile = fopen(argv[1],"r");
		if(myid==0 && !simfile) sflag = 0;
		MPI_Bcast(&sflag,1,MPI_INT,0,MPI_COMM_WORLD);
		if(sflag ==0){
		    if(myid==0) fprintf(stderr,"Can't open %s - exiting\n",argv[1]);
			MPI_Finalize();
			exit(0);
		}
   	}
	else {
		fprintf(stderr,"usage: namu.exe paramsfile \n");
		MPI_Finalize();
		exit(0);
	}

	if(myid ==0){
		tmppar = read_sim_parameter_file(simfile);
	}
	MPI_Bcast(&tmppar,sizeof(SimParameters),MPI_BYTE,0,MPI_COMM_WORLD);
	if(tmppar.IC == 1) Zeldovich(argc, argv);
	else SecondOrderLPT(argc, argv);
}
