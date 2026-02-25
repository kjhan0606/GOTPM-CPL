#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<mpi.h>
#include "Memory.h"

#include "srfftw_mpi.h"




int local_ny_after_transpose,local_y_start_after_transpose,total_local_size;

static rfftwnd_mpi_plan plan,iplan;

void fftwbroadcast_(int *, int *, int*, int *, int *, int *, int *, int *);

void fftwinit(int nx,int ny,int nz,int *local_nz, int *local_z_start){
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);


    plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,nx,ny,nz,FFTW_REAL_TO_COMPLEX,
                    FFTW_ESTIMATE);
    iplan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,nx,ny,nz,FFTW_COMPLEX_TO_REAL,
                    FFTW_ESTIMATE);
    rfftwnd_mpi_local_sizes(plan,local_nz,local_z_start,&local_ny_after_transpose,
                    &local_y_start_after_transpose,&total_local_size);
    printf("P%d has passed fftw %d\n",myid,total_local_size);
    if(*local_nz != local_ny_after_transpose)
            printf("P%d has strange parallel fftw domain %d : %d\n",myid,*local_nz,
                        local_ny_after_transpose);
    fftwbroadcast_(&nx,&ny,&nz,local_nz,local_z_start,&local_ny_after_transpose,
        &local_y_start_after_transpose,&total_local_size);
}
void fftwforward(float *den){
    int myid;
    fftw_real *work;
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    if(myid ==0) printf("CPU(F-FFT) starting\n");
    work = (fftw_real *) Malloc(sizeof(fftw_real)*total_local_size,PPTR(work));
    MPI_Barrier(MPI_COMM_WORLD);
    rfftwnd_mpi(plan,1,den,work,FFTW_TRANSPOSED_ORDER);
    Free(work);
}
void fftwbackward(float *den){
    int myid;
    fftw_real *work;
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    if(myid ==0) printf("CPU(B-FFT) starting\n");
    work = (fftw_real *) Malloc(sizeof(fftw_real)*total_local_size,PPTR(work));
    MPI_Barrier(MPI_COMM_WORLD);
    rfftwnd_mpi(iplan,1,den,work,FFTW_TRANSPOSED_ORDER);
    Free(work);
}
void fftwforward_(float *den){
    fftw_real *work;
    work = (fftw_real *) Malloc(sizeof(fftw_real)*total_local_size,PPTR(work));
    MPI_Barrier(MPI_COMM_WORLD);
    rfftwnd_mpi(plan,1,den,work,FFTW_NORMAL_ORDER);
    Free(work);
}
void fftwbackward_(float *den){
    fftw_real *work;
    work = (fftw_real *) Malloc(sizeof(fftw_real)*total_local_size,PPTR(work));
    MPI_Barrier(MPI_COMM_WORLD);
    rfftwnd_mpi(iplan,1,den,work,FFTW_NORMAL_ORDER);
    Free(work);
}
