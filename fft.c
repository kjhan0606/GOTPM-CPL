/* Modified for using the FFTW ver. 3 (Nov. 2023)*/
/* For historical reason in developing the GOTPM and Eunha codes 
   that are based on the Fortran language, the array for the fftw 
   is in column-major order. Therefore A(x,y,z) = A(x+mx*(y+ny*z)) 
   when expressed in the one-dimensional array. To summarize, the 
   order of dimensions in the code is reversed compared to the C 
   convention,  while we decided to use the Fortran convention.
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include<complex.h>
#include<omp.h>

#include "fftw3.h"
#include "fftw3-mpi.h"



static fftwf_plan plan, iplan, tplan, itplan;


static int fftwinitflag = 1;
int threads_ok;

void initfftw4MpiOpenMP(int argc, char **argv){
	if(fftwinitflag){
		int provided, ierr, resultlen;
		if( (ierr=MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,&provided)) )
		{
			char string[100];
			MPI_Error_string(ierr, string, &resultlen);
			printf("Something wrong in the MPI_Init_thread: %s\n", string);
		}
		int nthreads=1;
#ifdef _OPENMP
#pragma omp parallel
		{
#pragma omp master
			nthreads = omp_get_num_threads();
		}
#endif
		printf("Provided thread num %d  requested= %d & set nthreads= %d\n",provided, MPI_THREAD_FUNNELED,nthreads);
		threads_ok = provided >= MPI_THREAD_FUNNELED;
		if(threads_ok) threads_ok = fftwf_init_threads();
		fftwf_mpi_init();
		if(threads_ok) fftwf_plan_with_nthreads(nthreads);
		fftwinitflag = 0;
	}
}
void set3dFftwInfo(int nx, int ny, int nz, MPI_Comm mpicom, int *local_nz, int *local_z_start){
	ptrdiff_t Local_grid_size, Local_grid_size_after_transpose;
	ptrdiff_t Local_nz, Local_z_start;
	ptrdiff_t Local_ny_after_transpose, Local_y_start_after_transpose;
	Local_grid_size = fftwf_mpi_local_size_3d_transposed( nz, ny, nx/2+1, mpicom, 
			&Local_nz, &Local_z_start, 
			&Local_ny_after_transpose, &Local_y_start_after_transpose
			);
	*local_nz = Local_nz;
	*local_z_start = Local_z_start;

	{
		int local_y_after_transpose, local_y_start_after_transpose;
		local_y_after_transpose = Local_ny_after_transpose;
		local_y_start_after_transpose = Local_y_start_after_transpose;
		void fftwbroadcast_(int *, int *, int *, int *, int *, int *, int *);
		fftwbroadcast_(&nx,&ny,&nz,local_nz, local_z_start, 
				&local_y_after_transpose, &local_y_start_after_transpose);
	}
	int local_grid_size = Local_grid_size;

	float *rden = (float*)malloc(sizeof(fftwf_complex)*local_grid_size);
	fftwf_complex *cden = (fftwf_complex*)rden;
	printf("localsize = %d and starting point = %d\n", *local_nz, *local_z_start);
	fflush(stdout);

	/* transposed forward fftw */
	tplan = fftwf_mpi_plan_dft_r2c_3d(
			nz,ny,nx,
			rden, cden, mpicom,
			FFTW_MPI_TRANSPOSED_OUT | FFTW_ESTIMATE | FFTW_UNALIGNED 
			);
	/* transposed backward fftw */
	itplan = fftwf_mpi_plan_dft_c2r_3d(
			nz, ny, nx, 
			cden, rden, mpicom, 
			FFTW_MPI_TRANSPOSED_IN | FFTW_ESTIMATE | FFTW_UNALIGNED 
			);
	/* normal forward fftw */
	plan = fftwf_mpi_plan_dft_r2c_3d(
			nz,ny,nx,
			rden, cden, mpicom,
			FFTW_ESTIMATE | FFTW_UNALIGNED
			);
	/* normal backward fftw */
	iplan= fftwf_mpi_plan_dft_c2r_3d(
			nz,ny,nx,
			cden, rden, mpicom, 
			FFTW_ESTIMATE | FFTW_UNALIGNED
			);
	free(rden);
}

void fftwndnormalforward_(float *rden){
	fftwf_complex *cden = (fftwf_complex*) rden;
	fftwf_mpi_execute_dft_r2c(
			plan,
			rden,
			cden
			);
}
void fftwndnormalbackward_(fftwf_complex *cden){
	float *rden = (float *) cden;
	fftwf_mpi_execute_dft_c2r(
			iplan,
			cden,
			rden
			);
}
void fftwndtransforward_(float *rden){
    fftwf_complex *cden = (fftwf_complex*) rden;
    fftwf_mpi_execute_dft_r2c(
            tplan,
            rden,
            cden
            );
}
void fftwndtransbackward_(fftwf_complex *cden){
    float *rden = (float *) cden;
    fftwf_mpi_execute_dft_c2r(
            itplan,
            cden,
            rden
            );
}


/*
void malloc_fftw_den(FFTWGridInfo *fftwgrid, FFTW_INFO *fftwinfo){
	fftwgrid->den.rden = (DenType*) fftwf_malloc(sizeof(DenType)*fftwinfo->local_grid_size*2);
}
void free_fftw_den(FFTWGridInfo *fftwgrid, FFTW_INFO *fftwinfo){
	fftwf_free(fftwgrid->den.rden);
}
*/

void mpi_fftw_finalize(void){
	if(fftwinitflag == 0){
		fftwinitflag = 1;
		MPI_Finalize();
	}
}


/*
void set_mpi_fftw_3d_all_info(FFTWGridInfo *fftwgrid, FFTW_INFO *fftwinfo, MPI_Comm com){
	set_mpi_fftw_3d_info(fftwgrid,fftwinfo,com);
	malloc_fftw_den(fftwgrid,fftwinfo);
	set_mpi_fftw_3d_plans(fftwgrid,fftwinfo,com);
}
*/
