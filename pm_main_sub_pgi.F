c234567
      subroutine readpar(nspace,omep,omeplam,rth)
      integer nspace
      real omep,omeplam,rth

      print *, "nspace,omep,omeplam,rth ? "
      read(7,*) nspace,omep,omeplam,rth
      return
      end

c---  FFTW USE_WORK
      subroutine fftwinit(mx,my,mz,local_mz,local_zz_start)
      INCLUDE 'fftw_f77.i'
      INCLUDE 'mpif.h'
#include "fortran.i"
      integer nx,ny,nz,local_nz,local_z_start
      integer local_ny_after_transpose,local_y_start_after_transpose
      integer total_local_size
c     integer fftw_group_world, fftw_group,FFTW_COMM_WORLD,fftw_group_list
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      nx = mx
      ny = my
      nz = mz
#ifdef FFTW_DIVISIBLE_SIZE
#else
      call rfftw3d_f77_mpi_create_plan(plan,MPI_COMM_WORLD,nx,ny,nz,
     &       FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE+FFTW_IN_PLACE)
      call rfftw3d_f77_mpi_create_plan(iplan,MPI_COMM_WORLD,nx,ny,nz,
     &       FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE+FFTW_IN_PLACE)
#endif
cc***

c   Calculate and find the local array size for the argument of
c   "FFTW_TRANSPOSED_ORDER"
      call rfftwnd_f77_mpi_local_sizes(plan,local_nz,local_z_start,
     &          local_ny_after_transpose,local_y_start_after_transpose,
     &          total_local_size)

      local_mz = local_nz
      local_zz_start = local_z_start

      return
      end
c---  FFT with additional space for output.
      subroutine fftwinitaddspace(plan,iplan,nx,ny,nz,
     &          local_nz,local_z_start,
     &          local_ny_after_transpose,local_y_start_after_transpose,
     &          total_local_size)
      INCLUDE 'fftw_f77.i'
      INCLUDE 'mpif.h'
#include "fortran.i"
      integer nx,ny,nz,local_nz,local_z_start
      integer local_ny_after_transpose,local_y_start_after_transpose
      integer total_local_size

ccc**
c     rfftw3d_f77_mpi_create_plan : 3dimensional REAL-TO-COMPLEX FFT
c                               initialization with creating "plan" argument.
c                               This calculates FFT with in-place techniques.
      call rfftw3d_f77_mpi_create_plan(plan,MPI_COMM_WORLD,nx,ny,nz,
     &       FFTW_REAL_TO_COMPLEX,FFTW_MEASURE)
c     rfftw3d_f77_mpi_create_plan : 3dimensional COMPLEX-TO-REAL FFT
c                                initialization with creating "iplan" argument.
c                                This calculates FFT with in-place techniques.
      call rfftw3d_f77_mpi_create_plan(iplan,MPI_COMM_WORLD,nx,ny,nz,
     &       FFTW_COMPLEX_TO_REAL,FFTW_MEASURE)
cc***

c   Calculate and find the local array size for the argument of
c   "FFTW_TRANSPOSED_ORDER"
      call rfftwnd_f77_mpi_local_sizes(plan,local_nz,local_z_start,
     &          local_ny_after_transpose,local_y_start_after_transpose,
     &          total_local_size)

      return
      end
c--   FORWARD FFT with use_work
      subroutine fftforward(den)
      INCLUDE 'fftw_f77.i'
      INCLUDE 'mpif.h'
#include "fortran.i"
      real den(*)
      integer nx,ny,nz,local_nz,local_z_start
      integer local_ny_after_transpose,local_y_start_after_transpose
      integer total_local_size
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      real*4 CPU(2),stime
      real work
      integer myid,nid,ierror

      call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierror)
c     if(myid .EQ. 0) print *,' CPU(F-FFT) starting'

      if(myid.eq.0) stime = etime(CPU)
      call rfftwnd_f77_mpi(plan,1,den,work1,use_work,
     &    FFTW_TRANSPOSED_ORDER)
      if(myid .EQ. 0) print *,' CPU(F-FFT) ', etime(CPU)-stime
      return
      end

      subroutine fftbackward(den)
      INCLUDE 'fftw_f77.i'
      INCLUDE 'mpif.h'
#include "fortran.i"
      real den(*)
      integer nx,ny,nz,local_nz,local_z_start
      integer local_ny_after_transpose,local_y_start_after_transpose
      integer total_local_size
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      real*4 CPU(2),stime
      real work
      integer myid,nid,ierror

      call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierror)
      call MPI_Comm_size(MPI_COMM_WORLD,nid,ierror)
c     if(myid .EQ. 0) print *,' CPU(B-FFT) starting'

      if(myid.eq.0) stime = etime(CPU)
      call rfftwnd_f77_mpi(iplan,1,den,work1,use_work,
     &               FFTW_TRANSPOSED_ORDER)
      if(myid .EQ. 0) print *,' CPU(B-FFT) ', etime(CPU)-stime
      return
      end

c--------------------------------------------------------------------
c--   FORWARD FFT with additional work space
      subroutine fftforwardaddspace(plan,den,out,myid)
      INCLUDE 'fftw_f77.i'
#include "fortran.i"
c      integer*4 plan
c     integer*8 plan,iplan
      real*4 CPU(2),stime
      real work
      integer myid

      if(myid.eq.0) stime = etime(CPU)
      call rfftwnd_f77_mpi(plan,1,den,out,1,FFTW_TRANSPOSED_ORDER)
      if(myid .EQ. 0) print *,' CPU(F-FFT|out) ', etime(CPU)-stime
      return
      end
c--------------------------------------------------------------------
      subroutine fftbackwardaddspace(iplan,den,out,myid)
      INCLUDE 'fftw_f77.i'
c     integer*4 iplan
c     integer*8 plan,iplan
#include "fortran.i"
      real*4 CPU(2),stime
      real work
      integer myid

      if(myid.eq.0) stime = etime(CPU)
      call rfftwnd_f77_mpi(iplan,1,den,out,1,
     &               FFTW_TRANSPOSED_ORDER)
      if(myid .EQ. 0) print *,' CPU(B-FFT|out) ', etime(CPU)-stime
      return
      end
c--------------------------------------------------------------------
      subroutine fftwfinal()
      INCLUDE 'fftw_f77.i'
c     integer*4 plan,iplan
c     integer*8 plan,iplan
#include "fortran.i"
      integer nx,ny,nz,local_nz,local_z_start
      integer local_ny_after_transpose,local_y_start_after_transpose
      integer total_local_size
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      call rfftwnd_f77_mpi_destroy_plan(plan)
      call rfftwnd_f77_mpi_destroy_plan(iplan)

      return
      end
      subroutine fftwbroadcast(mx,my,mz,llocal_nz,llocal_z_start,
     &     llocal_ny_after_transpose,llocal_y_start_after_transpose,
     &     ttotal_local_size)
#include "fortran.i"
      integer nx,ny,nz,local_nz,local_z_start
      integer local_ny_after_transpose,local_y_start_after_transpose
      integer total_local_size
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      nx = mx
      ny = my
      nz = mz
      local_nz = llocal_nz
      local_z_start = llocal_z_start
      local_ny_after_transpose = llocal_ny_after_transpose
      local_y_start_after_transpose = llocal_y_start_after_transpose
      total_local_size = ttotal_local_size

      return
      end
