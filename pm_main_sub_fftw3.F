c234567
      subroutine readpar(nspace,omep,omeplam,rth)
      integer nspace
      real omep,omeplam,rth

      print *, "nspace,omep,omeplam,rth ? "
      read(7,*) nspace,omep,omeplam,rth
      return
      end

c---  FFTW USE_WORK
      subroutine fftwbroadcast(mx,my,mz,local_mz,local_zz_start,
     &        local_my, local_yy_start)
      INCLUDE 'fftw_f77.i'
      INCLUDE 'mpif.h'
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
      local_nz = local_mz
      local_z_start = local_zz_start
      local_ny_after_transpose = local_my
      local_y_start_after_transpose = local_yy_start

      return
      end
