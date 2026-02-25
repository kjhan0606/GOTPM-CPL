c234567
      subroutine psolver(den,poten,efold)
      IMPLICIT none
      include "mpif.h"
#include "fortran.i"
      real pi,fact,efold
      parameter(pi=3.1415926535,fact= -1./pi)
      integer nx,ny,nz,local_nz,local_z_start
      integer local_ny_after_transpose,local_y_start_after_transpose
      integer total_local_size
      real den(2*(nx/2+1),nz,local_ny_after_transpose)
      real  sfact
      real*8 rngc,rngm6,invrngc
      integer nxh,wy,wz,myid,nid,ierror
      real wave,factr,factor,delksq,rknyq2,poten,tpoten
      double precision poten1,poten2
      integer i,j,k,ii
      real*4 CPU(2),stime,etime
      integer ngx,ngy,ngz,nghx,nghy,nghz
      real xyz,xyy, pcorr
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierror)
      call MPI_Comm_size(MPI_COMM_WORLD,nid,ierror)


      rngc = dble(nx)*dble(ny)*dble(nz)
      invrngc = 1.d0/rngc
      rngm6 = 1.d0/rngc/rngc
      xyz = real(nx)/real(nz)
      xyy = real(nx)/real(ny)

      if(myid.eq.0) stime = etime(CPU)
      poten1 = 0.D0
      poten2 = 0.D0
      sfact = -(pi*efold/real(nx))**2
      nxh = nx/2
      rknyq2 = nxh**2
c$omp parallel do private(wy,wz,wave,factr,factor,delksq)
c$omp+ reduction(+:poten1,poten2)
      do j = 1, local_ny_after_transpose, +1
      wy = xyy*(j - 1 + local_y_start_after_transpose)
      if(wy.gt.nxh) wy = wy - nx
      do k = 1, nz, +1
      wz = xyz*(k - 1)
      if(wz.gt.nxh) wz = wz - nx
      do i = 0, nxh, +1  
         wave = i*i+wy*wy+wz*wz
         if(wave .EQ. 0. ) then
            den(2*i+1,k,j) = 0.
            den(2*i+2,k,j) = 0.
         else
            factr = fact/wave
            factor = factr*exp(wave*sfact)
c---  Division by n**3 is due to the length unit convention (dx^3=1/n**3)
            delksq = (den(2*i+1,k,j)**2+den(2*i+2,k,j)**2)*rngm6
            if(i.eq.0)then
               poten1 = poten1+factor*delksq
            else
               poten2 = poten2+factor*delksq
            end if
            den(2*i+1,k,j) = factor*den(2*i+1,k,j)*invrngc
            den(2*i+2,k,j) = factor*den(2*i+2,k,j)*invrngc
         endif
      enddo
      enddo
      enddo
      poten = 0.5*(poten1+2.*poten2)
      call MPI_REDUCE(poten,tpoten,1,MPI_REAL,MPI_SUM,0,
     &          MPI_COMM_WORLD,ierror)
      poten = tpoten
      call MPI_BCAST(poten,1,MPI_REAL,0,MPI_COMM_WORLD,ierror)

      if(myid.eq.0) print *,' CPU(Poisson) ',etime(cpu)-stime

      return
      end
