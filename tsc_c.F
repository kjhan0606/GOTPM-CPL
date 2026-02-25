c234567
      subroutine tsc(den,local_mz,nstart_z,
     &               particles,nn,nxpos,nend,pmas,myrank,nrank)
      IMPLICIT none    ! This statement is very important
      include 'mpif.h'
#include "fortran.i"
      integer NZWIDTH
      parameter(NZWIDTH=4) ! since 1 for lower, 2 for upper, and 1 for buffer for shfiting
      integer nx,ny,nz, nend,nstart_z,local_nz,local_mz
      integer nxpos
      real start_z
      real pmas, pmas0, p05
      integer nn
      real particles(nn,nend)
c     real px(*),py(*),pz(*)
      integer xsign,ysign,zsign
      integer myrank,nrank
      integer nz_per
      integer i,j,k,ii
      integer nearx,neary,nearz,i1,j1,k1,i2,j2,k2,i3,j3,k3
      real xp,yp,zp,xmin,ymin,zmin
      real wx1,wy1,wz1,wx2,wy2,wz2,wx3,wy3,wz3
      real xd1,yd1,zd1
      real rngc
      integer usndsize,dgetsize,dsndsize,ugetsize
      integer dest,src
      integer status(MPI_STATUS_SIZE),ierror,stag,rtag,tag
      real*4 CPU(2), etime,stime
      integer local_z_start, 
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      real den(2*(nx/2+1),ny,local_mz+NZWIDTH)

      
ccc   nz_per is the z size of local array to pretend periodic boundary conditions
ccc   in the calculation of local density.

      rngc = real(nx)*real(ny)*real(nz)
      if(myrank.eq.0) stime = etime(CPU)
      nz_per = local_mz+NZWIDTH
#ifdef VarPM
c$omp parallel do
      do 10 k = 1, local_mz+NZWIDTH
      do 10 j = 1, ny
      do 10 i = 1, 2*(nx/2+1)
         den(i,j,k) = 0.
10    continue
#else
c$omp parallel do
      do 10 k = 1, local_mz
      do 10 j = 1, ny
      do 10 i = 1, nx
         den(i,j,k) = -1.
10    continue
c$omp parallel do
      do 20 k = local_mz+1, local_mz+NZWIDTH, +1
      do 20 j = 1, ny
      do 20 i = 1, nx
         den(i,j,k) = 0.
20    continue
#endif
      start_z = nstart_z



      do i = 1, nend
#if defined(MULTIMASS)
         pmas0 = particles(nxpos,i)*rngc
#else
         pmas0 = pmas
#endif
         p05 = pmas0*0.5
         xp = particles(nxpos+1,i)
         yp = particles(nxpos+2,i)
         zp = particles(nxpos+3,i)-start_z
         nearx = nint(xp)
         neary = nint(yp)
         nearz = nint(zp)
         xmin = xp-nearx
         ymin = yp-neary
         zmin = zp-nearz
ccc Don`t use isign because this is the internal function of Fortran.
ccc So x,y,zsign were used instead.
         xsign = sign(1.,xmin)
         ysign = sign(1.,ymin)
         zsign = sign(1.,zmin)
         i1 = mod(nearx+nx,nx) + 1
         j1 = mod(neary+ny,ny) + 1
         k1 = mod(nearz+nz_per,nz_per) + 1
         i2 = mod(nearx+nx+xsign,nx) + 1
         j2 = mod(neary+ny+ysign,ny) + 1
         k2 = mod(nearz+nz_per+zsign,nz_per) + 1
         i3 = mod(nearx+nx-xsign,nx) + 1
         j3 = mod(neary+ny-ysign,ny) + 1
         k3 = mod(nearz+nz_per-zsign,nz_per) + 1
         xd1 = abs(xmin)
         yd1 = abs(ymin)
         zd1 = abs(zmin)
c        print *,k1,k2,k3,zd1
c        stop
         wx1 =(0.75-xd1**2)*pmas0
         wy1 = 0.75-yd1**2
         wz1 = 0.75-zd1**2
         wx2 = p05*(0.5+xd1)**2
         wy2 = 0.5*(0.5+yd1)**2
         wz2 = 0.5*(0.5+zd1)**2
         wx3 = p05*(0.5-xd1)**2
         wy3 = 0.5*(0.5-yd1)**2
         wz3 = 0.5*(0.5-zd1)**2
c---  Assign mass to the grid the particle overlaps
         den(i1,j1,k1) = den(i1,j1,k1) + wx1*wy1*wz1
         den(i2,j1,k1) = den(i2,j1,k1) + wx2*wy1*wz1
         den(i3,j1,k1) = den(i3,j1,k1) + wx3*wy1*wz1
         den(i1,j2,k1) = den(i1,j2,k1) + wx1*wy2*wz1
         den(i2,j2,k1) = den(i2,j2,k1) + wx2*wy2*wz1
         den(i3,j2,k1) = den(i3,j2,k1) + wx3*wy2*wz1
         den(i1,j3,k1) = den(i1,j3,k1) + wx1*wy3*wz1
         den(i2,j3,k1) = den(i2,j3,k1) + wx2*wy3*wz1
         den(i3,j3,k1) = den(i3,j3,k1) + wx3*wy3*wz1
         den(i1,j1,k2) = den(i1,j1,k2) + wx1*wy1*wz2
         den(i2,j1,k2) = den(i2,j1,k2) + wx2*wy1*wz2
         den(i3,j1,k2) = den(i3,j1,k2) + wx3*wy1*wz2
         den(i1,j2,k2) = den(i1,j2,k2) + wx1*wy2*wz2
         den(i2,j2,k2) = den(i2,j2,k2) + wx2*wy2*wz2
         den(i3,j2,k2) = den(i3,j2,k2) + wx3*wy2*wz2
         den(i1,j3,k2) = den(i1,j3,k2) + wx1*wy3*wz2
         den(i2,j3,k2) = den(i2,j3,k2) + wx2*wy3*wz2
         den(i3,j3,k2) = den(i3,j3,k2) + wx3*wy3*wz2
         den(i1,j1,k3) = den(i1,j1,k3) + wx1*wy1*wz3
         den(i2,j1,k3) = den(i2,j1,k3) + wx2*wy1*wz3
         den(i3,j1,k3) = den(i3,j1,k3) + wx3*wy1*wz3
         den(i1,j2,k3) = den(i1,j2,k3) + wx1*wy2*wz3
         den(i2,j2,k3) = den(i2,j2,k3) + wx2*wy2*wz3
         den(i3,j2,k3) = den(i3,j2,k3) + wx3*wy2*wz3
         den(i1,j3,k3) = den(i1,j3,k3) + wx1*wy3*wz3
         den(i2,j3,k3) = den(i2,j3,k3) + wx2*wy3*wz3
         den(i3,j3,k3) = den(i3,j3,k3) + wx3*wy3*wz3
      enddo
ccccc These lines send upper swap den(den(*,*,local_mz+1~2)) to a upper node.
ccccc Be careful! The small size of buffer resources makes the sending
ccccc process with one column! Not one slice!!!!! 
ccccc dest : send array to dest`th node
ccccc src : receivce array from src`th node
#ifdef VarPM
c---  rearrange the tsc mesh  in normal ordering for slicemigrate processes.
      do 30 k =  local_mz+2, 1 , -1
      do 30 j = 1, ny
      do 30 i = 1, nx
         den(i,j,k+1) = den(i,j,k)
30    continue
      k =  local_mz+NZWIDTH 
      do 40 j = 1, ny
      do 40 i = 1, nx
         den(i,j,1) = den(i,j,k)
40    continue
#else  
      dest = mod(myrank+1+nrank,nrank)
      src = mod(myrank-1+nrank,nrank)
      stag = 0
      rtag = 0
      usndsize = 2*(nx/2+1)*ny*2
      dgetsize = usndsize
      call MPI_SENDRECV_REPLACE(den(1,1,local_mz+1),usndsize,
     &   MPI_REAL,dest,stag,src,rtag,MPI_COMM_WORLD,status,ierror)
      do 30 k = 1, 2
      do 30 j = 1, ny
      do 30 i = 1, nx
         den(i,j,k) = den(i,j,k) + den(i,j,local_mz+k)
30    continue
ccccc These lines send lower swap(den(*,*,local_mz+NZWIDTH)) to a lower slice.
ccccc dest : send array to dest'th node
ccccc src : receivce array from src'th node
      dsndsize = 2*(nx/2+1)*ny
      ugetsize = dsndsize
      dest = mod(myrank-1+nrank,nrank)
      src = mod(myrank+1+nrank,nrank)
      stag = 0
      rtag = 0
      call MPI_SENDRECV_REPLACE(den(1,1,local_mz+NZWIDTH),dsndsize,
     &       MPI_REAL,dest,stag,src,rtag,MPI_COMM_WORLD,status,ierror)
      do 40 j = 1, ny
      do 40 i = 1, nx
         den(i,j,local_mz) = den(i,j,local_mz) + 
     &                               den(i,j,local_mz+NZWIDTH)
40    continue
#endif
c--------------------------------
c--------------------------------
      if(myrank.eq.0) print *,' CPU(TSC) ', etime(CPU)-stime
      return
      end
