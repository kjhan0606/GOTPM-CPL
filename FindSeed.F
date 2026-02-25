!234567
      subroutine findiseed(den,iseed)
!  only part of the power spectrum is needed to find the well-hehaving random number.
      parameter(NN=100)
      include 'mpif.h'
#include "fortran.i"
      integer nx,ny,nz,local_nz
      integer local_z_start,local_ny_after_transpose
      integer local_y_start_after_transpose,total_local_size,ierror
      integer myid,nid
      real den(2*(nx/2+1),ny,local_nz)
      real pk(NN+1),npk(NN+1),cor(NN+1),tpk(NN+1),tnpk(NN+1),tcor(NN+1)
      real corn1,a
      double precision rngm6,rngc
      save corn1
      real xyz, xyy
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size

      call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierror)
      call MPI_Comm_size(MPI_COMM_WORLD,nid,ierror)

      rngc = dble(nx)*dble(ny)*dble(nz)
      rngm6 = 1./rngc/rngc
      xyz = real(nx)/real(nz)
      pi = 3.141592654
      ng = nx
      rng = ng
      ngx = nx
      ngy = ny
      ngz = nz
      rngx = ngx
      rngy = ngy
      rngz = ngz
      xyz = rng/rngz
      xyy = rng/rngy
      ngxh = ngx/2
      ngh = ngxh
      ngyh = ngy/2
      ngzh = ngz/2
      do 90 i=1,NN+1
      pk(i) = 0.
      npk(i) = 0
      tpk(i) = 0.
      tnpk(i) = 0
   90 cor(i) = 0.

      do 100 k = 1, local_nz, +1
      wk = xyz*(k - 1 + local_z_start)
      if(wj.gt.ngh) wj = wj - ngx
      do 100 j = 1, ny, +1
      wj = xyy*(j - 1)
      if(wk.gt.ngh) wk = wk - ngx
      do 100 i = 0, ngh, +1  
      wave = sqrt(i*i+wj*wj+wk*wk)
c---  density array in k-space is transposed order: (i,j,k)-->(i,k,j)
      delksq = (den(2*i+1,j,k)**2+den(2*i+2,j,k)**2)/2.
      kbin1 = nint(wave)+1
      wt1 = wave-kbin1+1
      kbin2 = kbin1+sign(1.,wt1)
      wt1 = abs(wt1)
      if(kbin1 .le. NN+1) then
         pk(kbin1) = pk(kbin1)+(1.-wt1)*delksq
         npk(kbin1) = npk(kbin1)+1.-wt1
      endif
      if(kbin2 .le. NN+1) then
         pk(kbin2) = pk(kbin2)+wt1*delksq
         npk(kbin2) = npk(kbin2)+wt1
      endif
  100 continue

c--   Gathering pk in all processors and summing
      call MPI_REDUCE(pk,tpk,NN+1,MPI_REAL,MPI_SUM,0,
     &           MPI_COMM_WORLD,ierror)
      call MPI_REDUCE(npk,tnpk,NN+1,MPI_REAL,MPI_SUM,0,
     &           MPI_COMM_WORLD,ierror)
      if(myid.eq.0) then
         do i = 1, NN+1
            pk(i) = tpk(i)
            npk(i) = tnpk(i)
         enddo
      endif
      call MPI_BCAST(pk,NN+1,MPI_REAL,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(npk,NN+1,MPI_REAL,0,MPI_COMM_WORLD,ierror)
      
c---  only calculation on root(0) processor.
      if(myid.eq.0) then

      do 110 i=2,NN+1
      if(npk(i).gt.0)then
         pk(i)=pk(i)/npk(i)
      else
         pk(i) = 0.
      end if
  110 enddo

      do 130 i=2,NN+1
      seper = 2.*pi*(i-1.)/rngx
      do 120 k=2,ngx+1
  120 cor(i) = cor(i)+sin((k-1.)*seper)*(k-1.)*pk(k)
  130 cor(i) = 2.*rngx*cor(i)/(i-1.)


#ifdef AIX
      open(7,file='iseed.map',position='append')
#else
      open(7,file='iseed.map',access='append')
#endif
      write(7,*) iseed
      do i = 1, NN/2
         write(7,*) i,i-1,pk(i),cor(i)
      enddo
      close(7)
      end if
      return
      end
