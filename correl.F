c234567
      subroutine correl(den,a,nowflag,now,rth,nstep,pk,npk,cor,tpk,tnpk,
     &       tcor)
      include "mpif.h"
#include "fortran.i"
      integer nx,ny,nz,local_nz,nowflag
      integer local_z_start,local_ny_after_transpose
      integer local_y_start_after_transpose,total_local_size,ierror
      integer myid,nid
      real den(2*(nx/2+1),nz,local_ny_after_transpose)
c     real pk(NNX+1),npk(NNX+1),cor(NNX+1)
c     real tpk(NNX+1),tnpk(NNX+1),tcor(NNX+1)
      real pk(*),npk(*),cor(*)
      real tpk(*),tnpk(*),tcor(*)
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
      do 90 i=1,ngx+1
      pk(i) = 0.
      npk(i) = 0
      tpk(i) = 0.
      tnpk(i) = 0
   90 cor(i) = 0.

      do 100 j = 1, local_ny_after_transpose, +1
      wj = xyy*(j - 1 + local_y_start_after_transpose)
      if(wj.gt.ngh) wj = wj - ngx
	  if(wj .eq. 0) then
	  tscdeconvy = 1
	  else
      tscdeconvy = (sin(pi*wj/ng)/(pi*wj/ng))**3
	  endif
      do 100 k = 1, ngz, +1
      wk = xyz*(k - 1)
      if(wk.gt.ngh) wk = wk - ngx
	  if(wk.eq.0) then
	  tscdeconvz = 1
	  else
      tscdeconvz = (sin(pi*wk/ng)/(pi*wk/ng))**3
	  endif
      do 100 i = 0, ngh, +1  
      wave = sqrt(i*i+wj*wj+wk*wk)
c---  density array in k-space is transposed order: (i,j,k)-->(i,k,j)
      if(i.eq.0) then
	  tscdeconvx = 1
	  else
      tscdeconvx = (sin(pi*real(i)/ng)/(pi*real(i)/ng))**3
	  endif
      delksq = (den(2*i+1,k,j)**2+den(2*i+2,k,j)**2)*rngm6
	  tscdeconv = tscdeconvx * tscdeconvy * tscdeconvz
	  delksq = delksq /tscdeconv/tscdeconv
      kbin1 = nint(wave)+1
      wt1 = wave-kbin1+1
      kbin2 = kbin1+sign(1.,wt1)
      wt1 = abs(wt1)
      if(kbin1 .le. nx+1) then
         pk(kbin1) = pk(kbin1)+(1.-wt1)*delksq
         npk(kbin1) = npk(kbin1)+1.-wt1
      endif
      if(kbin2 .le. nx+1) then
         pk(kbin2) = pk(kbin2)+wt1*delksq
         npk(kbin2) = npk(kbin2)+wt1
      endif
  100 continue

c--   Gathering pk in all processors and summing
      call MPI_REDUCE(pk,tpk,nx+1,MPI_REAL,MPI_SUM,0,
     &           MPI_COMM_WORLD,ierror)
      call MPI_REDUCE(npk,tnpk,nx+1,MPI_REAL,MPI_SUM,0,
     &           MPI_COMM_WORLD,ierror)
      if(myid.eq.0) then
         do i = 1, nx+1
            pk(i) = tpk(i)
            npk(i) = tnpk(i)
         enddo
      endif
      call MPI_BCAST(pk,nx+1,MPI_REAL,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(npk,nx+1,MPI_REAL,0,MPI_COMM_WORLD,ierror)
      
c---  only calculation on root(0) processor.
      if(myid.eq.0) then
      std = 0.
      fact = 2.*pi*rth/rngx
      do 110 i=2,ngx+1
      if(npk(i).gt.0)then
         pk(i)=pk(i)/npk(i)
      else
         pk(i) = 0.
      end if
      thk = fact*(i-1.)
  110 std = std+pk(i)*(i-1.)**2*(sin(thk)-thk*cos(thk))**2/thk**6
      std = sqrt(4.*pi*std*9.)
      do 130 i=2,ngx+1
      seper = 2.*pi*(i-1.)/rngx
      do 120 k=2,ngx+1
  120 cor(i) = cor(i)+sin((k-1.)*seper)*(k-1.)*pk(k)
  130 cor(i) = 2.*rngx*cor(i)/(i-1.)
      if(cor(4).gt.0..and.cor(3).gt.0.)then
      corn2 = (alog(cor(4))-alog(cor(3)))/(alog(4.)-alog(3.))
      if(nowflag.eq.1)then
      if(a.ne.1..and.(corn1+1.8)*(corn2+1.8).le.0.)now=1
      if(now.eq.1) then
         print *,' Finishing at a=',a
         call MPI_FINALIZE(ierror)
         stop
      endif
      end if
      corn1 = corn2
      end if
#ifdef AIX
      open(7,file='pmevol.map',position='append')
#else
      open(7,file='pmevol.map',access='append')
#endif
      write(7,*)a,std,corn1
      do i = 1, ngh
         write(7,*) i,i-1,pk(i),cor(i)
      enddo
      close(7)
      endif
      return
      end

