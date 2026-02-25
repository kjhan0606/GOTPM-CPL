c234567
      subroutine fda4(particles,nn,nxpos,np,den,local_mz,
     &          a,astep,nextastep,amax,omei,omep,omeplam,pfact,
     &          nstart_z,fact1,fact2,ihflag1,ihflag2)
      IMPLICIT none
      include 'mpif.h'
#include "fortran.i"
      integer ihflag1,ihflag2
      real pi,rngx,rngy,rngz,rng,amax,omei,omeplam,omep
      parameter(pi=3.1415926535)
      integer np,nx,ny,nz,local_nz,nstart_z,local_mz
      integer nxpos
      real start_z
      integer nn
      real particles(nn,np)
c     real x(*),y(*),z(*),vx(*),vy(*),vz(*)
      real den(2*(nx/2+1),ny,local_nz+7)
      real astep,a,ainv,app,apsq,afact,bfact,fact1,fact2,fact3,pfact
      real fact1h,fact2h
      real nextastep
      real ffact1,ffact2
      real xp,yp,zp,xmin,ymin,zmin,xd1,yd1,zd1
      real wx1,wy1,wz1,wx2,wy2,wz2,wx3,wy3,wz3
      real wt1,wt2,wt3,wt4,wt5,wt6,wt7,wt8,wt9,wt10
      real wt11,wt12,wt13,wt14,wt15,wt16,wt17,wt18,wt19,wt20
      real wt21,wt22,wt23,wt24,wt25,wt26,wt27
      real fx,fy,fz
      real fx1,fx2,fx3,fx4,fx5,fx6,fx7,fx8,fx9,fx10
      real fx11,fx12,fx13,fx14,fx15,fx16,fx17,fx18,fx19,fx20
      real fx21,fx22,fx23,fx24,fx25,fx26,fx27
      real fy1,fy2,fy3,fy4,fy5,fy6,fy7,fy8,fy9,fy10
      real fy11,fy12,fy13,fy14,fy15,fy16,fy17,fy18,fy19,fy20
      real fy21,fy22,fy23,fy24,fy25,fy26,fy27
      real fz1,fz2,fz3,fz4,fz5,fz6,fz7,fz8,fz9,fz10
      real fz11,fz12,fz13,fz14,fz15,fz16,fz17,fz18,fz19,fz20
      real fz21,fz22,fz23,fz24,fz25,fz26,fz27
      integer nearx,neary,nearz,xsign,ysign,zsign
      integer i1,j1,k1,i2,j2,k2,i3,j3,k3
      integer i,ng,ngx,ngy,ngz

      integer xabp0(-10:NNX+10),xabp1(-10:NNX+10),xabm1(-10:NNX+10)
      integer yabp0(-10:NNY+10),yabp1(-10:NNY+10),yabm1(-10:NNY+10)
      integer zabp0(-10:LOCAL_NNZ+7+10),zabp1(-10:LOCAL_NNZ+7+10)
      integer zabm1(-10:LOCAL_NNZ+7+10)
      common /ptables/ xabp0,xabp1,xabm1,yabp0,yabp1,yabm1,zabp0,
     &       zabp1,zabm1

      integer i1m1,i1m2,i1p1,i1p2,i2m1,i2m2,i2p1,i2p2,i3m1,i3m2,i3p1
      integer i3p2,j1m1,j1m2,j1p1,j1p2,j2m1,j2m2,j2p1,j2p2,j3m1,j3m2
      integer j3p1,j3p2,k1m1,k1m2,k1p1,k1p2,k2m1,k2m2,k2p1,k2p2,k3m1
      integer k3m2,k3p1,k3p2
 
      integer nz_per,dest,src,sendsize,recvsize,stag,rtag
      integer status(MPI_STATUS_SIZE)
      integer myid,nid,ierror
      real*4 CPU(2),stime,etime
      integer ihalf
      integer local_z_start, 
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size

      call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierror)
      call MPI_Comm_size(MPI_COMM_WORLD,nid,ierror)



      if(myid .eq. 0) stime = etime(CPU)
ccc   nz_per is the z size of local array to pretend periodic 
ccc   boundary conditions in the calculation of local density.
      nz_per = local_mz+7
#ifndef VarPM
c** In the option of VarPM, the mesh is rearranged in the slicemigrate processes.
c***  Gather data of buffer den(*,*,local_mz+1:local_mz+7)
      sendsize = 2*(nx/2+1)*ny*3
      recvsize = sendsize
      dest = mod(myid+1+nid,nid)
      src = mod(myid-1+nid,nid)
      stag = 1
      rtag = 1
      call MPI_SENDRECV(den(1,1,local_mz-2),sendsize,MPI_REAL,dest,
     &        stag,den(1,1,nz_per-2),recvsize,MPI_REAL,src,rtag,
     &                 MPI_COMM_WORLD,status,ierror)
      sendsize = 2*(nx/2+1)*ny*4
      recvsize = sendsize
      dest = mod(myid-1+nid,nid)
      src = mod(myid+1+nid,nid)
      stag = 2
      rtag = 2
      call MPI_SENDRECV(den(1,1,1),sendsize,MPI_REAL,dest,stag,
     &          den(1,1,local_mz+1),recvsize,MPI_REAL,src,rtag,
     &                  MPI_COMM_WORLD,status,ierror)
#endif



c***  Solve the momentum equation by using Friedman equations
      ng = nx
      ngx = nx
      ngy = ny
      ngz = nz
      rng = ng
      rngx = ngx
      rngy = ngy
      rngz = ngz
      if(a.eq.1)astep = 0.5*astep
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ainv = 1./a
      app  = -4.*pi/3.*(ainv*ainv-2*omeplam/omep*a/amax**3)
      apsq = 8.*pi/3.*(ainv+1./omei-1+omeplam/omep*(a**2-1)/amax**3)
      afact = (2.+app*a/apsq)/2.*ainv
      bfact = ainv**3/apsq
      fact1 = (1.-afact*astep)/(1.+afact*astep)
      fact1h = (1.-afact*astep*0.5)/(1.+afact*astep*0.5)
      fact2 = bfact*astep/(1.+afact*astep)*rng
      fact2h = bfact*astep*0.5/(1.+afact*astep*0.5)*rng
      fact3 = fact2/fact1
c---  first time for the half step integration
      if(ihflag1 .eq. 1 .and. ihflag2 .eq.0) then
         fact1 = fact1h
         fact2 = fact2h
c---  second time for the half step integration
      else if(ihflag1 .eq.1 .and. ihflag2 .eq.1) then
         fact1 = fact1/fact1h ! Notice the order of fact1 and fact2 
         fact2 = (fact2-fact1*fact2h)
      endif


      if(a.eq.1)astep = 2.*astep
c---
      nextastep = astep
c---
      pfact = rng*astep

c     velocity update only accounting for the hubble expansion
c$omp parallel do
      do i = 1, np
         particles(nxpos+4,i) = fact1*particles(nxpos+4,i)
         particles(nxpos+5,i) = fact1*particles(nxpos+5,i)
         particles(nxpos+6,i) = fact1*particles(nxpos+6,i)
      enddo

c***  Force at grid points by 2-point FDA
      ffact1 = (4./3.)/2.
      ffact2 = (-1./3.)/4.

c***  Cell position run from 0 to n-1 and particles from 0 to n
c---  (xp,yp,zp) = particle position
c---  (ic,jc,kc) = major cell`s position the superparticle clould belongs to.
c---  (icc,jcc,kcc) = position of the minor cell
c---  Now find the velocity at a+astep/2
c---  CPU : Scalar (4.62), Vector (1.53,12.2,96) for ng = 32,64,128
      start_z = nstart_z
c$omp parallel do private(xp,yp,zp,nearx,neary,nearz,xmin,ymin,zmin)
c$omp+ private(xsign,ysign,zsign,i1,j1,k1,i2,j2,k2,i3,j3,k3)
c$omp+ private(xd1,yd1,zd1,wx1,wy1,wz1,wx2,wy2,wz2)
c$omp+ private(wx3,wy3,wz3,wt1,wt2,wt3,wt4,wt5,wt6,wt7,wt8)
c$omp+ private(wt9,wt10,wt11,wt12,wt13,wt14,wt15,wt16)
c$omp+ private(wt17,wt18,wt19,wt20,wt21,wt22,wt23,wt24)
c$omp+ private(wt25,wt26,wt27)
c$omp+ private(i1m1,i1m2,i1p1,i1p2,i2m1,i2m2,i2p1,i2p2)
c$omp+ private(i3m1,i3m2,i3p1,i3p2)
c$omp+ private(j1m1,j1m2,j1p1,j1p2,j2m1,j2m2,j2p1,j2p2)
c$omp+ private(j3m1,j3m2,j3p1,j3p2)
c$omp+ private(k1m1,k1m2,k1p1,k1p2,k2m1,k2m2,k2p1,k2p2)
c$omp+ private(k3m1,k3m2,k3p1,k3p2)
c$omp+ private(fx1,fx2,fx3,fx4,fx5,fx6,fx7,fx8,fx9,fx10)
c$omp+ private(fx11,fx12,fx13,fx14,fx15,fx16,fx17,fx18,fx19,fx20)
c$omp+ private(fx21,fx22,fx23,fx24,fx25,fx26,fx27)
c$omp+ private(fy1,fy2,fy3,fy4,fy5,fy6,fy7,fy8,fy9,fy10)
c$omp+ private(fy11,fy12,fy13,fy14,fy15,fy16,fy17,fy18,fy19,fy20)
c$omp+ private(fy21,fy22,fy23,fy24,fy25,fy26,fy27)
c$omp+ private(fz1,fz2,fz3,fz4,fz5,fz6,fz7,fz8,fz9,fz10)
c$omp+ private(fz11,fz12,fz13,fz14,fz15,fz16,fz17,fz18,fz19,fz20)
c$omp+ private(fz21,fz22,fz23,fz24,fz25,fz26,fz27,fx,fy,fz)
      do 160 i=1,np
      xp   = particles(nxpos+1,i)
      yp   = particles(nxpos+2,i)
      zp   = particles(nxpos+3,i)-start_z
      nearx= nint(xp)
      neary= nint(yp)
      nearz= nint(zp)
      xmin = xp-nearx
      ymin = yp-neary
      zmin = zp-nearz
      xsign= sign(1.,xmin)
      ysign= sign(1.,ymin)
      zsign= sign(1.,zmin)
      i1  = mod(nearx+ngx,ngx)+1
      j1  = mod(neary+ngy,ngy)+1
      k1  = mod(nearz+nz_per,nz_per)+1
      i2  = mod(nearx+ngx+xsign,ngx)+1
      j2  = mod(neary+ngy+ysign,ngy)+1
      k2  = mod(nearz+nz_per+zsign,nz_per)+1
      i3  = mod(nearx+ngx-xsign,ngx)+1
      j3  = mod(neary+ngy-ysign,ngy)+1
      k3  = mod(nearz+nz_per-zsign,nz_per)+1
      xd1 = abs(xmin)
      yd1 = abs(ymin)
      zd1 = abs(zmin)
      wx1 = 0.75-xd1**2
      wy1 = 0.75-yd1**2
      wz1 = 0.75-zd1**2
      wx2 = 0.5*(0.5+xd1)**2
      wy2 = 0.5*(0.5+yd1)**2
      wz2 = 0.5*(0.5+zd1)**2
      wx3 = 0.5*(0.5-xd1)**2
      wy3 = 0.5*(0.5-yd1)**2
      wz3 = 0.5*(0.5-zd1)**2
      wt1  = wx1*wy1*wz1
      wt2  = wx2*wy1*wz1
      wt3  = wx3*wy1*wz1
      wt4  = wx1*wy2*wz1
      wt5  = wx2*wy2*wz1
      wt6  = wx3*wy2*wz1
      wt7  = wx1*wy3*wz1
      wt8  = wx2*wy3*wz1
      wt9  = wx3*wy3*wz1
      wt10 = wx1*wy1*wz2
      wt11 = wx2*wy1*wz2
      wt12 = wx3*wy1*wz2
      wt13 = wx1*wy2*wz2
      wt14 = wx2*wy2*wz2
      wt15 = wx3*wy2*wz2
      wt16 = wx1*wy3*wz2
      wt17 = wx2*wy3*wz2
      wt18 = wx3*wy3*wz2
      wt19 = wx1*wy1*wz3
      wt20 = wx2*wy1*wz3
      wt21 = wx3*wy1*wz3
      wt22 = wx1*wy2*wz3
      wt23 = wx2*wy2*wz3
      wt24 = wx3*wy2*wz3
      wt25 = wx1*wy3*wz3
      wt26 = wx2*wy3*wz3
      wt27 = wx3*wy3*wz3
      i1m1 = xabm1(i1)
      i1m2 = xabm1(i1m1)
      i1p1 = xabp1(i1)
      i1p2 = xabp1(i1p1)
      i2m1 = xabm1(i2)
      i2m2 = xabm1(i2m1)
      i2p1 = xabp1(i2)
      i2p2 = xabp1(i2p1)
      i3m1 = xabm1(i3)
      i3m2 = xabm1(i3m1)
      i3p1 = xabp1(i3)
      i3p2 = xabp1(i3p1)

      j1m1 = yabm1(j1)
      j1m2 = yabm1(j1m1)
      j1p1 = yabp1(j1)
      j1p2 = yabp1(j1p1)
      j2m1 = yabm1(j2)
      j2m2 = yabm1(j2m1)
      j2p1 = yabp1(j2)
      j2p2 = yabp1(j2p1)
      j3m1 = yabm1(j3)
      j3m2 = yabm1(j3m1)
      j3p1 = yabp1(j3)
      j3p2 = yabp1(j3p1)

      k1m1 = zabm1(k1)
      k1m2 = zabm1(k1m1)
      k1p1 = zabp1(k1)
      k1p2 = zabp1(k1p1)
      k2m1 = zabm1(k2)
      k2m2 = zabm1(k2m1)
      k2p1 = zabp1(k2)
      k2p2 = zabp1(k2p1)
      k3m1 = zabm1(k3)
      k3m2 = zabm1(k3m1)
      k3p1 = zabp1(k3)
      k3p2 = zabp1(k3p1)


      fx1 = ffact1*(den(i1m1,j1,k1)-den(i1p1,j1,k1))+
     &      ffact2*(den(i1m2,j1,k1)-den(i1p2,j1,k1))
      fx2 = ffact1*(den(i2m1,j1,k1)-den(i2p1,j1,k1))+
     &      ffact2*(den(i2m2,j1,k1)-den(i2p2,j1,k1))
      fx3 = ffact1*(den(i3m1,j1,k1)-den(i3p1,j1,k1))+
     &      ffact2*(den(i3m2,j1,k1)-den(i3p2,j1,k1))
      fx4 = ffact1*(den(i1m1,j2,k1)-den(i1p1,j2,k1))+
     &      ffact2*(den(i1m2,j2,k1)-den(i1p2,j2,k1))
      fx5 = ffact1*(den(i2m1,j2,k1)-den(i2p1,j2,k1))+
     &      ffact2*(den(i2m2,j2,k1)-den(i2p2,j2,k1))
      fx6 = ffact1*(den(i3m1,j2,k1)-den(i3p1,j2,k1))+
     &      ffact2*(den(i3m2,j2,k1)-den(i3p2,j2,k1))
      fx7 = ffact1*(den(i1m1,j3,k1)-den(i1p1,j3,k1))+
     &      ffact2*(den(i1m2,j3,k1)-den(i1p2,j3,k1))
      fx8 = ffact1*(den(i2m1,j3,k1)-den(i2p1,j3,k1))+
     &      ffact2*(den(i2m2,j3,k1)-den(i2p2,j3,k1))
      fx9 = ffact1*(den(i3m1,j3,k1)-den(i3p1,j3,k1))+
     &      ffact2*(den(i3m2,j3,k1)-den(i3p2,j3,k1))

      fx10 = ffact1*(den(i1m1,j1,k2)-den(i1p1,j1,k2))+
     &       ffact2*(den(i1m2,j1,k2)-den(i1p2,j1,k2))
      fx11 = ffact1*(den(i2m1,j1,k2)-den(i2p1,j1,k2))+
     &       ffact2*(den(i2m2,j1,k2)-den(i2p2,j1,k2))
      fx12 = ffact1*(den(i3m1,j1,k2)-den(i3p1,j1,k2))+
     &       ffact2*(den(i3m2,j1,k2)-den(i3p2,j1,k2))
      fx13 = ffact1*(den(i1m1,j2,k2)-den(i1p1,j2,k2))+
     &       ffact2*(den(i1m2,j2,k2)-den(i1p2,j2,k2))
      fx14 = ffact1*(den(i2m1,j2,k2)-den(i2p1,j2,k2))+
     &       ffact2*(den(i2m2,j2,k2)-den(i2p2,j2,k2))
      fx15 = ffact1*(den(i3m1,j2,k2)-den(i3p1,j2,k2))+
     &       ffact2*(den(i3m2,j2,k2)-den(i3p2,j2,k2))
      fx16 = ffact1*(den(i1m1,j3,k2)-den(i1p1,j3,k2))+
     &       ffact2*(den(i1m2,j3,k2)-den(i1p2,j3,k2))
      fx17 = ffact1*(den(i2m1,j3,k2)-den(i2p1,j3,k2))+
     &       ffact2*(den(i2m2,j3,k2)-den(i2p2,j3,k2))
      fx18 = ffact1*(den(i3m1,j3,k2)-den(i3p1,j3,k2))+
     &       ffact2*(den(i3m2,j3,k2)-den(i3p2,j3,k2))

      fx19 = ffact1*(den(i1m1,j1,k3)-den(i1p1,j1,k3))+
     &       ffact2*(den(i1m2,j1,k3)-den(i1p2,j1,k3))
      fx20 = ffact1*(den(i2m1,j1,k3)-den(i2p1,j1,k3))+
     &       ffact2*(den(i2m2,j1,k3)-den(i2p2,j1,k3))
      fx21 = ffact1*(den(i3m1,j1,k3)-den(i3p1,j1,k3))+
     &       ffact2*(den(i3m2,j1,k3)-den(i3p2,j1,k3))
      fx22 = ffact1*(den(i1m1,j2,k3)-den(i1p1,j2,k3))+
     &       ffact2*(den(i1m2,j2,k3)-den(i1p2,j2,k3))
      fx23 = ffact1*(den(i2m1,j2,k3)-den(i2p1,j2,k3))+
     &       ffact2*(den(i2m2,j2,k3)-den(i2p2,j2,k3))
      fx24 = ffact1*(den(i3m1,j2,k3)-den(i3p1,j2,k3))+
     &       ffact2*(den(i3m2,j2,k3)-den(i3p2,j2,k3))
      fx25 = ffact1*(den(i1m1,j3,k3)-den(i1p1,j3,k3))+
     &       ffact2*(den(i1m2,j3,k3)-den(i1p2,j3,k3))
      fx26 = ffact1*(den(i2m1,j3,k3)-den(i2p1,j3,k3))+
     &       ffact2*(den(i2m2,j3,k3)-den(i2p2,j3,k3))
      fx27 = ffact1*(den(i3m1,j3,k3)-den(i3p1,j3,k3))+
     &       ffact2*(den(i3m2,j3,k3)-den(i3p2,j3,k3))

      fy1 = ffact1*(den(i1,j1m1,k1)-den(i1,j1p1,k1))+
     &      ffact2*(den(i1,j1m2,k1)-den(i1,j1p2,k1))
      fy2 = ffact1*(den(i2,j1m1,k1)-den(i2,j1p1,k1))+
     &      ffact2*(den(i2,j1m2,k1)-den(i2,j1p2,k1))
      fy3 = ffact1*(den(i3,j1m1,k1)-den(i3,j1p1,k1))+
     &      ffact2*(den(i3,j1m2,k1)-den(i3,j1p2,k1))
      fy4 = ffact1*(den(i1,j2m1,k1)-den(i1,j2p1,k1))+
     &      ffact2*(den(i1,j2m2,k1)-den(i1,j2p2,k1))
      fy5 = ffact1*(den(i2,j2m1,k1)-den(i2,j2p1,k1))+
     &      ffact2*(den(i2,j2m2,k1)-den(i2,j2p2,k1))
      fy6 = ffact1*(den(i3,j2m1,k1)-den(i3,j2p1,k1))+
     &      ffact2*(den(i3,j2m2,k1)-den(i3,j2p2,k1))
      fy7 = ffact1*(den(i1,j3m1,k1)-den(i1,j3p1,k1))+
     &      ffact2*(den(i1,j3m2,k1)-den(i1,j3p2,k1))
      fy8 = ffact1*(den(i2,j3m1,k1)-den(i2,j3p1,k1))+
     &      ffact2*(den(i2,j3m2,k1)-den(i2,j3p2,k1))
      fy9 = ffact1*(den(i3,j3m1,k1)-den(i3,j3p1,k1))+
     &      ffact2*(den(i3,j3m2,k1)-den(i3,j3p2,k1))
      fy10 = ffact1*(den(i1,j1m1,k2)-den(i1,j1p1,k2))+
     &       ffact2*(den(i1,j1m2,k2)-den(i1,j1p2,k2))
      fy11 = ffact1*(den(i2,j1m1,k2)-den(i2,j1p1,k2))+
     &       ffact2*(den(i2,j1m2,k2)-den(i2,j1p2,k2))
      fy12 = ffact1*(den(i3,j1m1,k2)-den(i3,j1p1,k2))+
     &       ffact2*(den(i3,j1m2,k2)-den(i3,j1p2,k2))
      fy13 = ffact1*(den(i1,j2m1,k2)-den(i1,j2p1,k2))+
     &       ffact2*(den(i1,j2m2,k2)-den(i1,j2p2,k2))
      fy14 = ffact1*(den(i2,j2m1,k2)-den(i2,j2p1,k2))+
     &       ffact2*(den(i2,j2m2,k2)-den(i2,j2p2,k2))
      fy15 = ffact1*(den(i3,j2m1,k2)-den(i3,j2p1,k2))+
     &       ffact2*(den(i3,j2m2,k2)-den(i3,j2p2,k2))
      fy16 = ffact1*(den(i1,j3m1,k2)-den(i1,j3p1,k2))+
     &       ffact2*(den(i1,j3m2,k2)-den(i1,j3p2,k2))
      fy17 = ffact1*(den(i2,j3m1,k2)-den(i2,j3p1,k2))+
     &       ffact2*(den(i2,j3m2,k2)-den(i2,j3p2,k2))
      fy18 = ffact1*(den(i3,j3m1,k2)-den(i3,j3p1,k2))+
     &       ffact2*(den(i3,j3m2,k2)-den(i3,j3p2,k2))

      fy19 = ffact1*(den(i1,j1m1,k3)-den(i1,j1p1,k3))+
     &       ffact2*(den(i1,j1m2,k3)-den(i1,j1p2,k3))
      fy20 = ffact1*(den(i2,j1m1,k3)-den(i2,j1p1,k3))+
     &       ffact2*(den(i2,j1m2,k3)-den(i2,j1p2,k3))
      fy21 = ffact1*(den(i3,j1m1,k3)-den(i3,j1p1,k3))+
     &       ffact2*(den(i3,j1m2,k3)-den(i3,j1p2,k3))
      fy22 = ffact1*(den(i1,j2m1,k3)-den(i1,j2p1,k3))+
     &       ffact2*(den(i1,j2m2,k3)-den(i1,j2p2,k3))
      fy23 = ffact1*(den(i2,j2m1,k3)-den(i2,j2p1,k3))+
     &       ffact2*(den(i2,j2m2,k3)-den(i2,j2p2,k3))
      fy24 = ffact1*(den(i3,j2m1,k3)-den(i3,j2p1,k3))+
     &       ffact2*(den(i3,j2m2,k3)-den(i3,j2p2,k3))
      fy25 = ffact1*(den(i1,j3m1,k3)-den(i1,j3p1,k3))+
     &       ffact2*(den(i1,j3m2,k3)-den(i1,j3p2,k3))
      fy26 = ffact1*(den(i2,j3m1,k3)-den(i2,j3p1,k3))+
     &       ffact2*(den(i2,j3m2,k3)-den(i2,j3p2,k3))
      fy27 = ffact1*(den(i3,j3m1,k3)-den(i3,j3p1,k3))+
     &       ffact2*(den(i3,j3m2,k3)-den(i3,j3p2,k3))

      fz1 = ffact1*(den(i1,j1,k1m1)-den(i1,j1,k1p1))+
     &      ffact2*(den(i1,j1,k1m2)-den(i1,j1,k1p2))
      fz2 = ffact1*(den(i2,j1,k1m1)-den(i2,j1,k1p1))+
     &      ffact2*(den(i2,j1,k1m2)-den(i2,j1,k1p2))
      fz3 = ffact1*(den(i3,j1,k1m1)-den(i3,j1,k1p1))+
     &      ffact2*(den(i3,j1,k1m2)-den(i3,j1,k1p2))
      fz4 = ffact1*(den(i1,j2,k1m1)-den(i1,j2,k1p1))+
     &      ffact2*(den(i1,j2,k1m2)-den(i1,j2,k1p2))
      fz5 = ffact1*(den(i2,j2,k1m1)-den(i2,j2,k1p1))+
     &      ffact2*(den(i2,j2,k1m2)-den(i2,j2,k1p2))
      fz6 = ffact1*(den(i3,j2,k1m1)-den(i3,j2,k1p1))+
     &      ffact2*(den(i3,j2,k1m2)-den(i3,j2,k1p2))
      fz7 = ffact1*(den(i1,j3,k1m1)-den(i1,j3,k1p1))+
     &      ffact2*(den(i1,j3,k1m2)-den(i1,j3,k1p2))
      fz8 = ffact1*(den(i2,j3,k1m1)-den(i2,j3,k1p1))+
     &      ffact2*(den(i2,j3,k1m2)-den(i2,j3,k1p2))
      fz9 = ffact1*(den(i3,j3,k1m1)-den(i3,j3,k1p1))+
     &      ffact2*(den(i3,j3,k1m2)-den(i3,j3,k1p2))
      fz10 = ffact1*(den(i1,j1,k2m1)-den(i1,j1,k2p1))+
     &       ffact2*(den(i1,j1,k2m2)-den(i1,j1,k2p2))
      fz11 = ffact1*(den(i2,j1,k2m1)-den(i2,j1,k2p1))+
     &       ffact2*(den(i2,j1,k2m2)-den(i2,j1,k2p2))
      fz12 = ffact1*(den(i3,j1,k2m1)-den(i3,j1,k2p1))+
     &       ffact2*(den(i3,j1,k2m2)-den(i3,j1,k2p2))
      fz13 = ffact1*(den(i1,j2,k2m1)-den(i1,j2,k2p1))+
     &       ffact2*(den(i1,j2,k2m2)-den(i1,j2,k2p2))
      fz14 = ffact1*(den(i2,j2,k2m1)-den(i2,j2,k2p1))+
     &       ffact2*(den(i2,j2,k2m2)-den(i2,j2,k2p2))
      fz15 = ffact1*(den(i3,j2,k2m1)-den(i3,j2,k2p1))+
     &       ffact2*(den(i3,j2,k2m2)-den(i3,j2,k2p2))
      fz16 = ffact1*(den(i1,j3,k2m1)-den(i1,j3,k2p1))+
     &       ffact2*(den(i1,j3,k2m2)-den(i1,j3,k2p2))
      fz17 = ffact1*(den(i2,j3,k2m1)-den(i2,j3,k2p1))+
     &       ffact2*(den(i2,j3,k2m2)-den(i2,j3,k2p2))
      fz18 = ffact1*(den(i3,j3,k2m1)-den(i3,j3,k2p1))+
     &       ffact2*(den(i3,j3,k2m2)-den(i3,j3,k2p2))

      fz19 = ffact1*(den(i1,j1,k3m1)-den(i1,j1,k3p1))+
     &       ffact2*(den(i1,j1,k3m2)-den(i1,j1,k3p2))
      fz20 = ffact1*(den(i2,j1,k3m1)-den(i2,j1,k3p1))+
     &       ffact2*(den(i2,j1,k3m2)-den(i2,j1,k3p2))
      fz21 = ffact1*(den(i3,j1,k3m1)-den(i3,j1,k3p1))+
     &       ffact2*(den(i3,j1,k3m2)-den(i3,j1,k3p2))
      fz22 = ffact1*(den(i1,j2,k3m1)-den(i1,j2,k3p1))+
     &       ffact2*(den(i1,j2,k3m2)-den(i1,j2,k3p2))
      fz23 = ffact1*(den(i2,j2,k3m1)-den(i2,j2,k3p1))+
     &       ffact2*(den(i2,j2,k3m2)-den(i2,j2,k3p2))
      fz24 = ffact1*(den(i3,j2,k3m1)-den(i3,j2,k3p1))+
     &       ffact2*(den(i3,j2,k3m2)-den(i3,j2,k3p2))
      fz25 = ffact1*(den(i1,j3,k3m1)-den(i1,j3,k3p1))+
     &       ffact2*(den(i1,j3,k3m2)-den(i1,j3,k3p2))
      fz26 = ffact1*(den(i2,j3,k3m1)-den(i2,j3,k3p1))+
     &       ffact2*(den(i2,j3,k3m2)-den(i2,j3,k3p2))
      fz27 = ffact1*(den(i3,j3,k3m1)-den(i3,j3,k3p1))+
     &       ffact2*(den(i3,j3,k3m2)-den(i3,j3,k3p2))
c     print *,den(i1,j1,k1m1),i1,j1,k1m1
c     print *,den(i1,j1,k1p1),i1,j1,k1p1
c     print *,fz1,fz2,fz3,fz4,fz5,fz6,fz7,fz8
c     return

      fx=(wt1*fx1+wt2*fx2+wt3*fx3+wt4*fx4  +wt5 *fx5 +wt6 *fx6+
     & wt7 *fx7 +wt8 *fx8 +wt9 *fx9 +wt10*fx10+wt11*fx11+wt12*fx12+
     & wt13*fx13+wt14*fx14+wt15*fx15+wt16*fx16+wt17*fx17+wt18*fx18+
     & wt19*fx19+wt20*fx20+wt21*fx21+wt22*fx22+wt23*fx23+wt24*fx24+
     & wt25*fx25+wt26*fx26+wt27*fx27)
      fy=(wt1*fy1+wt2*fy2+wt3*fy3+wt4*fy4  +wt5 *fy5 +wt6 *fy6+
     & wt7 *fy7 +wt8 *fy8 +wt9 *fy9 +wt10*fy10+wt11*fy11+wt12*fy12+
     & wt13*fy13+wt14*fy14+wt15*fy15+wt16*fy16+wt17*fy17+wt18*fy18+
     & wt19*fy19+wt20*fy20+wt21*fy21+wt22*fy22+wt23*fy23+wt24*fy24+
     & wt25*fy25+wt26*fy26+wt27*fy27)
      fz=(wt1*fz1+wt2*fz2+wt3*fz3+wt4*fz4  +wt5 *fz5 +wt6 *fz6+
     & wt7 *fz7 +wt8 *fz8 +wt9 *fz9 +wt10*fz10+wt11*fz11+wt12*fz12+
     & wt13*fz13+wt14*fz14+wt15*fz15+wt16*fz16+wt17*fz17+wt18*fz18+
     & wt19*fz19+wt20*fz20+wt21*fz21+wt22*fz22+wt23*fz23+wt24*fz24+
     & wt25*fz25+wt26*fz26+wt27*fz27)

C#ifdef INCLUDE_TREE_FORCE
C       particles(nforcepos+1,i) = fx
C       particles(nforcepos+2,i) = fy
C       particles(nforcepos+3,i) = fz
C#else
c     fx = fx * fact2
c     fy = fy * fact2
c     fz = fz * fact2
      particles(nxpos+4,i) = particles(nxpos+4,i) + fact2*fx
      particles(nxpos+5,i) = particles(nxpos+5,i) + fact2*fy
      particles(nxpos+6,i) = particles(nxpos+6,i) + fact2*fz
C The position of partcles should be updated after tree correction 
C      particles(nxpos+1,i)=amod(particles(nxpos+1,i)+
C     &                     particles(nxpos+4,i)*pfact+rngx,rngx)
C      particles(nxpos+2,i)=amod(particles(nxpos+2,i)+
C     &                     particles(nxpos+5,i)*pfact+rngy,rngy)
C      particles(nxpos+3,i)=amod(particles(nxpos+3,i)+
C     &                     particles(nxpos+6,i)*pfact+rngz,rngz)
C#endif

c     x(i)=amod(x(i)+vx(i)*pfact+rng,rng)
c     y(i)=amod(y(i)+vy(i)*pfact+rng,rng)
c     z(i)=amod(z(i)+vz(i)*pfact+rngz,rngz)
  160 continue
      if(myid .EQ. 0) print *,'CPU(force) =', etime(CPU)-stime

      return
      end
ccc****************************************************************
      subroutine onestepforwardposition(particles,np,nn,nxpos,pfact,
     &          nx,ny,nz,nstep)
      include 'mpif.h'
      integer nn,np,nxpos,nstep
      integer nx,ny,nz,ierror
      real rngx,rngy,rngz
      real particles(nn,np)
      real dmax,xd,yd,zd,rd,tdmax
      integer myid

      call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierror)
      rngx = nx
      rngy = ny
      rngz = nz
      dmax = -1.e20
c$omp parallel do private(xd,yd,zd,rd)
c$omp+ reduction(MAX:dmax)
      do i = 1, np
         xd = particles(nxpos+4,i)*pfact
         yd = particles(nxpos+5,i)*pfact
         zd = particles(nxpos+6,i)*pfact
         rd = xd**2+yd**2+zd**2
         dmax = amax1(dmax,rd)
         particles(nxpos+1,i)=amod(particles(nxpos+1,i)+xd+rngx,rngx)
         particles(nxpos+2,i)=amod(particles(nxpos+2,i)+yd+rngy,rngy)
         particles(nxpos+3,i)=amod(particles(nxpos+3,i)+zd+rngz,rngz)
      enddo
      dmax = sqrt(dmax)
      call MPI_Reduce(dmax,tdmax,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,
     &        ierror)
      if(myid .eq. 0) print *,'Maximum displacement = ',nstep,tdmax
      
      return
      end
      subroutine fda4sup(nx,ny,nz,local_nz)
      IMPLICIT none
      integer nx,ny,nz,local_nz
      integer xabp0(-10:NNX+10),xabp1(-10:NNX+10),xabm1(-10:NNX+10)
      integer yabp0(-10:NNY+10),yabp1(-10:NNY+10),yabm1(-10:NNY+10)
      integer zabp0(-10:LOCAL_NNZ+7+10),zabp1(-10:LOCAL_NNZ+7+10)
      integer zabm1(-10:LOCAL_NNZ+7+10)
      common /ptables/ xabp0,xabp1,xabm1,yabp0,yabp1,yabm1,zabp0,
     &       zabp1,zabm1
      integer ngx,ngy,ngz,i
      print *,nx,ny,nz,local_nz
      ngx = nx
      ngy = ny
      ngz = local_nz+7
c---  Making periodic xables
      do 10 i=-10,0
      xabm1(i)=ngx+i-1
      xabp0(i)=ngx+i
      xabp1(i)=ngx+i+1
      yabm1(i)=ngy+i-1
      yabp0(i)=ngy+i
      yabp1(i)=ngy+i+1
      zabm1(i)=ngz+i-1
      zabp0(i)=ngz+i
   10 zabp1(i)=ngz+i+1
      xabp1(0)=1
      yabp1(0)=1
      zabp1(0)=1
      do 20 i=1,ngx
      xabm1(i)=i-1
      xabp0(i)=i
   20 xabp1(i)=i+1
      xabm1(1)=ngx
      xabp1(ngx)=1
      do 21 i=1,ngy
      yabm1(i)=i-1
      yabp0(i)=i
   21 yabp1(i)=i+1
      yabm1(1)=ngy
      yabp1(ngy)=1

      do 30 i=1,ngz
      zabm1(i) = i-1
      zabp0(i) = i
   30 zabp1(i) = i+1
      zabm1(1) = ngz
      zabp1(ngz)=1
      do 40 i=ngx+1,ngx+10
      xabm1(i)=i-ngx-1
      xabp0(i)=i-ngx
   40 xabp1(i)=i+1-ngx
      xabm1(ngx+1)=ngx
      do 41 i=ngy+1,ngy+10
      yabm1(i)=i-ngy-1
      yabp0(i)=i-ngy
   41 yabp1(i)=i+1-ngy
      yabm1(ngx+1)=ngy

      do 50 i=ngz+1,ngz+10
      zabm1(i)=i-ngz-1
      zabp0(i)=i-ngz
   50 zabp1(i)=i+1-ngz
      zabm1(ngz+1) = ngz


      return
      end
