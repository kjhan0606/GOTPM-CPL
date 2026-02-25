CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Parallel version of PMSEED using density array with local_nz+7 zwidth
C     It does not allow inhomogeneous widths of local domains.  Feb. 16, 2006 by Juhan Kim
C     Now under recording to admit the inhomogeneity of domain decompositions.
C     local_nz : local z width in fft
C     local_z_start : starting of z in fft
C     zstart : starting point in z for load balance
C     zwidth : width in z for load balance
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c234567
      subroutine pmseedforce(np,particles,nn,nxpos,den,omep,omeplam,ken,
     &     denpad,zstart,zwidth,vamp,omei,fx,fy,fz)
      IMPLICIT none
      include 'mpif.h'
      include "fftw_f77.i"
#include "fortran.i"
      real iceil
      external iceil
      integer nzsend,noffset
      real zstart,zwidth,zfinal
      integer izwidth,izstart
      integer m2z,p2z
      character*1 dir1,dir2,dir3,dir4,dir5
      integer TID,NTHREADS,OMP_GET_THREAD_NUM
      integer OMP_GET_NUM_THREADS,myistart,myiend
      real bpowernow

      integer np,mbias,nxpos,NZBUFFER
      parameter(mbias=4,NZBUFFER=4)
      integer nx,ny,nz,local_nz,local_z_start
      integer myid,nid,dest,src,stag,rtag
      real den(2*(nx/2+1),ny,local_nz+NZBUFFER),denpad(2,ny,nz)
      real fx(2*(nx/2+1),ny,*)
      real fy(2*(nx/2+1),ny,*)
      real fz(2*(nx/2+1),ny,*)
      real displ,mdispl
      integer status(MPI_STATUS_SIZE),root

      integer nn
      real particles(nn,np)
      integer*8 lnp,tlnp
      integer mlocalstartz,mlocalnz

      real w1, w2

      external powint,pow2int,bpowint
      real ken
      double precision tmpken,tken
      real ktot,cpu(2),stime,etime
      character*3 dark
      integer iseed,nspace,iiseed
      real size
      real smooth,h,rth,bias,zi
      integer knyq,ngz,nnp,nznp,ng,ngx,ngy
      real npow
      real rngx,rngy
      real xyz,xyy,omei,powh,pi,omehs,rknyq,rknyq2,rkinv
      real twopi,rngz,rngzc,ran3,xxx,rng,rngc
      real pmas,pmass,ai,anow,wavemax,pamp,pamp1,pamp2,vamp
      real xnow, xini,bnow,bini
      real vfact,vfact2,pfact,amax,astep,omep,omeplam
      real efold,sfac
      real const,poten,growth,red
      integer iindex

      external gasdev
      external Pgasdev0,Pgasdev1,Pgasdev2,Pgasdev3,Pgasdev4
      external Pgasdev5,Pgasdev6,Pgasdev7
      real Pgasdev0,Pgasdev1,Pgasdev2,Pgasdev3,Pgasdev4
      real Pgasdev5,Pgasdev6,Pgasdev7
      
      integer ngh,i,j,k,ky,kz,nghx,nghy,nghz
      real rk,velpot,gasdev,Pgasdev,xksq,power,bpower,bpower1
      integer sendsize,recvsize,ierror
  
      integer ii,jj,kk,iim1,jjm1,kkm1,iim2,jjm2,kkm2
      integer iip1,jjp1,kkp1,iip2,jjp2,kkp2
      integer nadd,nzparticles,nzprefixel

      real work
      real*8 rincube
      common/powsp/rknyq,rkinv,powh,smooth,omehs,dark
      common/tophat/rth
      real com_omep,com_omepb,com_omeplam,com_h,com_size,omepb
      common/cospar/com_omep,com_omepb,com_omeplam,com_h,com_size

      integer local_ny_after_transpose,
     &     local_y_start_after_transpose,total_local_size
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
C      common/seed2/ai,ken,size,iseed,amax,astep
      character*190 powfilename,asciipkfile
      real asciirk(8192), asciipk(8192)
      integer nascii
      integer powreadflag
      common /ASCIIPOWER/ powreadflag,nascii,asciirk,asciipk

      integer seedflag


      call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierror)
      call MPI_Comm_size(MPI_COMM_WORLD,nid,ierror)

      if(myid .eq. 0) stime = etime(CPU)

      print *, myid,zstart,zwidth

C This is a hook into the C functions
      call getpowerparams(nspace,size,amax,astep,
     &                        npow,smooth,h,rth,bias,iseed,omepb,
     &                    powreadflag,powfilename,asciipkfile)
      dark = 'cdm'
      com_omep=omep
      com_omepb = omepb
      com_omeplam=omeplam
      com_h = h
      com_size = size
      if(iseed .ge. 0) then
          seedflag = 1
          iseed = -1
      else
          seedflag = 0
      endif


      rincube = 1.d0/dble(nx)/dble(ny)/dble(nz)
      ngx = nx
      ngy = ny
      ngz = nz
      ng = nx
      nghx = ngx/2
      ngh = nghx
      nghy = ngy/2
      nghz = ngz/2
      rng = real(ng)
      rngx = real(ngx)
      rngy = real(ngy)
      rngz = real(ngz)
      ngz = nz
      rngc = rngx*rngy*rngz
      rngz = real(ngz)
      rngzc = rngx*rngy*rngz
      xyz  = rngx/rngz
      xyy  = rngx/rngy
      pi = 3.141592654
      twopi = 2.*pi
      nnp  = ng/nspace
      nznp = ngz/nspace





      knyq = ngx/2
      rknyq2 = real(knyq)**2
      rknyq = twopi*real(knyq)
      rkinv = 1./rknyq
      powh = 0.5*real(npow)
      ai = 1.
      zi = amax-1
      anow = ai*(1.+zi)
c     omei = omep*(1.+zi)/(1.+omep*zi)
      omei = omep*(1+zi)**3
     &     /(omep*(1+zi)**3+omeplam+(1-omep-omeplam)*(1+zi)**2)
      omehs = 1./(omep*h*size)
      smooth = (smooth/sqrt(2.))/rng/sqrt(2.)
      rth = rth/rng
c---  Total mass is ngz/ng so that the mean density = 1
      pmas  = real(nspace)**3
      pmass = pmas/rngc


      wavemax = 100000
c     if(size.lt.33)wavemax = 300
c     from wavemin=0.0001 to wavemax=100
#ifdef EH98
c     THIS IS THE POWER FIT BY EISENSTEIN & HU 1998
      pamp = 0
      w1 = 0.001*size
      do i = 1, 5
         w2 = w1*10.
         call qromb(bpowint,w1,w2,pamp1)
         pamp = pamp + pamp1
         w1 = w2
      enddo
c     call qromb(bpowint,0.0001,1000.,pamp1)
c     call qromb(bpowint,1000.,wavemax,pamp2)
c     pamp = pamp1+pamp2
#else
c     THIS IS THE POWER CALCULATED BY THE CAMB
      if(powreadflag .ne. 0) then
         do i=0,nid-1
            if(myid.eq.i) then
               call rInitPower(omep,omepb,omeplam,h,npow,zi,bias,
     &                pamp,size,rng,powreadflag,powfilename,myid)
               powfilename = trim(powfilename)//char(0)
            endif
            call MPI_Barrier(MPI_COMM_WORLD,ierror)
         enddo
         if(powreadflag .eq. 2) then
         if(myid.eq.0) then
            open(91,file=asciipkfile)
            i = 1
899         read(91,*,end=900) asciirk(i), asciipk(i)
            if(i .gt. 8192) then
               print *,'Exit from reading external ascii powerspectrum '
               stop
            endif
            i = i + 1
            goto 899
900         close(91)
            nascii = i - 1
         endif
         call MPI_Bcast(nascii,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
         call MPI_Bcast(asciirk,nascii,MPI_REAL,0,MPI_COMM_WORLD,ierror)
         call MPI_Bcast(asciipk,nascii,MPI_REAL,0,MPI_COMM_WORLD,ierror)
         ! This is to measure the amplitude of the input power which is modified.
         pamp = 0
         w1 = 0.001*size
         do i = 1, 5 ! integration from k=0.0001 to 1000.
            w2 = w1*10.
            call qromb(bpowint,w1,w2,pamp1)
            pamp = pamp + pamp1
            w1 = w2
         enddo
         endif
      else 
         call InitPower(omep,omepb,omeplam,h,npow,zi,bias,pamp,size,rng)
      endif

#endif

c     call qsimp(pow2int,0.00001,wavemax,pamp2)
      if(myid .EQ. 0) print *,' <1>=',pamp
c---  Here pamp = initial ms of top-hat smoothed density fluctuation
c---  pamp = 50.06, 800.95, 12815 for 64, 128, 256 h^-1 Mpc
c---  Now pamp*sqrt(P(k)) give rms fluctuation equal to thamp initially.
c---  After the expansion factor grows by the factor anow, the linear
c---  amplitude of fluctuation  will be (1/bias) over the scale rth.
c---  I. Einstein-de Sitter model

      if(omep.eq.1.)then
        pamp = (1./bias)/sqrt(pamp)*ai/anow
        vamp = 1
c---  II. Open model
      else if(omep.lt.1.and.omeplam.eq.0)then
        xnow = (1.-omep)/omep
        xini = xnow/(1.+zi)
        bini = 1.+3./xini+3.*sqrt(1.+xini)/xini**1.5*alog(sqrt(1.+xini)
     &         -sqrt(xini))
        bnow = 1.+3./xnow+3.*sqrt(1.+xnow)/xnow**1.5*alog(sqrt(1.+xnow)
     &         -sqrt(xnow))
        pamp = (1./bias)/sqrt(pamp)*bini/bnow
        vamp = omei**0.6
      else if(omep.lt.1.and.omeplam.gt.0)then
c---  III. Zero curvature with cosmological constant model
        bini = growth(1.,  amax,omep)
        bnow = growth(amax,amax,omep)
        pamp = (1./bias)/sqrt(pamp)*bini/bnow
        vamp = (growth(1.+0.005,amax,omep)-growth(1.-0.005,amax,omep))
     &        /(0.01)*1./bini
        if(myid.eq.0) then
           print *,omep,' bini=',bini,' bnow=',bnow,vamp
        endif
      end if



      ngz = nz

      efold = 0.585*sqrt(2.)
      sfac   = -(pi*efold/size)**2
c***  Velocity Potential
      pamp = pamp*rngzc*sqrt(0.5)
ccc   P(k) = (pamp*bpower1(rk))**2 * boxsize**3; here rk is k and rk in the
ccc   code is 2pi*u


901   continue


      if(seedflag .ne. 0) then
#ifdef _OPENMP
!$OMP PARALLEL private(i,j,k,ky,kz,xksq,rk,velpot,TID,NTHREADS)
!$OMP+ private (myistart,myiend,xxx,iiseed )
#endif
      TID = OMP_GET_THREAD_NUM()
      NTHREADS = OMP_GET_NUM_THREADS()
      myistart = TID*((nghx+NTHREADS-1)/NTHREADS)
      if(TID .eq. NTHREADS-1) then
         myiend = nghx
      else
         myiend   = MIN(myistart + (nghx+NTHREADS-1)/NTHREADS,nghx)-1
      endif

      iiseed = iseed*(TID+NTHREADS*myid+1)
      xxx = ran3(iiseed)
#ifndef _OPENMP
      do i=0,nghx
#else
      do i=myistart,myiend
#endif
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.nghx)ky=ky-ngx
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.nghx)kz=kz-ngx
            xksq = i**2+ky**2+kz**2
            if(xksq.eq.0..or.xksq.gt.rknyq2)then
              den(2*i+1,j,k)=0.
              den(2*i+2,j,k)=0.
            else
              rk = twopi*sqrt(xksq)
              den(2*i+1,j,k) = gasdev(iiseed)
              den(2*i+2,j,k) = gasdev(iiseed)
            end if
         enddo
      enddo
      enddo
C$OMP END PARALLEL
      else
#ifdef _OPENMP
C$OMP PARALLEL private(i,j,k,ky,kz,xksq,rk,velpot,TID,NTHREADS)
C$OMP+ private (myistart,myiend,xxx,iiseed)
      TID = OMP_GET_THREAD_NUM()
      NTHREADS = OMP_GET_NUM_THREADS()
      myistart = TID*((nghx+NTHREADS-1)/NTHREADS)
      if(TID .eq. NTHREADS-1) then
         myiend = nghx
      else
         myiend   = MIN(myistart+(nghx+NTHREADS-1)/NTHREADS,nghx)-1
      endif

      iiseed = iseed*(TID+NTHREADS*myid+1)
c     xxx = ran3(iiseed)
      do i=myistart,myiend
#else
      iiseed = iseed*(myid+1)
      xxx = ran3(iiseed)
      do i=0,nghx
#endif
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.nghx)ky=ky-ngx
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.nghx)kz=kz-ngx
            xksq = i**2+ky**2+kz**2
            if(xksq.eq.0..or.xksq.gt.rknyq2)then
              den(2*i+1,j,k)=0.
              den(2*i+2,j,k)=0.
            else
              rk = twopi*sqrt(xksq)
#ifdef _OPENMP
c  This line should not include bpower function because it is not thread safe.
c  But gasdev seems to be thread safe at first sight.
              velpot = pamp/rk**2
              if(TID .eq. 0) then
              den(2*i+1,j,k) = velpot*Pgasdev0(iiseed)
              den(2*i+2,j,k) = velpot*Pgasdev0(iiseed)
              else if(TID .eq.1) then
              den(2*i+1,j,k) = velpot*Pgasdev1(iiseed)
              den(2*i+2,j,k) = velpot*Pgasdev1(iiseed)
              else if(TID .eq.2) then
              den(2*i+1,j,k) = velpot*Pgasdev2(iiseed)
              den(2*i+2,j,k) = velpot*Pgasdev2(iiseed)
              else if(TID .eq.3) then
              den(2*i+1,j,k) = velpot*Pgasdev3(iiseed)
              den(2*i+2,j,k) = velpot*Pgasdev3(iiseed)
              else if(TID .eq.4) then
              den(2*i+1,j,k) = velpot*Pgasdev4(iiseed)
              den(2*i+2,j,k) = velpot*Pgasdev4(iiseed)
              else if(TID .eq.5) then
              den(2*i+1,j,k) = velpot*Pgasdev5(iiseed)
              den(2*i+2,j,k) = velpot*Pgasdev5(iiseed)
              else if(TID .eq.6) then
              den(2*i+1,j,k) = velpot*Pgasdev6(iiseed)
              den(2*i+2,j,k) = velpot*Pgasdev6(iiseed)
              else if(TID .eq.7) then
              den(2*i+1,j,k) = velpot*Pgasdev7(iiseed)
              den(2*i+2,j,k) = velpot*Pgasdev7(iiseed)
              endif
#else
              velpot = pamp*bpower(rk)/rk**2
              den(2*i+1,j,k) = velpot*gasdev(iiseed)
              den(2*i+2,j,k) = velpot*gasdev(iiseed)
#endif
            end if
         enddo
      enddo
      enddo
#ifdef _OPENMP
C$OMP SINGLE
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.nghx)ky=ky-ngx
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.nghx)kz=kz-ngx
            xksq = i**2+ky**2+kz**2
            if(xksq.eq.0..or.xksq.gt.rknyq2)then
              den(2*i+1,j,k)=0.
              den(2*i+2,j,k)=0.
            else
              rk = twopi*sqrt(xksq)
              velpot = bpower(rk)
              den(2*i+1,j,k) = den(2*i+1,j,k)*velpot
              den(2*i+2,j,k) = den(2*i+2,j,k)*velpot
            end if
         enddo
      enddo
      enddo
C$OMP END SINGLE
#endif
#ifdef _OPENMP
C$OMP END PARALLEL
#endif
      endif

      if(myid .EQ. 0) print *,'last #=',ran3(iiseed)
cccc
c***  Make an array, den(*,*,local_nz+1:local_nz+2) which is needed to
ccc   follow the complex conjugate relation for a real array.
ccc   mpi_gatherv is gathering values with displacement of local_nz.
ccc   0 is the root node number.
ccc   Gathering all the den(*,*,local_nz+1:local_nz+2) to 
ccc     den(*,*,local_nz+1:local_nz+2) of 0'th processor.
ccc   rcounts is the size of sending buffer and also displs is displacement of
ccc   receive buffer saved to 0'th node.

c---  Satisfy complex conjugation relation for a real array
ccc   Only on rank 0 processor.
      do k = 1, local_nz
      do j = 1, ny
      do i = 1, 2
          denpad(i,j,local_z_start+k) = den(i,j,k)
      enddo
      enddo
      enddo
      if(myid .eq. 0) then
         do i = 1, nid-1
           call MPI_RECV(mlocalstartz,1,MPI_INTEGER,i,0,MPI_COMM_WORLD,
     &              status,ierror)
           call MPI_RECV(mlocalnz,1,MPI_INTEGER,i,0,MPI_COMM_WORLD,
     &              status,ierror)
           call MPI_RECV(denpad(1,1,mlocalstartz+1),2*ny*mlocalnz,
     &              MPI_REAL,i,0,MPI_COMM_WORLD,status,ierror)
#ifdef DEBUG
           print *,'P0 received data from',i
#endif
         enddo
      else
         call MPI_SEND(local_z_start,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,
     &             ierror)
         call MPI_SEND(local_nz,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,
     &             ierror)
         call MPI_SEND(denpad(1,1,local_z_start+1),2*ny*local_nz,
     &     MPI_REAL,0,0,MPI_COMM_WORLD,ierror)
      endif
c---  Satisfy complex conjugation relation for a real array
      if(myid .eq. 0) then
         do 12 k = 2, nghz
            denpad(1,1,ngz+2-k) = denpad(1,1,k)
            denpad(2,1,ngz+2-k) =-denpad(2,1,k)
            do 12 j=2,nghy
            denpad(1,ngy+2-j,ngz+2-k) = denpad(1,j,k)
12          denpad(2,ngy+2-j,ngz+2-k) =-denpad(2,j,k)

         do 13 j=2,nghy
            denpad(1,ngy+2-j,1) = denpad(1,j,1)
13          denpad(2,ngy+2-j,1) =-denpad(2,j,1)

         do 14 k=2,nghz
         do 14 j=nghy+2,ngy
            denpad(1,ngy+2-j,ngz+2-k) = denpad(1,j,k)
14          denpad(2,ngy+2-j,ngz+2-k) =-denpad(2,j,k)
      endif
      sendsize = 2*ny*nz
      call MPI_BCAST(denpad,sendsize,MPI_REAL,0,MPI_COMM_WORLD,ierror)

      if(myid.eq.0) print *,'just after broadcasting the conjugate'
      do k = 1, local_nz
      do j = 1, ny
      do i = 1, 2
          den(i,j,k) = denpad(i,j,local_z_start+k)
      enddo
      enddo
      enddo

C   If this is to find the seed number then goto 901 or finish it.
      if(seedflag .ne.0) then
         call findiseed(den,iseed)
         iseed = iseed -1
         if(iseed .gt. -50) then
             goto 901
         else
            call MPI_Finalize()
            stop
         endif
      endif

      if(myid.eq.0) print *,'just before making force array'
c$omp parallel do private(i,j,k,ky,kz)
      do i=0,nghx
      do j=1,ngy
      ky=xyy*(j-1)
      if(ky.gt.ngh)  ky=ky-ng
      do k=1,local_nz
         kz=xyz*((k-1)+local_z_start)
         if(kz.gt.ngh)kz=kz-ng
         fx(2*i+1,j,k) = -twopi*i*den(2*i+2,j,k)
         fx(2*i+2,j,k) =  twopi*i*den(2*i+1,j,k)
         fy(2*i+1,j,k) = -twopi*ky*den(2*i+2,j,k)
         fy(2*i+2,j,k) =  twopi*ky*den(2*i+1,j,k)
         fz(2*i+1,j,k) = -twopi*kz*den(2*i+2,j,k)
         fz(2*i+2,j,k) =  twopi*kz*den(2*i+1,j,k)
      enddo
      enddo
      enddo

ccc   FFT COMPLEX-TO-REAL 
c$omp parallel 
c$omp single
#ifdef FFTW_DIVISIBLE_SIZE
      if(myid.eq.0) print *,'passing fft for fx '
      if(myid.lt. FFTW_DIVISIBLE_SIZE) then
      call rfftwnd_f77_mpi(iplan,1,fx,work,USE_WORK,FFTW_NORMAL_ORDER)
      endif
      if(myid.eq.0) print *,'passing fft for fy '
      if(myid.lt. FFTW_DIVISIBLE_SIZE) then
      call rfftwnd_f77_mpi(iplan,1,fy,work,USE_WORK,FFTW_NORMAL_ORDER)
      endif
      if(myid.eq.0) print *,'passing fft for fz '
      if(myid.lt. FFTW_DIVISIBLE_SIZE) then
      call rfftwnd_f77_mpi(iplan,1,fz,work,USE_WORK,FFTW_NORMAL_ORDER)
      endif
#else
      if(myid.eq.0) print *,'passing fft for fx '
      call rfftwnd_f77_mpi(iplan,1,fx,work,USE_WORK,FFTW_NORMAL_ORDER)
      if(myid.eq.0) print *,'passing fft for fy '
      call rfftwnd_f77_mpi(iplan,1,fy,work,USE_WORK,FFTW_NORMAL_ORDER)
      if(myid.eq.0) print *,'passing fft for fz '
      call rfftwnd_f77_mpi(iplan,1,fz,work,USE_WORK,FFTW_NORMAL_ORDER)
#endif
c$omp end single
c$omp end parallel 

       do k = 1, local_nz
c$omp parallel do
       do j = 1, ny
       do i = 1, nx
           fx(i,j,k) = fx(i,j,k)*rincube
           fy(i,j,k) = fy(i,j,k)*rincube
           fz(i,j,k) = fz(i,j,k)*rincube
       enddo
       enddo
       enddo
c      dir1 = char(48+myid/10000)
c      dir2 = char(48+mod(myid,10000)/1000)
c      dir3 = char(48+mod(myid,1000)/100)
c      dir4 = char(48+mod(myid,100)/10)
c      dir5 = char(48+mod(myid,10))
c      do ii = 0 , nid-1
c         if(myid.eq.ii) then
c             open(33,file='X_IFA'//dir1//dir2//dir3//dir4//dir5,
c                      form='unformatted')
c             write(33) nx,ny,local_nz
c             write(33) ((fx(i,j,k),i=1,nx),j=1,ny,k=1,local_nz)
c             close(33)
c             open(33,file='Y_IFA'//dir1//dir2//dir3//dir4//dir5,
c                      form='unformatted')
c             write(33) nx,ny,local_nz
c             write(33) ((fy(i,j,k),i=1,nx),j=1,ny,k=1,local_nz)
c             close(33)
c             open(33,file='Z_IFA'//dir1//dir2//dir3//dir4//dir5,
c                      form='unformatted')
c             write(33) nx,ny,local_nz
c             write(33) ((fz(i,j,k),i=1,nx),j=1,ny,k=1,local_nz)
c             close(33)
c         endif
c         call MPI_Barrier(MPI_COMM_WORLD,ierror)
c      enddo
       if(myid.eq.0) print *,'Ending of saving initial force array'

       return
       end

c---------------------------------------------------------------------
       subroutine getinnerparticlesforce(np,pvr,nn,nxpos,
     &     nzparticlewidth,zstartpos,
     &     nzprefixel,omep,omeplam,ken,zstart,zwidth,vamp,nspace,
     &     omei,amax,fx,fy,fz)
c      IMPLICIT none
      include 'mpif.h'
#include "fortran.i"
      include "fftw_f77.i"
      real omei,amax
      integer*8 lnp,tlnp
      double precision delx,tdelx,tmpken,tken
       integer np,nn,nxpos,nspace,mp
      integer nx,ny,nz,nzper
       real omep,omeplam,ken,zstart,zwidth,vamp
       real pvr(3,*)
       real fx(2*(nx/2+1),ny,*)
       real fy(2*(nx/2+1),ny,*)
       real fz(2*(nx/2+1),ny,*)
       real vfact,vfact2
       integer lindx,i,j,k
       real zstartpos,pfact,ai,pi
      real displ,mdispl,tmdispl
       integer nzparticlewidth,nzprefixel,ngz
       integer NZBUFFER,rng
      parameter(NZBUFFER=4)
c     integer xabp0(-10:NNX+10),xabp1(-10:NNX+10),xabm1(-10:NNX+10)
c     integer yabp0(-10:NNY+10),yabp1(-10:NNY+10),yabm1(-10:NNY+10)
c     integer zabp0(-10:LOCAL_NNZ+NZBUFFER+10)
c     integer zabp1(-10:LOCAL_NNZ+NZBUFFER+10)
c     integer zabm1(-10:LOCAL_NNZ+NZBUFFER+10)
c     common /ptables/ xabp0,xabp1,xabm1,yabp0,yabp1,yabm1,zabp0,
c    &       zabp1,zabm1
      real rnx,rny,rnz
      integer local_ny_after_transpose,local_nz,local_z_start,
     &     local_y_start_after_transpose,total_local_size
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      integer myid,nid,ierror
      integer ii,jj,kk,iim1,jjm1,kkm1,iim2,jjm2,kkm2,iip1,iip2,jjp1,jjp2
      integer kkp1,kkp2
      integer local_mz
      real zfinal

      rng = rnz

c     vfact = (4./3.)*rng/2.*vamp
      vfact = vamp
      np = 0
      do 20 k=2,nzparticlewidth-1
      do 20 j=1,ny/nspace
      do 20 i=1,nx/nspace
      np = np + 1
      ii = nspace*(i-1)+1
      jj = nspace*(j-1)+1
      kk = nspace*(k-1)+1+nzprefixel
c--------------------------------------------------------------
      pvr(1,np) = vfact*fx(ii,jj,kk)
      pvr(2,np) = vfact*fy(ii,jj,kk)
      pvr(3,np) = vfact*fz(ii,jj,kk)
c--------------------------------------------------------------
   20 continue
      return
      end
c---------------------------------------------------------------------
       subroutine setinnerparticlesforce(np,particles,nn,nxpos,
     &     nzparticlewidth,zstartpos,
     &     nzprefixel,omep,omeplam,ken,zstart,zwidth,vamp,nspace,
     &     omei,amax)
c      IMPLICIT none
      include 'mpif.h'
#include "fortran.i"
      include "fftw_f77.i"
      real omei,amax
      integer*8 lnp,tlnp
      double precision delx,tdelx,tmpken,tken
       integer np,nn,nxpos,nspace,mp
      integer nx,ny,nz,nzper
       real omep,omeplam,ken,zstart,zwidth,vamp
       real particles(nn,np)
       real vfact,vfact2
       integer lindx,i,j,k
       real zstartpos,pfact,ai,pi
      real displ,mdispl,tmdispl
       integer nzparticlewidth,nzprefixel,ngz
       integer NZBUFFER,rng
      parameter(NZBUFFER=4)
c     integer xabp0(-10:NNX+10),xabp1(-10:NNX+10),xabm1(-10:NNX+10)
c     integer yabp0(-10:NNY+10),yabp1(-10:NNY+10),yabm1(-10:NNY+10)
c     integer zabp0(-10:LOCAL_NNZ+NZBUFFER+10)
c     integer zabp1(-10:LOCAL_NNZ+NZBUFFER+10)
c     integer zabm1(-10:LOCAL_NNZ+NZBUFFER+10)
c     common /ptables/ xabp0,xabp1,xabm1,yabp0,yabp1,yabm1,zabp0,
c    &       zabp1,zabm1
      real rnx,rny,rnz
      integer local_ny_after_transpose,local_nz,local_z_start,
     &     local_y_start_after_transpose,total_local_size
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      integer myid,nid,ierror
      integer status(MPI_STATUS_SIZE),root
      real ktot,cpu(2),stime,etime
      real const,poten
      integer ii,jj,kk,iim1,jjm1,kkm1,iim2,jjm2,kkm2,iip1,iip2,jjp1,jjp2
      integer kkp1,kkp2
      integer local_mz
      real zfinal
      pi = 3.141592654
      root = 0
      zfinal = zstart+zwidth
      local_mz = int(zfinal) - int(zstart) + 1
      if(zfinal .eq. real(int(zfinal))) local_mz = local_mz - 1
      call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierror)
      call MPI_Comm_size(MPI_COMM_WORLD,nid,ierror)
      ai = 1.
      rnx = nx
      rny = ny
      rnz = nz
      rng = rnz
c---  Periodic table for mesh
ccc   local_mz+NZBUFFER is the local periodic condition in the current processor
      nzper = local_mz+NZBUFFER

c---  For positions we use the pixel size as the length unit
      mp = np
c     vfact = (4./3.)*rng/2.*vamp
c     vfact2 = -(1./3.)*rng/4.*vamp
      vfact = vamp

c---  For (n,size,seed)=(64,64,17), vx(1,1,1)=-0.0051 & vy(1,1,1)=-0.0049

c***  Push : position at a=ai (dx/da is constant from a=0 to ai)
c---  Displacement from ai*dx/da at particle position
c     pfact = rng*ai/vamp
      pfact = rng*ai
      do k=2,nzparticlewidth-1
      do j=1,ny/nspace
      do i=1,nx/nspace
         np = np + 1
#ifdef XYZDBL
         particles(nxpos+1,np) = 0
         particles(nxpos+2,np) = 0
         particles(nxpos+3,np) = 0
#else
         particles(nxpos+1,np) = nspace*(i-1.)
         particles(nxpos+2,np) = nspace*(j-1.)
         particles(nxpos+3,np) = nspace*(k-1.)+zstartpos
#endif
      enddo
      enddo
      enddo

#ifdef INDEX
      call indexingip(particles,np,mp,nzparticlewidth,nspace,nx,ny,nz)
#endif
c--------------------------------------------------------------------
c$omp parallel do private(i)
      do i = mp+1, np
#ifdef XYZDBL
         particles(nxpos+1,i) = particles(nxpos+1,i)+
     &                             particles(nxpos+4,i)*pfact
         particles(nxpos+2,i) = particles(nxpos+2,i)+
     &                             particles(nxpos+5,i)*pfact
         particles(nxpos+3,i) = particles(nxpos+3,i)+
     &                             particles(nxpos+6,i)*pfact
#else
         particles(nxpos+1,i) = amod(particles(nxpos+1,i)+
     &         particles(nxpos+4,i)*pfact+rnx,rnx)
         particles(nxpos+2,i) = amod(particles(nxpos+2,i)+
     &         particles(nxpos+5,i)*pfact+rny,rny)
         particles(nxpos+3,i) = amod(particles(nxpos+3,i)+
     &         particles(nxpos+6,i)*pfact+rnz,rnz)
#endif
      enddo
      delx = 0.D0
      displ = 0
      mdispl = 0.
c$omp parallel do private(i,displ)
c$omp+ reduction(+:delx)
c$omp+ reduction(MAX:mdispl)
      do i = 1, np
         delx = delx + particles(nxpos+4,i)**2 + 
     &         particles(nxpos+5,i)**2+particles(nxpos+6,i)**2
         displ = sqrt(particles(nxpos+4,i)**2+particles(nxpos+5,i)**2+
     &                particles(nxpos+6,i)**2)*pfact
         mdispl = amax1(displ,mdispl)
      enddo
c--------------------------------------------------------------------
      call MPI_REDUCE(delx,tdelx,1,MPI_DOUBLE_PRECISION,MPI_SUM,root,
     &                MPI_COMM_WORLD,ierror)
      delx = tdelx
      call MPI_BCAST(delx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &          ierror)
      lnp = np
      call MPI_REDUCE(lnp,tlnp,1,MPI_INTEGER8,MPI_SUM,root,
     &                MPI_COMM_WORLD,ierror)
      call MPI_BCAST(tlnp,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierror)
      if(myid.eq.0)print *, 'P',myid,delx,pfact, tlnp
      if(myid.eq.0)print *,'rms displacement=',dsqrt(delx)*pfact
     &           /sqrt(real(tlnp))
      call MPI_REDUCE(mdispl,tmdispl,1,MPI_REAL,MPI_MAX,0,
     &                MPI_COMM_WORLD,ierror)
      if(myid.eq.0)print *, 'Maximum displacement =',tmdispl


c***  Energy
      tmpken=0.D0
c$omp parallel do reduction(+:tmpken)
      do 70 i=1,np
          tmpken=tmpken+ 
     &        (particles(nxpos+4,i)**2.+particles(nxpos+5,i)**2.
     &        +particles(nxpos+6,i)**2.)
70    continue
      call MPI_REDUCE(tmpken,tken,1,MPI_DOUBLE_PRECISION,MPI_SUM,root,
     &                MPI_COMM_WORLD,ierror)
      tmpken = tken
      call MPI_BCAST(tmpken,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &                ierror)
      if(myid.eq.0)then
      print *,'sum of dx/da**2=',tmpken
      tmpken=0.5D0*8.D0*pi/3.D0*(1.D0/ai+1.D0/omei-1.)*tmpken
C
C possible error?
C
      ktot = 0.
      poten = 0.
      const = ai**3.*tmpken+poten+ktot
      print *,' poten,ken,const=',poten,tmpken,const
      endif
      ken = tmpken
 
    
      if(myid .eq. 0) 
     &     print *,'Generating Initial Condition ',etime(CPU)-stime

      
      return
      end
c----- First initialize boundary particles and migrate them
c----- to neighbor domains.
      subroutine setboundaryparticlesforce(np,particles,nn,nxpos,
     &     nzparticlewidth,zstartpos,
     &   nzprefixel,omep,omeplam,ken,zstart,zwidth,vamp,nspace,
     &     omei,amax,fx,fy,fz)
c      IMPLICIT none
      include 'mpif.h'
#include "fortran.i"
      include "fftw_f77.i"
      real omei,amax
      integer*8 lnp,tlnp
      double precision delx,tdelx,tmpken,tken
       integer np,nn,nxpos,nspace
      integer nx,ny,nz,nzper
       real omep,omeplam,ken,zstart,zwidth,vamp
       real particles(nn,np)
       real fx(2*(nx/2+1),ny,*)
       real fy(2*(nx/2+1),ny,*)
       real fz(2*(nx/2+1),ny,*)
       real vfact,vfact2
       integer lindx,i,j,k
       real zstartpos,pfact,ai,pi
      real displ,mdispl,tmdispl
       integer nzparticlewidth,nzprefixel,ngz
       integer NZBUFFER,rng
      parameter(NZBUFFER=4)
c     integer xabp0(-10:NNX+10),xabp1(-10:NNX+10),xabm1(-10:NNX+10)
c     integer yabp0(-10:NNY+10),yabp1(-10:NNY+10),yabm1(-10:NNY+10)
c     integer zabp0(-10:LOCAL_NNZ+NZBUFFER+10)
c     integer zabp1(-10:LOCAL_NNZ+NZBUFFER+10)
c     integer zabm1(-10:LOCAL_NNZ+NZBUFFER+10)
c     common /ptables/ xabp0,xabp1,xabm1,yabp0,yabp1,yabm1,zabp0,
c    &       zabp1,zabm1
      real rnx,rny,rnz
      integer local_ny_after_transpose,local_nz,local_z_start,
     &     local_y_start_after_transpose,total_local_size
      common /pmcommon/plan,iplan,nx,ny,nz,local_nz,local_z_start,
     &    local_ny_after_transpose,local_y_start_after_transpose,
     &    total_local_size
      integer myid,nid,ierror
      integer status(MPI_STATUS_SIZE),root
      real ktot,cpu(2),stime,etime
      real const,poten
      integer ii,jj,kk,iim1,jjm1,kkm1,iim2,jjm2,kkm2,iip1,iip2,jjp1,jjp2
      integer kkp1,kkp2
      integer local_mz
      real zfinal

      pi = 3.141592654
      root = 0
      zfinal = zstart+zwidth
      local_mz = int(zfinal) - int(zstart) + 1
      if(zfinal .eq. real(int(zfinal))) local_mz = local_mz - 1



      call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierror)
      call MPI_Comm_size(MPI_COMM_WORLD,nid,ierror)

      ai = 1.
      
      rnx = nx
      rny = ny
      rnz = nz
      rng = rnz
c---  Periodic table for mesh
ccc   local_mz+NZBUFFER is the local periodic condition in the current processor
      nzper = local_mz+NZBUFFER
c     do 4 i=-10,0
c     xabp0(i)=nx+i
c     xabp1(i)=nx+i+1
c     yabp0(i)=ny+i
c     yabp1(i)=ny+i+1
c     zabp0(i)=nzper+i
c   4 zabp1(i)=nzper+i+1
c     xabp1(0)=1
c     yabp1(0)=1
c     zabp1(0)=1

c     do 5 i=1,nx
c     xabp0(i)=i
c   5 xabp1(i)=i+1
c     xabp1(nx)=1
c     do 51 i=1,ny
c     yabp0(i)=i
c  51 yabp1(i)=i+1
c     yabp1(ny)=1

c     do 6 i=1,nzper
c     zabp0(i) = i
c   6 zabp1(i) = i+1
c     zabp1(nzper)=1
c     do 7 i=nx+1,nx+10
c     xabp0(i)=i-nx
c   7 xabp1(i)=i+1-nx
c     do 71 i=ny+1,ny+10
c     yabp0(i)=i-ny
c  71 yabp1(i)=i+1-ny

c     do 8 i=nzper+1,nzper+10
c     zabp0(i)=i-nzper
c   8 zabp1(i)=i+1-nzper
c------------------------------------------------------------------------
c---  For positions we use the pixel size as the length unit
c     vfact = (4./3.)*rng/2.*vamp
c     vfact2 = -(1./3.)*rng/4.*vamp
      vfact = vamp
      np = 0
      nzstride = max(nzparticlewidth-1,1)
      do k=1,nzparticlewidth,nzstride
      do j=1,ny/nspace
      do i=1,nx/nspace
      np = np + 1
c--------------------------------------------------------------
#ifdef XYZDBL
      particles(nxpos+1,np) = 0
      particles(nxpos+2,np) = 0
      particles(nxpos+3,np) = 0
#else
      particles(nxpos+1,np) = nspace*(i-1.)
      particles(nxpos+2,np) = nspace*(j-1.)
      particles(nxpos+3,np) = nspace*(k-1.)+zstartpos
#endif
c--------------------------------------------------------------
      ii = nspace*(i-1)+1
      jj = nspace*(j-1)+1
      kk = nspace*(k-1)+1+nzprefixel
c--------------------------------------------------------------
      particles(nxpos+4,np) = vfact*fx(ii,jj,kk)
      particles(nxpos+5,np) = vfact*fy(ii,jj,kk)
      particles(nxpos+6,np) = vfact*fz(ii,jj,kk)
c--------------------------------------------------------------
      enddo
      enddo
      enddo
c--------------------------------------------------------------------
#ifdef INDEX
      call indexingbp(particles,np,nzparticlewidth,nspace,nx,ny,nz)
#endif
c--------------------------------------------------------------------
c     pfact = rng*ai/vamp
      pfact = rng*ai
c$omp parallel do private(i)
      do i = 1, np
#ifdef XYZDBL
         particles(nxpos+1,i) = particles(nxpos+1,i) +
     &                             particles(nxpos+4,i)*pfact
         particles(nxpos+2,i) = particles(nxpos+2,i) +
     &                             particles(nxpos+5,i)*pfact
         particles(nxpos+3,i) = particles(nxpos+3,i) +
     &                             particles(nxpos+6,i)*pfact
#else
         particles(nxpos+1,i) = amod(particles(nxpos+1,i)+
     &         particles(nxpos+4,i)*pfact+rnx,rnx)
         particles(nxpos+2,i) = amod(particles(nxpos+2,i)+
     &         particles(nxpos+5,i)*pfact+rny,rny)
         particles(nxpos+3,i) = amod(particles(nxpos+3,i)+
     &         particles(nxpos+6,i)*pfact+rnz,rnz)
#endif
      enddo

      return
      end

c---  Growth factor in the universe with zero curvature
c---  omep= present matter density parameter
c---  omeplam = 1-omep is assumed (zero curvature inflation model)
      function growth(a,amax,omep)
      external coscon
      common/pass/x

      x0 = (1./omep-1)**(1./3.)
      red = amax/a-1
      x = (1./omep-1)**(1./3.)/(1+red)
      call qsimp(coscon,0.,1.,gfac)
      growth = (x/x0)*sqrt(x**3+1)*gfac
      return
      end

      function coscon(y)
      common/pass/x
      coscon = 1./(sqrt(1+x**3*y**1.2))**3
      return
      end
