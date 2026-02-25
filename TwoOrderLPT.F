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
      subroutine TwoLPTpmseedforce(den,omep,
     &     omeplam,wlam0,wlam1, Cs2,ken,
     &     denpad,zstart,zwidth,damp1,vamp1,damp2,vamp2,omei,fx,fy,den2,
     &     fNL,gNL)
      IMPLICIT none
      include 'mpif.h'
      include "fftw_f77.i"
#include "fortran.i"
      real iceil,fNL,gNL
      external iceil
      integer nzsend,noffset
      real zstart,zwidth,zfinal
      integer izwidth,izstart
      integer m2z,p2z
      character*1 dir1,dir2,dir3,dir4,dir5
      integer TID,NTHREADS,OMP_GET_THREAD_NUM
      integer OMP_GET_NUM_THREADS,myistart,myiend
      real bpowernow

      integer*8 np
      integer mbias,NZBUFFER
      parameter(mbias=4,NZBUFFER=4)
      integer nx,ny,nz,local_nz,local_z_start
      integer myid,nid,dest,src,stag,rtag
      real den(2*(nx/2+1),ny,local_nz+NZBUFFER),denpad(2,ny,nz)
      real fx(2*(nx/2+1),ny,*)
      real fy(2*(nx/2+1),ny,*)
      real den2(2*(nx/2+1),ny,*)
      real displ,mdispl
      integer status(MPI_STATUS_SIZE),root

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
      real smooth,h,rth,bias,zi, As
      integer knyq,ngz,nnp,nznp,ng,ngx,ngy
      real npow
      real rngx,rngy
      real xyz,xyy,omei,powh,pi,omehs,rknyq,rknyq2,rkinv
      real twopi,rngz,rngzc,ran3,xxx,rng,rngc
      real pmas,pmass,ai,anow,wavemax,pamp,damp1,damp2,vamp1,vamp2
      real tmppamp
      real xnow, xini,bnow,bini
      real vfact,vfact2,pfact,amax,astep,omep,omeplam,wlam0,wlam1
      real Cs2
      real efold,sfac
      real const,poten,growth,red,growth2,growthgen, growthcpl
      integer iindex

      external gasdev
      integer ngh,i,j,k,ky,kz,nghx,nghy,nghz
      real rk,velpot,gasdev,xksq,power,bpower,bpower1,growthfactor
      integer sendsize,recvsize,ierror
  
      integer ii,jj,kk,iim1,jjm1,kkm1,iim2,jjm2,kkm2
      integer iip1,jjp1,kkp1,iip2,jjp2,kkp2

      integer nadd,nzparticles,nzprefixel
      real work
      real*8 rincube
      real omepk
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
      character*80 powfilename,asciipkfile
      real asciirk(8192), asciipk(8192)
      integer nascii
      integer powreadflag
      common /ASCIIPOWER/ powreadflag,nascii,asciirk,asciipk

      integer seedflag, PkRescalingFlag
      real*8 cubicsize


      call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierror)
      call MPI_Comm_size(MPI_COMM_WORLD,nid,ierror)

      if(myid .eq. 0) stime = etime(CPU)


C This is a hook into the C functions
      call getpowerparams(nspace,size,amax,astep,
     &                        npow,smooth,h,rth,bias,As,iseed,omepb,
     &                    powreadflag,powfilename,asciipkfile,
     &                        PkRescalingFlag)
      print *, myid,zstart,zwidth

      dark = 'cdm'
      com_omep=omep
      com_omepb = omepb
      com_omeplam=omeplam
      com_h = h
      com_size = size
      cubicsize = size**3
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
      omei = omep*(1+zi)**3
     &     /(omep*(1+zi)**3+omeplam*(1+zi)**(3+3*wlam0)+
     &         (1-omep-omeplam)*(1+zi)**2)
      omehs = 1./(omep*h*size)
      smooth = (smooth/sqrt(2.))/rng/sqrt(2.)
      rth = rth/rng
c---  Total mass is ngz/ng so that the mean density = 1
      pmas  = real(nspace)**3
      pmass = pmas/rngc


      wavemax = 100000


#ifdef EH98
c     THIS IS THE POWER FIT BY EISENSTEIN & HU 1998
      pamp = 0
      w1 = 0.001*size
      do i = 1, 5
         w2 = w1*10.
         call qromb(bpowint,w1,w2,tmppamp)
         pamp = pamp + tmppamp
         w1 = w2
      enddo
#else
c     THIS IS THE POWER CALCULATED BY THE CAMB
      if(powreadflag .eq. 1) then
         do i=0,nid-1
            if(myid.eq.i) then
               call rInitPower(omep,omepb,omeplam, wlam0,wlam1, Cs2,h,
     &            npow,zi,  As,pamp,size,rng,
     &            powreadflag,powfilename,myid)
            endif
            call MPI_Barrier(MPI_COMM_WORLD,ierror)
         enddo
       else  if(powreadflag .eq. 2) then
            if(myid.eq.0) then
               open(91,file=asciipkfile)
               read(91,*)
               read(91,*)
               read(91,*)
               read(91,*)
               i = 1
899            read(91,*,end=900) asciirk(i), asciipk(i)
               if(i .gt. 8192) then
                  print *,'overflow when reading ascii P(k) ',i
                  stop
               endif
               i = i + 1
               goto 899
900            close(91)
               nascii = i - 1
               print *, 'Total Ascii lines are ', nascii
            endif
            call MPI_Bcast(nascii,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
            call MPI_Bcast(asciirk,nascii,MPI_REAL,0,MPI_COMM_WORLD,
     &           ierror)
            call MPI_Bcast(asciipk,nascii,MPI_REAL,0,MPI_COMM_WORLD,
     &           ierror)
            ! This is to measure the amplitude of the input power which is modified.
            !  0.0001*boxsize <w<1000*boxsize where w = u. Therefore, 2pi*0.0001< k < 2pi*1000.
            pamp = 0
            w1 = 0.0001*size
            do i = 1, 6 ! integration from k=0.0001*2pi to 1000*2pi.
               w2 = w1*10.
               call qromb(bpowint,w1,w2,tmppamp)
               pamp = pamp + tmppamp
               w1 = w2
            enddo
            ! k/2pi = u/Lbox ==> k^2dk/(2pi)^3 = u^2du/Lbox^3

      else if(powreadflag .eq.0) then
         call InitPower(omep,omepb,omeplam,wlam0,wlam1,Cs2,h,npow,zi,
     &            As, pamp,size,rng)
      else 
         print *,'error in powreadflag ', powreadflag
         call MPI_Finalize()
         stop
      endif


#endif



c     call qsimp(pow2int,0.00001,wavemax,pamptmp)
      if(myid .EQ. 0) print *,' <1>=',pamp
c---  Here pamp = initial ms of top-hat smoothed density fluctuation
c---  pamp = 50.06, 800.95, 12815 for 64, 128, 256 h^-1 Mpc
c---  Now pamp*sqrt(P(k)) give rms fluctuation equal to thamp initially.
c---  After the expansion factor grows by the factor anow, the linear
c---  amplitude of fluctuation  will be (1/bias) over the scale rth.
c---  I. Einstein-de Sitter model
      omepk = 1-omep-omeplam

       if(wlam0 .ne.-1 .and. omepk .ne.0) then
          print *, 'Error! Beyond the parameter ranges!'
          call MPI_Finalize()
          stop
c      CPL model
       else if(wlam1 .ne. 0) then
          bini = growthcpl(1., amax, omep, omeplam, wlam0, wlam1,Cs2)
          bnow = growthcpl(amax, amax, omep, omeplam, wlam0, wlam1,Cs2)
          pamp = (1./bias)/sqrt(pamp)*bini/bnow
          damp1 = 1
          vamp1=(growthcpl(1.+0.005,amax,omep,omeplam,wlam0, wlam1,Cs2)-
     &          growthcpl(1.-0.005,amax,omep, omeplam,wlam0, wlam1,Cs2))
     &        /(0.01)*1./bini
          if(myid.eq.0) then
             print *, 'wlambda0,1 case: ',omep,' bini=', bini, 
     &             ' bnow=',bnow, vamp1
          endif
c      Quintessence model
       else if(wlam0 .ne.-1) then
         bini = growthgen(1., amax,omep,omeplam,wlam0, Cs2)
         bnow = growthgen(amax,amax,omep,omeplam, wlam0, Cs2)
         pamp = (1./bias)/sqrt(pamp)*bini/bnow
         damp1 = 1
         vamp1 = (growthgen(1.+0.005, amax,omep,omeplam,wlam0, Cs2)-
     &        growthgen(1.-0.005, amax,omep,omeplam,wlam0, Cs2))
     &       /(0.01)*1./bini
         if(myid.eq.0) then
            print *,'wlambda0 case: ',omep,' bini=', bini, 
     &            ' bnow=', bnow,vamp1
         endif

       else if(omep.eq.1.)then
         pamp = (1./bias)/sqrt(pamp) *ai/anow
         damp1 = 1.
         vamp1 = 1.
c---  II. Open model
       else if(omep.lt.1.and.omeplam.eq.0)then
         xnow = (1.-omep)/omep
         xini = xnow/(1.+zi)
         bini = 1.+3./xini+3.*sqrt(1.+xini)/xini**1.5*alog(sqrt(1.+xini)
     &         -sqrt(xini))
         bnow = 1.+3./xnow+3.*sqrt(1.+xnow)/xnow**1.5*alog(sqrt(1.+xnow)
     &         -sqrt(xnow))
         pamp = (1./bias)/sqrt(pamp) *bini/bnow
         damp1 = 1.
         vamp1 = omei**0.6
       else if(omep.lt.1.and.omeplam.gt.0.and.abs(omepk) .le. 1.E-5)then
c---  III. Zero curvature with cosmological constant model
        bini = growth(1.,  amax,omep)
        bnow = growth(amax,amax,omep)
        pamp = (1./bias)/sqrt(pamp)*bini/bnow
        damp1 = 1.
        vamp1 = (growth(1.+0.005,amax,omep)-growth(1.-0.005,amax,omep))
     &        /(0.01)*1./bini
        if(myid.eq.0) then
           print *,omep,' bini=',bini,' bnow=',bnow,vamp1
           print *, 'pamp = ', pamp
        endif
       else if(omeplam .gt. 0 .and. abs(omepk) .gt.1.E-5) then
c---  IV. Non-zero curvature with cosmological constant model
        bini = growth2(1.,  amax,omep,omeplam)
        bnow = growth2(amax,amax,omep,omeplam)
        pamp = (1./bias)/sqrt(pamp)*bini/bnow
        damp1 = 1.
        vamp1=(growth2(1.+0.005,amax,omep)-growth2(1.-0.005,amax,omep))
     &        /(0.01)*1./bini
        if(myid.eq.0) then
           print *,omep,' bini=',bini,' bnow=',bnow,vamp1
        endif
       end if

      if(PkRescalingFlag .eq. 0) then
         pamp = 1./dsqrt(cubicsize) ! from k space to u space conversion factor
         if(myid.eq.0) print *, 'normalized by pamp = ', pamp
      endif



      ngz = nz

      efold = 0.585*sqrt(2.)
      sfac   = -(pi*efold/size)**2
c***  power spectrum amplitude
      pamp = pamp*rngzc*sqrt(0.5)

c***  Velocity potential coefficients
      damp2 = -3./7.
      vamp2 = damp2 * 2 * omei**(4./7.)

ccc   P(k) = (pamp*bpower1(rk))**2 * boxsize**3; here rk is k and rk in the
ccc   code is 2pi*u

901   continue

!     if(myid.eq.0) then
!        velpot = bpower(twopi*10)
!     endif
!     stop


      if(seedflag .ne. 0) then
      xxx = ran3(iiseed)
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
              den(2*i+1,j,k) = gasdev(iiseed)
              den(2*i+2,j,k) = gasdev(iiseed)
            end if
         enddo
      enddo
      enddo
      else
      iiseed = iseed*(myid+1)
      xxx = ran3(iiseed)
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
              velpot = pamp*bpower(rk)/rk**2
              den(2*i+1,j,k) = velpot*gasdev(iiseed)
              den(2*i+2,j,k) = velpot*gasdev(iiseed)
            end if
         enddo
      enddo
      enddo

      endif

      if(myid .EQ. 0) print *,'last #=',ran3(iiseed), bpower(twopi)

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


c--   TwoLPT has an input array of 1'st order displacement potential in Fourier space
c     and returns an array of 2'nd order displacement potential in Fourier space.
      growthfactor = bini/bnow
      call TwoLPT(den,fx,fy,den2,myid,nid,fNL,gNL,size,amax,omep,
     &     Growthfactor)

       if(myid.eq.0) print *,'Ending of saving initial force array'



       return
       end

