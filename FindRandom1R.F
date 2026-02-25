c234567
      parameter(IMAX=1000)
      parameter(mx=10,my=5)
      integer iseed,iiseed
      integer nx,ny,nz
      character*80 instring

      real xxx
      real den(mx,my,5)
      real pk(512)
      real npk(512)

      call getarg(1,instring)
      read(instring,*) nx
      ny = nx





      do kk = 1, IMAX
         iseed = -kk
         xxx = ran3(iseed)
         do i = 1, 10
         do j = 1, ny
         do k = 1, 5
            if(i .le. mx .and. j .le. my) then
                den(i,j,k) = gasdev(iseed)
            else
                xxx = gasdev(iseed)
            endif
         enddo
         enddo
         enddo

         do i = 1, 512
            pk(i)  =0
            npk(i)  =0
         enddo

         do k = 1, 5
         wk = k-1
         do j = 1, my
         wj = j-1
         do i = 0, mx/2-1
         delksq = den(2*i+1,j,k)**2 + den(2*i+2,j,k)**2
         delksq = delksq /2.
         wave = sqrt(i*i+wj*wj+wk*wk)
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

         enddo
         enddo
         enddo

         do i=2,6
         if(npk(i).gt.0)then
            pk(i)=pk(i)/npk(i)
         else
            pk(i) = 0.
         end if
         enddo
         write(6,'i4,5(1x,g14.5)') -kk,(pk(i),i=2,6)
      enddo
      stop
      end

      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran3(idum)-1.
        v2=2.*ran3(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
