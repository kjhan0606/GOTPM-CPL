c---  Reads position or position & velocity of particles
c---  npkid = number of particles in the subset
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
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


      function power(rk)
c---  Calculates sqrt(Power spectrum) with gaussian smoothing and with
c---  biquadratic cutoff at Nyquist frequency
c---  Power spectrum consists of primordial scale-free power law and 
c---  transfer function
      character*3 dark
      common/powsp/rknyq,rkinv,powh,smooth,omehs,dark

c***  Primordial scale-free power spectrum : power = sqrt(P(k))
c---  rk = 2*pi*(wave number) & rkinv = 1./(2*pi*knyq)
c---  smooth = (smoothing length in pixels /ng) /sqrt(2)
      if(rk.gt.rknyq.or.rk.eq.0.)then
      power=0.
      return
      else
c      power=rk**powh*(1.-(rk*rkinv)**4.)*exp(-(rk*smooth)**2.)
      power=rk**powh
      end if



c***  Transfer function (BBKS 1986, p60) 
c---  omehs = 1/(omega*h*v[100])
      qq = rk*omehs
c---  Cold dark matter, adiabatic fluctuations with 3 flavors of 
c---  relativistic neutrinos plus photons and negligible baryon
      if(dark.eq.'cdm')then
          power=power*alog(1+2.34*qq)/(2.34*qq)*(1+3.89*qq
     &     +(16.1*qq)**2+(5.46*qq)**3+(6.71*qq)**4)**(-0.25)

c---  Hot dark matter, adiabatic fluctuations with 2 species of
c---  relativiatic neutrinos
      else if(dark.eq.'hdm')then 
          damp=2.6*qq
          power=power*exp(-0.16*damp-0.5*damp**2)/(1+1.6*qq
     &     +(4.*qq)**1.5+(0.92*qq)**2)
      end if
ccccccccccccccccccccc a cutoff at k=8
c     if(rk/2./3.141592654.gt.8)power = 0
      return
      end


c---  Integrand for calculating (dM/M)**2 in top hat window
      function powint(xk)
      common/tophat/rth
      data twopi/6.28318531/,fopi/12.566371/

c---  xk = wave number
      rk = twopi*xk
c---  Power spectrum(without knyquist cutoff) times volumn factor
      powint = fopi*xk*xk*(power1(rk))**2
c---  Fourier transform of top-hat window function
      xk1   = rth*rk
      if(xk1.gt.0.02)then
      powint = powint*9*((sin(xk1)-xk1*cos(xk1))/xk1**3)**2
      else
      powint = powint*(1-xk1*xk1*0.1)**2
      end if
      return
      end

c---  Integrand for calculating (dM/M)**2 in top hat window
      function bpowint(xk)
      common/tophat/rth
      data twopi/6.28318531/,fopi/12.566371/


c---  xk = wave number (i**2+j**2+k**2)
      rk = twopi*xk
c---  Power spectrum(without knyquist cutoff) times volumn factor
      powint = fopi*xk*xk*(bpower1(rk))**2
c---  Fourier transform of top-hat window function
      xk1   = rth*rk
      if(xk1.gt.0.02)then
      bpowint = powint*9*((sin(xk1)-xk1*cos(xk1))/xk1**3)**2
      else
      bpowint = powint*(1-xk1*xk1*0.1)**2
      end if
      return
      end
c---  Integrand for calculating integral(P(k)*k^2) in top hat window
      function bpow2int(xk)
      common/tophat/rth
      data twopi/6.28318531/,fopi/12.566371/

c---  xk = wave number
      rk = twopi*xk
c---  Power spectrum times volumn factor
      powint = fopi*xk**2*rk**2*(bpower1(rk))**2
c---  Fourier transform of top-hat window function
      xk1   = rth*rk
      if(xk1.gt.0.02)then
      bpow2int = powint*9*((sin(xk1)-xk1*cos(xk1))/xk1**3)**2
      else
      bpow2int = powint*(1-xk1*xk1*0.1)**2
      end if
      return
      end


c---  Integrand for calculating integral(P(k)*k^2) in top hat window
      function pow2int(xk)
      common/tophat/rth
      data twopi/6.28318531/,fopi/12.566371/

c---  xk = wave number
      rk = twopi*xk
c---  Power spectrum times volumn factor
      powint = fopi*xk**2*rk**2*(power1(rk))**2
c---  Fourier transform of top-hat window function
      xk1   = rth*rk
      if(xk1.gt.0.02)then
      pow2int = powint*9*((sin(xk1)-xk1*cos(xk1))/xk1**3)**2
      else
      pow2int = powint*(1-xk1*xk1*0.1)**2
      end if
      return
      end

      function power1(rk)
c---  Calculates sqrt(Power spectrum) with gaussian smoothing and with
c---  biquadratic cutoff at Nyquist frequency
c---  Power spectrum consists of primordial scale-free power law and 
c---  transfer function
      character*3 dark
      common/powsp/rknyq,rkinv,powh,smooth,omehs,dark

c***  Primordial scale-free power spectrum : power = sqrt(P(k))
c---  rk = 2*pi*(wave number) & rkinv = 1./(2*pi*knyq)
c---  smooth = (smoothing length in pixels /ng) /sqrt(2)
      power=rk**powh*exp(-(rk*smooth)**2)

c***  Transfer function (BBKS 1986, p60) 
c---  omehs = 1/(omega*h*v[100])
      qq = rk*omehs
c---  Cold dark matter, adiabatic fluctuations with 3 flavors of 
c---  relativistic neutrinos plus photons and negligible baryon
      if(dark.eq.'cdm')then
      if(qq.lt.0.002)then
          power=power*(1-1.17*qq+(1.351*qq)**2)*(1+3.89*qq
     &     +(16.1*qq)**2+(5.46*qq)**3+(6.71*qq)**4)**(-0.25)
      else
          power=power*alog(1+2.34*qq)/(2.34*qq)*(1+3.89*qq
     &     +(16.1*qq)**2+(5.46*qq)**3+(6.71*qq)**4)**(-0.25)
      end if

c---  Hot dark matter, adiabatic fluctuations with 2 species of
c---  relativiatic neutrinos
      else if(dark.eq.'hdm')then 
          damp=2.6*qq
          power=power*exp(-0.16*damp-0.5*damp**2)/(1+1.6*qq
     &     +(4*qq)**1.5+(0.92*qq)**2)
      end if
      power1 = power
ccccccccccccccccccccccccccc
c     if(rk/2./3.141592654.gt.8)power = 0
      return
      end

#ifdef EH98
      function bpower(rk)
c---  Calculates sqrt(Power spectrum) with gaussian smoothing and with
c---  cutoff at Nyquist frequency
c---  Power spectrum consists of primordial scale-free power law and
c---  transfer function
c---  Ref. Eisenstein and Hu, 1998, ApJ, 496,605
      character*3 dark
      common/powsp/rknyq,rkinv,powh,smooth,omehs,dark
      common/cospar/omep,omepb,omeplam,h,ssize
      save in1,cmbt27,o,ob,oc,zeq,rkeq,b1,b2,zd,Rd,Req,sound,rksk
      save y,gy,alphab,a11,a22,alphac,b11,b22,betac,betab,betan
      data in1/0/

c---  Important cosmological scales and epochs
c---  Tcmb = 2.728+-0.0004K (Fixen et al. 1996) adopted
c---  rkeq = scale of particle horizon at equality epoch in unit of [h Mpc^-1]
c---  zd = drag epoch (baryon released from the drag of photons)
c---  sound = sound horizon at drag epoch in [h^-1 Mpc]
c---  rksk = Silk damping scale in unit of [h Mpc^-1]
      if(in1.eq.0)then
      cmbt27 = 2.728/2.7
      o   = omep*h**2
      ob  = omepb*h**2
      oc  = (omep-omepb)*h**2

      zeq = 2.50e4*o/cmbt27**4
      rkeq = 7.46e-2*omep*h/cmbt27**2
      b1  = 0.313/o**0.419*(1+0.607*o**0.674)
      b2  = 0.238*o**0.223
      zd  = 1291*o**0.251/(1+0.659*o**0.828)*(1+b1*ob**b2)
      Rd  = 31.5*ob/cmbt27**4/(zd/1000.)
      Req = 31.5*ob/cmbt27**4/(zeq/1000.)
      sound=2./3./rkeq*sqrt(6./Req)
     &      *alog((sqrt(1+Rd)+sqrt(Rd+Req))/(1+sqrt(Req)))
      rksk = 1.6*ob**0.52*o**0.73*(1+(10.4*o)**(-0.95))/h

      y = (1+zeq)/(1+zd)
      gy = y*(-6*sqrt(1+y)+(2+3*y)*alog((sqrt(1+y)+1)/(sqrt(1+y)-1)))
      alphab = 2.07*rkeq*sound/(1+Rd)**(3./4.)*gy

      a11  = (46.9*o)**0.670*(1+(32.1*o)**(-0.532))
      a22  = (12.0*o)**0.424*(1+(45.0*o)**(-0.582))
      alphac = a11**(-ob/o)*a22**(-ob/o*ob/o*ob/o)
      b11  = 0.944/(1+(458*o)**(-0.708))
      b22  = 1./(0.395*o)**0.0266
      betac = 1./(1+b11*((oc/o)**b22-1))
      betab = 0.5+ob/o+(3-2*ob/o)*sqrt((17.2*o)**2+1)
      betan = 8.41*o**0.435

      in1 = 1
      end if

c***  Primordial scale-free power spectrum : power = sqrt(P(k))
c---  rk = 2*pi*(wave number) & rkinv = 1./(2*pi*knyq)
c---  smooth = (smoothing length in pixels /ng) /sqrt(2)
      if(rk.gt.rknyq.or.rk.eq.0.)then
      bpower=0.
      return
      else
c     power=rk**powh*(1.-(rk*rkinv)**4.)*exp(-(rk*smooth)**2.)
      power=rk**powh
      end if

c***  Baryonic effects included
      q   = rk*omehs*cmbt27**2
      rks = rk/ssize*sound
      f   = 1./(1+(rks/5.4)**4)
      transc = f*tnot(q,1.0,betac)+(1-f)*tnot(q,alphac,betac)

      soundt = sound/(1+(betan/rks)**3)**(1./3.)
      rkstilde = rk/ssize*soundt
      rkrksk = rk/ssize/rksk
      transb = (tnot(q,1.0,1.0)/(1+(rks/5.2)**2)
     & +alphab/(1.+(betab/rks)**3)*exp(-rkrksk**1.4))
     & *sin(rkstilde)/rkstilde
      
      bpower = power*abs(ob*transb+oc*transc)/o

c***  Transfer function (BBKS 1986, p60)
c---  omehs = 1/(omega*h*v[100])
      qq = rk*omehs
c---  Cold dark matter, adiabatic fluctuations with 3 flavors of
c---  relativistic neutrinos plus photons and negligible baryon
c         power=power*alog(1+2.34*qq)/(2.34*qq)*(1+3.89*qq
c    &     +(16.1*qq)**2+(5.46*qq)**3+(6.71*qq)**4)**(-0.25)

      return
      end

      function tnot(q,alphac,betac)
c     data e/2.718281828/
      e = 2.718281828
      c = 14.2/alphac+386./(1+69.9*q**1.08)
      tnot = alog(e+1.8*betac*q)/(alog(e+1.8*betac*q)+c*q**2)
      return
      end

      function bpower1(rk)
c---  Calculates sqrt(Power spectrum) with gaussian smoothing and with
c---  cutoff at Nyquist frequency
c---  Power spectrum consists of primordial scale-free power law and
c---  transfer function
c---  Ref. Eisenstein and Hu, 1998, ApJ, 496,605
      character*3 dark
      common/powsp/rknyq,rkinv,powh,smooth,omehs,dark
      common/cospar/omep,omepb,omeplam,h,ssize
      save in1,cmbt27,o,ob,oc,zeq,rkeq,b1,b2,zd,Rd,Req,sound,rksk
      save y,gy,alphab,a11,a22,alphac,b11,b22,betac,betab,betan
      data in1/0/

c---  Important cosmological scales and epochs
c---  Tcmb = 2.728+-0.0004K (Fixen et al. 1996) adopted
c---  rkeq = scale of particle horizon at equality epoch in unit of [h Mpc^-1]
c---  zd = drag epoch (baryon released from the drag of photons)
c---  sound = sound horizon at drag epoch in [h^-1 Mpc]
c---  rksk = Silk damping scale in unit of [h Mpc^-1]
      if(in1.eq.0)then
      cmbt27 = 2.728/2.7
      o   = omep*h**2
      ob  = omepb*h**2
      oc  = (omep-omepb)*h**2

      zeq = 2.50e4*o/cmbt27**4
      rkeq = 7.46e-2*omep*h/cmbt27**2
      b1  = 0.313/o**0.419*(1+0.607*o**0.674)
      b2  = 0.238*o**0.223
      zd  = 1291*o**0.251/(1+0.659*o**0.828)*(1+b1*ob**b2)
      Rd  = 31.5*ob/cmbt27**4/(zd/1000.)
      Req = 31.5*ob/cmbt27**4/(zeq/1000.)
      sound=2./3./rkeq*sqrt(6./Req)
     &      *alog((sqrt(1+Rd)+sqrt(Rd+Req))/(1+sqrt(Req)))
      rksk = 1.6*ob**0.52*o**0.73*(1+(10.4*o)**(-0.95))/h

      y = (1+zeq)/(1+zd)
      gy = y*(-6*sqrt(1+y)+(2+3*y)*alog((sqrt(1+y)+1)/(sqrt(1+y)-1)))
      alphab = 2.07*rkeq*sound/(1+Rd)**(3./4.)*gy

      a11  = (46.9*o)**0.670*(1+(32.1*o)**(-0.532))
      a22  = (12.0*o)**0.424*(1+(45.0*o)**(-0.582))
      alphac = a11**(-ob/o)*a22**(-ob/o*ob/o*ob/o)
      b11  = 0.944/(1+(458*o)**(-0.708))
      b22  = 1./(0.395*o)**0.0266
      betac = 1./(1+b11*((oc/o)**b22-1))
      betab = 0.5+ob/o+(3-2*ob/o)*sqrt((17.2*o)**2+1)
      betan = 8.41*o**0.435

      in1 = 1
      end if

c***  Primordial scale-free power spectrum : power = sqrt(P(k))
c---  rk = 2*pi*(wave number) & rkinv = 1./(2*pi*knyq)
c---  smooth = (smoothing length in pixels /ng) /sqrt(2)
c     if(rk.gt.rknyq.or.rk.eq.0.)then
      if(rk.eq.0.)then
      bpower1=0.
      return
      else 
c     power=rk**powh*(1.-(rk*rkinv)**4.)*exp(-(rk*smooth)**2.)
      power=rk**powh
      end if

c***  Baryonic effects included
      q   = rk*omehs*cmbt27**2
      rks = rk/ssize*sound
      f   = 1./(1+(rks/5.4)**4)
      transc = f*tnot(q,1.0,betac)+(1-f)*tnot(q,alphac,betac)

      soundt = sound/(1+(betan/rks)**3)**(1./3.)
      rkstilde = rk/ssize*soundt
      rkrksk = rk/ssize/rksk
      transb = (tnot(q,1.0,1.0)/(1+(rks/5.2)**2)
     &         +alphab/(1.+(betab/rks)**3)*exp(-rkrksk**1.4))
     &         *sin(rkstilde)/rkstilde
      
      bpower1 = power*abs(ob*transb+oc*transc)/o
      return
      end
#else
c---------------
c     rk is the 2pi * u (or xk) where u is the wave number
c     rks is the k = 2pi * u /boxsize, real wave vector
      function bpower1(rk)
      character*3 dark
      common/powsp/rknyq,rkinv,powh,smooth,omehs,dark
      common/cospar/omep,omepb,omeplam,h,ssize
      data in1/0/
      real matterpk(8192),rkarray(8192),Trans(8192),minkh,dlnkh
      integer point,nmatterpk,nn,nmin,nmax
      save in1
      common /cambpower/ matterpk,rkarray,nmatterpk,minkh,dlnkh
      real asciirk(8192), asciipk(8192)
      integer nascii,powreadflag
      common /ASCIIPOWER/powreadflag,nascii,asciirk,asciipk
      rks = rk/ssize
      if(in1 .eq.0) then
         if(powreadflag .eq.0 .or. powreadflag .eq. 1) then
            nmatterpk = 8192 ! this will be modified properly in GetPower()
            call GetPower(minkh,dlnkh,matterpk,nmatterpk,Trans)
            do i = 1, nmatterpk
               rkarray(i) = minkh*exp((i-1)*dlnkh)
               matterpk(i) = sqrt(matterpk(i))
            enddo
         else if(powreadflag .eq. 2) then
            minkh = 5e-5
            dlnkh = 0.02
            nmatterpk = 900
            do i = 1, nmatterpk
               rkarray(i) = minkh*exp((i-1)*dlnkh)
            enddo
            call getpkfromascii(rkarray,matterpk,nmatterpk,asciirk,
     &               asciipk,nascii,Trans)
            do i = 1, nmatterpk
               matterpk(i) = sqrt(matterpk(i))
            enddo
         endif
         in1 = 1
      endif


      if(rk .eq. 0 ) then
          bpower1 = 0
          return 
      endif 
      if(rks .lt. rkarray(1)) then
         print *,'############################################'
         print *,'############################################'
         print *,'Less than lowest boundary ', rks, rkarray(1)
         print *,'############################################'
         print *,'############################################'
         stop
      else if(rks .gt. rkarray(nmatterpk)) then
         print *,'############################################'
         print *,'############################################'
         print *,'Larger than highest boundary ',rks,rkarray(nmatterpk)
         print *,'############################################'
         print *,'############################################'
         stop
      endif
      point = (alog(rks/minkh)/dlnkh+1)
      nmin = max(1,point-2)
      nmax = min(nmatterpk,point+2)
      nn = nmax-nmin+1
      call polint(rkarray(nmin),matterpk(nmin),nn,rks,y,dy)
      bpower1 =y
      return
      end
c---------------
c     bpower = sqrt(PK)
c     rk = 2pi*u 
c     rks = 2pi*u /boxsize (real k)
c---------------
      function bpower(rk)
      character*3 dark
      common/powsp/rknyq,rkinv,powh,smooth,omehs,dark
      common/cospar/omep,omepb,omeplam,h,ssize
      data in1/0/
      real matterpk(8192),rkarray(8192),Trans(8192),minkh,dlnkh
      integer point,nmatterpk,nn,nmin,nmax
      save in1
      common /cambpower/ matterpk,rkarray,nmatterpk,minkh,dlnkh
      real asciirk(8192), asciipk(8192)
      integer nascii,powreadflag
      common /ASCIIPOWER/powreadflag,nascii,asciirk,asciipk
      rks = rk/ssize
      if(in1 .eq.0) then
         if(powreadflag .eq.0 .or. powreadflag .eq. 1) then
            nmatterpk = 8192 ! this will be modified in GetPower() properly.
            call GetPower(minkh,dlnkh,matterpk,nmatterpk,Trans)
            do i = 1, nmatterpk
               rkarray(i) = minkh*exp((i-1)*dlnkh)
               matterpk(i) = sqrt(matterpk(i))
            enddo
         else if(powreadflag .eq. 2) then
            nmatterpk = 900
            minkh = 5e-5
            dlnkh = 0.02
            do i = 1, nmatterpk
               rkarray(i) = minkh*exp((i-1)*dlnkh)
            enddo
            call getpkfromascii(rkarray,matterpk,nmatterpk,asciirk,
     &               asciipk,nascii,Trans)
            do i = 1, nmatterpk
               rkarray(i) = minkh*exp((i-1)*dlnkh) ! reassure the value
               matterpk(i) = sqrt(matterpk(i))
            enddo
         endif
         in1 = 1
      endif


      if(rk .eq. 0 .or. rk .gt. rknyq) then
          bpower = 0
          return 
      endif 
      if(rks .lt. rkarray(1)) then
         print *,'############################################'
         print *,'############################################'
         print *,'Less than lowest boundary ', rks, rkarray(1)
         print *,'############################################'
         print *,'############################################'
         stop
      else if(rks .gt. rkarray(nmatterpk)) then
         print *,'############################################'
         print *,'############################################'
         print *,'Larger than highest boundary ',rks,rkarray(nmatterpk)
         print *,'############################################'
         print *,'############################################'
         stop
      endif
      point = (alog(rks/minkh)/dlnkh+1)
!     nmin = max(1,point-2)
!     nmax = min(nmatterpk,point+2)
!     nn = nmax-nmin+1
!     call polint(rkarray(nmin),matterpk(nmin),nn,rks,y,dy)
!     bpower = y
      bpower = (matterpk(point+1)-matterpk(point))
     &        /(rkarray(point+1)-rkarray(point))
     &        *(rks-rkarray(point))
     &        + matterpk(point)
      return
      end
c###############
      subroutine getpkfromascii(rk,pk,npk,ark,apk,napk,Ts2)
      integer npk,napk,n1,n2,nn
      real rk(npk),pk(npk),Ts2(npk),ark(napk),apk(napk)
      real y
c --- change scaling: note that the number of effective digits might be changed.
      do i = 1, npk
         rk(i) = alog10(rk(i))
      enddo
      do i = 1, napk
         ark(i) = alog10(ark(i))
         apk(i) = alog10(apk(i))
      enddo
      do i = 1, npk
         if(rk(i) .lt. ark(1)) then
            y = (apk(2)-apk(1))/(ark(2)-ark(1))*(rk(i)-ark(1)) + apk(1)
         else if(rk(i) .gt. ark(napk)) then
            y = (apk(napk)-apk(napk-1))/(ark(napk)-ark(napk-1))*
     &             (rk(i)-ark(napk)) + apk(napk)
         else
            do j = 2, napk-1
!              if(rk(i).ge. ark(j) .and. rk(i).le.ark(j+1)) goto 10
               if(rk(i).ge. ark(j) .and. rk(i).le.ark(j+1)) then
                  y = (apk(j+1)-apk(j))/(ark(j+1)-ark(j))*(rk(i)-ark(j))
     &                    + apk(j)
                  goto 10
               endif
            enddo
10          continue
!           n1 = j
!           n2 = j+1
!           n1 = max(n1,1)
!           n2 = min(n2,napk)
!           nn = n2 - n1 + 1
!           call polint(ark(n1),apk(n1),nn,rk(i),y,dy)
         endif
         pk(i) = 10**y
c        pk(i) = 10.**(y+rk(i))*Ts2(i) ! Prava's request
c        print *, 10**rk(i), pk(i),y,Ts2(i)
      enddo
c     stop
c --- change scaling: note that the number of effective digits might be changed.
      do i = 1, npk
         rk(i) = 10.**(rk(i))
      enddo
      do i = 1, napk
         ark(i) = 10.**(ark(i))
         apk(i) = 10.**(apk(i))
      enddo
  
      return
      end
#endif
