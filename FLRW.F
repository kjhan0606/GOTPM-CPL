c----
c     Note: 
c       (1) For the whole part of this code, we are using the scale factor which
c            changes from 0 to amax
c       (2) For the quintessence and CPL model, one needs to initialize the cosmological models
c            using the common block of "cosmoblock".
c234567 
      function getOmei(omep, omeplam, wlam0, wlam1, zi)
      real omep, omeplam, wlam0, wlam1, zi
      real amax, anow
      anow =1
      
      if(wlam1.eq.0) then
      getOmei = quintOmep(anow)
      else 
      getOmei = cplOmegaMatter(anow)
      endif
      return
      end

      subroutine appapsq(romep, romeplam, rwlam0, rwlam1, ramax,anow, 
     &      app, apsq)
      real omep, omeplam, wlam0, wlam1, amax, anow, app, apsq
      real CsDE
      common /cosmoblock/ amax, omep, omeplam, wlam0,wlam1,CsDE
      real ainv,omei,zi,omemnow,omeLnow
      real PI,ww
      external cplOmegaMatter, getOmei, scaleDE
      real expfDE,red, fi

      omep = romep
      omeplam = romeplam
      wlam0 = rwlam0
      wlam1 = rwlam1
      amax = ramax

      PI = 3.141592653
      ainv = 1./anow
      zi = amax - 1
      red = amax/anow-1
      ww = wlam0 + wlam1*red/(1+red)

      if(wlam1.eq.0) then
         omei = getOmei(omep, omeplam, wlam0, wlam1, zi)
         app = -4.*PI/3.*(ainv*ainv+(1.+3*wlam0)*omeplam/omep/
     &      anow**(2.+3*wlam0)*amax**(3*wlam0));
         apsq = 8.*PI/3.*(ainv+1./omei-1.+omeplam/omep*
     &      (anow**(-1.-3*wlam0)-1.)*amax**(3*wlam0))
c         apsq = 8.*PI/3.*(ainv+omeplam/omep*
c     &      (anow**(-1.-3*wlam0))*amax**(3*wlam0))
      else 
         omei = getOmei(omep, omeplam, wlam0, wlam1, zi)
         expfDE = scaleDE(amax,anow, wlam0, wlam1)
         fi = scaleDE(amax, 1.,wlam0, wlam1)
         app = -4.*PI/3.*(ainv*ainv+omeplam/omep*anow/amax**3
     &       *(1+3*ww)*expfDE)
         apsq = 8.*PI/3.*(ainv+1./omei-1.+(omeplam/omep)*
     &      ((expfDE*anow**2 - fi)/amax**3))
c         apsq = 8.*PI/3.*(ainv+(omeplam/omep)*
c     &      (anow**2/amax**3*expfDE))         
      endif
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


c--- growth factor for non-flat and non-zero cosmological constant
      function growth2(a,amax,omep,omeplam)
      external coscon2

      omepk = 1-omep-omeplam
      aa = a/amax
      growth2 = sqrt(omeplam*aa**3+omepk*aa+omep)/aa**1.5

      call qsimp(coscon2,0.,aa,gfac)
      growth2 = growth2*gfac

      return
      end

      function coscon2(y)
      real com_omep,com_omepb,com_omeplam,com_h,com_size,omepb
      common/cospar/com_omep,com_omepb,com_omeplam,com_h,com_size
      com_omepk = 1-com_omep-com_omeplam
      coscon2 = y**1.5/(com_omeplam*y**3+com_omepk*y+com_omep)**1.5
      return
      end

c---  Omega Matter for Quintessence model
      function quintOmep(anow)
      real amax, omep, omeplam, wlam0, wlam1
      real CsDE
      common /cosmoblock/ amax,omep,omeplam, wlam0, wlam1,CsDE
      real redp1
      if(anow.eq.0) then
         if(wlam0 .eq.0) then
            quintOmep = omep/(omep+omeplam)
         else if(wlam0 .gt.0) then
            quintOmep = 0
         else
            quintOmep = 1
         endif
      else
        redp1 = amax/anow
        quintOmep=redp1**3*omep/
     &        (redp1**3*omep + redp1**(3*(1+wlam0))*omeplam+
     &        redp1**2*(1-omep-omeplam))
      endif
      return
      end
      
c---  Omega DE for Quintessence model
      function quintOmeplam(anow)
      real amax, omep, omeplam, wlam0, wlam1
      real CsDE
      common /cosmoblock/ amax,omep,omeplam, wlam0, wlam1,CsDE
      real redp1
      if(anow.eq.0) then
c---  This is only for w<0
         if(wlam0 .eq. 0) then
            quintOmeplam = omeplam/(omep+omeplam)
         else if(wlam0 .gt.0) then
            quintOmeplam = 1
         else 
            quintOmeplam = 0
         endif
      else
      redp1 = amax/anow
      quintOmeplam=redp1**(3*(1+wlam0))*omeplam /
     &        (redp1**3*omep + redp1**(3*(1+wlam0))*omeplam+
     &        redp1**2*(1-omep-omeplam))
      endif
      return
      end


c--- CPL growth factor for matter
      function cplIntMG(anow)
      real amax, omep, omeplam, w0, w1, cplIntMG
      real CofW
      real CsDE
      common /cosmoblock/ amax,omep,omeplam, w0, w1,CsDE
      real cplClusterDpAll, cplOmegaMatter
      real cplSmoothDpAll,cplDotSmoothDpAll
      if(CsDE.eq.0) then
         if(anow.eq.0) then
            cplIntMG = cplOmegaMatter(anow)/amax
         else
            cplIntMG =(2.5-1.5*cplClusterDpAll(anow)*amax/anow)
     &        *cplOmegaMatter(anow)/amax
         endif
      else
         CofW = Cw(amax,anow,omep, omeplam,w0,w1)
         if(anow.eq.0) then
            cplIntMG = 1./CofW
         else
            cplIntMG = cplDotSmoothDpAll(anow)/anow/CofW
         endif
         print *, anow, cplDotSmoothDpAll(anow), CofW,cplIntMG
      endif
      return
      end

c--- CPL growth factor for DE
      function cplIntQG(anow)
      real amax, omep, omeplam, wlam0, wlam1
      real CsDE
      common /cosmoblock/ amax,omep,omeplam, wlam0, wlam1,CsDE
      external cplClusterDpAll, cplOmegaDE
      if(anow.eq.0) then 
         cplIntQG = cplOmegaDE(anow)/amax
      else
         cplIntQG = (2.5-1.5*cplClusterDpAll(anow)*amax/anow)*
     &           cplOmegaDE(anow)/amax
      endif
      return
      end


c---  Integrand of the Quintessence growth factor for matter
      function quintMGrowth(anow)
      real amax, omep, omeplam, wlam0, wlam1
      real CsDE
      common /cosmoblock/ amax,omep,omeplam, wlam0, wlam1,CsDE
      external quintGrowth,quintOmep
      if(anow.eq.0) then
         quintMGrowth = quintOmep(anow)/amax
      else
         quintMGrowth = (2.5-1.5*quintGrowth(anow)*amax/anow)*
     &         quintOmep(anow)/amax
      endif
      return
      end

c---  Integrand of the Quintessence growth factor for DE
      function quintQGrowth(anow)
      real amax, omep, omeplam, wlam0, wlam1
      real CsDE
      common /cosmoblock/ amax,omep,omeplam, wlam0, wlam1,CsDE
      external quintGrowth
      if(anow.eq.0) then 
         quintQGrowth = quintOmeplam(anow)/amax
      else
         quintQGrowth = (2.5-1.5*quintGrowth(anow)*amax/anow)*
     &                  quintOmeplam(anow)/amax
      endif
      return
      end

c--- Poisson correction factor for Quintessence & CPL models.
      function PoissonCor(anow,ramax,romep, romeplam, rwlam0,
     &    rwlam1,iCsDE)
      real amax, omep, omeplam, wlam0, wlam1, w0, w1
      real CsDE,iCsDE
c---  initialization block for the cosmological model
      common /cosmoblock/ amax,omep,omeplam, wlam0, wlam1,CsDE, w0, w1
      external quintQGrowth, quintMGrowth
      external cplIntQG, cplIntMG, midsql, midsqu
      real red, cplOmegaMatter, cplOmegaDE
      omep = romep
      omeplam = romeplam
      wlam0 = rwlam0
      wlam1 = rwlam1
      w0 = rwlam0
      w1 = rwlam1   
      CsDE = iCsDE
      amax = ramax
      red = amax/anow-1


      if(CsDE .eq.0) then
      if(rwlam1.eq.1) then
      call qromo(quintQGrowth,0., anow, growthXQ, midsql)
      call qromo(quintMGrowth,0., anow, growthXM, midsql)
      PoissonCor = (1 + (1+wlam0)*growthXQ/growthXM)
      else
c      call qromo(cplIntQG,0., anow, growthXQ, midsql)
c      call qromo(cplIntMG,0., anow, growthXM, midsql)
      Cofw=stCw(amax,anow,omep,omeplam,w0,w1,CsDE)
      growthXQ = DEgrowthcpl(anow)
      call qsimp(cplIntMG,0.,anow,growthXM)
      PoissonCor = (1 + cplOmegaDE(anow) * growthXQ /
     & cplOmegaMatter(anow) / growthXM)
      endif
      else if(CsDE .eq.1) then
      PoissonCor = 1
      endif

      return
      end

c---  w for the Quintessence & CPL
      function wlambda(red, w0, w1)
      real wlambda, red, w0, w1
      wlambda = w0 + w1*red/(1+red)
      return
      end

c     y = ln(1+red)
      function exponentDE(y)
      real amax, omep, omeplam, w0, w1,y
      integer CsDE
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE
      real exponentDE, wlambda
      real red
      external wlambda
      red = exp(y)-1
      exponentDE = 3*(1+wlambda(red, w0, w1))

      return
      end


      function scaleDE(amax,anow, w0, w1)
      real scaleDE, anow, res, exponentDE, red,y
      external exponentDE,midsqu,midexp,midinf, midpnt,midsql3
      if(anow.le.0) then
        y = 10
      else
        red = amax/anow - 1
        y = alog(1.+red)
      endif
      if(w0 .eq. -1 .and. w1 .eq.0) then
         scaleDE = 1
         return 
      else if (w1 .eq.0 .or. y .ge. 10) then
         scaleDE = (amax/anow)**(3*(1+w0))
         return
      else if (anow .eq. amax) then
         scaleDE = 1
         return
      else
          scaleDE=(anow/amax)**(-3*(1+w0+w1))*exp(-3*w1*(1-anow/amax))
c          call qromo(exponentDE,  0., y, res, midsql3)
c          scaleDE = exp(res)
          return
      endif
      end

      function HofEz(amax, anow, omep, omeplam, w0, w1)
      real HofEz, amax, anow, omep, omeplam, w0, w1
      real red, zp1
      external scaleDE
      red = amax/anow -1
      zp1 = red + 1
      HofEz = sqrt(zp1**3*omep+zp1*zp1*(1-omep-omeplam)+omeplam * 
     &      scaleDE(amax, anow, w0, w1) )
      return
      end


c---  Omega matter for various DM models including CPL
      function cplOmegaMatter(anow)
      real amax, omep, omeplam, w0,w1
      integer CsDE
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE
      real cplOmegaMatter, anow
      real zp1, hoverz
      real red,ww
      external HofEz
      red = amax/anow-1
      if(anow .le.0) then
         ww = w0 + w1
         if(ww .eq. 0) then
            cplOmegaMatter = omep/(omep+omeplam)
         else if(ww .gt.0) then
            cplOmegaMatter = 0
         else 
            cplOmegaMatter = 1
         endif
      else 
         zp1 = red + 1
         hoverez = HofEz(amax, anow, omep, omeplam, w0, w1)
         cplOmegaMatter = zp1**3*omep/hoverez**2
      endif
      return
      end

c---  Omega DE for various DM models including CPL
      function cplOmegaDE(anow)
      real amax, omep, omeplam, w0,w1
      integer CsDE
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE
      real cplOmegaDE,  anow
      real red
      real hoverez
      external HofEz, scaleDE
      if(anow.le.0) then
      cplOmegaDE = 0
      return
      endif
      hoverez = HofEz(amax, anow, omep, omeplam, w0, w1)
      cplOmegaDE = scaleDE(amax, anow, w0,w1)*omeplam/hoverez**2
      return
      end

      function stCw(ramax, anow, romep, romeplam, rw0,rw1, iCsDE)
      real Cw, amax, anow, omep, omeplam, w0,w1, red
      real CsDE,iCsDE
c---  initialization block for the cosmological model
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE
      external wlambda, cplOmegaDE, cplOmegaMatter
      amax = ramax
      omep = romep
      omeplam = romeplam
      CsDe = iCsDE
      w0 = rw0
      w1 = rw1

      red = amax/anow-1
      stCw = 1. + (1+wlambda(red, w0,w1))*
     &    cplOmegaDE(anow)/cplOmegaMatter(anow)
      return
      end

      function Cw(amax, anow, omep, omeplam, w0,w1)
      real Cw, amax, anow, omep, omeplam, w0,w1, red
      external wlambda, cplOmegaDE, cplOmegaMatter
      if(anow.eq.0) then
      ww = w0+w1
      else
      red = amax/anow-1
      ww = wlambda(red,w0,w1)
      endif
      Cw = 1. + (1+ww)*cplOmegaDE(anow)/cplOmegaMatter(anow)
      return
      end
c--   Sefusatti et al. 2011
      function growthQSmooth(anow,ramax, romep, romeplam, rw0,rw1)
      real amax, anow, omep, omeplam, w0,w1,Cw
      integer CsDE
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE
      real res,amaxlog,anowlog
      external fPsmooth,midsqu
      amax = ramax
      omep = romep
      romeplam = omeplam
      w0 = rw0
      w1 = rw1

      anowlog = alog(anow)
      amaxlog = alog(amax)
c--- ***
      call qromb(fPsmooth,-5., anowlog,  res)
      call qromb(fPsmooth,-5., amaxlog,  g0)
      growthQSmooth = exp(-res)/exp(-g0)
c--- ***

      return
      end

c--   Sefusatti et al. 2011
      function fpSmooth(anowlog)
      real amax, anow, omep, omeplam, w0,w1,Cw
      real smoothGamma,red,Hsqr
      real omegam, omegalam,wlambda,expfDE
      integer CsDE
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE

      anow = exp(anowlog)

      red = amax/anow-1
      ww = wlambda(red, w0,w1)

      omegalam = cplOmegaDE(anow)
      smoothGamma = 0.545 +0.007*omegalam+(1+ww)*(0.025+0.005*omegalam)

      fpSmooth = cplOmegaMatter(anow)**smoothGamma

      return
      end

      function CwOveraH3(anow)
      real CwOveraH3, amax, anow, omep, omeplam, w0,w1,Cw
      integer CsDE
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE
      external Cw,HofEz
      real red
      if(anow.le.0) then
         CwOveraH3 = 0
      else 
         if(CsDE.eq.0) then
         CwOveraH3 = Cw(amax,anow,omep, omeplam,w0,w1)/
     &      (anow*HofEz(amax,anow,omep, omeplam,w0,w1))**3
         else
         CwOveraH3 = 1./(anow*HofEz(amax,anow,omep, omeplam,w0,w1))**3
         endif
      endif
      return
      end


      function cplClusterDpAll(anow)
      real cplClusterDpAll, amax, anow, omep, omeplam, w0,w1
      integer CsDE
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE
      real red, res
      external midsqu, CwOveraH3,midpnt,midsql
      if(anow .le.0) then
         cplClusterDpAll = 0
      else
         if(anow/amax .le. 1.E-3) then
                 cplClusterDpAll = anow/amax
         else
         call qromb(CwOveraH3, 0., anow, res)
         cplClusterDpAll = 5./2.*omep*
     &  HofEz(amax, anow, omep, omeplam,w0,w1)*res*amax*amax
         endif
      endif
      return
      end

      function Dmslope(anow)
      real Dmslope, amax, anow, omep, omeplam, w0,w1
      integer CsDE
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE
      real a1,a2,a3, red
      external cplClusterDpAll, cplOmegaMatter
      if(anow .gt. 0) then
      red = amax/anow-1
      Dmslope = (2.5-1.5*cplClusterDpAll(anow)/(anow/amax))
     &      *cplOmegaMatter(anow)/amax
      else 
      Dmslope = 1
      endif
      return
      end
c-----
      function cplSmoothDpAll(anow)
      INTEGER KMAXX,NMAX,NVAR
      PARAMETER (KMAXX=200,NMAX=50,NVAR=2)
      INTEGER i,kmax,kount,nbad,nok,nrhs
      REAL dxsav,eps,h1,hmin,x1,x2,x,y,ystart(NVAR)
      COMMON /path/ kmax,kount,dxsav,x(KMAXX),y(NMAX,KMAXX)
      COMMON nrhs
      real cplSmoothDpAll, amax, anow, omep, omeplam, w0,w1
      integer CsDE
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE
      external derivs, rkqs
      real red, res

      if(anow .le.0) then
         cplSmoothDpAll = 0
      else
         nrhs = 0
         x1 = alog(0.01/amax)
         x2 = alog(anow/amax)
         eps = 1.e-5
         h1 = .1
         hmin = 0
         hmax = 100
         dxsav = (x2-x1)/100
         ystart(1) = 0.01/amax
         ystart(2) = 0.01/amax
         call odeint(ystart, 2, x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)
         cplSmoothDpAll = ystart(1)
      endif
      return
      end
c-----
      function cplDotSmoothDpAll(anow)
      INTEGER KMAXX,NMAX,NVAR
      PARAMETER (KMAXX=200,NMAX=50,NVAR=2)
      INTEGER i,kmax,kount,nbad,nok,nrhs
      REAL dxsav,eps,h1,hmin,x1,x2,x,y,ystart(NVAR)
      COMMON /path/ kmax,kount,dxsav,x(KMAXX),y(NMAX,KMAXX)
      COMMON nrhs
      real cplSmoothDpAlldot, amax, anow, omep, omeplam, w0,w1
      integer CsDE
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE
      external derivs, rkqs
      real red, res

      if(anow .le.0) then
         cplDotSmoothDpAll = 0
      else
         nrhs = 0
         x1 = alog(0.01/amax)
         x2 = alog(anow/amax)
         eps = 1.e-5
         h1 = .1
         hmin = 0
         hmax = 100
         dxsav = (x2-x1)/100
         ystart(1) = 0.01/amax
         ystart(2) = 0.01/amax
         call odeint(ystart, 2, x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)
         cplDotSmoothDpAll = ystart(2)
      endif
      return
      end



c---  for the smooth case of CPL
      subroutine derivs(x,y,dydx)
      integer nrhs
      real x, y(*), dydx(*)
      real Omegam, CofW, dlnHdlna
      integer CsDE
      real amax, anow, omep, omeplam, w0,w1
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE
      COMMON nrhs
      nrhs = nrhs + 1
      anow = exp(x)*amax
      Omegam = cplOmegaMatter(anow)
      CofW = Cw(amax,anow,omep, omeplam,w0,w1)
      dlnHdlna = -1.5*Omegam*CofW
      dydx(1) = y(2)
      dydx(2) = -(2+dlnHdlna)*y(2) - dlnHdlna*y(1)/CofW
      return
      end

c---  for cplDEgrowth, dD+/ dlna  2024.2
      subroutine cplDEderivs(x,y,dydx)
      integer nrhs, CsDE
      real x, y(*), dydx(*)
      real hoverez, Omegam, CofW, dDpClusterAlloverda
      real dDpClusterAlloverdla, res, aa
      real amax, anow, omep, omeplam, w0, w1, wlam
      common /cosmoblock/ amax, omep, omeplam, w0, w1, CSDE
      COMMON nrhs
      external CwOveraH3
      nrhs = nrhs + 1
      anow = exp(x) * amax 
      hoverez = HofEz(amax, anow, omep, omeplam, w0, w1)
      Omegam = cplOmegaMatter(anow)
      CofW = Cw(amax,anow,omep,omeplam,w0,w1)
      call qromb(CwOveraH3, 0., anow, res)
      aa = -1.5 * hoverez * Omegam * CofW / anow
      dDpClusterAlloverda = 2.5 * omep * aa*res 
     & + 2.5 * omep * CofW / (anow**3) / (hoverez**2)
      dDpClusterAlloverdla = anow * amax * amax * dDpClusterAlloverda
      wlamb = w0 + w1 * (1-anow/amax)
      dydx(1) = 3 * wlamb * y(1) + (1+wlamb)/CofW * dDpClusterAlloverdla 
      return
      end  

c--- cpl DE growth ode solver 2024.2
      function DEgrowthcpl(anow)
      INTEGER KMAXX,NMAX,NVAR
      PARAMETER (KMAXX=200,NMAX=50,NVAR=1)
      INTEGER i,kmax,kount,nbad,nok,nrhs
      REAL dxsav,eps,h1,hmin,x1,x2,x,y,ystart(NVAR)
      COMMON /path/ kmax,kount,dxsav,x(KMAXX),y(NMAX,KMAXX)
      COMMON nrhs
      real DEgrowthcpl, amax, anow, omep, omeplam, w0,w1
      integer CsDE
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE
      external cplDEderivs, rkqs
      real red, res

      if(anow .le.0) then
         DEgrowthcpl = 0
      else
         nrhs = 0
         x1 = alog(0.01/amax)
         x2 = alog(anow/amax)
         eps = 1.e-5
         h1 = .1
         hmin = 0
         hmax = 100
         dxsav = (x2-x1)/100
         ystart(1) = 0.01/amax * (1+(w0+w1)) / (1-3*(w0+w1))
         call odeint(ystart, 1, x1,x2,eps,h1,hmin,nok,nbad
     &,cplDEderivs,rkqs)
         DEgrowthcpl = ystart(1)
      endif
      return
      end

      
c---  matter growth factor in CPL model + smooth Quintessence
      function growthcpl(anow, ramax, romep, romeplam, rw0,rw1, iCsDE)
      real CsDE,iCsDE
      real growthcpl, ramax, romep, romeplam, rw0,rw1, Cofw
      real amax, anow, omep, omeplam, w0,w1, cplIntMG
c---  initialization block for the cosmological model
      common /cosmoblock/ amax, omep, omeplam, w0,w1,CsDE
      external midsqu, Dmslope, midsql,cplIntMG,midpnt
      amax = ramax
      omep = romep
      omeplam = romeplam
      w0 = rw0
      w1 = rw1
      CsDE = iCsDE

      if(CsDE.eq.0.) then
         call qsimp(cplIntMG, 0., anow, growthcpl)
c         call qromb(cplIntMG, 0., anow, growthcpl)
      else if (CsDE.eq.1.) then
         Cofw=stCw(amax,anow,omep,omeplam,w0,w1,CsDE)
         growthcpl=cplSmoothDpall(anow)
c        call qromo(cplIntMG, 0., anow, growthcpl, midsql)
c         call qsimp(cplIntMG, 0., anow, growthcpl)
      else
         print *, 'Choose an integer between 0 and 1'
         stop
      endif
      return
      end




c---  growth factor of matter for clustering quintessence
      function growthgen(anow,ramax,romep,romeplam, rwlam0,iCsDE)
      real anow
      external quintMGrowth,midsql
      real amax, omep, omeplam, wlam0, wlam1
      real CsDE,iCsDE
c---  initialization block for the cosmological model
      common /cosmoblock/ amax,omep,omeplam, wlam0, wlam1,CsDE
      amax = ramax
      omep = romep
      omeplam = romeplam
      wlam0 = rwlam0
      CsDE = iCsDE
c---- Matter growth factor: Sefusatti & Vernizzi 2011 JCAP 047
      call qromo(quintMGrowth,0.,anow,growthgen, midsql)
      growthgen = growthgen/amax
      return
      end

c---  growth factor of DM for quintessence
c---- Total growth factor: Sefusatti & Vernizzi 2011 JCAP 047
      function quintGrowth(anow)
      complex a1,b1,c1,z,a2,b2,c2,hypgeo
      complex geo1,geo2
      real amax, omep, omeplam, wlam0, wlam1
      real CsDE
      common /cosmoblock/ amax,omep,omeplam, wlam0, wlam1,CsDE
      a = anow/amax

      x = (1-omep)/omep * (a)**(-3*wlam0)

      a1 = cmplx(1.5,0)
      b1 = cmplx(-5./6./wlam0,0)
      c1 = cmplx(1.-5./6./wlam0,0)
      z = cmplx(-x,0)
      a2 = cmplx(1.5,0)
      b2 = cmplx(1.-5./6./wlam0,0)
      c2 = cmplx(2.-5./6./wlam0,0)

      geo1 = hypgeo(a1,b1,c1,z)
      geo2 = hypgeo(a2,b2,c2,z)

      quintGrowth = a*sqrt(1+x)*(real(geo1)
     &     + x*5.*(1+wlam0)/(5.-6.*wlam0)*real(geo2) )

      return
      end
