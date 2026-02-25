c234567
      real amax, a, omep,omeplam, wlam
      external apprXgrowth,Xgrowthm, XgrowthQ

      amax = 100
      omep = 0.26
      omeplam = 0.74
      wlam = -0.5

      bini = growth(1., amax, omep)
      bnow = growth(amax, amax, omep)

      print *, 'cosmological constant:  ',bini,bnow,bini/bnow




c     a = growthXComp(amax, amax, omep, omeplam, wlam, apprXgrowthm)
c     b = growthXComp(1., amax, omep, omeplam, wlam, apprXgrowthm)
c     c = ApprGrowthXtotal(1./amax)
c     d = ApprGrowthXtotal(1.)
c     write(6,'(6f10.6)') b,a,b/a, c,d,c/d

      a = growthXComp(amax, amax, omep, omeplam, wlam, Xgrowthm)
      b = growthXComp(1., amax, omep, omeplam, wlam, Xgrowthm)
      c = GrowthXtotal(1./amax)
      d = GrowthXtotal(1.)
      write(6,'(6e16.6)') b,a,b/a, c,d,c/d

      a = growthXComp(amax, amax, omep, omeplam, wlam, XgrowthQ)
      b = growthXComp(1., amax, omep, omeplam, wlam, XgrowthQ)
      c = GrowthXtotal(1./amax)
      d = GrowthXtotal(1.)
      write(6,'(6e16.6)') b,a,b/a, c,d,c/d

      stop

      do i = 1, amax
        anow = i
        a = growth(anow,amax,omep)
        b = growth(amax,amax,omep)
        c = growthXComp(anow, amax, omep, omeplam, wlam, Xgrowthm)
        d = growthXComp(amax, amax, omep, omeplam, wlam, Xgrowthm)
        print *, amax/anow-1, a/b, c/d, c/a
      enddo

      stop
      end

      function omenowQ(a)
      common /wgrowth/ omep,omeplam, wlam
      real redp1
      if(a.eq.0) then
         omenowQ = 0
      else
         redp1 = 1./a
         omenowQ=redp1**(3*(1+wlam))*omeplam/
     &        (redp1**3*omep + redp1**(3*(1+wlam))*omeplam)
      endif
      return
      end

      function XgrowthQ(a)
      real GrowthXtotal
      external GrowthXtotal
      if(a.eq.0) then
         XgrowthQ = omenowQ(a)
      else
         XgrowthQ = (2.5-1.5*GrowthXtotal(a)/a)*omenowQ(a)
      endif
      return
      end

      function omenow(a)
      common /wgrowth/ omep,omeplam, wlam
      real redp1
      if(a.eq.0) then
         omenow = 1
      else
         redp1 = 1./a
         omenow=redp1**3*omep/
     &        (redp1**3*omep + redp1**(3*(1+wlam))*omeplam)
      endif
      return
      end

      function Xgrowthm(a)
      real GrowthXtotal
      external GrowthXtotal
      if(a.eq.0) then
         Xgrowthm = omenow(a)
      else
         Xgrowthm = (2.5-1.5*GrowthXtotal(a)/a)*omenow(a)
      endif
      return
      end

      function apprXgrowthm(a)
      real ApprGrowthXtotal
      external ApprGrowthXtotal
      if(a.eq.0) then
         apprXgrowthm = omenow(a)
      else
         apprXgrowthm = (2.5-1.5*ApprGrowthXtotal(a)/a)*omenow(a)
      endif
      return
      end


      function growthXComp(anow,amax,romep,romeplam, rwlam, growthit)
      real anow,growthit
      external growthit,midsql
      common /wgrowth/ omep,omeplam, wlam
      omep = romep
      omeplam = romeplam
      wlam = rwlam
      a = anow/amax
      call qsimp(growthit,0.,a,growthXComp)
c     call qromb(growthit,0.,a,growthXComp)
      return
      end
c---  Sefusatti & Vernizzi 2011 JCAP 047
      function ApprGrowthXtotal(a)
      real amax,a,redp1,omepz,omeplamz,ApprGrowthXtotal
      common /wgrowth/ omep,omeplam, wlam
      if(a.eq.0) then
      omepz = 1
      omeplamz = 0
      else
      redp1 = 1./a
      omepz = omep*(redp1)**3/
     &   (omep*(redp1)**3+omeplam*(redp1)**(3*(1+wlam)))
      omeplamz = omeplam*(redp1)**(3*(1+wlam))/
     &    (omep*(redp1)**3+omeplam*(redp1)**(3*(1+wlam)))
      endif

      ApprGrowthXtotal = a*2.5*omepz/(omepz**(4./7.) + 1.5*omepz
     &      + (1/70.-(1+wlam)/4.)*omeplamz*(1+omepz/2.))
      return
      end

c---  Sefusatti & Vernizzi 2011 JCAP 047
      function GrowthXtotal(a)
      complex a1,b1,c1,z,a2,b2,c2, hypgeo
      real a,amax,omep,omeplam,wlam,x,GrowthXtotal
      complex geo1, geo2
      common /wgrowth/ omep,omeplam, wlam


      x = (1-omep)/omep * (a)**(-3*wlam)


      a1 = cmplx(1.5,0)
      b1 = cmplx(-5./6./wlam,0)
      c1 = cmplx(1.-5./6./wlam,0)

      z = cmplx(-x,0)

      a2 = cmplx(1.5,0)
      b2 = cmplx(1.-5./6./wlam,0)
      c2 = cmplx(2.-5./6./wlam,0)

      geo1 = hypgeo(a1,b1,c1,z)
      geo2 = hypgeo(a2,b2,c2,z)


      GrowthXtotal = a*sqrt(1.+x)*(real(geo1)
     &     + x*5.*(1.+wlam)/(5.-6.*wlam)*real(geo2) )

      return
      end
c---  Growth factor in the universe with zero curvature
c---  omep= present matter density parameter
c---  omeplam = 1-omep is assumed (zero curvature inflation model)
      function growth(a,amax,omep)
      external coscon
      common/pass/x

      x0 = (1./omep-1)**(1./3.)
      redp1 = amax/a
      x = (1./omep-1)**(1./3.)/redp1
      call qsimp(coscon,0.,1.,gfac)
      growth = (x/x0)*sqrt(x**3+1)*gfac
      return
      end


      function coscon(y)
      common/pass/x
      coscon = 1./(sqrt(1+x**3*y**1.2))**3
      return
      end

C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      FUNCTION hypgeo(a,b,c,z)
      COMPLEX hypgeo,a,b,c,z
      REAL EPS
      PARAMETER (EPS=1.e-6)
CU    USES bsstep,hypdrv,hypser,odeint
      INTEGER kmax,nbad,nok
      EXTERNAL bsstep,hypdrv
      COMPLEX z0,dz,aa,bb,cc,y(2)
      COMMON /hypg/ aa,bb,cc,z0,dz
      COMMON /path/ kmax
      kmax=0
      if (real(z)**2+aimag(z)**2.le.0.25) then
        call hypser(a,b,c,z,hypgeo,y(2))
        return
      else if (real(z).lt.0.) then
        z0=cmplx(-0.5,0.)
      else if (real(z).le.1.0) then
        z0=cmplx(0.5,0.)
      else
        z0=cmplx(0.,sign(0.5,aimag(z)))
      endif
      aa=a
      bb=b
      cc=c
      dz=z-z0
      call hypser(aa,bb,cc,z0,y(1),y(2))
      call odeint(y,4,0.,1.,EPS,.1,.0001,nok,nbad,hypdrv,bsstep)
      hypgeo=y(1)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE hypser(a,b,c,z,series,deriv)
      INTEGER n
      COMPLEX a,b,c,z,series,deriv,aa,bb,cc,fac,temp
      deriv=cmplx(0.,0.)
      fac=cmplx(1.,0.)
      temp=fac
      aa=a
      bb=b
      cc=c
      do 11 n=1,1000
        fac=fac*aa*bb/cc
        deriv=deriv+fac
        fac=fac*z/n
        series=temp+fac
        if (series.eq.temp) return
        temp=series
        aa=aa+1.
        bb=bb+1.
        cc=cc+1.
11    continue
      pause 'convergence failure in hypser'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *rkqs)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.e-30)
      INTEGER i,kmax,kount,nstp
      REAL dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) pause
     *'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
      pause 'too many steps in odeint'
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER nv,NMAX,KMAXX,IMAX
      REAL eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv),SAFE1,SAFE2,
     *REDMAX,REDMIN,TINY,SCALMX
      PARAMETER (NMAX=50,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25,SAFE2=.7,
     *REDMAX=1.e-5,REDMIN=.7,TINY=1.e-30,SCALMX=.1)
CU    USES derivs,mmid,pzextr
      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      REAL eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest,xnew,
     *a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(NMAX),ysav(NMAX),
     *yseq(NMAX)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      EXTERNAL derivs
      DATA first/.true./,epsold/-1./
      DATA nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then
        hnext=-1.e29
        xnew=-1.e29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do 11 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,KMAXX
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+
     *1)))
12        continue
13      continue
        epsold=eps
        do 14 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif
      h=htry
      do 15 i=1,nv
        ysav(i)=y(i)
15    continue
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=x+h
        if(xnew.eq.x)pause 'step size underflow in bsstep'
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 16 i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1./(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1./err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.e35
      do 18 kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
18    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE hypdrv(s,y,dyds)
      REAL s
      COMPLEX y(2),dyds(2),aa,bb,cc,z0,dz,z
      COMMON /hypg/ aa,bb,cc,z0,dz
      z=z0+s*dz
      dyds(1)=y(2)*dz
      dyds(2)=(aa*bb*y(1)-(cc-(aa+bb+1.)*z)*y(2))*dz/(z*(1.-z))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
      INTEGER nstep,nvar,NMAX
      REAL htot,xs,dydx(nvar),y(nvar),yout(nvar)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
      INTEGER i,n
      REAL h,h2,swap,x,ym(NMAX),yn(NMAX)
      h=htot/nstep
      do 11 i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      call derivs(x,yn,yout)
      h2=2.*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        call derivs(x,yn,yout)
13    continue
      do 14 i=1,nvar
        yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
      INTEGER iest,nv,IMAX,NMAX
      REAL xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=50)
      INTEGER j,k1
      REAL delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),x(IMAX)
      SAVE qcol,x
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1./(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE qsimp(func,a,b,s)
      INTEGER JMAX
      REAL a,b,func,s,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20)
CU    USES trapzd
      INTEGER j
      REAL os,ost,st
      ost=-1.e30
      os= -1.e30
      do 11 j=1,JMAX
        call trapzd(func,a,b,st,j)
        s=(4.*st-ost)/3.
        if (abs(s-os).lt.EPS*abs(os)) return
        os=s
        ost=st
11    continue
      pause 'too many steps in qsimp'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE midinf(funk,aa,bb,s,n)
      INTEGER n
      REAL aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL a,b,ddel,del,sum,tnm,func,x
      func(x)=funk(1./x)/x**2
      b=1./aa
      a=1./bb
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE qromo(func,a,b,ss,choose)
      INTEGER JMAX,JMAXP,K,KM
      REAL a,b,func,ss,EPS
      EXTERNAL func,choose
      PARAMETER (EPS=1.e-6, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint
      INTEGER j
      REAL dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call choose(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.
11    continue
      pause 'too many steps in qromo'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE midsql(funk,aa,bb,s,n)
      INTEGER n
      REAL aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL ddel,del,sum,tnm,x,func,a,b
      func(x)=2.*x*funk(aa+x**2)
      b=sqrt(bb-aa)
      a=0.
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE midexp(funk,aa,bb,s,n)
      INTEGER n
      REAL aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL ddel,del,sum,tnm,x,func,a,b
      func(x)=funk(-log(x))/x
      b=exp(-aa)
      a=0.
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE midsqu(funk,aa,bb,s,n)
      INTEGER n
      REAL aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL ddel,del,sum,tnm,x,func,a,b
      func(x)=2.*x*funk(bb-x**2)
      b=sqrt(bb-aa)
      a=0.
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      pause 'too many steps in qromb'
      END
