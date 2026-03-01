      function wlambda(red, w0, w1)
      real wlambda, red, w0, w1
      wlambda = w0 + w1*red/(1+red)
      return
      end

      function DEfunc(tau)
      real amax, omep, omeplam, w0, w1
      common /cosmoblock/ amax, omep, omeplam, w0,w1
      real DEfunc, tau, wlambda
      real red
      if(tau.le.0) then
      red = 1.e20
      else
      red = 1./tau-1
      endif
      DEfunc = 3*(1+wlambda(red, w0, w1))*(1+red)

      return
      end


      function gzDE(amax,anow, w0, w1)
      real gzDE, anow, res, DEfunc, tau,red
      external DEfunc,midsqu,midexp,midinf, midpnt,midsql3
      tau = anow/amax
      if(w0 .eq. -1 .and. w1 .eq.0) then
         gzDE = 1
         return
      else if (w1 .eq.0 .or. tau .le. 0.001) then
         gzDE = (amax/anow)**(3*(1+w0))
         return
      else if (anow .eq. amax) then
         gzDE = 1
         return
      else
          call qromo3(DEfunc,  tau, 1., res, midsql3)
          gzDE = exp(res)
          return
      endif
      end

      function HofEz(amax, anow, omep, omeplam, w0, w1)
      real HofEz, amax, anow, omep, omeplam, w0, w1
      real red, zp1
      red = amax/anow -1
      zp1 = red + 1
      HofEz = sqrt(zp1**3*omep+zp1*zp1*(1-omep-omeplam)+omeplam *
     &      gzDE(amax, anow, w0, w1) )
      return
      end
      SUBROUTINE qromo3(func,a,b,ss,choose)
      INTEGER JMAX,JMAXP,K,KM
      REAL a,b,func,ss,EPS
      EXTERNAL func,choose
      PARAMETER (EPS=1.e-5, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
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
      pause 'too many steps in qromo3'
      END
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
      SUBROUTINE midsql3(funk,aa,bb,s,n)
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
