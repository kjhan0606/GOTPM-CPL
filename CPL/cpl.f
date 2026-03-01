c234567
      real app, apsq
      real amax, anow
      real wlam0, wlam1
      real omep, omeplam
      character*88 arg

      omep = 0.26
      omeplam = 0.74
      amax = 100.



      call getarg(1, arg)
      read(arg, *) wlam0
      call getarg(2, arg)
      read(arg, *) wlam1


      print *, 'w0= ', wlam0, 'w1= ',wlam1
      do i = 1, 100
      anow = i
c     wlam0 = -1.
c     wlam1 = 1.5
      bini = growthall(anow, amax, omep, omeplam, wlam0, wlam1)
      bnow = growthall(amax, amax, omep, omeplam, wlam0, wlam1)
      a1 =  bini/bnow

      w0 = -1.
      w1 = 0.0
      bini = growthall(anow, amax, omep, omeplam, w0, w1)
      bnow = growthall(amax, amax, omep, omeplam, w0, w1)
      a2 =  bini/bnow
      ratio = a1/a2
      red = amax/anow - 1
      print *, anow, red,a1, a2, ratio
      enddo
      stop


c     do i = 1, 100
c     anow = i
c     call appapsq(omep, omeplam, wlam0, wlam1, amax,anow, 
c    &      app, apsq)
c     print *, anow,app,apsq
c     enddo

c     stop



      vamp1 = (growthall(1.+0.005,amax,omep,omeplam,wlam0, wlam1)-
     &          growthall(1.-0.005,amax,omep, omeplam,wlam0, wlam1))
     &    /(0.01)*1./bini
      print *, 'wlambda0,1 case: ',omep,' bini=', bini,
     &        ' bnow=',bnow, vamp1

      stop


      bini = growthall(1., amax, omep, omeplam, wlam0, wlam1)
      print *, 'bini= ', bini
      bnow = growthall(amax, amax, omep, omeplam, wlam0, wlam1)
      print *,'bnow',  bnow




      vamp1 = growthall( 1-0.005,amax,omep,omeplam,wlam0, wlam1)
      vamp2 = growthall(1+0.005,amax,omep,omeplam,wlam0, wlam1)
      print *, vamp1, vamp2

      anow = 10.
      a=growthgen(anow, amax, omep, omeplam, wlam0)
      print *, a, cplDpAll(anow)
      a=growthall(anow, amax, omep, omeplam, wlam0, wlam1)
      print *,a, cplDpAll(anow)






      stop
      end
