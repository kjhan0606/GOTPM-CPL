c234567
      real app, apsq
      real amax, anow
      real wlam0, wlam1
      real omep, omeplam

      omep = 0.26
      omeplam = 0.74
      amax = 100.





      do i = 1, 100
      anow = i
      wlam0 = -1.
      wlam1 = 1.
      a1 = PoissonCor(anow, amax, omep, omeplam, wlam0, wlam1)

      wlam0 = -1.
      wlam1 = -1
      a2 = PoissonCor(anow, amax, omep, omeplam, wlam0, wlam1)
      print *, anow, a1, a2
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
	  print *, a, DplusAll(anow)
      a=growthall(anow, amax, omep, omeplam, wlam0, wlam1)
      print *,a, DplusAll(anow)






      stop
      end
