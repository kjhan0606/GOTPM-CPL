c234567
      real app, apsq
      real app0, apsq0
	  real amax, anow
	  real wlam0, wlam1
	  real omep, omeplam

      omep = 0.3
	  omeplam = 0.7

	  amax = 100
	  anow = 1
	  do i = 1, amax
	     anow = i
	     wlam0 = -1.1
		 wlam1 = 0.1
	     call appapsq(omep, omeplam, wlam0, wlam1, amax, anow,
     &        app, apsq)
	     aaa= PoissonCor(anow, amax, omep, omeplam, wlam0, wlam1)
	     wlam0 = -1
	     wlam1 = 0
	     call appapsq(omep, omeplam, wlam0, wlam1, amax, anow,
     &        app0, apsq0)
	     aaa0= PoissonCor(anow, amax, omep, omeplam, wlam0, wlam1)
	     print *, anow, app/app0, apsq/apsq0,aaa
	  enddo



	  stop
	  end
