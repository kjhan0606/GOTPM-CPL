c ifort -o testCw.exe testCw.f ../FLRW.o  ../NR.o  -g

c234567
      real app, apsq
      real app0, apsq0
      real amax, anow
      real w0, w1,Cw
      real omep, omeplam, growthXQ, growthXM
      external cplQGrowth, cplMgrowth,midsql,midsqu
      external quintQGrowth, quintMgrowth

      omep = 0.279
      omeplam = 0.721

      amax = 100
      anow = 1
      do i = 1, amax
         anow = i
         w0 = -0.5
         w1 = 0.5
         Cofw = stCw(amax, anow, omep, omeplam, w0,w1)
         aa = cplDpAll(anow)
         call qromo(cplQGrowth,0., anow, growthXQ, midsql)
         call qromo(cplMGrowth,0., anow, growthXM, midsql)
         call appapsq(omep,omeplam, w0, w1, amax, anow, app, apsq)
         w0 = -1
         w1 = 0
         Cofw0 = stCw(amax, anow, omep, omeplam, w0,w1)
         aa0 = cplDpAll(anow)
         call qromo(quintQGrowth,0., anow, growthXQ0,midsql)
         call qromo(quintMGrowth,0., anow, growthXM0,midsql)
c        print *, amax/anow-1, Cofw, aa/aa0, growthXM/growthXM0,
c    &     growthXQ/growthXQ0
         call appapsq(omep,omeplam, w0, w1, amax, anow, app0, apsq0)
         print *, amax/anow-1, app/app0, app, app0, apsq/apsq0
      enddo



      stop
      end
