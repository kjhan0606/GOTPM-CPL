c ifort -o testpcorr.exe testpcorr.f ../FLRW.o  ../NR.o  -g

c234567
      real app, apsq
      real app0, apsq0
      real amax, anow
      real w0, w1,Cw
      real omep, omeplam, growthXQ, growthXM
      external cplQGrowth, cplMgrowth,midsql,midsqu,midpnt, midinf
      real growthall,growthall0
      real aa0, aa
      external quintQGrowth, quintMgrowth

      omep = 0.26
      omeplam = 0.74

      amax = 100
      anow = 1
      do i = 1, amax
         anow = i
         w0 = -0.5
         w1 =  0.5
         Cofw = stCw(amax, anow, omep, omeplam, w0,w1)
         growthall = cplMGrowth(anow)
         call qromo(cplMGrowth, 0., anow, growthall,midsql)
         call qromo(cplMGrowth, 0., amax, aa,midsql)
         growthall = growthall/aa
         w0 = -1
         w1 = 0
         Cofw = stCw(amax, anow, omep, omeplam, w0,w1)
         growthall0 = cplMGrowth(anow)
         call qromo(cplMGrowth, 0., anow, growthall0,midsql)
         call qromo(cplMGrowth, 0., amax, aa0,midsql)
         growthall0 = growthall0/aa0
         print *, i,amax/anow-1, growthall, growthall0,
     &        growthall/growthall0
      enddo



      stop
      end
