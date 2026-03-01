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

      amax = 1000
      anow = 1
      do i = 1, amax
         anow = i
         w0 = -0.5
         w1 =  0.5
         Cofw = stCw(amax, anow, omep, omeplam, w0,w1)
         growthall = cplDpAll(anow)
         aa = cplOmep(anow)
         w0 = -1
         w1 = 0
         Cofw = stCw(amax, anow, omep, omeplam, w0,w1)
         growthall0 = cplDpAll(anow)
         aa0 = cplOmep(anow)
         print *, i,amax/anow-1, aa, aa0
      enddo



      stop
      end
