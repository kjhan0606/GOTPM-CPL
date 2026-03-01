c ifort -o testD.exe testD.f ../FLRW.o  ../NR.o  -g

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

      amax = 10000
      anow = 1
      do i = 0, amax
         anow = i
         w0 = -0.5
         w1 =  0.5
         Cofw = stCw(amax, anow, omep, omeplam, w0,w1)
         aa = cplMGrowth(anow)
         bb = cplDpAll(anow)
         cc = cplOmep(anow)
         w0 = -0.1
         w1 = 0
         Cofw0 = stCw(amax, anow, omep, omeplam, w0,w1)
         aa0 = cplMGrowth(anow)
         bb0 = cplDpAll(anow)
         cc0 = cplOmep(anow)
         print *, amax/anow-1, aa,aa0,cc,cc0
      enddo



      stop
      end
