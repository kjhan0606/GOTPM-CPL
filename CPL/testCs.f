c ifort -o testCs.exe testCs.f ../FLRW.o  ../NR.o  -g

c234567
      real app, apsq,growthcpl
      real app0, apsq0
      real amax, anow
      real w0, w1,Cw
      real omep, omeplam, growthXQ, growthXM
      real CsDE
      external cplQGrowth, cplMgrowth,midsql3 ,midsql,midsqu,midinf
      external cplIntMG, cplIntQG
      external quintQGrowth, quintMgrowth,growthcpl
      external PoissonCor,ScaleDE, exponentDE, appapsq
      external interp1d
      external getOmei
      integer Fig2_1,Fig2_2,Fig2_3,Fig2_4,testw1
c---
c      integer, parameter :: max_lines = 1000
c      integer :: i, num_lines
c      character(100) :: filename
c      real, dimension(max_lines) :: x
c      integer ll
c      real xarr(2), yarr(2)
c      real ddx
c      integer n
c--- 
      omep = 0.26
      omeplam = 0.74

      amax = 48
c      anow = 1
      Fig2_1 = 0
      Fig2_2 = 0
      Fig2_3 = 0
      Fig2_4 = 0
      testw1 = 1
      if(Fig2_1) then
      do i = 1, amax
         anow = i
         w0 = -1.1
         w1 =  0.0
         Cofw = stCw(amax, anow, omep, omeplam, -1.2,w1,0)
         aa = cplClusterDpall(anow)
         Cofw = stCw(amax, anow, omep, omeplam, -1.1,w1,0)
         bb = cplClusterDpall(anow)
         Cofw = stCw(amax, anow, omep, omeplam, -0.9,w1,0)
         cc = cplClusterDpall(anow)
         Cofw = stCw(amax, anow, omep, omeplam, -0.8,w1,0)
         d = cplClusterDpall(anow)

         Cofw = stCw(amax, anow, omep, omeplam, -1.2,w1,1)
         aa0 = cplSmoothDpAll(anow)
         Cofw = stCw(amax, anow, omep, omeplam, -1.1,w1,1)
         bb0 = cplSmoothDpAll(anow)
         Cofw = stCw(amax, anow, omep, omeplam, -0.9,w1,1)
         cc0 = cplSmoothDpAll(anow)
         Cofw = stCw(amax, anow, omep, omeplam, -0.8,w1,1)
         d0 = cplSmoothDpAll(anow)

         write(*, '(9(F12.7,1X))'),amax/anow-1,aa,aa0,bb,bb0,cc,cc0,d,d0
      enddo
      else if(Fig2_3) then
      do i = 1, amax
         anow = i
         w0 = -0.9
         w1 =  0.0
         Cofw = stCw(amax, anow, omep, omeplam, w0,w1,0)
         aa = cplClusterDpall(anow)
         Cofw = stCw(amax, anow, omep, omeplam, -1.1,w1,0)
         bb = cplClusterDpall(anow)
         Cofw = stCw(amax, anow, omep, omeplam, w0,w1,1)
         cc = cplSmoothDpall(anow)
         Cofw = stCw(amax, anow, omep, omeplam, -1.1,w1,1)
         dd = cplSmoothDpall(anow)
         
c        Cofw = stCw(amax, anow, omep, omeplam, w0,w1,1)
c        aa = cplSmoothDpall(anow)

c--      LCDM case
         w0 = -1
         w1 =  0.0
         Cofw = stCw(amax, anow, omep, omeplam, w0,w1,1)
         aa0 = cplSmoothDpAll(anow)
         write(*,'(6(F12.7,1X))'), amax/anow-1,aa,bb,cc,dd,aa0
      enddo

      else if(Fig2_2) then
      w0 = -0.8
      w1 =  0
      do i = 1, amax
          anow = i
          aa = growthcpl(anow,amax,omep,omeplam,w0,w1,0)
          aa0 = growthcpl(anow,amax,omep,omeplam,w0,w1,1)
          print *, amax/anow-1,aa/aa0,aa,aa0
c         aa = growthcpl(anow, amax, omep, omeplam, w0,w1, 0)
c         aa0 = growthcpl(anow, amax, omep, omeplam, w0,w1, 1)
c         print *, amax/anow-1, aa/aa0,aa,aa0
c         if(i.eq.1) stop
      enddo

      else if (Fig2_4) then
      w0 = -1.0
      w1 = 0.5
      do i = 1, amax
         anow = i
         Cofw = stCw(amax, anow, omep, omeplam, w0,w1,0)
c         aa   = cplSmoothDpAll(anow)
c         aa0 = growthcpl(amax,amax,omep,omeplam,w0,w1,0)
         aa = growthcpl(anow,amax,omep,omeplam,w0,w1,1)
c         aa= PoissonCor(anow,amax,omep,omeplam,w0,w1,0)
c          aa = cplIntMG(anow)
c         Cofw = stCw(amax, anow, omep, omeplam, w0,w1,0)
c         aa = cplClusterDpAll(anow)
         print *, amax/anow-1, aa
      enddo
      else if (testw1) then
      w0 = -0.8
      w1 =  0
      CsDE = 0
      do i=1, amax
        anow = i
        call appapsq(omep,omeplam,w0,w1,amax,anow,app,apsq)
        print *, anow, app, apsq
      enddo

c      aa = PoissonCor(anow,amax,omep,omeplam,w0,w1,0)
c      aa = growthcpl(anow,amax,omep,omeplam,w0,w1,1)
c      aa = cplsmoothDpAll(anow)
c      aa = DEgrowthcpl(anow)
      print *, aa
c      filename = 'Dm_clustering_-1_0.5.txt'
c
c      open(unit=10, file=trim(filename), status='old')
c
c      num_lines = 0
c      do
c        read(10, *, iostat=i) x(num_lines+1)
c        if (i /= 0) exit
c        num_lines = num_lines + 1
c        if (num_lines >= max_lines) exit
c      enddo
c      print *, x(100)
c      close(10)
c      
c      ll=4.5 / 0.25+0.00001
c      print *, ll 

c      w0 = -0.5
c      w1 =  0
c      CsDE = 0
c      anow = 1 


c      aa=HofEz(amax,anow,omep,omeplam,w0,w1)
c      print *, aa
c      bb=ScaleDE(amax,anow,w0,w1)
c      red= amax/anow-1
c      print *,bb,  aa, res, cc, dd
c      anow = 50
c      do i=1, amax
c        anow = i
c        aa = PoissonCor(anow,amax,omep,omeplam,w0,w1,0)
c        print *, anow/amax, aa
c       enddo 
c      Cofw = stCw(amax, anow, omep, omeplam, w0,w1,0)
c         dd = CwOveraH3(anow)
c       print *, Cw(amax,anow,omep, omeplam,w0,w1)
c       print *, HofEz(amax,anow,omep, omeplam,w0,w1)
c      dd = cplClusterDpall(anow)
c      ee = cplOmegaMatter(anow)
c       call qromo(cplIntMG, 0., anow, dd, midsql)
c       call qromo(cplIntMG, 0., anow, ff, midsqu)
c       call qsimp(cplIntMG, 0., anow, gg)
c       call qromb(cplIntMG,0,anow,ee)
c      bb = cplIntMG(anow)
c      aa = PoissonCor(anow,amax,omep,omeplam,w0,w1,0)
      
      endif
      stop
      end
