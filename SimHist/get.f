c234567
      real ymin,ymax
	  real tickmin,tickmax

	  tickmin = 2
	  ymin=1067
	  tickmax = 10
	  ymax=139
	  yscale = (tickmax-tickmin)/(ymax-ymin)

	  open(1,file='data.txt')
10    read(1,*,end=20) x, y, year
      y = (y-ymin)*yscale + tickmin
	  print *,year, y
	  goto 10
20    continue
      close(1)
	  stop
	  end
