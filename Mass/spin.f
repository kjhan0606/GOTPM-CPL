c234567
      subroutine triaxialshape(nend,x,y,z,q,s,d,v,cx,cy,cz)
      real x(nend),y(nend),z(nend)
      real d(3),v(3,3)
      real fimom(3,3)
      real q,s
      real*8 cx,cy,cz
      integer n, np

      n = 3
      np = 3
      q = 1.
      s = 1.
c     do i = 1, 3
c     do j = 1, 3
c        v(j,i) = 0.
c     enddo
c     enddo
c     do i = 1, 3
c        v(i,i) = 1.
c     enddo
ccccc Iteration procedure ccccccccccccccccccccccccccccccccccc
      do 20 k = 1, 1
      do 10 j = 1, 3
      do 10 i = 1, 3
         fimom(i,j) = 0.
 10   continue
      do i = 1, nend
         distx = x(i) - cx
         disty = y(i) - cy
         distz = z(i) - cz
         fimom(1,1) = fimom(1,1) + distx*distx
         fimom(1,2) = fimom(1,2) + distx*disty
         fimom(1,3) = fimom(1,3) + distx*distz
         fimom(2,2) = fimom(2,2) + disty*disty
         fimom(2,3) = fimom(2,3) + disty*distz
         fimom(3,3) = fimom(3,3) + distz*distz
      enddo
      fimom(2,1) = fimom(1,2)
      fimom(3,1) = fimom(1,3)
      fimom(3,2) = fimom(2,3)
      call jacobi(fimom,np,n,d,v,nrot)
      call eigsrt(d,v,np,n)
      xxm = d(1)
      yym = d(2)
      zzm = d(3)
      q = sqrt(yym/xxm)
      s = sqrt(zzm/xxm)
 20   continue
      tmp = v(2,1)*v(3,2) - v(3,1)*v(2,2)
c---  conter-clockwise rotation.
      if(tmp/v(1,3) .le. -0.) then
          v(1,3) = -v(1,3)
          v(2,3) = -v(2,3)
          v(3,3) = -v(3,3)
      endif
ccccc Iteration procedure ccccccccccccccccccccccccccccccccccc
      return
      end

      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      INTEGER n,np,nrot,NMAX
      REAL a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     *g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      pause 'too many iterations in jacobi'
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
      SUBROUTINE eigsrt(d,v,n,np)
      INTEGER n,np
      REAL d(np),v(np,np)
      INTEGER i,j,k
      REAL p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).ge.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +)-*1a311.
