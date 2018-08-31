*
*     programa para el calculo de integrales (AA|BB)
*     a partir de la subrutina coul0sim
*
c***********************************************************************
c     subroutine heq                                                   *
c                                                                      *
c     subroutine for calculating the basic matrices h and hg in case   *
c     of equal exponents:   exa .eq. exb.                              *
c     it uses the one-range translation formula of the stos.           *
c     see reference  (5).                                              *
c                                                                      *
c***********************************************************************
      subroutine heq (exa,a,lmax,nmax,h,hg)
      implicit real * 8 (a-h,o-z)
c     parameter (ldim = 25, ndim = 25)
      parameter (ldim = 40, ndim = 40)
      include 'comun.inc'
      dimension  h(0:ldim,0:ndim),hg(0:ldim,0:ndim)
      parameter (pi = 3.1415926535897932384626433832795059d0)
      x = exa * a
      dx=(pi+pi) * dexp(-x)
c
c   equal exponents
c
      zl=dx/x
      ax=.5d0/exa
      do 10 l=0,lmax
         l1=l+1
         dln=ax* real(l1+l1)
         zl=zl*a
         zln=zl*a
         h(l,0)=zl * ri(l1)
         h(l,1)=(zln+dln*h(l,0)) * ri(l1+1)
         hg(l,0)=h(l,1) * ri(l1)
         do 20 n=2, nmax-l
            nm1=n-1
            zln=zln*a
            cln=dln+ax
            h(l,n)=(zln+ real(n)*cln*h(l,nm1)) * ri(l1+n)
            hg(l,nm1)=(h(l,n)+ real(nm1)*dln*hg(l,n-2)) * ri(l+n)
            dln=cln
   20    continue
   10 continue
      return
      end
c***********************************************************************
c     subroutine hpeq                                                  *
c                                                                      *
c     subroutine for calculating the matrix  h                         *
c         h(l,n) = < (n+l) l 0 ] 0 0 0 >                               *
c     see eq. (d1) of reference (3)                                    *
c                                                                      *
c***********************************************************************
      subroutine hpeq (lm0p,nmax,a,fl,gl,ak,ai,h)
      implicit real*8 (a-h,o-z)
c     parameter (ldim = 25, ldim2 = 30, ndim = 25)
      parameter (ldim = 40, ldim2 = 40, ndim = 40)
      parameter (pi = 3.1415926535897932384626433832795059d0)
      dimension fl(0:ldim2), gl(0:ldim2), ai(0:ldim,0:ndim)
      dimension ak(0:ldim,0:ndim)
      dimension h(0:ldim,0:ndim)
      include 'comun.inc'
      an4pi = 4.d0 * pi * a
      lsup = lm0p
      do 150 n = 0, nmax
         bux = 1.d0
         aux = an4pi
         l21 = 1
         do 160 l = 0, lsup
            h(l,n) = bux * ( aux * ri(l21) * ( fl(l)*ai(l,n)
     +               + gl(l)*ak(l,n) ) )
            l21 = l21 + 2
            aux = aux * a
  160    continue
         lsup = lsup - 1
         an4pi = an4pi * a
  150 continue

      return
      end
c***********************************************************************
c     subroutine hgran                                                 *
c                                                                      *
c     subroutine for calculating the matrix  hg                        *
c         hg(l,n) = < (n+l) l 0 ] 1 0 0 >                              *
c                                                                      *
c***********************************************************************
      subroutine hgran (lm0p, nmax, x, y, a, expny, gl, ai, ak,
     ,                  h, hg)
      implicit real*8 (a-h,o-z)
c     parameter (ldim = 25, ldim2 = 30, ndim = 25)
      parameter (ldim = 40, ldim2 = 40, ndim = 40)
      parameter (pi = 3.1415926535897932384626433832795059d0)
      dimension  gl(0:ldim2), ai(0:ldim,0:ndim)
      dimension ak(0:ldim,0:ndim)
      dimension h(0:ldim,0:ndim), hg(0:ldim,0:ndim)
      include 'comun.inc'
      lsup = lm0p

c     rows with   l > 0: relation (26)  of reference 5

      do 10 n = 0, nmax-2
         l21 = 3
         n2 = n + 2
         do 20 l = 1, lsup
            hg(l,n) = ri(l21) * y * ( h(l-1,n2) - h(l+1,n) )
            l21 = l21 + 2
   20    continue
         lsup = lsup - 1
   10 continue

c     row with  l = 0: equation  (a.11)  of reference 2.

      y3 = y * ri(3)
      aux = 4.d0*pi*a*a
      bux = gl(0)+y*y3 * gl(1)
      yi = 1.d0 / y
      hg(0,0) = aux * ((ai(0,0) - y3 * ai(1,0)) * expny
     &  + ( gl(0) * ak(1,0) - bux * ak(0,0)) * yi)

      if ((x-y) .gt. 1.d-10) then

c     different exponents: equation (17) of reference 7

         do 263 n = 1, nmax-1
            aux = aux * a
            hg(0,n) = aux * ((ai(0,n) - y3 * ai(1,n)) * expny
     &              + ( gl(0) * ak(1,n) - bux * ak(0,n)) * yi)
  263    continue


      else

c     equal exponents: equation (34) of reference 5

         axmed = .5d0 * a / x
         do 40 n = 1, nmax
            n1 = n + 1
            hg(0,n) =  axmed * n * hg(0,n-1) + ri(n1) * h(0,n1)
   40    continue

      endif
      return
      end
c***********************************************************************
c     subroutine kpos                                                  *
c                                                                      *
c     subroutine for calculating the matrix  ak                        *
c        ak(l,n) = exp(-b) * k(l,n) / b**(n+1)                         *
c                                                                      *
c     where b = (exa+exb)*a                                            *
c     and k(l,n) is defined in equation (9) of reference 2             *
c                                                                      *
c***********************************************************************
      subroutine kpos (lmax0, y2, bi, ypb, expnb, ak)
      implicit real * 8 (a-h,o-z)
c     parameter (ldim = 25, ndim = 25)
      parameter (ldim = 40, ndim = 40)
      dimension ak(0:ldim,0:ndim)
      include 'comun.inc'

      ak(0,0)= expnb * bi
      do 10 n=0,lmax0-1
         aux=( (n+1)*ak(0,n) + expnb ) * bi
         ak(0,n+1)=aux
   10 continue
c
c     calculation of  ak(1,n)  by means of the relation:
c         ak(1,n) =  y*expnb/(x+y) + (1 + (n+1)*y/(x+y)) ak(0,n)
c
      aux = ypb * expnb
      bux = 1.d0 + ypb
      do 50 n = 0, lmax0
         ak(1,n) = aux + bux * ak(0,n)
         bux = bux + ypb
   50 continue
c
c     calation of the remaining  ak(l,n)  by means of
c        ak(l+1,n) =  ak(l,n) + y2 * ak(l-1,n+2) / ((2l+1)(2l-1))
c     equation (25) of reference (2)
c
      nmaxl = lmax0 - 1
      l2m1 = 1
      l2p1 = 3
      do 60 l = 1, lmax0-1
         lm1 = l - 1
         lp1 = l + 1
         auxl = y2 * ri(l2m1) * ri(l2p1)
         nmaxl = nmaxl - 1
         do 70 n = 0, nmaxl
            ak(lp1,n) = ak(l,n) + auxl * ak(lm1,n+2)
   70    continue
         l2m1 = l2p1
         l2p1 = l2p1 + 2
   60 continue
      return
      end
c***********************************************************************
c     subroutine ipos                                                  *
c                                                                      *
c     subroutine for calculating the matrix  ai                        *
c        ai(l,n) = exp(-x) * i(l,n)                                    *
c                                                                      *
c     i(l,n) is defined in equation (8) of reference 2                 *
c                                                                      *
c***********************************************************************
      subroutine ipos (lm, lmax0, nmax, x, y2, expnx, x2my2, gl, f, ai,
     ,                 ierr)
      implicit real*8 (a-h,o-z)
      include 'comun.inc'
c     parameter (ldim = 25, ldim2 = 30, ndim = 25, ndim2 = 150)
      parameter (ldim = 40, ldim2 = 40, ndim = 40, ndim2 = 500)
      dimension gl(0:ldim2), ai(0:ldim,0:ndim), f(0:ndim2)

c     i(lm,0)  by means of equation (8) of reference (2)
      ierr = 0
      y2med = .5d0 * y2
      sum1 = 0.d0
      aux1 = 1.d0
      indi = 2*lm+3
      jfin = ndim2/2 - indi
      do 10 j = 0, jfin
         s1 = aux1 * f(indi)
         sum1 = sum1 + s1
         if (s1 .le. 1.d-15*sum1) go to 500
            j1 = j + 1
            aux1 = aux1 * y2med / (j1*indi)
            indi = indi + 2
   10 continue
      ierr = 1
  500 continue
c     ai(lm,0) = sum1 * expnx   ; ahora la exp(-x) va en las f's
      ai(lm,0) = sum1 

c     ai(l,0)  and  ai(l,1)  by means of equations (17) and (18)
c         of reference (2)

      l2 = lm + lm
      do 20 l = lm, 1, -1
         lm1 = l - 1
         l2p1 = l2 + 1
         ai(lm1,1) = ( expnx * gl(l) + x * ai(l,0) ) * ri(l2p1)
         ai(lm1,0) = ( expnx * gl(lm1) + ( expnx * x * gl(l) +
     +     x2my2 * ai(l,0) ) * ri(l2p1) ) * ri(l2)
         l2 = l2 - 2
   20 continue

c     remaining elements of the  ai  matrix: equation (19) of ref. (2)

      lsup = lmax0 - 1
      do 30 n = 0, nmax - 2
         n2 = n + 2
         l21 = 3
         do 40 l = 1, lsup
            l23 = l21 + 2
            ai(l-1,n2) = ai(l,n) + y2 * ai(l+1,n) * ri(l21) * ri(l23)
            l21 = l23
   40    continue
         lsup = lsup - 1
   30 continue
      return
      end
c***********************************************************************
c     subroutine funfl                                                 *
c                                                                      *
c     subroutine for calculating the auxiliary functions:              *
c        fl(l)   = exp(-y) * phi(l)    (capital phi)                   *
c     where  phi  is defined in equation (6) of reference 2            *
c                                                                      *
c***********************************************************************
      subroutine funfl (expny, lm0p, y, y2, fl)
      implicit real*8(a-h,o-z)
c     parameter (ldim2 = 30)
      parameter (ldim2 = 40)
      include 'comun.inc'
      dimension fl(0:ldim2)
      fl(0) = expny
      fl(1) = expny * (1.d0 + y)
      l2m1 = 1
      l2p1 = 3
      do 10 l = 1, lm0p-1
         fl(l+1) = fl(l) + y2 * fl(l-1) * ri(l2m1) * ri(l2p1)
         l2m1 = l2p1
         l2p1 = l2p1 + 2
  10  continue
      return
      end
c***********************************************************************
c     subroutine fn3                                                   *
c   OJO    MODIFICACION EN CURSO: meto la exponencial                  *
c     subroutine for calculating the auxiliary functions:              *
c        f(n) = exp(-x) 1f1 (1;n;x) / (n-1) = gamma(n-1,x)/x**(n-1)    *
c     where  1f1(1;n;x)  is defined in equation (20) reference 2.      *
c                                                                      *
c***********************************************************************
      subroutine fn3 (n, n0, x, f, ierr)
      implicit real*8(a-h,o-z)
      logical t1 , t2
      parameter (ndim2 = 500)
      dimension f(0:ndim2)
      include 'comun.inc'
c   si el argumento o n son grandes no se usan pades.
      ierr = 0
      ex = dexp(-x)
      if (x.gt.1.d2 .or. n.gt.200) then
         m = n
         xn = real(n)
         t1 = x.gt.7.d2
         t2 = x.gt.(-35.d0+1.41d0*n+5.d-4*n**2) .and. n.gt.100
         if (t1.or.t2) then
c           print*,'fn3. recurre hacia arriba'
            xi = 1.d0 / x
            f(2) = ( 1.d0 - ex ) * xi
            do 22 i = 3 , n
               f(i) = xi * ( f(i-1)*real(i-2) - ex )
   22       continue
            return
         else
c           print*,'fn3. serie + rec. hacia abajo'
            c = 1.d0
            y = 1.d0
            aux0 = xn-1.d0
            aux = aux0
            do 1 j = 1 , 1000
               aux = aux + 1.d0
               c = c * x / aux
               y = y + c
               if ( c.lt.1.d-15*y ) go to 2
    1       continue
            print*,'no converge la serie de fn'
            ierr = 1
    2       continue
            f(n) = y * ex / aux0
            go to 1500
         endif
      endif

      if(n.gt.105) goto 150
        if(n.gt.50) then
           call fn105(x,fm)
           m=105
        else
           if(n.gt.25) then
              call fn50(x,fm)
              m=50
           else
              call fn25(x,fm)
              m=25
           endif
        endif
      goto 1000
  150 if(n.gt.150) then
         call fn200(x,fm)
         m=200
      else
         call fn150(x,fm)
         m=150
      endif
 1000 f(m)=fm * ex / real(m-1)
 1500 rj=real(m-1)
      do 2000 j=m-1,max(n0,2),-1
      rj=rj-1.d0
 2000 f(j) = (ex + x*f(j+1) ) / rj
      return
      end
c***********************************************************************
c     subroutine fn25                                                  *
c                                                                      *
c     subroutine for calculating  1f1(1;25;x)                          *
c                                                                      *
c***********************************************************************
      subroutine fn25(x,f25)
      implicit real*8(a-h,o-z)
c coefficients for the polynomial f(25)  with 0<x<2
      data a0/  1.000000000000000d+00/, a1/  4.000000000000000d-02/,
     &a2/  1.538461538461538d-03/, a3/  5.698005698005698d-05/,
     &a4/  2.035002035002035d-06/, a5/  7.017248396558740d-08/,
     &a6/  2.339082798852913d-09/, a7/  7.545428383396495d-11/,
     &a8/  2.357946369811405d-12/, a9/  7.145292029731528d-14/,
     &a10/  2.101556479332802d-15/, a11/  6.004447083808007d-17/
c rational approximations for  f(25)  with 2<x<20
      data b1/  0.999999999999860d+00/, c1/  0.100000000000000d+01/,
     &b2/ -0.104714737640631d+00/, c2/ -0.144714737640909d+00/,
     &b3/  0.465734042761888d-02/, c3/  0.890746839503547d-02/,
     &b4/ -0.107890542447286d-03/, c4/ -0.298531277438548d-03/,
     &b5/  0.129752008992129d-05/, c5/  0.574582565486514d-05/,
     &b6/ -0.642051572571287d-08/, c6/ -0.602004159824384d-07/,
     &                             c7/  0.267808535531380d-09/,
     &                             c8/  0.417000629047170d-14/,
     &                             c9/  0.562537647721706d-16/
c rational approximations for  f(25)  with 20<x<30
      data d1/  0.997838292575181d+00/, e1/  0.100000000000000d+01/,
     &d2/ -0.404704646482292d-01/, e2/ -0.814870574903038d-01/,
     &d3/  0.586936497614565d-03/, e3/  0.252333823052715d-02/,
     &d4/  0.244097392117162d-04/, e4/ -0.351886121051830d-04/,
     &d5/ -0.217369731011646d-05/, e5/  0.186436972724487d-06/,
     &d6/  0.126577273940265d-06/,
     &d7/ -0.479533856185303d-08/,
     &d8/  0.125942566294675d-09/,
     &d9/ -0.198132854998548d-11/,
     &d10/  0.164072629016418d-13/
c coefficients for the truncated exponential of  f(25)  with 30<x<86
      data f0/  0.100000000000000d+01/, f1/  0.100000000000000d+01/,
     &f2/  0.500000000000000d+00/, f3/  0.166666666666667d+00/,
     &f4/  0.416666666666667d-01/, f5/  0.833333333333333d-02/,
     &f6/  0.138888888888889d-02/, f7/  0.198412698412698d-03/,
     &f8/  0.248015873015873d-04/, f9/  0.275573192239859d-05/,
     &f10/  0.275573192239859d-06/, f11/  0.250521083854417d-07/,
     &f12/  0.208767569878681d-08/, f13/  0.160590438368216d-09/,
     &f14/  0.114707455977297d-10/, f15/  0.764716373181981d-12/,
     &f16/  0.477947733238738d-13/, f17/  0.281145725434552d-14/,
     &f18/  0.156192069685862d-15/, f19/  0.822063524662433d-17/,
     &f20/  0.411031762331216d-18/, f21/  0.195729410633913d-19/,
     &f22/  0.889679139245057d-21/, f23/  0.386817017063068d-22/
      if(x.gt.2.d0) goto 100
        if(x.gt.1.d0) goto 20
          if(x.gt.1.d-1) goto 10
      f25= a0+x*(a1+x*(a2+x*(a3+x*(a4+x*a5))))
      return
   10 f25= a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*(a6+x*(a7+x*(a8+x*a9))))))
     &     ))
      return
   20 f25= a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*(a6+x*(a7+x*(a8+x*(a9+x*
     &     (a10+x*a11))))))))))
      return
  100 if(x.gt.3.d1) goto 150
        if(x.gt.2.d1) goto 130
c--------------
c chebyshev approximations 6/9
c--------------
      fn=b1+x*(b2+x*(b3+x*(b4+x*(b5+x*b6))))
      fd=c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*(c8+x*c9)))))))
      f25=fn/fd
      return
c--------------
c chebyshev approximations 10/5
c--------------
  130 fn=d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*(d8+x*(d9+x*d10
     &   ))))))))
      fd=e1+x*(e2+x*(e3+x*(e4+x*e5)))
      f25=fn/fd
      return
  150   if(x.gt.86.d0) goto 170
c--------------
c     closed formula
c--------------
      f25=f0+x*(f1+x*(f2+x*(f3+x*(f4+x*(f5+x*(f6+x*(f7+x*(f8+x*(f9+x*(
     &    f10+x*(f11+x*(f12+x*(f13+x*(f14+x*(f15+x*(f16+x*(f17+x*(f18+
     &    x*(f19+x*(f20+x*(f21+x*(f22+x*f23))))))))))))))))))))))
      f25=(dexp(x)-f25)/x**24
      f25=f25*6.20448401733239d23
      return
c--------------
c     assymptotic formula  24"*exp(x)/x**24
c--------------
  170 f25=dexp(x)/x**24
      f25=f25*6.20448401733239d23
      return
      end
c***********************************************************************
c     subroutine fn50                                                  *
c                                                                      *
c     subroutine for calculating  1f1(1;50;x)                          *
c                                                                      *
c***********************************************************************
      subroutine fn50(x,f50)
      implicit real*8(a-h,o-z)
c coefficients for the polynomial f(50)  with 0<x<1
      data a0/  1.000000000000000d+00/, a1/  2.000000000000000d-02/,
     &a2/  3.921568627450980d-04/, a3/  7.541478129713423d-06/,
     &a4/  1.422920401832721d-07/, a5/  2.635037781171706d-09/,
     &a6/  4.790977783948557d-11/, a7/  8.555317471336708d-13/,
     &a8/  1.500932889708194d-14/
c rational approximations for  f(50)  with 1<x<25
      data b1/  0.999999999999990d+00/, c1/  0.100000000000000d+01/,
     &b2/ -0.606384113439031d-01/, c2/ -0.806384113439206d-01/,
     &b3/  0.142011960962117d-02/, c3/  0.264073097376706d-02/,
     &b4/ -0.151740904702204d-04/, c4/ -0.439072816703196d-04/,
     &b5/  0.623262071180348d-07/, c5/  0.370731842737299d-06/,
     &b6/  0.666616564411615d-14/, c6/ -0.127193705182704d-08/
c rational approximations for  f(50)  with 25<x<40
      data d1/  0.999997802698584d+00/, e1/  0.100000000000000d+01/,
     &d2/ -0.376177104812775d-01/, e2/ -0.576183817242153d-01/,
     &d3/  0.498893376540307d-03/, e3/  0.125919438921593d-02/,
     &d4/ -0.223360940137912d-05/, e4/ -0.123705033636332d-04/,
     &d5/ -0.862623601247297d-10/, e5/  0.460958759313318d-07/,
     &d6/  0.138129726359552d-10/,
     &d7/ -0.158664944212928d-12/,
     &d8/  0.250464951989587d-14/
c rational approximations for  f(50)  with 40<x<50
      data f1/  0.997721752279347d+00/, g1/  0.100000000000000d+01/,
     &f2/ -0.175113983565115d-01/, g2/ -0.379665571317971d-01/,
     &f3/ -0.115276526520580d-03/, g3/  0.290872071261466d-03/,
     &f4/  0.354986033428588d-05/, g4/  0.325220662333212d-05/,
     &f5/ -0.749562821931817d-08/, g5/  0.588213247464554d-08/,
     &                             g6/ -0.141105303546913d-08/,
     &                             g7/  0.887136240917339d-11/,
     &                             g8/  0.284609486633007d-12/,
     &                             g9/ -0.730925594571724d-14/,
     &                             g10/  0.731598302588279d-16/,
     &                             g11/ -0.275645588904962d-18/
c rational approximations for  f(50)  with 50<x<60
      data fb1/  0.159050723434975d+01/, gb1/  0.100000000000000d+01/,
     &fb2/ -0.875646510575019d-01/, gb2/ -0.112428987895828d-01/,
     &fb3/  0.650475606259438d-02/, gb3/ -0.209710987123228d-03/,
     &fb4/ -0.229548632811819d-03/, gb4/  0.147669473956403d-05/,
     &fb5/  0.386038262239120d-05/, gb5/  0.966410592843145d-08/,
     &fb6/  0.580251601022677d-09/, gb6/ -0.589575837899969d-09/,
     &fb7/ -0.118417896276571d-08/, gb7/  0.577211111237690d-11/,
     &fb8/  0.168440391139320d-10/, gb8/  0.182155907646624d-12/,
     &fb9/  0.161368056865524d-13/, gb9/  0.470631278067002d-14/,
     &fb10/ -0.216442002377797d-14/, gb10/ -0.163358709513760d-15/,
     &fb11/  0.151342326832851d-16/, gb11/  0.364357985656402d-18/,
     &                               gb12/  0.152687431201180d-19/,
     &                               gb13/ -0.918169699087600d-22/
c rational approximations for  f(50)  with 60<x<80
      data h1/  0.360937927130168d-24/, o1/  0.100000000000000d+01/,
     &h2/ -0.398845543839481d-25/, o2/ -0.848187362634474d-01/,
     &h3/  0.184024482155724d-26/, o3/  0.297530185782044d-02/,
     &h4/ -0.464877350556152d-28/, o4/ -0.521205192398682d-04/,
     &h5/  0.714282522507662d-30/, o5/  0.391424066764965d-06/,
     &h6/ -0.715833421288194d-32/, o6/  0.943580379783575d-09/,
     &h7/  0.516696371565316d-34/, o7/ -0.351067318449105d-10/,
     &h8/ -0.296640354285062d-36/, o8/  0.170527824281788d-12/,
     &h9/  0.117982638246899d-38/,
     &h10/ -0.163123552468383d-41/
c rational approximations for  f(50)  with 80<x<100
      data p1/  0.125331451794748d-30/, q1/  0.100000000000000d+01/,
     &p2/ -0.199745401977752d-32/, q2/ -0.614056661792172d-01/,
     &p3/ -0.252171922458197d-33/, q3/  0.143321129914182d-02/,
     &p4/  0.117958788425682d-34/, q4/ -0.150798106936193d-04/,
     &p5/ -0.228725953865918d-36/, q5/  0.606588922853750d-07/,
     &p6/  0.245287682214720d-38/
     &p7/ -0.161107855626302d-40/
     &p8/  0.681416557134086d-43/
     &p9/ -0.168399467728711d-45/
     &p10/  0.185258048923584d-48/
      if(x.gt.25.d0) goto 100
        if(x.gt.1.d0) goto 20
          if(x.gt.1.d-1) goto 10
      f50= a0+x*(a1+x*(a2+x*(a3+x*(a4+x*a5))))
      return
   10 f50= a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*(a6+x*(a7+x*a8)))))))
      return
c-------------
c chebyshev approximations 6/6
c-------------
   20 fn=b1+x*(b2+x*(b3+x*(b4+x*(b5+x*b6))))
      fd=c1+x*(c2+x*(c3+x*(c4+x*(c5+x*c6))))
      f50=fn/fd
      return
  100 if(x.gt.6.d1) goto 150
        if(x.gt.4.d1) goto 130
c--------------
c chebyshev approximations 8/5
c--------------
      fn=d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8))))))
      fd=e1+x*(e2+x*(e3+x*(e4+x*e5)))
      f50=fn/fd
      return
  130 if(x.gt.5.d1) goto 135
c--------------
c chebyshev approximations 5/11
c--------------
      fn=f1+x*(f2+x*(f3+x*(f4+x*f5)))
      fd=g1+x*(g2+x*(g3+x*(g4+x*(g5+x*(g6+x*(g7+x*(g8+x*(g9+x*(g10+x*g11
     &   )))))))))
      f50=fn/fd
      return
c--------------
c chebyshev approximations 11/13
c--------------
  135 fn=fb1+x*(fb2+x*(fb3+x*(fb4+x*(fb5+x*(fb6+x*(fb7+x*(fb8+x*(fb9+x*
     &   (fb10+x*fb11)))))))))
      fd=gb1+x*(gb2+x*(gb3+x*(gb4+x*(gb5+x*(gb6+x*(gb7+x*(gb8+x*(gb9+x*
     &   (gb10+x*(gb11+x*(gb12+x*gb13)))))))))))
      f50=fn/fd
      return
  150   if(x.gt.8.d1) goto 170
c--------------
c chebyshev approximations 10/8   maximum relative error = .57e-12
c--------------
      fn=h1+x*(h2+x*(h3+x*(h4+x*(h5+x*(h6+x*(h7+x*(h8+x*(h9+x*h10))))))
     &   ))
      fd=o1+x*(o2+x*(o3+x*(o4+x*(o5+x*(o6+x*(o7+x*o8))))))
      f50=(fn/fd)*(6.d1/x)**49
      f50=f50*dexp(x)
      return
c--------------
c chebyshev approximations 10/5   maximum relative error = .17e-12
c--------------
  170 fn=p1+x*(p2+x*(p3+x*(p4+x*(p5+x*(p6+x*(p7+x*(p8+x*(p9+x*p10))))))
     &   ))
      fd=q1+x*(q2+x*(q3+x*(q4+x*q5)))
      f50=(fn/fd)*(8.d1/x)**49
      f50=f50*dexp(x)
      return
      end
c***********************************************************************
c     subroutine fn105                                                 *
c                                                                      *
c     subroutine for calculating  1f1(1;105;x)                         *
c                                                                      *
c***********************************************************************
      subroutine fn105(x,f105)
      implicit real*8(a-h,o-z)
c coefficients for the polynomial f(105)  with 0<x<1
      data a0/  1.000000000000000d+00/, a1/  9.523809523809524d-03/,
     &a2/  8.984725965858041d-05/, a3/  8.396940155007515d-07/,
     &a4/  7.774944587969921d-09/, a5/  7.132976686210937d-11/,
     &a6/  6.484524260191760d-13/, a7/  5.841913747920505d-15/
c rational approximations for  f(105)  with 1<x<40
      data b1/  0.1000000000000011d+01/,c1/  0.1000000000000000d+01/,
     &b2/ -0.2550001591594343d-01/,c2/ -0.3502382543973739d-01/,
     &b3/  0.2202419723388067d-03/,c3/  0.4639549549574303d-03/,
     &b4/ -0.6436784034595878d-06/,c4/ -0.2755196298058278d-05/,
     &b5/  0.8216485139813122d-14/,c5/  0.6189243597971857d-08/
c rational approximations for  f(105)  with 40<x<70
       data d1/  0.999999951254663d+00/, e1/  0.100000000000000d+01/,
     &d2/ -0.234866121409395d-01/, e2/ -0.330104286816774d-01/,
     &d3/  0.187088497313817d-03/, e3/  0.411626699054356d-03/,
     &d4/ -0.503961309069964d-06/, e4/ -0.229802705914890d-05/,
     &d5/  0.115267974932120d-12/, e5/  0.484653372049368d-08/,
     &d6/  0.463219099562710d-14/
c rational approximations for  f(105)  with 70<x<100
      data f1/  0.100560610512446d+01/, g1/  0.100000000000000d+01/,
     &f2/  0.323165497961713d-02/, g2/ -0.549539509516998d-02/,
     &f3/  0.332246465874758d-04/, g3/ -0.556620479544504d-04/,
     &f4/ -0.361560071773655d-05/, g4/ -0.144968104971354d-05/,
     &f5/  0.506444358322703d-07/, g5/  0.155746773529450d-07/,
     &f6/ -0.879093755346411d-09/, g6/  0.224932045729219d-10/,
     &f7/  0.127280408793787d-10/, g7/  0.120685885754946d-11/,
     &f8/ -0.105087439383639d-12/, g8/ -0.197264163522160d-13/,
     &f9/  0.542398316691009d-15/, g9/  0.184815677040377d-16/,
     &f10/ -0.189978437550371d-17/, g10/  0.588036182278456d-18/,
     &f11/  0.359750457178732d-20/, g11/ -0.190998574782540d-20/
      if(x.gt.1.d0) goto 100
         if(x.gt.4.d-1) goto 10
      f105= a0+x*(a1+x*(a2+x*(a3+x*(a4+x*a5))))
       return
   10 f105= a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*(a6+x*a7))))))
       return
  100    if(x.gt.4.d1) goto 110
c---------------
c chebyshev approximations 5/5
c---------------
       fn= b1+x*(b2+x*(b3+x*(b4+x*b5)))
       fd= c1+x*(c2+x*(c3+x*(c4+x*c5)))
       f105=fn/fd
       return
  110       if(x.gt.7.d1) goto 120
c--------------
c chebyshev approximations 6/5
c--------------
       fn= d1+x*(d2+x*(d3+x*(d4+x*(d5+x*d6))))
       fd= e1+x*(e2+x*(e3+x*(e4+x*e5)))
       f105=fn/fd
       return
c--------------
c chebyshev approximations 11/11
c--------------
  120  fn= f1+x*(f2+x*(f3+x*(f4+x*(f5+x*(f6+x*(f7+x*(f8+x*(f9+x*(f10
     &    +x*f11)))))))))
       fd= g1+x*(g2+x*(g3+x*(g4+x*(g5+x*(g6+x*(g7+x*(g8+x*(g9+x*(g10
     &    +x*g11)))))))))
       f105=fn/fd
       return
       end
c***********************************************************************
c     subroutine fn150                                                 *
c                                                                      *
c     subroutine for calculating  1f1(1;150;x)                         *
c                                                                      *
c***********************************************************************
      subroutine fn150(x,f150)
      implicit real*8(a-h,o-z)
c coefficients for the polynomial f(150)  with 0<x<2
      data d0/  1.000000000000000d+00/, d1/  6.666666666666666d-03/,
     &d2/  4.415011037527593d-05/, d3/  2.904612524689206d-07/,
     &d4/  1.898439558620396d-09/, d5/  1.232752960143114d-11/,
     &d6/  7.953244904149124d-14/, d7/  5.098233912916105d-16/
c rational approximations for  f(150)  with 2<x<90
      data e1/  0.999999999999925d+00/, f1/  0.100000000000000d+01/,
     &e2/ -0.173950678287829d-01/, f2/ -0.240617344955064d-01/,
     &e3/  0.102079397417238d-03/, f3/  0.218340850358374d-03/,
     &e4/ -0.201827108800031d-06/, f4/ -0.885565797914679d-06/,
     &e5/  0.310273864891050d-13/, f5/  0.135459243125970d-08/,
     &e6/  0.129121194327191d-15/,
     &e7/  0.409580383814601d-18/,
     &e8/  0.215314858882139d-20/
      if(x.gt.2.d0) goto 10
      f150=d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*d7))))))
      return
c------------
c chebyshev approximations 8/5
c------------
   10 fn=e1+x*(e2+x*(e3+x*(e4+x*(e5+x*(e6+x*(e7+x*e8))))))
      fd=f1+x*(f2+x*(f3+x*(f4+x*f5)))
      f150=fn/fd
      return
      end
c***********************************************************************
c     subroutine fn200                                                 *
c                                                                      *
c     subroutine for calculating  1f1(1;200;x)                         *
c                                                                      *
c***********************************************************************
      subroutine fn200(x,f200)
      implicit real*8(a-h,o-z)
c coefficients for the polynomial f(200)  with 0<x<2
      data a0/  1.000000000000000d+00/, a1/  5.000000000000000d-03/,
     &a2/  2.487562189054726d-05/, a3/  1.231466430225112d-07/,
     &a4/  6.066337094704984d-10/, a5/  2.973694654267149d-12/,
     &a6/  1.450582758179097d-14/
c rational approximations for  f(200)  with 2<x<100
      data b1/  0.999999999999988d+00/, c1/  0.100000000000000d+01/,
     &b2/ -0.138653653980072d-01/, c2/ -0.188653653980136d-01/,
     &b3/  0.646251749791174d-04/, c3/  0.134076380079788d-03/,
     &b4/ -0.101216076731367d-06/, c4/ -0.425456923774057d-06/,
     &b5/  0.455225367253129d-15/, c5/  0.508624450975601d-09/,
     &b6/  0.939486000203512d-18/
      if(x.gt.2.d0) goto 10
      f200=a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*a6)))))
      return
c------------
c chebyshev approximations 6/5
c------------
   10 fn=b1+x*(b2+x*(b3+x*(b4+x*(b5+x*b6))))
      fd=c1+x*(c2+x*(c3+x*(c4+x*c5)))
      f200=fn/fd
      return
      end
c***********************************************************************
c     subroutine fn2                                                   *
c                                                                      *
c     subroutine for calculating the auxiliary functions:              *
c        f(n) = 1f1 (1;n;b) / (n-1)                                    *
c     where  1f1(1;n;b)  is defined in equation (20) reference 2.      *
c                                                                      *
c***********************************************************************
      subroutine fn2 (n, n0, x, f, ierr)
      implicit real*8(a-h,o-z)
      parameter (ndim2 = 150)
      dimension f(0:ndim2)
      include 'comun.inc'
c   if the argument is too large, the general formula (infinite
c   expansion) is applied for obtaining the element with highest n
      ierr = 0
      if (x.gt.1.d2 .or. n.gt.200) then
         m = max( n, int(x) )
         if (m .gt. ndim2) m = ndim2
         bux = 1.d0
         sum = bux
         aux = m
         do 1 j=1,1000000
            bux = bux * x / aux
            sum = sum + bux
            aux = aux + 1.d0
            if ( bux.lt.1.d-15*sum ) goto 10
    1    continue
         ierr = 1
   10    f(m) = sum / dfloat(m-1)
         goto 1500
      endif
      if(n.gt.105) goto 150
        if(n.gt.50) then
           call fn105(x,fm)
           m=105
        else
           if(n.gt.25) then
              call fn50(x,fm)
              m=50
           else
              call fn25(x,fm)
              m=25
           endif
        endif
      goto 1000
  150 if(n.gt.150) then
         call fn200(x,fm)
         m=200
      else
         call fn150(x,fm)
         m=150
      endif
 1000 f(m)=fm / dfloat(m-1)
 1500 rj=dfloat(m-1)
      jn=max(25,n0)
      do 2000 j=m-1,jn,-1
      rj=rj-1.d0
 2000 f(j) = (1.d0 + x*f(j+1) ) / rj
      if (n0 .ge. 25) return
      f(24) = (1.d0 + x*f(25) ) * ri(23)
      f(23) = (1.d0 + x*f(24) ) * ri(22)
      f(22) = (1.d0 + x*f(23) ) * ri(21)
      f(21) = (1.d0 + x*f(22) ) * ri(20)
      f(20) = (1.d0 + x*f(21) ) * ri(19)
      f(19) = (1.d0 + x*f(20) ) * ri(18)
      f(18) = (1.d0 + x*f(19) ) * ri(17)
      f(17) = (1.d0 + x*f(18) ) * ri(16)
      f(16) = (1.d0 + x*f(17) ) * ri(15)
      f(15) = (1.d0 + x*f(16) ) * ri(14)
      f(14) = (1.d0 + x*f(15) ) * ri(13)
      f(13) = (1.d0 + x*f(14) ) * ri(12)
      f(12) = (1.d0 + x*f(13) ) * ri(11)
      f(11) = (1.d0 + x*f(12) ) * ri(10)
      f(10) = (1.d0 + x*f(11) ) * ri(9)
      f(9)  = (1.d0 + x*f(10) ) * ri(8)
      f(8)  = (1.d0 + x*f(9)  ) * ri(7)
      f(7)  = (1.d0 + x*f(8)  ) * ri(6)
      f(6)  = (1.d0 + x*f(7)  ) * ri(5)
      f(5)  = (1.d0 + x*f(6)  ) * ri(4)
      f(4)  = (1.d0 + x*f(5)  ) * ri(3)
      f(3)  = (1.d0 + x*f(4)  ) * ri(2)
      f(2)  = (1.d0 + x*f(3)  ) * ri(1)
      f(1)  = (1.d0 + x*f(2)  )
      return
      end
