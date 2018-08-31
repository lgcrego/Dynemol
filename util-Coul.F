c  ********************************************************************
c
      subroutine mete(vbielx,ibielx,i,j,k,l,x)
      implicit real*8 (a-h,o-z)
      include 'comun.inc'
      integer*1 ibielx
      common / imete / imete
      common /kincom/ kint, kint1, kint2, isABCC, isABAC, isABCD

      dimension vbielx(mxbuff), ibielx(5*mxbuff)

      kint = kint + 1
      imete = imete + 1

c conviene ordenar los indices para estar seguros

      if (dabs(x) .gt. 1.d3) then
         print*, 'Error in the integral of indices = ', i, j, k, l
         print*, 'Value of the integral = ', x
         print*, 'Stop'
         stop
      endif

      if(i.gt.j) then
        in=i
        jn=j
      else
        in=j
        jn=i
      endif
      if(k.gt.l) then
        kn=k
        ln=l
      else
        kn=l
        ln=k
      endif
      ij = ind(in-1) + jn
      kl = ind(kn-1) + ln
      if (kl.gt.ij) then
         m=in
         in=kn
         kn=m
         m=jn
         jn=ln
         ln=m
      endif
      ijkl=1
      if(in.eq.jn) ijkl=ijkl+5
      if(in.eq.kn) ijkl=ijkl+1
      if(in.eq.ln) ijkl=ijkl-7
      if(jn.eq.kn) ijkl=ijkl+2
      if(jn.eq.ln) ijkl=ijkl+3
      if(kn.eq.ln) ijkl=ijkl+6
      vbielx(imete) = x

*     Usa los cuatro bits mas significativos del byte ibielx(5+ishft)
*     para ampliar los indices  in, jn, kn, ln  hasta   511
*     Esta forma de almacenar los indices condiciona su lectura,
*     por ello es necesario que el programa que lea las integrales
*     tenga la estructura apropiada (version rico2000 y ricoop2000
*     y posteriores; las anteriores no sirven)

      ishft = 5*(imete-1)
      ijkl = ijkl-128
      if (in .gt. 255) then
         ibielx(1+ishft) = in - 384
         ijkl = ijkl + 128
      else
         ibielx(1+ishft) = in - 128
      endif
      if (jn .gt. 255) then
         ibielx(2+ishft) = jn - 384
         ijkl = ijkl + 64
      else
         ibielx(2+ishft) = jn - 128
      endif
      if (kn .gt. 255) then
         ibielx(3+ishft) = kn - 384
         ijkl = ijkl + 32
      else
         ibielx(3+ishft) = kn - 128
      endif
      if (ln .gt. 255) then
         ibielx(4+ishft) = ln - 384
         ijkl = ijkl + 16
      else
         ibielx(4+ishft) = ln - 128
      endif
      ibielx(5+ishft) = ijkl
      if (imete.eq.mxbuff) then
         imete = 0
      endif
      return
      end
c
c   ********************************************************************
c
      subroutine zzzhv
c
c  subrutina para generar el data de las h�bridas
c
      implicit real*8 (a-h,o-z)
      include 'comun.inc'
      common / zzzh / z(izzzh)
      common / zzzdat / lmax , mxlb
      data pi/0.31415926535897932384626433832795059d+01/
 1000 format(4i4,d25.15)
      lmax = 2 * mxn
      mxlb = mxl
      do 1 m = 0 , 2*mxl
         if (m.eq.0) then
            z( indzzz(0,0,0,0) ) = 2.d0 * pi
         else
            do 2 k = 0, 2*m , 2
               aux = 0.d0
               if (k.ge.2) aux = aux - k*(k-1)/dfloat(2*k-1)*
     &                     z( indzzz(m-1,m-1,m-1,k-2) )
               if (k.le.2*m-4) aux = aux - (k+1)*(k+2)/dfloat(2*k+3)*
     &                     z( indzzz(m-1,m-1,m-1,k+2) )
               if (k.le.2*m-2) aux = aux + ( 2*k+1
     &           - (k+1)**2 / dfloat(2*k+3)
     &           - k**2 / dfloat(2*k-1) ) *
     &                     z( indzzz(m-1,m-1,m-1,k) )
               z( indzzz(m,m,m,k) ) = aux * (2*m-1)**2 / dfloat(2*k+1)
   2        continue
         endif
         do 3 l2 = m , 2*mxl
         do 3 l1 = l2, 2*mxn
         if(l1.gt.m) then
         do 4 k = l1-l2 , l1+l2 , 2
          aux = 0.d0
          if(l1-1.ge.m .and. k.le.l1+l2-2 .and. k.ge.iabs(l1-1-l2)-1)
     &    aux = aux + (k+1)* z( indzzz(m,l1-1,l2,k+1) )
          if(l1-1.ge.m .and. k.le.l1+l2 .and. k.ge.iabs(l1-1-l2)+1)
     &    aux = aux + k* z( indzzz(m,l1-1,l2,k-1) )
          aux = aux * (2*l1-1)/dfloat( (l1-m)*(2*k+1) )
          if(l1-2.ge.m .and. k.le.l1+l2-2 .and. k.ge.iabs(l1-2-l2))
     &    aux = aux - (l1+m-1)/dfloat(l1-m)*z( indzzz(m,l1-2,l2,k) )
          z( indzzz(m,l1,l2,k) ) = aux
    4    continue
         endif
    3    continue
   1  continue
      m = 0
      i = 0
      do 5 l2 = m , 2*mxl
      do 5 l1 = l2, 2*mxn
      do 5 k = l1-l2 , l1+l2 , 2
        i = i + 1
        z(i) = z(i)*2.d0
    5 continue
      return
      end

c***********************************************************************
c     subroutine dup                                                   *
c                                                                      *
c     subroutine for copying matrices. copies a  matrix 'a'  into a    *
c     matrix 'b' from row  0 to  lmax, and from column  0  to  nmax.   *
c                                                                      *
c***********************************************************************
      subroutine dup ( lmax, nmax, b, a )
      implicit real * 8 ( a-h,o-z )
c     parameter (ldim = 25, ndim = 25)
      parameter (ldim = 40, ndim = 40)
      dimension a(0:ldim,0:ndim), b(0:ldim,0:ndim)
      do 1 n = 0, nmax
      do 1 l = 0, lmax
         b(l,n) = a(l,n)
    1 continue
      return
      end
c***********************************************************************
c     subroutine reqmb                                                 *
c                                                                      *
c     recurrence relation for increasing the m quantum number of       *
c     the stos in the overlap integrals.                               *
c     see eq. (34) of reference (3)                                    *
c                                                                      *
c***********************************************************************
      subroutine reqmb (lmax, nmin, nmax, m, mb, h)
      implicit real * 8 (a-h,o-z)
c     parameter (ldim = 25, ndim = 25)
      parameter (ldim = 40, ndim = 40)
      dimension h(0:ldim,0:ndim)
      include 'comun.inc'
      aux = dosl1(mb)
      if (nmax.eq.lmax) then
         nmax = nmax - 1
      else
         nmax = nmax - 2
      endif
      lmax = lmax - 1
      lsup = lmax
      do 10 n = nmin, nmax
         h(m,n) = 0.d0
         do 20 l = m, lsup
            q = rll1(l-m) * h(l+1,n)
            if (l.gt.0) q = q - rll1(l+m) * h(l-1,n+2)
            h(l,n) = -aux * q * dosl1p(l)
   20    continue
      if (n.gt.0) lsup = lsup - 1
   10 continue
      return
      end
c***********************************************************************
c     subroutine reqlb                                                 *
c                                                                      *
c     recurrence relation for increasing the lb quantum number of      *
c     the sto centered at point 'b' in the overlap integrals.          *
c     see eq. (35) of reference (3)                                    *
c                                                                      *
c***********************************************************************
      subroutine reqlb (lmax, nmin, nmax, m, lb, mb, a, h, h2)
      implicit real * 8 (a-h,o-z)
c     parameter (ldim = 25, ndim = 25)
      parameter (ldim = 40, ndim = 40)
      dimension h(0:ldim,0:ndim), h2(0:ldim,0:ndim)
      include 'comun.inc'
      den = 1.d0/real(lb-mb+1)
      aux = dosl1(lb) * den
      bux = real(lb+mb) * den
      if (lb.gt.mb) then
         lmaux = lmax + 1
c ojo. willy 11-4-94
c----------------------
         nmaux = nmax+1
c----------------------
         call reqnb (lmaux, nmin, nmaux, m, a, h2)
      endif
      if (nmax.eq.lmax) then
         nmax = nmax - 1
      else
         nmax = nmax - 2
      endif
      lmax = lmax - 1
      lsup = lmax
      m1 = m - 1
      do 30 n = nmin, nmax
         do 40 l = m, lsup
            q = real(l-m1) * h(l+1,n)
            if (l.gt.m) q = q + real(l+m) * h(l-1,n+2)
            q = aux * (q * dosl1p(l) - a * h(l,n))
            if (lb.gt.mb) q = q - bux * h2(l,n)
            h2(l,n) = h(l,n)
            h(l,n) = q
   40    continue
         h2(lsup+1,n) = h(lsup+1,n)
         if (n.gt.0) lsup = lsup - 1
   30 continue
      do 50 n = nmax+1, nmax+2
      do 50 l = 0, lmax+1-n
         h2(l,n) = h(l,n)
   50 continue
      return
      end
c***********************************************************************
c     subroutine reqnb                                                 *
c                                                                      *
c     recurrence relation for increasing the nb quantum number of      *
c     the sto centered at point 'b' in the overlap integrals.          *
c     see eq. (36) of reference (3)                                    *
c                                                                      *
c***********************************************************************
      subroutine reqnb (lmax, nmin, nmax, m, a, h)
      implicit real * 8 (a-h,o-z)
      include 'comun.inc'
c     parameter (ldim = 25, ndim = 25)
      parameter (ldim = 40, ndim = 40)
      dimension h(0:ldim,0:ndim)
      a2 = a * a
      dosa = a + a
      lmax = lmax - 1
      nmax = nmax - 2
      lsup = lmax
      m1 = m - 1
      do 10 n = nmin, nmax
         do 20 l = m, lsup
            q = real(l-m1) * h(l+1,n)
            if (l.gt.m) q = q + real(l+m) * h(l-1,n+2)
            h(l,n) = h(l,n+2) + a2 * h(l,n) - dosa * q * dosl1p(l)
   20    continue
         if (n.gt.0) lsup = lsup - 1
   10 continue
      return
      end
c***********************************************************************
c     subroutine h2r                                                   *
c                                                                      *
c     subroutine for calculating the basic matrices h and hg in case   *
c     of different exponents:     exa .ne. exb.                        *
c     it uses the two-range translation formula ot the stos and some   *
c     relations derived from the one-range expansion.                  *
c     see references (2), (3) and (5)                                  *
c                                                                      *
c***********************************************************************
      subroutine h2r (exa, exb, a, lm0p, nmax, fl, gl, ai, ak, f, h, hg)
      implicit real*8 (a-h,o-z)
      parameter (ldim = 40, ldim2 = 40, ndim = 40, ndim2 = 250)
      dimension h(0:ldim,0:ndim), hg(0:ldim,0:ndim)
      dimension fl(0:ldim2), gl(0:ldim2), ai(0:ldim,0:ndim)
      dimension ak(0:ldim,0:ndim)
      dimension f(0:ndim2)
      data a0/4.819808871d0/,a1/4.069856966d0/,a2/-0.408033291d0/
      data b0/10.278911565d0/,b1/1.069387755d0/,b2/-2.923836201d-3/
      x = exa*a
      y = exb*a
      x2 = x * x
      y2 = y * y
      y2p4 = .25d0 * y2
      b = x + y
      bi= 1.d0 / b
      ypb=y * bi
      ypb2=ypb*ypb
      xmypb = (x - y) * bi
      x2my2 = x2 - y2
      expnx = dexp(-x)
      expny = dexp(-y)
      expnb = expnx * expny
c fl function
      call funfl (expny, lm0p, y, y2, fl)
c gl function
      lm = lm0p + 1
      rm = lm
      call fungl(lm, rm, y, y2p4, gl, ierr)
      if (ierr .eq. 1) then
         write (6,*) 'the series in subr. fungl has not converged'
      endif
c f  function
      if(y.le.5.d0) then
         jmax=a0+y*(a1+y*a2)
      elseif(y .lt. 180.d0) then
         jmax=b0+y*(b1+y*b2)
      else
         jmax = 110
      endif
      n0 = 0
      n1 = 2 * (lm+jmax) + 3
      call fn3 (n1, n0, x, f, ierr)
      if (ierr .eq. 1) then
         write (6,*) 'the series in subr. fn3 has not converged'
      endif
c ak  matrix
      call kpos (lm0p, y2, bi, ypb, expnb, ak)
c ai matrix
      call ipos (lm, lm0p, nmax, x, y2, expnx, x2my2, gl, f, ai, ierr)
      if (ierr .eq. 1) then
         write (6,*) 'the series in subr. ipos  has not converged'
      endif
c h  matrix
      call hpeq (lm0p,nmax,a,fl,gl,ak,ai,h)
c hg  matrix
      call hgran (lm0p, nmax, x, y, a, expny, gl, ai, ak, h, hg)
      return
      end
c***********************************************************************
c     subroutine fungl                                                 *
c                                                                      *
c     subroutine for calculating the auxiliary functions:              *
c        gl(l)   = phi(l)              (lower case phi)                *
c     defined in equation (7) of reference 2                           *
c                                                                      *
c     mayo de 1993  (version modificada por willy)                     *
c     junio de 1993  (version modificada por willy)                    *
c                                                                      *
c***********************************************************************
      subroutine fungl(lm,rm,yy,yp,gl,ierr)
      implicit real*8(a-h,o-z)
c     parameter (ldim2 = 30)
      parameter (ldim2 = 40)
      dimension gl(0:ldim2)
c coefficients for the power series of gl  with  l=25 and  l=26
      data a1/0.377358490566038d-01/,b1/  0.363636363636364d-01/,
     @a2/   0.686106346483705d-03/ ,b2/   0.637958532695375d-03/,
     @a3/   0.802463563138836d-05/ ,b3/   0.720857099090819d-05/,
     @a4/   0.680053867066810d-07/ ,b4/   0.590866474664606d-07/,
     @a5/   0.445936962011023d-09/ ,b5/   0.375153317247369d-09/,
     @a6/   0.235945482545515d-11/ ,b6/   0.192386316537112d-11/,
     @a7/   0.103712300020007d-13/ ,b7/   0.820410731501544d-14/,
     @a8/   0.386986194104502d-16/ ,b8/   0.297250265036791d-16/,
     @a9/   0.124633234816265d-18/ ,b9/   0.930360766938314d-19/,
     @a10/  0.351079534693703d-21/ ,b10/  0.254893360805018d-21/,
     @a11/  0.874419762624417d-24/ ,b11/  0.617923298921255d-24/,
     @a12/  0.194315502805426d-26/ ,b12/  0.133749631801137d-26/,
     @a13/  0.388242762848004d-29/ ,b13/  0.260466663682838d-29/,
     @a14/  0.702066478929483d-32/ ,b14/  0.459376831892131d-32/
c chebyshev approximations for g(25,x) with 10<x<20
      data c1/  0.1000019256738100d+01/,d1/  0.1000000000000000d+01/,
     &c2/ -0.5632114880223869d-01/,d2/ -0.5630187635784173d-01/,
     &c3/  0.7032777640545537d-02/,d3/ -0.2410084263005416d-02/,
     &c4/ -0.2980914600658174d-03/,d4/  0.2355712204692262d-03/,
     &c5/  0.1841087904412635d-04/,d5/ -0.2218578589302566d-05/,
     &c6/ -0.5391960985554504d-06/,d6/ -0.2798947073751880d-06/,
     &c7/  0.1945068785133128d-07/,d7/  0.1141149654664462d-07/,
     &c8/ -0.3064202784569949d-09/,d8/ -0.1788819910312548d-09/,
     &c9/  0.5792725201121395d-11/,d9/  0.1059497283455104d-11/
c chebyshev approximations for g(26,x) with 10<x<20
      data q1/  0.9999117414938084d+00/,r1/  0.1000000000000000d+01/,
     &q2/ -0.8563185974723568d-02/,r2/ -0.8658289107394246d-02/,
     &q3/  0.2508054681791693d-02/,r3/ -0.6534990617504875d-02/,
     &q4/  0.1224213739475198d-03/,r4/  0.1861670875525153d-03/,
     &q5/ -0.5027876758310956d-05/,r5/  0.1777369906431536d-04/,
     &q6/  0.5337472123731052d-06/,r6/ -0.1343219855843926d-05/,
     &q7/ -0.1463812946516166d-07/,r7/  0.3839355677180214d-07/,
     &q8/  0.4330229336365580d-09/,r8/ -0.5273627010607001d-09/,
     &q9/ -0.3554476231970785d-11/,r9/  0.2912170916302903d-11/
c chebyshev approximations for g(25,x) with 20<x<30
      data e1/ 0.1713704300995508d+01/,f1/  0.1000000000000000d+01/,
     &e2/ -0.4369193234204283d+00/,f2/ -0.4069429289956433d-01/,
     &e3/  0.1095016969291855d+00/,f3/  0.4970953844296807d-03/,
     &e4/ -0.1528769750693318d-01/,f4/  0.5933925195289208d-05/,
     &e5/  0.1515352349550369d-02/,f5/ -0.3158014291243555d-05/,
     &e6/ -0.1010159690337459d-03/,f6/  0.1495251600134320d-06/,
     &e7/  0.4733387733217231d-05/,f7/ -0.1385971819539630d-08/,
     &e8/ -0.1440304726779274d-06/,f8/  0.1318206661233724d-10/,
     &e9/  0.2540557815697403d-08/,f9/ -0.3627897780405481d-11/,
     &e10/ -0.5697393935599115d-11/,f10/  0.1364528709797339d-12/,
     &e11/ -0.6406263981556108d-12/,f11/ -0.1956014199498570d-14/,
     &e12/  0.1168725459428224d-13/,f12/  0.1016999158485181d-16/
c chebyshev approximations for g(26,x) with 20<x<30
      data s1/  0.2258844040686752d+00/,t1/  0.1000000000000000d+01/,
     &s2/  0.4340495000232318d+00/,t2/ -0.1547063271983727d-01/,
     &s3/ -0.1043188443389071d+00/,t3/  0.2408033444087682d-02/,
     &s4/  0.1700605606177594d-01/,t4/ -0.1917275089806831d-03/,
     &s5/ -0.1582804650575731d-02/,t5/  0.3916283431166363d-05/,
     &s6/  0.9850305867366136d-04/,t6/ -0.2800529075030164d-07/,
     &s7/ -0.3463444504775444d-05/,t7/  0.5324925868849985d-09/,
     &s8/  0.4795054211435089d-07/,t8/ -0.2532808035362497d-10/,
     &s9/  0.1411019717896250d-08/,t9/  0.1937562077071697d-11/,
     &s10/ -0.2259079293383414d-10/,t10/ -0.3014368727017189d-13/,
     &s11/ -0.3486968887693469d-11/,t11/ -0.2136312805855105d-14/,
     &s12/  0.2118929018907108d-12/,t12/  0.8783411172328588d-16/,
     &s13/ -0.4943913216654270d-14/,t13/ -0.1206194411603270d-17/,
     &s14/  0.5315234221021871d-16/,t14/  0.5803857645788197d-20/
c chebyshev approximations for g(25,x) with 30<x<50
      data g1/  0.2136784974859081d+27/,h1/  0.1000000000000000d+01/,
     &g2/ -0.5307537207822698d+26/,h2/ -0.4597687465820516d-01/,
     &g3/  0.3991958284862451d+25/,h3/  0.2777228416870579d-02/,
     &g4/  0.1108803616436256d+24/,h4/ -0.5953181072505075d-04/,
     &g5/ -0.3711221382420830d+23/,h5/  0.1683449752266743d-05/,
     &g6/  0.2454556302050156d+22/,h6/ -0.1661569902915213d-07/,
     &g7/ -0.6599012875510564d+20/,h7/  0.2804061533093740d-09/,
     &g8/  0.1844239669280550d+18/,h8/ -0.8725730161861490d-12/,
     &g9/  0.1787345641585511d+17/,h9/  0.1030105909414584d-13/
c chebyshev approximations for g(26,x) with 30<x<50
      data u1/ -0.9227788016377853d+28/,v1/  0.1000000000000000d+01/,
     &u2/  0.3899095143731592d+28/,v2/ -0.4490924683341691d-01/,
     &u3/ -0.7205768107226046d+27/,v3/  0.2608880306113027d-02/,
     &u4/  0.7541005203892215d+26/,v4/ -0.5416251351112540d-04/,
     &u5/ -0.4795501295532321d+25/,v5/  0.1465695683459566d-05/,
     &u6/  0.1811824968593083d+24/,v6/ -0.1389119812874038d-07/,
     &u7/ -0.3420863788902199d+22/,v7/  0.2237773180622521d-09/,
     &u8/  0.5075988131226513d+19/,v8/ -0.6574103042963774d-12/,
     &u9/  0.6454603661090363d+18/,v9/  0.7455552778285355d-14/
c chebyshev approximations for g(25,x) with 50<x<70
      data o1/ -0.6141811427602216d+28/,p1/  0.1000000000000000d+01/,
     &o2/  0.1751360571301923d+28/,p2/ -0.1699296252380966d-01/,
     &o3/ -0.1887369638589229d+27/,p3/  0.8842430633309536d-03/,
     &o4/  0.1005551236891022d+26/,p4/ -0.4331172810245365d-05/,
     &o5/ -0.2691211987905570d+24/,p5/  0.1563370383720422d-06/,
     &o6/  0.2811307601877234d+22/,p6/  0.4895038102689216d-10/,
     &o7/  0.4871744832582661d+19/,p7/  0.4970536534837672d-11/
c chebyshev approximations for g(26,x) with 50<x<70
      data w1/ -0.7075147701304132d+30/,z1/  0.1000000000000000d+01/,
     &w2/  0.1437432458976818d+30/,z2/ -0.1635529244186009d-01/,
     &w3/ -0.1214456312182020d+29/,z3/  0.8305650334772285d-03/,
     &w4/  0.5333574224217863d+27/,z4/ -0.4140710290149234d-05/,
     &w5/ -0.1208488093637420d+26/,z5/  0.1380808827439971d-06/,
     &w6/  0.1067974314983611d+24/,z6/  0.7170201185335186d-12/,
     &w7/  0.2231586612184770d+21/,z7/  0.4115534362275280d-11/
c  general formula for large lm
      y = dabs(yy)
      ierr = 0
c     print*,'fungl. lm, y ',lm,y
      if(lm.gt.25.or.y.gt.7.d1) goto 2000
         if(y.gt.10.d0) goto 100
           if(y.gt.2.d0) goto 30
             if(y.gt.2.d-1) goto 20
       g25=1.d0+yp*(a1+yp*(a2+yp*(a3+a4*yp)))
       g26=1.d0+yp*(b1+yp*(b2+yp*(b3+b4*yp)))
       goto 1000
   20  g25=1.d0+yp*(a1+yp*(a2+yp*(a3+yp*(a4+yp*(a5+yp*(a6+yp*a7))))))
       g26=1.d0+yp*(b1+yp*(b2+yp*(b3+yp*(b4+yp*(b5+yp*(b6+yp*b7))))))
       goto 1000
   30        if(y.gt.5.d0) goto 50
       g25=1.d0+yp*(a1+yp*(a2+yp*(a3+yp*(a4+yp*(a5+yp*(a6+yp*
     &    (a7+yp*(a8+yp*(a9+yp*a10)))))))))
       g26=1.d0+yp*(b1+yp*(b2+yp*(b3+yp*(b4+yp*(b5+yp*(b6+yp*
     &    (b7+yp*(b8+yp*(b9+yp*b10)))))))))
       goto 1000
   50  g25=1.d0+yp*(a1+yp*(a2+yp*(a3+yp*(a4+yp*(a5+yp*(a6+yp*
     &    (a7+yp*(a8+yp*(a9+yp*(a10+yp*(a11+yp*(a12+yp*(a13+yp*a14
     &    )))))))))))))
       g26=1.d0+yp*(b1+yp*(b2+yp*(b3+yp*(b4+yp*(b5+yp*(b6+yp*
     &    (b7+yp*(b8+yp*(b9+yp*(b10+yp*(b11+yp*(b12+yp*(b13+yp*b14
     &    )))))))))))))
       goto 1000
  100      if(y.gt.3.d1) goto 150
             if(y.gt.2.d1) goto 120
c--------------------
c chebyshev approximations 10-20
c--------------------
       pn=c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*
     &    (c8+y*c9)))))))
       pd=d1+y*(d2+y*(d3+y*(d4+y*(d5+y*(d6+y*(d7+y*
     &    (d8+y*d9)))))))
       g25=pn/pd
       pn=q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*
     &    (q8+y*q9)))))))
       pd=r1+y*(r2+y*(r3+y*(r4+y*(r5+y*(r6+y*(r7+y*
     &    (r8+y*r9)))))))
       g26=pn/pd
       goto 1000
c--------------------
c chebyshev approximations 20-30
c--------------------
  120  pn=e1+y*(e2+y*(e3+y*(e4+y*(e5+y*(e6+y*(e7+y*
     &    (e8+y*(e9+y*(e10+y*(e11+y*e12))))))))))
       pd=f1+y*(f2+y*(f3+y*(f4+y*(f5+y*(f6+y*(f7+y*
     &    (f8+y*(f9+y*(f10+y*(f11+y*f12))))))))))
       g25=pn/pd
       pn=s1+y*(s2+y*(s3+y*(s4+y*(s5+y*(s6+y*(s7+y*
     &    (s8+y*(s9+y*(s10+y*(s11+y*(s12+y*(s13+y*s14))))))))))))
       pd=t1+y*(t2+y*(t3+y*(t4+y*(t5+y*(t6+y*(t7+y*
     &    (t8+y*(t9+y*(t10+y*(t11+y*(t12+y*(t13+y*t14))))))))))))
       g26=pn/pd
       goto 1000
  150    if(y.gt.5.d1) goto 170
c--------------------
c chebyshev approximations 30-50
c--------------------
       pn=g1+y*(g2+y*(g3+y*(g4+y*(g5+y*(g6+y*(g7+y*
     &    (g8+y*g9)))))))
       pd=h1+y*(h2+y*(h3+y*(h4+y*(h5+y*(h6+y*(h7+y*
     &    (h8+y*h9)))))))
       g25=pn/pd
       pn=u1+y*(u2+y*(u3+y*(u4+y*(u5+y*(u6+y*(u7+y*
     &    (u8+y*u9)))))))
       pd=v1+y*(v2+y*(v3+y*(v4+y*(v5+y*(v6+y*(v7+y*
     &    (v8+y*v9)))))))
       g26=pn/pd
       goto 999
c--------------------
c chebyshev approximations 50-70
c--------------------
  170  pn=o1+y*(o2+y*(o3+y*(o4+y*(o5+y*(o6+y*o7)))))
       pd=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
       g25=pn/pd
       pn=w1+y*(w2+y*(w3+y*(w4+y*(w5+y*(w6+y*w7)))))
       pd=z1+y*(z2+y*(z3+y*(z4+y*(z5+y*(z6+y*z7)))))
       g26=pn/pd
  999 fy=dexp(y)/y**25
      g25=g25*fy
      g26=g26*fy/y
 1000 continue
c llega aqui despues de calcular g25 , g26 con los pades
      gl(26)=g26
      gl(25)=g25
      lrec = 24
      go to 9000
 2000 aux1=1.d0
      aux2=1.d0
      sum1=1.d0
      sum2=1.d0
      xl1=rm+.5d0
      rj=0.d0
      ierr = 0
c
c  aumento el numero maximo de terminos en la serie
c
      do 2300 j=1,1000
        rj=rj+1.d0
        aj=(xl1+rj)*rj
        aux1=aux1*yp/aj
        if(aux1.lt.1.d-15*sum1) goto 2400
        aux2=aux2*yp/(aj+rj)
        sum1=sum1+aux1
        sum2=sum2+aux2
 2300 continue
c     the series has not converged
c     print*,'err rel ',aux1/sum1
      ierr = 1
 2400 continue
      gl(lm)=sum1
      gl(lm+1)=sum2
      lrec = lm-1
 9000 continue
      do 2500 l = lrec , 0 , -1
        gl(l)=gl(l+1) + yp * gl(l+2)/( (l+2.5)*(l+1.5) )
 2500 continue
      return
      end
