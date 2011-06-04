*     Version de la subroutina  coul0  adaptada para el calculo
*     aprovechando la simetria local. Difiere de la original en
*     que no se rotan las integrales al sistema de ejes molecular
*
*     Version de Rafa de septiembre de 1999
*
*
*=============================================================      
      subroutine coul0sim(nn1,ll1,xx1,nn2,ll2,xx2,nn3,ll3,xx3
     &   ,nn4,ll4,xx4,rn1,rn2,rn3,rn4,a, coul)
*=============================================================      
      implicit real * 8 (a-h,o-z)
      include 'comun.inc'
      include 'm2cdat.inc'
      dimension coul(-mxl:mxl,-mxl:mxl,-mxl:mxl,-mxl:mxl)
c
c calculo las integrales en el sistema de ejes z alineados.
c

      call  coulom(nn1,ll1,xx1,nn2,ll2,xx2,nn3,ll3,xx3,
     &             nn4,ll4,xx4,a,coul)
c
c normalizo todo
c
      rn = rn1*rn2*rn3*rn4
      do 101 m1 = -ll1,ll1
      ang1 = ang(ll1,iabs(m1))
      do 101 m2 = -ll2,ll2
      ang2 = ang(ll2,iabs(m2))
      do 101 m3 = -ll3,ll3
      ang3 = ang(ll3,iabs(m3))
      do 101 m4 = -ll4,ll4
         coul(m1,m2,m3,m4) = coul(m1,m2,m3,m4) * ang1*ang2*ang3*
     &                      ang(ll4,iabs(m4)) * rn

  101 continue

      return
      end
c
c 
c
*================================================================      
      subroutine coulom(nn1,ll1,xx1,nn2,ll2,xx2,nn3,ll3,xx3,
     &                  nn4,ll4,xx4,a,coul)
*================================================================      
c
c subrutina para calcular integrales coulombianas bicentricas.
c -------------------- version mas reciente ---------------------------
c lee los datos en la forma (aa,bb). devuelve los datos en coul
c
      implicit real * 8 (a-h,o-z)
      logical test1 , test2 , pab
      include 'comun.inc'
      dimension cj(0:2*mxl,0:2*mxl,0:2*mxl)
      dimension coul(-mxl:mxl,-mxl:mxl,-mxl:mxl,-mxl:mxl)
c
c utilizo el siguiente truco sucio para poder pasar alm y blm como
c argumentos a la subrutina mdes.
c
      dimension alm(0:mxkcof,0:2*mxl+1), blm(0:mxkcof,0:2*mxl+1)
      equivalence ( alm(0,0) , app(0,0) )
      equivalence ( blm(0,0) , bpp(0,0) )

*     escribe en coulom

      pab = xx3+xx4 .gt. xx1+xx2
      if ( pab ) then
         fase = (-1)**(ll1+ll2+ll3+ll4)
         n1 = nn3
         n2 = nn4
         n3 = nn1
         n4 = nn2
         l1 = ll3
         l2 = ll4
         l3 = ll1
         l4 = ll2
         x1 = xx3
         x2 = xx4
         x3 = xx1
         x4 = xx2
      else
         n1 = nn1
         n2 = nn2
         n3 = nn3
         n4 = nn4
         l1 = ll1
         l2 = ll2
         l3 = ll3
         l4 = ll4
         x1 = xx1
         x2 = xx2
         x3 = xx3
         x4 = xx4
      endif
      lamin = iabs(l1-l2)
      lamax = l1+l2
      lbmin = iabs(l3-l4)
      lbmax = l3+l4
c
c  calculo las integrales coulombianas basicas
c
      call jbas (n1,l1,n2,l2,n3,l3,n4,l4,x1,x2,x3,x4,a,cj)
c
c  calculo las integrales coulombianas finales
c
      l1l1 = ind(l1)
      l2l2 = ind(l2)
      l3l3 = ind(l3)
      l4l4 = ind(l4)
      test1 = l1 .eq. l2
      test2 = l3 .eq. l4

      do 1 m1 = -l1 , l1
      do 1 m2 = -l2 , l2
        call emes ( m1, m2, ms, md, ss, sd )
        call indab( l1 , m1 , l2 , m2 , l1l1 , l2l2 , k12 )
        do 1 m3 = -l3 , l3
        do 1 m4 = -l4 , l4
          call emes ( m3, m4, msb, mdb, ssb, sdb )
          call indab( l3 , m3 , l4 , m4 , l3l3 , l4l4 , k34 )

          aux = 0.d0
          if (ms.eq.msb) then
           msa = iabs(ms)
           call mdes(k12,k34,lamin,lamax,lbmin,lbmax,msa,cj,alm,alm,bux)
           aux = aux + bux * ss * ssb
          endif
          if (ms.eq.mdb.and.sdb.ne.0.d0) then
           msa = iabs(ms)
           call mdes(k12,k34,lamin,lamax,lbmin,lbmax,msa,cj,alm,blm,bux)
           aux = aux + bux * ss * sdb
          endif
          if (md.eq.msb.and.sd.ne.0.d0) then
           mda = iabs(md)
           call mdes(k12,k34,lamin,lamax,lbmin,lbmax,mda,cj,blm,alm,bux)
           aux = aux + bux * sd * ssb
          endif
          if (md.eq.mdb.and.sd*sdb.ne.0.d0) then
           mda = iabs(md)
           call mdes(k12,k34,lamin,lamax,lbmin,lbmax,mda,cj,blm,blm,bux)
           aux = aux + bux * sd * sdb
          endif
          coul(m1,m2,m3,m4) = aux

    1 continue
c
c se deshace la permutacion si la hubo
c
      if ( pab ) then
         lsup = max( ll1 , ll2 , ll3 , ll4 )
         do 50 m2 = -lsup , lsup
         do 50 m4 = -lsup , lsup
         do 50 m1 = -lsup , lsup
         do 50 m3 = -lsup , m1
            aux = coul(m1,m2,m3,m4)
            coul(m1,m2,m3,m4) = coul(m3,m2,m1,m4)
            coul(m3,m2,m1,m4) = aux
   50    continue
         do 51 m1 = -lsup , lsup
         do 51 m3 = -lsup , lsup
         do 51 m2 = -lsup , lsup
         do 51 m4 = -lsup , m2
            aux = coul(m1,m2,m3,m4)
            coul(m1,m2,m3,m4) = coul(m1,m4,m3,m2) * fase
            coul(m1,m4,m3,m2) = aux * fase
   51    continue
      endif

c     integrales finales sin normalizar : coul(m1,m2,m3,m4)

      return
      end
c
c 
c
*======================================================================      
      subroutine mdes (ka,kb,lamin,lamax,lbmin,lbmax,m,cj,co1,co2,bux)
*======================================================================      
      implicit real * 8 (a-h,o-z)
      include 'comun.inc'
      dimension co1(0:mxkcof,0:2*mxl+1), co2(0:mxkcof,0:2*mxl+1)
     &         ,cj(0:2*mxl,0:2*mxl,0:2*mxl)

      bux = 0.d0

      ma = iabs(m)
      do 1 l = lamin , lamax , 2
      do 1 lp= lbmin , lbmax , 2
         if (l .ge. ma .and. lp .ge. ma) then
            bux = bux + co1(ka,l) * co2(kb,lp)  * cj(l,lp,m)
         endif
    1 continue

      return
      end
c
c
c
*================================================================
      subroutine jbas
     & (n1,l1,n2,l2,n3,l3,n4,l4,x1,x2,x3,x4,a,cj)
*================================================================
c
c  subrutina para el calculo de las integrales coulombianas basicas
c
      implicit real * 8 ( a-h,o-z )
      include 'comun.inc'
      parameter (ldim = 40, ldim2 = 40, ndim = 40, ndim2 = 500)
      dimension g(0:ldim,0:ndim), g1(0:ldim,0:ndim)
     &         ,g2(0:ldim,0:ndim), g3(0:ldim,0:ndim)
     &        , f(0:ldim,0:ndim), f1(0:ldim,0:ndim)
     &        ,f2(0:ldim,0:ndim), f3(0:ldim,0:ndim)
     &        , h(0:ldim,0:ndim), h1(0:ldim,0:ndim)
     &        ,h2(0:ldim,0:ndim), h3(0:ldim,0:ndim)
     &        , cj(0:2*mxl,0:2*mxl,0:2*mxl) , xaf(0:ldim) , xafi(0:ldim)
      dimension fl(0:ldim2),gl(0:ldim2)
     &         ,ai(0:ldim,0:ndim),ak(0:ldim,0:ndim),fn2(0:ndim2)

*     'escribe en jbas'

      do 4321 n = 0 , ndim
      do 4321 l = 0 , ldim
        h (l,n) = 0.d0
 4321   h2(l,n) = 0.d0
      na = n1 + n2 - 1
      nb = n3 + n4 - 1
      lamax = l1 + l2
      lamin = abs ( l1 - l2 )
      lbmax = l3 + l4
      lbmin = abs ( l3 - l4 )
c
c  lbpar es 0 o 1 segun sea lbmax par o impar
c
      lbpar =  lbmax - 2*(lbmax/2)
      xa = x1 + x2
      xb = x3 + x4
      tope = 40.d0
      if ( xb*a .gt. tope ) then
         call jmult ( lamax, lbmax, lamin, lbmin , na , nb ,
     &                xa , xb , a , cj )
         return
      endif
      xai = 1.d0 / xa
c
c   vectores auxiliares
c
      xaf(0) = 1.d0
      xafi(0) = 1.d0
      do 9 i = 1 , na + lamax + 1
         xaf(i) = xaf(i-1) * xa * ri(i)
         xafi(i) = 1.d0 / xaf(i)
    9 continue
c
c   calcula las matrices h, f, y g
c
      lnb = ( nb + lbmax ) / 2
      lg = lamax
      ng = lnb
      lf = lamax
      nf = lnb+1
      lh = nb + max (lamax  + lnb + 1 , na )
      nh = lh
      naux = ( nb - lbmax ) / 2
      if ( 2*naux.eq.(nb-lbmax) ) then
          call h2r  ( xa, xb, a, lh, nh, fl,gl,ai,ak,fn2, h, h2)
          call gmat ( xa, xb, a, lg, ng, h, h2, g, g2)
          call fmat ( xb, a, lf, nf, f, f2)
      else
          call h2r  ( xa, xb, a, lh, nh, fl,gl,ai,ak,fn2, h2, h)
          call gmat ( xa, xb, a, lg, ng, h2, h, g2, g)
          call fmat ( xb, a, lf, nf, f2, f)
      endif
c     call matout ( h , lh , nh )
      do 20 m = 0, min(lamax,lbmax)
      do 20 lb= 0, lbmax
      do 20 la= 0, lamax
         cj(la,lb,m) = 0.d0
   20 continue
c
c   ahora se recurre sobre las matrices f, g y h para aumentar los
c   indices. se van acumulando los resultados en la matriz cj con la
c   subrutina jota
c
 1111 ninf = 0
      nmin = 0
      mpar = 1
      do 1 m = 0, min ( lamax, lbmax )
         lpar = mpar
         mpar = 1 - mpar
         if ( m.gt.0 ) then
            m1 = m - 1
            call reqmg ( lg, ng, m1, xai, h, g )
            call reqmf ( lf, nf, m1, f )
            call reqmb (      lh, nmin, nh, m1, m1, h )
c           call reqmb (ninf, lh, nmin, nh, m1, m1, h )
         endif
         lh1 = lh
         nh1 = nh
         ng1 = ng
         nf1 = nf
         call dup ( lg+ng , ng, g1, g )
         call dup ( lf+nf , nf, f1, f )
         call dup ( lh, nh, h1, h )
         do 2 l = m, lbmax
            lpar = 1 - lpar
            if ( l.gt.m ) then
               lb1 = l - 1
               call reqlg ( lg, ng1, m, lb1, xai, a, h1, h2, g1, g2 )
               call reqlf ( lf, nf1, m, lb1, a, f1, f2 )
               call reqlb (      lh1, nmin, nh1, m, lb1, m, a, h1, h2)
c              call reqlb (ninf, lh1, nmin, nh1, m, lb1, m, a, h1, h2)
            endif
            if ( (lpar .eq. lbpar) .and. (l.ge.lbmin) ) then
               lh3 = lh1
               nh3 = nh1
               ng3 = ng1
               nf3 = nf1
               call dup ( lg+ng1 , ng1, g3, g1)
               call dup ( lf+nf1 , nf1, f3, f1)
               call dup ( lh1, nh1, h3, h1)
               do 3 n = l+2, nb, 2
                  call reqng ( lg, ng3, m, xai, a, h3, g3 )
                  call reqnf ( lf, nf3, m, a, f3 )
                  call reqnb ( lh3, nmin, nh3, m, a, h3)
c                 call reqnb (ninf, lh3, nmin, nh3, m, a, h3)
    3          continue
               call jota ( na, l, m, lamin, lamax, xai, xaf, xafi
     &                   , f3, g3, h3, cj )
            endif
    2    continue
    1 continue
 7777 continue

      return
      end
c
c
c
*===============================================================
      subroutine reqmg ( lmax, nmax, m, xai, h, g )
*===============================================================
      implicit real * 8 ( a-h,o-z )
      include 'comun.inc'
      parameter (ldim = 40, ldim2 = 40, ndim = 40, ndim2 = 500)
      dimension h(0:ldim,0:ndim), g(0:ldim,0:ndim)
      nmax = nmax - 1
      aux = - dosl1(m)
      m1 = m + 1
      do 1 n = 0, nmax
      n1 = n + 1
      do 1 l = lmax + n, max(m1,n), -1
         bux = rll1( l-m ) * g(l+1,n1)
         l1 = l - 1
         cux = h(l1,2)
         if ( l.gt.n ) then
             p = 2.d0 * real(l-n) * xai
             q = p - xai
             cux = cux + p * ( h(l1,1) + q * g(l1,n) )
         endif
         g(l,n) =( bux - rll1( l+m ) * cux ) * aux * dosl1p(l)
    1  continue
       return
       end
c
c 
c
*==============================================================
      subroutine reqng ( lmax, nmax, m, xai, a, h, g )
*==============================================================
      implicit real * 8 ( a-h,o-z )
      logical test
      include 'comun.inc'
      parameter (ldim = 40, ldim2 = 40, ndim = 40, ndim2 = 500)
      dimension h(0:ldim,0:ndim), g(0:ldim,0:ndim)
      nmax = nmax - 1
      a2 = a * a
      aa = a + a
      do 1 n = 0, nmax
      n1 = n + 1
      do 1 l = lmax+n, max(m,n), -1
         l1 = l - 1
         test = ( l.gt.m )
         if ( test ) then
             aux = h(l1,2)
         else
             aux = 0.d0
         endif
         bux = h(l,2)
         if ( l.gt.n ) then
             p = 2.d0 * real(l-n) * xai
             q = p - xai
             bux = bux + p * ( h(l,1) + q * g(l,n1) )
             if ( test ) aux = aux + p * ( h(l1,1) + q * g(l1,n) )
         endif
         aux = real( l+m ) * aux
         cux = ( aux + real(l-m+1)*g(l+1,n1) ) * aa * dosl1p(l)
         g(l,n) = a2 * g(l,n) -  cux + bux
    1  continue
       return
       end
c
c
c
*==================================================================
      subroutine reqlg ( lmax, nmax, m, lb, xai, a, h1, h2, g1, g2)
*==================================================================
      implicit real * 8 ( a-h,o-z )
      logical test
      include 'comun.inc'
      parameter (ldim = 40, ldim2 = 40, ndim = 40, ndim2 = 500)
      dimension h1(0:ldim,0:ndim), h2(0:ldim,0:ndim)
     &         ,g1(0:ldim,0:ndim), g2(0:ldim,0:ndim)
      test = ( lb.gt.m )
      if ( test ) then
         nmax2 = nmax + 1
         call reqng ( lmax, nmax2, m, xai, a, h2, g2 )
      endif
      nmax = nmax - 1
      auxl =  dosl1(lb)
      buxl =  real(lb+m)
      cuxl =  1.d0 / real( lb-m+1 )
      do 1 n = 0, nmax
      n1 = n + 1
      do 1 l = lmax+n, max(m,n), -1
         l1 = l - 1
         if ( l.gt.m ) then
             aux = h1(l1,2)
             if ( l.gt.n ) then
                p = 2.d0 * real( l-n) * xai
                q = p - xai
                aux = aux + p * ( h1(l1,1) + q * g1(l1,n) )
             endif
         else
             aux = 0.d0
         endif
         aux = ( real(l-m+1)*g1(l+1,n1) + real(l+m)*aux ) * dosl1p(l)
         aux = ( aux - a * g1(l,n) ) * auxl
         if ( test ) aux = aux - buxl * g2(l,n)
         g2(l,n) = g1(l,n)
         g1(l,n) = aux * cuxl
    1 continue
      n1 = nmax + 1
      do 3 l = lmax+n1, n1, -1
         g2(l,n1) = g1(l,n1)
    3 continue
      return
      end
c
c 
c
*===============================================================
       subroutine gmat ( xa, xb, a, lmax, nmax, h, hg, g, g2)
*===============================================================

************************************************************************
*                                                                      *
*  subrutina para calcular las matrices basicas g y g2 definidas asi:  *
*                                                                      *
*    g(l,n) = (2l-2n)³ / xa**(2l-2n)  *                                *
*           suma de j=0 a 2l-2n  de ( h(l,j-2l+2n) * xa**j / j³ )      *
*                                                                      *
*   g2(l,n) = (2l-2n)³ / xa**(2l-2n)  *                                *
*           suma de j=0 a 2l-2n  de (hg(l,j-2l+2n) * xa**j / j³ )      *
*                                                                      *
************************************************************************

      implicit real * 8 ( a-h,o-z )
      include 'comun.inc'
      parameter (ldim = 40, ldim2 = 40, ndim = 40, ndim2 = 500)
      dimension g(0:ldim,0:ndim), g2(0:ldim,0:ndim)
     &         ,h(0:ldim,0:ndim), hg(0:ldim,0:ndim)
     &        , en(0:2*ldim)
      data pi/0.31415926535897932384626433832795059d+01/
      ln = lmax + nmax
      x = xa * a
      ex = dexp(-x)
      xi = 1.d0 / x
      en(0) = 4.d0 * pi * ex * xi * xi
      do 10 i = 1, 2*ln
         en(i) = en(0) + real(i) * xi * en(i-1)
   10 continue
      xa2 = 1.d0 / ( xa * xa )
      ba =  xb * xb * xa2
      ba2= ( xb + xb ) * xa2
c
c   diagonales  g(l,l), g2(l,l)
c
      do 1 l = 0 , ln
         g(l,l) = h(l,0)
         g2(l,l) = hg(l,0)
    1 continue
c
c   resto de g y g2
c
      al = a
      do 2 l = 1, ln
      l1 = l + 1
      ln2 = 0
      al = al * a
      do 2 n = l, max(1,l-lmax), -1
            rn = 1.d0 / real(n)
            n1 = n - 1
            bux = real(l1-n)
            aux = bux * ( en(ln2)* al - ba * g(l,n) )
            aux = aux + real(l1)*h(l,0)
            g(l,n1) = aux * rn
            aux = bux * ( ba2 * g(l,n) - ba * g2(l,n) )
            aux = aux + real(l1)*hg(l,0)
            g2(l,n1) = aux * rn
            ln2 = ln2 + 2
    2    continue
      return
      end
c
c 
c
*===============================================================
      subroutine jota ( na, lb, m, lamin, lamax, xai, xaf, xafi
     &                , f, g, h, cj )
*===============================================================

************************************************************************
*                                                                      *
*  subrutina que a¥ade las integrales de solapamiento g y h a la matriz*
*  j nlm n'l'm para todos los l, l' y m.                               *
*                                                                      *
************************************************************************

      implicit real * 8 ( a-h,o-z )
      logical test
      include 'comun.inc'
      parameter (ldim = 40, ldim2 = 40, ndim = 40, ndim2 = 500)
      dimension cj(0:2*mxl,0:2*mxl,0:2*mxl)
     &        , h(0:ldim,0:ndim), g(0:ldim,0:ndim)
     &        , f(0:ldim,0:ndim) , xaf(0:ldim) , xafi(0:ldim)
      data pi4/ 0.125663706143592d+02 /
      test = ( m.gt.0 )
c
c  contribucion de las sumas de integrales de solapamiento con
c  indices negativos ( matriz g ) y de la matriz f.
c
      do 1 l = lamin , lamax, 2
         if ( l.ge.m ) then
            nl1 = na + l + 1
            aux = ( f(l,0) - xaf(l+l) * g(l,0) ) * xafi(nl1)
            nl = na - l
            sum = 0.d0
            t = fact(nl1) / fact(nl)
            ll1 = l + l + 1
            do 3 i = 0, nl-1
               b = 1.d0 - fact(i) * t / fact(ll1+i)
               sum = sum + b * xaf(i) * h(l,i+1)
    3       continue
            sum = sum * xafi(nl)
 1234 format(3i4,3d15.5)
            aux = aux + sum
            if ( test ) aux = aux * 0.5d0
            cj(l,lb,m)  = aux * pi4 * dosl1p(l) * xai
         endif
    1 continue
      return
      end
c
c 
c
*==============================================
      subroutine reqmf ( lmax, nmax, m, f )
*==============================================
      implicit real * 8 ( a-h,o-z )
      parameter (ldim = 40, ndim = 40)
      include 'comun.inc'
      dimension f(0:ldim,0:ndim)
      nmax = nmax - 1
      aux = - dosl1(m)
      m1 = m + 1
      do 1 n = 0, nmax
      n1 = n + 1
      do 1 l = lmax + n, m1, -1
         bux = rll1( l-m ) * f(l+1,n1)
         cux = rll1( l+m ) * f(l-1,n)
         f(l,n) = aux * ( bux - cux ) * dosl1p(l)
    1 continue
      return
      end
c
c
c
*==============================================
      subroutine reqnf ( lmax, nmax, m, a, f )
*==============================================
      implicit real * 8 ( a-h,o-z )
      parameter (ldim = 40, ndim = 40)
      include 'comun.inc'
      dimension f(0:ldim,0:ndim)
      nmax = nmax - 1
      a2 = a * a
      aa = a + a
      do 1 n = 0, nmax
      n1 = n + 1
      do 1 l = lmax+n, m, -1
         aux = real(l-m+1) * f(l+1,n1)
         if ( l.gt.m ) aux = aux + real(l+m) * f(l-1,n)
         f(l,n) = a2 * f(l,n) + f(l,n1) - aa * aux * dosl1p(l)
    1 continue
      return
      end
c
c
c
*======================================================
      subroutine reqlf ( lmax, nmax, m, lb, a, f1, f2 )
*======================================================
      implicit real * 8 ( a-h,o-z )
      parameter (ldim = 40, ndim = 40)
      logical test
      include 'comun.inc'
      dimension f1(0:ldim,0:ndim), f2(0:ldim,0:ndim)
      test = ( lb.gt.m )
      if ( test ) then
         nmax2 = nmax + 1
         call reqnf ( lmax, nmax2, m, a, f2 )
      endif
      nmax = nmax - 1
      auxl = dosl1(lb)
      buxl = real( lb+m )
      cuxl =  1.d0 / real( lb-m+1 )
      do 1 n = 0, nmax
      n1 = n + 1
      do 1 l = lmax+n, m, -1
         aux = real(l-m+1) * f1(l+1,n1)
         if ( l.gt.m ) aux = aux + real(l+m) * f1(l-1,n)
         aux = auxl * ( aux * dosl1p(l) - a * f1(l,n) )
         if ( test ) aux = aux - buxl * f2(l,n)
         f2(l,n) = f1(l,n)
         f1(l,n) = aux * cuxl
    1 continue
      n1 = nmax + 1
      do 3 l = lmax+n1, m, -1
         f2(l,n1) = f1(l,n1)
    3 continue
      return
      end
c
c
c
*======================================================
      subroutine fmat ( xb, a, lmax, nmax, f1, f2)
*======================================================

***********************************************************************
*    subrutina para calcular las matrices basicas f1 y f2 :           *
*                                                                     *
*    f1(l,n) = hp (l, 2n-2l) ( exa = 0, exb )                         *
*    f2(l,n) = hg (l, 2n-2l) ( exa = 0, exb )                         *
***********************************************************************
      implicit real * 8 ( a-h,o-z )
      include 'comun.inc'
      parameter (ldim = 40, ldim2 = 40, ndim = 40, ndim2 = 500)
      dimension f1(0:ldim,0:ndim), f2(0:ldim,0:ndim), cl(0:ldim,0:ndim)
     &        , dl(0:ldim,0:ndim), fl(0:ldim2), gl(0:ldim2)
      data pi/0.31415926535897932384626433832795059d+01/
      y = xb * a
      y2 = y * y
      y2p4 = y2 * 0.25d0
      yi = 1.d0 / y
      y2i = yi * yi
      ll = lmax + nmax
c
c  funciones fl(l)
c
      fl(0) = 1.d0
      fl(1) = 1.d0 + y
      do 11 l = 1, ll
         l1 = l - 1
         fl(l+1) = fl(l) + y2 * fl(l1) * dosl1p(l1) * dosl1p(l)
   11 continue
c
c   funciones gl
c
      lm = ll + 1
      rm = lm
      call fungl(lm,rm,y,y2p4,gl,ierr)
      if (ierr .eq. 1) then
         write (6,*) 'fmat. the series in subr. fungl has not converged'
      endif
c     call fungl(lm,rm,y,y2p4,gl)
c
c   matriz dl
c
      do 1 l = 0, lmax
         if ( l.eq.0 ) then
             dl(0,0) = yi
         else
             l1 = l - 1
             dl(l,0) = fl(l1) * dosl1p(l1)
         endif
         do 1 n = 0, nmax-1
            ln = l + n
            dl(ln+1,n+1) = ( fl(ln) + 2.d0 * real(n+1) * dl(ln,n) )
     &                   * dosl1p(ln)
    1 continue
      do 12 n = 1, nmax
         nn = n + n
         dl(0,n) = yi * ( 1.d0 + real(nn) * yi *
     &                  ( 1.d0 + real(nn-1) * dl(0,n-1) ) )
         do 12 l = 1, nmax-n
            ln = l + n
            l1 = l - 1
            dl(l,ln)=( fl(l1) + real(ln+ln)*dl(l1,ln-1) )*dosl1p(l1)
   12 continue
c
c   matriz cl
c
c   calculo de la diagonal superior
c
      if ( y.ge.15.d0 ) then
         l1 = lmax
         n1 = 0
         call clln ( l1, n1, y2p4, sum )
         cl(l1,n1) = sum
         do 2 n = 1, nmax
            aux = gl(l1) - cl(l1,n-1) * real(n+n)
            l1 = l1 + 1
            cl(l1,n) = dosl1(l1) * aux * y2i
    2    continue
      else
         l1 = lmax + nmax
         n1 = nmax
         call clln ( l1, n1, y2p4, sum )
         cl(l1,n1) = sum
         do 3 n = nmax-1, 0, -1
            n1 = n + 1
            l2 = l1
            l1 = l1 - 1
            aux = gl(l1) - y2 * cl(l2,n1) * dosl1p(l2)
            cl(l1,n) = aux * 0.5d0 / real(n1)
    3    continue
      endif
c
c  el resto, por recurrencia
c
      do 44 n = 0, nmax
      aux = dosl1(lmax)
      do 44 l = lmax+n, 1, -1
         aux = aux - 2.d0
         cl(l-1,n) = ( gl(l) + aux * cl(l,n) ) * dosl1p(l)
   44 continue
c
c  matriz f1
c
      py = 4.d0 * pi * dexp(-y)
      an = a
      a2 = a * a
      ai = 1.d0 / a
      do 4 n = 0, nmax
         pyn = py * an
         an = an * a2
         al = 1.d0
         do 4 l = 0, lmax+n
            f1(l,n) = pyn * ( fl(l)*cl(l,n) + gl(l)*dl(l,n) )
     &               * dosl1p(l) * al
            al =al * ai
    4 continue
c
c   matriz f2
c
      xb2 = 1.d0 / ( xb * xb )
      ay = 2.d0 * a * yi
      f2(0,0) = ay * f1(0,0) - py * xb2
      do 50 n = 1, nmax
         nn = n + n - 1
         f2(0,n) = ay*f1(0,n) + rll1(nn) * f2(0,n-1) * xb2
   50 continue
      do 5 n = 0, nmax - 1
      n1 = n + 1
      do 5 l = 1, lmax+n
         f2(l,n) = y * ( f1(l-1,n) - f1(l+1,n1) ) * dosl1p(l)
    5 continue
      return
      end
c
c
c
*============================================
      subroutine clln ( ll, nn, y2p4, sum )
*============================================
c
c  subrutina que calcula el elemento cl (ll, nn) con la formula directa
c
      implicit real * 8 ( a-h,o-z )
      c = 1.d0
      rn = nn + 1.d0
      sum = 1.d0 / rn
      rj = 0.d0
      rlj = ll + 0.5d0
      do 4 j = 1, 100
         rj = rj + 1.d0
         rlj = rlj + 1.d0
         c = c * y2p4 / ( rj * rlj )
         cc = c / ( rn + rj )
         sum = sum + cc
         if ( cc.lt.1.d-15*sum ) goto 5
    4 continue
      write (6,*) ' no converge cl '
    5 sum = sum * 0.5d0
      return
      end
c
c
c
*=================================================================
      subroutine jmult ( lamax, lbmax, lamin, lbmin , na , nb ,
     &                   xa , xb , a , cj )
*=================================================================
      implicit real * 8 ( a-h,o-z )
      include 'comun.inc'
      dimension cj(0:2*mxl,0:2*mxl,0:2*mxl)
      data pi/0.31415926535897932384626433832795059d+01/
      opi = 8.d0 * pi * pi
      do 1 m = 0 , min( lamax , lbmax )
      do 1 lp = lbmin , lbmax , 2
         if ( lp.ge.m ) then
            qb = fact(nb+1+lp)/xb**(nb+2+lp)
            do 2 l = lamin , lamax , 2
               if ( l.ge.m ) then
                  qa = fact(na+1+l)/xa**(na+2+l)
                  cj(l,lp,m) =opi*qa * qb * (-1)**(lp+m) * fact(l+lp) /
     &            ( a**(l+lp+1) * (2*l+1) * (2*lp+1)
     &            * fact(l-m) * fact(lp-m) )
                  if (m.eq.0) cj(l,lp,m) = cj(l,lp,m) * 2.d0
               endif
    2       continue
         endif
    1 continue

      return
      end
